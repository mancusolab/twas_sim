#! /usr/bin/env python
import argparse as ap
import importlib
import re
import sys
import time

import numpy as np
import pandas as pd
import scipy.linalg as linalg
from numpy.linalg import multi_dot as mdot
from pandas_plink import read_plink
from scipy import stats
from scipy.stats import invgamma
from sklearn import linear_model as lm
from glimix_core.lmm import LMM
from numpy_sugar.linalg import economic_qs


class NumCausalSNPs:
    """Helper class to keep track of the number of causal variants from the command line.
    Seems like overkill, but generalizes the internal logic to handle either:
        1) finite number of causal SNPs
        2) percentage of total SNPs at locus
        3) a random value sampled from truncated Poisson, with lower-bound of 1
    """

    def __init__(self, value):
        self._is_pct = False

        match = re.match("^([0-9]+)(pct|avg)?$", value, flags=re.IGNORECASE)
        if match:
            num_tmp = float(match.group(1))  # num
            num_mod = match.group(2)  # modifier
            if num_mod:
                num_mod = num_mod.upper()
                if num_mod == "PCT":
                    if not (0 < num_tmp <= 100):
                        raise ValueError("Percentage of causal SNPs must be in (0, 1].")
                    num_tmp /= 100.0
                    self._value = num_tmp
                    self._is_pct = True
                elif num_mod == "AVG":
                    if num_tmp == 0:
                        raise ValueError(
                            "Average number of causal SNPs must be at least 1"
                        )
                    num_smpl = 0
                    while num_smpl == 0:
                        # sample from poisson using average num of causals truncated below by 1
                        num_smpl = np.random.poisson(num_tmp)
                    self._value = int(num_smpl)
            else:
                if num_tmp < 1:
                    raise ValueError("Number of causal SNPs must be at least 1")
                self._value = int(num_tmp)
        else:
            raise ValueError("Invalid number of causal SNPs")

        return

    def __repr__(self):
        return str(self)

    def __str__(self):
        if self._is_pct:
            return f"NumCausals := {self._value * 100}% of observed SNPs"
        else:
            return f"NumCausals := {self._value} SNPs"

    def get(self, n_snps):
        if self._is_pct:
            return int(np.ceil(self._value * n_snps))
        else:
            # handle case where either avg number of causals or explicit number of causals
            # is greater than observed number
            return min(int(self._value), n_snps)


class NumCausalSNPsAction(ap.Action):
    """Custom action to parse the number of causal SNPs with percentage and average modifiers"""

    def __call__(self, parser, namespace, values, option_string=None):
        ncausals = NumCausalSNPs(values)
        setattr(namespace, self.dest, ncausals)
        return


def compute_s2g(L, beta):
    """
    Compute genetic variance given betas and LD cholesky factor

    s2g := beta' V beta = beta' L L' beta

    :param L: numpy.ndarray lower cholesky factor of the p x p LD matrix for the population
    :param beta: numpy.ndarray genotype effects

    :return: float s2g
    """
    Ltb = np.dot(L.T, beta)
    s2g = np.dot(Ltb.T, Ltb)

    return s2g


def fit_lasso(Z, y, h2g, b_qtls, args=None):
    """
    Infer eqtl coefficients using LASSO regression. Uses the PLINK-style coordinate descent algorithm
    that is bootstrapped by the current h2g estimate.

    :param Z: numpy.ndarray n x p genotype matrix
    :param y: numpy.ndarray gene expression for n individuals
    :param h2g: float the -estimated- h2g from reference panel
    :param b_qtls: numpy.ndarray latent eQTL effects for the causal gene

    :return: (numpy.ndarray, float, float) tuple of the LASSO coefficients, the r-squared score, and log-likelihood
    """
    return _fit_sparse_penalized_model(Z, y, h2g, lm.Lasso)


def fit_enet(Z, y, h2g, b_qtls, args=None):
    """
    Infer eqtl coefficients using ElasticNet regression. Uses the PLINK-style coordinate descent algorithm
    that is bootstrapped by the current h2g estimate.

    :param Z:  numpy.ndarray n x p genotype matrix
    :param y: numpy.ndarray gene expression for n individuals
    :param h2g: float the -estimated- h2g from reference panel
    :param b_qtls: numpy.ndarray latent eQTL effects for the causal gene

    :return: (numpy.ndarray, float, float) tuple of the ElasticNet coefficients, the r-squared score, and log-likelihood
    """
    return _fit_sparse_penalized_model(Z, y, h2g, lm.ElasticNet)


def fit_ridge(Z, y, h2g, b_qtls, args=None):
    """
    Infer eqtl coefficients using Ridge regression. Uses the optimal ridge penality defined from REML.

    :param Z:  numpy.ndarray n x p genotype matrix
    :param y: numpy.ndarray gene expression for n individuals
    :param h2g: float the -estimated- h2g from reference panel
    :param b_qtls: numpy.ndarray latent eQTL effects for the causal gene

    :return: (numpy.ndarray, float, float) tuple of the Ridge coefficients, the r-squared score, and log-likelihood
    """
    n, p = Z.shape
    lambda_r = (1 - h2g) / (h2g / p)

    model = lm.Ridge(alpha=lambda_r)
    model.fit(Z, y)
    coef, r2, logl = _get_model_info(model, Z, y)

    return coef, r2, logl


def _fit_sparse_penalized_model(Z, y, h2g, model_cls=lm.Lasso, args=None):
    """
    Infer eqtl coefficients using L1/L2 penalized regression. Uses the PLINK-style coordinate descent algorithm
    that is bootstrapped by the current h2g estimate.

    :param Z: numpy.ndarray n x p genotype matrix
    :param y: numpy.ndarray gene expression for n individuals
    :param h2g: float the -estimated- h2g from reference panel
    :param model_cls: linear_model from sklearn. Must be either Lasso or ElasticNet

    :return: (numpy.ndarray, float, float) tuple of the LASSO or ElasticNet coefficients, the r-squared score, and log-likelihood
    """
    if model_cls not in [lm.Lasso, lm.ElasticNet]:
        raise ValueError("penalized model must be either Lasso or ElasticNet")

    n, p = Z.shape

    def _gen_e():
        e = np.random.normal(size=n)
        return np.linalg.norm(Z.T.dot(e), np.inf)

    # PLINK-style LASSO
    lambda_max = np.linalg.norm(Z.T.dot(y), np.inf) / float(n)

    min_tmp = np.median([_gen_e() for _ in range(1000)])
    sige = np.sqrt(1.0 - h2g + (1.0 / float(n)))
    lambda_min = (sige / n) * min_tmp

    # 100 values spaced logarithmically from lambda-min to lambda-max
    alphas = np.exp(np.linspace(np.log(lambda_min), np.log(lambda_max), 100))

    # fit solution using coordinate descent, updating with consecutively smaller penalties
    model = model_cls(fit_intercept=True, warm_start=True)
    for penalty in reversed(alphas):
        model.set_params(alpha=penalty)
        model.fit(Z, y)

    coef, r2, logl = _get_model_info(model, Z, y)

    return coef, r2, logl


def fit_external(Z, y, h2g, b_qtls=None, args=None):
    try:
        load_mod = importlib.import_module(args.external_module)
        coef, r2, logl = load_mod.fit(Z, y, h2g, b_qtls, args)
    except ModuleNotFoundError as ex:
        raise ModuleNotFoundError(
            f"ERROR! Cannot find external module {args.external_module}!"
        ) from ex
    except AttributeError as ex:
        raise AttributeError(
            f"ERROR! Cannot find 'fit' function in module {args.external_module}!"
        ) from ex
    except TypeError as ex:
        raise TypeError(
            (
                f"ERROR! Argument signature for 'fit' in module {args.external_module} ",
                "should be 'fit(Z, y, h2g, b_qtls, args)'!",
            )
        ) from ex

    return coef, r2, logl


def fit_trueqtl(Z, y, h2g, b_qtls, args=None):
    """
    Return true latent eQTL effects for the causal gene.

    :param Z:  numpy.ndarray n x p genotype matrix
    :param y: numpy.ndarray gene expression for n individuals
    :param h2g: float the -estimated- h2g from reference panel
    :param b_qtls: numpy.ndarray latent eQTL effects for the causal gene
    """

    return b_qtls, None, None


def _get_model_info(model, Z, y):
    """
    Helper function to get fitted coefficients, R2, and log-likelihood
    """
    n, p = Z.shape
    coef = model.coef_

    r2 = model.score(Z, y)
    ystar = model.predict(Z)
    s2e = sum((y - ystar) ** 2) / (n - 1)

    logl = sum(stats.norm.logpdf(y, loc=ystar, scale=np.sqrt(s2e)))

    return coef, r2, logl


def sim_beta(L, ncausal, eqtl_h2, rescale=True):
    """
    Sample qtl effects under a specified architecture.

    :param L: numpy.ndarray lower cholesky factor of the p x p LD matrix for the population
    :param ncausal: NumCausalSNPs class containing number of causal SNPs to select.
        Encapsulates diff logic for integer, percentage, and averages.
    :param eqtl_h2: float the heritability of gene expression
    :param rescale: bool whether to rescale effects such that var(b V b) = h2 (default=True)

    :return: numpy.ndarray of causal effects
    """

    n_snps = L.shape[0]
    
    if eqtl_h2 != 0:
        # choose how many eQTLs
        n_qtls = ncausal.get(n_snps)

        # select which SNPs are causal
        c_qtls = np.random.choice(range(int(n_snps)), n_qtls)
        b_qtls = np.zeros(int(n_snps))

        # sample effects from normal prior
        b_qtls[c_qtls] = np.random.normal(
            loc=0, scale=np.sqrt(eqtl_h2 / n_qtls), size=n_qtls
        )
        if rescale:
            s2g = compute_s2g(L, b_qtls)
            b_qtls *= np.sqrt(eqtl_h2 / s2g)
    else:
        b_qtls = np.zeros(int(n_snps))

    return b_qtls


def sim_trait(g, h2g):
    """
    Simulate a complex trait as a function of latent genetic value and env noise.
    :param g: numpy.ndarray of latent genetic values
    :param h2g: float the heritability of the trait in the population

    :return: (numpy.ndarray, float) simulated phenotype, sd of Y
    """
    n = len(g)

    if h2g > 0:
        s2g = np.var(g, ddof=1)
        s2e = s2g * ((1.0 / h2g) - 1)
        e = np.random.normal(0, np.sqrt(s2e), n)
        y = g + e
    else:
        e = np.random.normal(0, 1, n)
        y = e

    # standardize
    y -= np.mean(y)

    return y


def sim_geno(L, n):
    """
    Sample genotypes from an MVN approximation.

    :param L: numpy.ndarray lower cholesky factor of the p x p LD matrix for the population
    :param n: int the number of genotypes to sample

    :return: numpy.ndarray n x p centered/scaled genotype matrix
    """
    p, p = L.shape

    Z = L.dot(np.random.normal(size=(n, p)).T).T
    Z -= np.mean(Z, axis=0)
    Z /= np.std(Z, axis=0)

    return Z


def regress(Z, pheno):
    """
    Perform a marginal linear regression for each snp on the phenotype.

    :param Z: numpy.ndarray n x p genotype matrix to regress over
    :param pheno: numpy.ndarray phenotype vector

    :return: pandas.DataFrame containing estimated beta and standard error
    """
    betas = []
    ses = []
    pvals = []
    for snp in Z.T:
        beta, inter, rval, pval, se = stats.linregress(snp, pheno)
        betas.append(beta)
        ses.append(se)
        pvals.append(pval)

    gwas = pd.DataFrame({"beta": betas, "se": ses, "pval": pvals})

    return gwas

def estimate_her(Z, y):
    """
    Calculate proportion of expression variation explained by genotypes (cis-heritability; :math:`h_g^2`).
    """
    n, p = Z.shape

    Z -= np.mean(Z, axis=0)
    Z /= np.std(Z, axis=0)

    y -= np.mean(y)
    y /= np.std(y)

    covar = np.ones(n)

    GRM = np.dot(Z, Z.T) / p
    # normalize the covariance matrix as suggested by Limix
    # https://horta-limix.readthedocs.io/en/api/_modules/limix/her/_estimate.html#estimate
    # and https://horta-limix.readthedocs.io/en/api/_modules/limix/qc/kinship.html#normalise_covariance
    # here, we calculate GRM using p, instead of p-1, so jnp.diag.mean should be equivalent to jnp.trace/(p-1)
    GRM /= jnp.diag(GRM).mean()
    QS = economic_qs(GRM)
    method = LMM(y, covar, QS, restricted=True)
    method.fit(verbose=False)

    g = method.scale * (1 - method.delta)
    e = method.scale * method.delta
    v = np.var(method.mean())
    h2g = g / (v + g + e)

    return h2g

def sim_gwasfast(L, ngwas, beta, h2ge):
    """
    Simulate a GWAS using `ngwas` individuals such that genetics explain `h2ge` of phenotype.
        This function differs from `sim_gwas` in that it samples GWAS summary statistics directly
        Using an MVN approximation, rather than generating genotype, phenotype, and performing
        a marginal regression for each simulated variant. Runtime should be `O(p^2)` where
        `p` is number of variants.

    :param L: numpy.ndarray lower cholesky factor of the p x p LD matrix for the population
    :param ngwas: int the number of GWAS genotypes to sample
    :param b_qtls: numpy.ndarray latent eQTL effects for the causal gene
    :param h2ge: float the amount of phenotypic variance explained by genetic component of gene expression

    :return: (pandas.DataFrame) estimated GWAS beta and standard error

    """

    n_snps = L.shape[0]

    s2g = compute_s2g(L, beta)

    if h2ge > 0:
        s2e = s2g * ((1.0 / h2ge) - 1)
    else:
        s2e = 1.0  # var[y]; could be diff from 1, but here we assume 1

    dof = ngwas - 1
    tau2 = s2e / ngwas
    se_gwas = np.sqrt(invgamma.rvs(a=0.5 * dof, scale=0.5 * dof * tau2, size=n_snps))
    DL = se_gwas[:, np.newaxis] * L

    beta_adj = mdot([DL, L.T, (beta / se_gwas)])  # D @ L @ Lt @ inv(D) @ beta

    # b_gwas ~ N(D @ L @ L.t inv(D) @ beta, D @ L @ Lt @ D), but fast
    b_gwas = beta_adj + np.dot(DL, np.random.normal(size=(n_snps,)))

    Z = b_gwas / se_gwas
    pvals = 2 * stats.norm.sf(abs(Z))

    gwas = pd.DataFrame({"beta": b_gwas, "se": se_gwas, "pval": pvals})

    return gwas


def sim_gwas(L, ngwas, beta, h2ge):
    """
    Simulate a GWAS using `ngwas` individuals such that genetics explain `h2ge` of phenotype.
        This function approximates genotypes under an LD structure using an MVN model. Generating
        genotype for `ngwas` individuals takes `O(np^2)` time. Simulating a phenotype and
        performing marginal regression for each variant takes `O(np)` time. If `n > p`, which is
        typically the case for GWAS (at a fixed region), this function should be slower than `sim_gwasfast`,
        which requires `O(p^2)` time to simulate GWAS summary results.

    :param L: numpy.ndarray lower cholesky factor of the p x p LD matrix for the population
    :param ngwas: int the number of GWAS genotypes to sample
    :param b_qtls: numpy.ndarray latent eQTL effects for the causal gene
    :param h2ge: float the amount of phenotypic variance explained by genetic component of gene expression

    :return: (pandas.DataFrame) estimated GWAS beta and standard error
    """
    Z_gwas = sim_geno(L, ngwas)

    # h2ge should only reflect that due to genetics
    g = np.dot(Z_gwas, beta)
    y = sim_trait(g, h2ge)

    gwas = regress(Z_gwas, y)

    return gwas


def sim_eqtl(L, b_qtls, args):
    """
    Simulate an eQLT study using `nqtl` individuals.

    :param L: numpy.ndarray lower cholesky factor of the p x p LD matrix for the population
    :param nqtl: int the number of eQTL-panel genotypes to sample
    :param b_qtls: numpy.ndarray latent eQTL effects for the causal gene
    :param eqtl_h2: float the amount of expression variance explained by linear model of SNPs
    :param linear_model: str the name of linear model to fit gene expression on genotype

    :return:  (pandas.DataFrame, numpy.ndarray, float) DataFrame of eQTL scan, vector of fitted eQTL coefficients,
        and estimated h2g.
    """
    nqtl = args.nqtl
    eqtl_h2 = args.eqtl_h2
    linear_model = args.linear_model

    Z_qtl = sim_geno(L, nqtl)
    n, p = [float(x) for x in Z_qtl.shape]

    # simulate gene expression
    gexpr = sim_trait(np.dot(Z_qtl, b_qtls), eqtl_h2)

    # get marginal eQTLs for reporting
    eqtl = regress(Z_qtl, gexpr)

    # fit predictive model
    h2g = estimate_her(Z_qtl, gexpr)

    # sample eQTL reference pop genotypes from MVN approx and perform eQTL scan + fit penalized linear model
    if linear_model == "lasso":
        pred_func = fit_lasso
    elif linear_model == "enet":
        pred_func = fit_enet
    elif linear_model == "ridge":
        pred_func = fit_ridge
    elif linear_model == "trueqtl":
        pred_func = fit_trueqtl
    elif linear_model == "external":
        pred_func = fit_external
    else:
        raise ValueError("Invalid linear model")

    # fit penalized to get predictive weights
    coef, r2, logl = pred_func(Z_qtl, gexpr, h2g, b_qtls, args)

    return gexpr, eqtl, coef, h2g, Z_qtl


def compute_twas(gwas, coef, LD):
    """
    Compute the TWAS test statistics.

    :param gwas: pandas.DataFrame containing estimated GWAS beta and standard error
    :param coef: numpy.ndarray LASSO eQTL coefficients
    :param LD:  numpy.ndarray p x p LD matrix

    :return: (float, float) the TWAS test statistics and p-value
    """
    # compute Z scores
    Z = gwas.beta.values / gwas.se.values

    # score and variance
    score = np.dot(coef, Z)
    within_var = mdot([coef, LD, coef])

    if within_var > 0:
        z_twas = score / np.sqrt(within_var)
        p_twas = 2 * stats.norm.sf(np.abs(z_twas))
    else:
        # on underpowered/low-h2g genes LASSO can set all weights to 0 and effectively break the variance estimate
        z_twas = 0
        p_twas = 1

    return z_twas, p_twas


def main(args):
    argp = ap.ArgumentParser(
        description="Simulate TWAS using real genotype data",
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )
    argp.add_argument(
        "prefix", help="Prefix to PLINK-formatted data for GWAS LD information"
    )
    argp.add_argument(
        "--eqtl-prefix",
        default=None,
        help="Optional prefix to PLINK-formatted data for eQTL LD information. Otherwise use GWAS LD.",
    )
    argp.add_argument(
        "--test-prefix",
        default=None,
        help="Optional prefix to PLINK-formatted data for LD information in TWAS test statistic. Otherwise use GWAS LD.",
    )
    argp.add_argument(
        "--fast-gwas-sim",
        default=False,
        action="store_true",
        help="If set then simulate GWAS summary data directly from LD",
    )
    argp.add_argument(
        "--ngwas", default=100000, type=int, help="Sample size for GWAS panel"
    )
    argp.add_argument(
        "--nqtl", default=500, type=int, help="Sample size for eQTL panel"
    )
    argp.add_argument("--IDX", type=int, help="Simulation index")
    argp.add_argument(
        "--ncausal",
        default="1",
        action=NumCausalSNPsAction,
        help=(
            "Number of causal SNPs for gene expression/trait. Can represent explicit number (e.g., 1, 10),"
            " a percentage using the 'pct' modifier (e.g., '1pct', '10pct'),"
            " or an average under a truncated Poisson model (e.g., '1avg', '10avg')."
        ),
    )
    argp.add_argument(
        "--ld-ridge", default=0.1, type=float, help="Offset to add to LD Diagonal"
    )
    argp.add_argument(
        "--linear-model",
        choices=["lasso", "enet", "ridge", "trueqtl", "external"],
        default="lasso",
        help=(
            "Linear model to predict gene expression from genotype. "
            "Use external to indicate an external module should be loaded."
        ),
    )
    argp.add_argument(
        "--external-module",
        help=(
            "Path to external Python file with custom `fit` function. Only used if `--linear-module=external`. "
            "E.g., if `my_module.py` contains `fit function then pass in `my_module`."
        ),
    )
    argp.add_argument(
        "--eqtl-h2",
        default=0.1,
        type=float,
        help="The SNP heritability of gene expression.",
    )
    argp.add_argument(
        "--h2ge",
        default=0.01,
        type=float,
        help="Phenotypic variance explained by genetic component of gene expression,",
    )
    argp.add_argument(
        "--indep-gwas",
        default=False,
        action="store_true",
        help="Generate GWAS effect-sizes independently from eQTLs.",
    )
    argp.add_argument(
        "--h2g-gwas",
        default=0.01,
        type=float,
        help="The SNP heritability of downstream phenotype. Only used when `--indep-gwas` is set.",
    )
    argp.add_argument("-o", "--output", help="Output prefix")
    argp.add_argument(
        "-c",
        "--compress",
        action="store_true",
        default=False,
        help="Compress output (gzip)",
    )

    argp.add_argument("--seed", type=int, help="Seed for random number generation")

    args = argp.parse_args(args)

    def get_ld(prefix):
        # return cholesky L and ldscs
        bim, fam, G = read_plink(prefix, verbose=False)
        G = G.T

        # estimate LD for population from PLINK data
        n, p = [float(x) for x in G.shape]
        p_int = int(p)
        mafs = (np.mean(G, axis=0) / 2).compute()
        G -= mafs * 2
        G /= np.std(G, axis=0)

        # regularize so that LD is PSD
        LD = np.dot(G.T, G) / n + np.eye(p_int) * args.ld_ridge

        # re-adjust to get proper correlation matrix
        LD = LD / (1 + args.ld_ridge)

        # compute cholesky decomp for faster sampling/simulation
        L = linalg.cholesky(LD, lower=True)

        # compute LD-scores for reports
        # weird dask issues require us to call compute here
        ldscs = np.sum(LD**2, axis=0).compute()

        return (L, mafs, ldscs, bim)

    np.random.seed(args.seed)
    real_time_start = time.time()
    cpu_time_start = time.process_time()

    # compute GWAS LD information from reference genotype data
    L_pop, mafs, ldscs, bim = get_ld(args.prefix)
    pop_p = len(ldscs)

    if args.h2ge > 0 and not args.indep_gwas:
        # we don't need to sample since alpha is determined by h2 and h2ge
        # and we've already normalized b_qtls to be on h2g scale [rescale=True above]
        sign = np.random.choice([-1, 1])
        alpha = np.sqrt(args.h2ge / args.eqtl_h2) * sign
    else:
        alpha = 0.0

    # simulate eQTL effects
    b_qtls = sim_beta(L_pop, args.ncausal, args.eqtl_h2, rescale=True)

    # simulate downstream trait/pheno effects
    if args.indep_gwas:
        # downstream trait/pheno effects are independent from eQTL
        beta = sim_beta(L_pop, args.ncausal, args.h2g_gwas, rescale=True)
    else:
        # downstream trait/pheno effects are mediated by eQTL
        beta = b_qtls * alpha

    # determine LD to use for eQTL ref panel
    if args.eqtl_prefix is not None:
        L_eqtl, mafs_eqtl, ldscs_eqtl, bim_eqtl = get_ld(args.eqtl_prefix)
        eqtl_p = L_eqtl.shape[0]
        if eqtl_p != pop_p:
            raise ValueError("The number of SNPs in eQTL LD must match GWAS LD")
    else:
        L_eqtl = L_pop

    # fit prediction model on simulated eQTL ref panel and grab eQTL coefficients
    gexpr, eqtl, coef, eqtl_h2g_hat, Z_qtl = sim_eqtl(L_eqtl, b_qtls, args)

    # determine LD to use for TWAS test
    if args.test_prefix is not None:
        L_test, mafs_test, ldscs_test, bim_test = get_ld(args.test_prefix)
        test_p = L_test.shape[0]
        if test_p != pop_p:
            raise ValueError(
                "The number of SNPs in TWAS testing ref panel must match GWAS LD"
            )
    else:
        L_test = L_pop

    if args.fast_gwas_sim:
        sim_mode, sim_func = ("fast", sim_gwasfast)
    else:
        sim_mode, sim_func = ("std", sim_gwas)

    gwas = sim_func(L_pop, args.ngwas, beta, args.h2ge)

    # compute some global GWAS summary statistics
    min_p_val = np.min(gwas.pval.values)
    mean_chi2 = np.mean((gwas.beta.values / gwas.se.values) ** 2)
    med_chi2 = np.median((gwas.beta.values / gwas.se.values) ** 2)

    # compute LD for the TWAS test
    LD_test = np.dot(L_test, L_test.T)

    # compute ldscore -at- the causals
    causals = b_qtls != 0
    if np.sum(causals) > 0:
        ldsc_causals = np.sum(LD_test[:, causals] ** 2, axis=1)
    else:
        ldsc_causals = np.zeros(pop_p)

    # compute TWAS statistics
    z_twas, p_twas = compute_twas(gwas, coef, LD_test)

    # compute real time and cpu time
    real_time_end = time.time()
    real_time = round(real_time_end - real_time_start, 2)
    cpu_time_end = time.process_time()
    cpu_time = round(cpu_time_end - cpu_time_start, 2)

    # output the GWAS, eQTL, and LASSO estimates
    output = bim.drop(columns=["cm", "i"])
    output["maf"] = mafs
    output["ld.score"] = ldscs
    output["ld.score.causal"] = ldsc_causals
    output["gwas.sim"] = [sim_mode] * len(mafs)
    output["gwas.true"] = b_qtls * alpha
    output["gwas.beta"] = gwas.beta
    output["gwas.se"] = gwas.se
    output["eqtl.true"] = b_qtls
    output["eqtl.beta"] = eqtl.beta
    output["eqtl.se"] = eqtl.se
    output["eqtl.model"] = [args.linear_model] * len(mafs)
    output["eqtl.model.beta"] = coef

    # output a summary that contains the actual TWAS test statistic
    df = pd.DataFrame(
        {
            "gwas.sim": [sim_mode],
            "real.time": [real_time],
            "cpu.time": [cpu_time],
            "linear_model": [args.linear_model],
            "h2ge": [args.h2ge],
            "snp_model": [args.ncausal],
            "nsnps": [int(pop_p)],
            "ngwas": [args.ngwas],
            "nqtl": [args.nqtl],
            "h2g": [args.eqtl_h2],
            "h2g.hat": [eqtl_h2g_hat],
            "avg.ldsc": [np.mean(ldscs)],
            "min.gwas.p": [min_p_val],
            "mean.gwas.chi2": [mean_chi2],
            "median.gwas.chi2": [med_chi2],
            "twas.z": [z_twas],
            "twas.p": [p_twas],
            "alpha": [alpha],
        }
    )

    scan_out = f"{args.output}.scan.tsv"
    summary_out = f"{args.output}.summary.tsv"

    if args.compress:
        scan_out += ".gz"
        summary_out += ".gz"

    output.to_csv(scan_out, sep="\t", index=False)
    df.to_csv(summary_out, sep="\t", index=False)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
                                        