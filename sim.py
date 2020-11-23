#! /usr/bin/env python
import argparse as ap
import sys

import numpy as np
import limix.her as her
import pandas as pd
import scipy.linalg as linalg

from numpy.linalg import multi_dot as mdot
from pandas_plink import read_plink
from scipy import stats
from sklearn import linear_model as lm

mvn = stats.multivariate_normal


def fit_lasso(Z, y, h2g):
    """
    Infer eqtl coefficients using LASSO regression. Uses the PLINK-style coordinate descent algorithm
    that is bootstrapped by the current h2g estimate.

    :param Z:  numpy.ndarray n x p genotype matrix
    :param y: numpy.ndarray gene expression for n individuals
    :param h2g: float the -estimated- h2g from reference panel

    :return: (numpy.ndarray, float, float) tuple of the LASSO coefficients, the r-squared score, and log-likelihood
    """
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

    # fit LASSO solution using coordinate descent, updating with consecutively smaller penalties
    lasso = lm.Lasso(fit_intercept=True, warm_start=True)
    for penalty in reversed(alphas):
        lasso.set_params(alpha=penalty)
        lasso.fit(Z, y)

    coef = lasso.coef_

    r2 = lasso.score(Z, y)
    ystar = lasso.predict(Z)
    s2e = sum((y - ystar) ** 2) / (n - 1)

    logl = sum(stats.norm.logpdf(y, loc=ystar, scale=np.sqrt(s2e)))

    return coef, r2, logl


def sim_beta(model, eqtl_h2, n_snps):
    """
    Sample qtl effects under a specified architecture.

    :param model: str the model to simulate under. choices="10pct", "1pct", "1"
    :param eqtl_h2: float the heritability of gene expression
    :param n_snps: the total number of snps at the region

    :return: numpy.ndarray of causal effects
    """
    # simulate trait architecture
    mapper = {"10pct": 0.1 * n_snps, "1pct": 0.01 * n_snps, "1snp": 1}
    n_qtls = int(mapper[model])

    # select which SNPs are causal
    c_qtls = np.random.choice(range(int(n_snps)), n_qtls)
    b_qtls = np.zeros(int(n_snps))

    # sample effects from normal prior
    b_qtls[c_qtls] = np.random.normal(loc=0, scale=np.sqrt(eqtl_h2 / n_qtls), size=n_qtls)

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
        s2e = s2g * ( (1.0 / h2g ) - 1 )
        e = np.random.normal(0, np.sqrt(s2e), n)
        y = g + e
    else:
        e = np.random.normal(0, 1, n)
        y = e

    # standardize
    y -= np.mean(y)
    y_std = np.std(y)
    y /= y_std
    y_std = y_std.item()

    return y, y_std


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

    gwas = pd.DataFrame({"beta":betas, "se":ses, "pval":pvals})

    return gwas


def sim_gwas(L, ngwas, b_qtls, var_explained):
    """
    Simulate a GWAS using `ngwas` individuals such that genetics explain `var_explained` of phenotype.

    :param L: numpy.ndarray lower cholesky factor of the p x p LD matrix for the population
    :param ngwas: int the number of GWAS genotypes to sample
    :param b_qtls: numpy.ndarray latent eQTL effects for the causal gene
    :param var_explained: float the amount of phenotypic variance explained by genetic component of gene expression

    :return: (pandas.DataFrame, float) estimated GWAS beta and standard error, causal GE effect
    """
    Z_gwas = sim_geno(L, ngwas)

    # var_explained should only reflect that due to genetics
    gwas_expr = np.dot(Z_gwas, b_qtls)
    if var_explained > 0:
        alpha = np.random.normal(loc=0, scale=1)
    else:
        alpha = 0
    y, y_std = sim_trait(gwas_expr * alpha, var_explained)

    gwas = regress(Z_gwas, y)

    # correct alpha for original SD of y
    alpha /= y_std

    return (gwas, alpha)


def sim_eqtl(L, nqtl, b_qtls, eqtl_h2):
    """
    Simulate an eQLT study using `nqtl` individuals.

    :param L: numpy.ndarray lower cholesky factor of the p x p LD matrix for the population
    :param nqtl: int the number of eQTL-panel genotypes to sample
    :param b_qtls: numpy.ndarray latent eQTL effects for the causal gene
    :param eqtl_h2: float the amount of expression variance explained by linear model of SNPs

    :return:  (pandas.DataFrame, numpy.ndarray, numpy.ndarray, float) DataFrame of eQTL scan, vector of LASSO eQTL coefficients,
        LD estimated from eQTL reference panel, and original gene expression SD.
    """
    Z_qtl = sim_geno(L, nqtl)
    n, p = [float(x) for x in  Z_qtl.shape]

    # GRM and LD
    A = np.dot(Z_qtl, Z_qtl.T) / p
    LD_qtl = np.dot(Z_qtl.T, Z_qtl) / n

    # simulate gene expression
    gexpr, gexpr_std = sim_trait(np.dot(Z_qtl, b_qtls), eqtl_h2)

    # get marginal eQTLs for reporting
    eqtl = regress(Z_qtl, gexpr)

    # fit predictive model using LASSO
    h2g = her.estimate(gexpr, "normal", A, verbose=False)

    # fit LASSO to get predictive weights
    coef, r2, logl = fit_lasso(Z_qtl, gexpr, h2g)

    return (eqtl, coef, LD_qtl, gexpr_std)


def compute_twas(gwas, coef, LD):
    """
    Compute the TWAS test statistics.

    :param gwas: pandas.DataFrame containing estimated GWAS beta and standard error
    :param coef: numpy.ndarray LASSO eQTL coefficients
    :param LD:  numpy.ndarray p x p LD matrix

    :return: (float, float) the TWAS score and variance estimates
    """
    # compute Z scores
    Z = gwas.beta.values / gwas.se.values

    # score and variance
    score = np.dot(coef, Z)
    within_var = mdot([coef, LD, coef])

    return score, within_var


def main(args):
    argp = ap.ArgumentParser(description="Simulate TWAS using real genotype data",
                             formatter_class=ap.ArgumentDefaultsHelpFormatter)
    argp.add_argument("prefix",
                      help="Prefix to PLINK-formatted data")

    argp.add_argument("--ngwas", default=100000, type=int, help="Sample size for GWAS panel")
    argp.add_argument("--nqtl", default=500, type=int, help="Sample size for eQTL panel")
    argp.add_argument("--model", choices=["10pct", "1pct", "1snp"], default="10pct",
                      help="SNP model for generating gene expression. 10pct = 10%% of SNPs, 1pct = 1%% of SNPs, 1snp = 1 SNP")
    argp.add_argument("--eqtl-h2", default=0.1, type=float, help="The narrow-sense heritability of gene expression")
    argp.add_argument("--var-explained", default=0.01, type=float,
                      help="Variance explained in complex trait by gene expression")
    argp.add_argument("-o", "--output", help="Output prefix")
    argp.add_argument("--seed", type=int, help="Seed for random number generation")

    args = argp.parse_args(args)

    # read in plink data
    bim, fam, G = read_plink(args.prefix, verbose=False)
    G = G.T 

    # estimate LD for population from PLINK data
    n, p = [float(x) for x in G.shape]
    p_int = int(p)
    mafs = np.mean(G, axis=0) / 2
    G -= mafs * 2
    G /= np.std(G, axis=0)

    # regularize so that LD is PSD
    LD = np.dot(G.T, G) / n + np.eye(p_int) * 0.1

    # compute cholesky decomp for faster sampling/simulation
    L = linalg.cholesky(LD, lower=True)

    # compute LD-scores for reports
    ldscs = np.sum(LD ** 2, axis=0)

    # set random seed from argument for simulation reproducibity
    np.random.seed(args.seed)
    
    b_qtls = sim_beta(args.model, args.eqtl_h2, p)

    # simulate GWAS under assumption that expression => downstream trait
    gwas, alpha = sim_gwas(L, args.ngwas, b_qtls, args.var_explained)

    # sample eQTL reference pop genotypes from MVN approx and perform eQTL scan + fit LASSO
    eqtl, coef, LD_qtl, gexpr_std = sim_eqtl(L, args.nqtl, b_qtls, args.eqtl_h2)

    # compute TWAS statistics
    score, within_var = compute_twas(gwas, coef, LD)

    min_p_val = np.min(gwas.pval.values)
    mean_chi2 = np.mean((gwas.beta.values / gwas.se.values) ** 2)
    med_chi2 = np.median((gwas.beta.values / gwas.se.values) ** 2)

    if within_var > 0:
        z_twas = score / np.sqrt(within_var)
        p_twas = 2 * stats.norm.sf(np.abs(z_twas))
    else:
        # on underpowered/low-h2g genes LASSO can set all weights to 0 and effectively break the variance estimate 
        z_twas = 0
        p_twas = 1

    # output the GWAS, eQTL, and LASSO estimates
    output = bim.drop(columns=["cm", "i"])
    output["maf"] = mafs
    output["ld.score"] = ldscs
    output["gwas.beta"] = gwas.beta
    output["gwas.se"] = gwas.se
    output["gwas.true"] = b_qtls * alpha
    output["eqtl.beta"] = eqtl.beta
    output["eqtl.se"] = eqtl.se
    output["eqtl.true"] = b_qtls / gexpr_std
    output["eqtl.lasso"] = coef
    output.to_csv("{}.scan.tsv".format(args.output), sep="\t", index=False)

    # correct alpha for original SD of gene expression 
    # (this needs to come after construction of gwas.true above)
    alpha *= gexpr_std

    # output a summary that contains the actual TWAS test statistic
    df = pd.DataFrame({"stat": ["ngwas", "nqtl", "nsnps", "h2ge", "h2g", "avg.ldsc",
                                "min.gwas.p", "mean.gwas.chi2", "median.gwas.chi2", "twas.z", "twas.p",
                                "alpha"],
                       "values": [args.ngwas, args.nqtl, int(p), args.var_explained, args.eqtl_h2, np.mean(ldscs),
                                  min_p_val, mean_chi2, med_chi2, z_twas, p_twas, alpha]})
    df.to_csv("{}.summary.tsv".format(args.output), sep="\t", index=False)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
