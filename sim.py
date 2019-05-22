#! /usr/bin/env python
import argparse as ap
import os
import sys

import numpy as np
import limix.her as her
import pandas as pd

from numpy.linalg import multi_dot as mdot
from pandas_plink import read_plink
from scipy import stats
from sklearn import linear_model as lm

mvn = stats.multivariate_normal


def sim_trait(g, h2g):
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
    y /= np.std(y)

    return y


def sim_gwas(LD, ngwas, b_qtls, eqtl_h2, var_explained):
    Z_gwas = mvn.rvs(cov=LD, size=ngwas)
    Z_gwas -= np.mean(Z_gwas, axis=0)
    Z_gwas /= np.std(Z_gwas, axis=0)

    # var_explained should only reflect that due to genetics
    #gwas_expr = sim_trait(np.dot(Z_gwas, b_qtls), eqtl_h2)
    gwas_expr = np.dot(Z_gwas, b_qtls) 
    alpha = np.random.normal(loc=0, scale=1)
    y = sim_trait(gwas_expr * alpha, var_explained)

    betas = []
    ses = []
    pvals = []
    for snp in Z_gwas.T:
        beta, inter, rval, pval, se = stats.linregress(snp, y)
        betas.append(beta)
        ses.append(se)
        pvals.append(pval)

    gwas = pd.DataFrame({"beta":betas, "se":ses, "pval":pvals})

    return gwas


def fit_lasso(Z, y, h2g):
    n, p = Z.shape

    # PLINK-style LASSO
    lambda_max = np.linalg.norm(Z.T.dot(y), np.inf) / float(n)
    def _gen_e():
        e = np.random.normal(size=n)
        return np.linalg.norm(Z.T.dot(e), np.inf)

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

    return (coef, r2, logl)


def compute_twas(Z_qtl, y, A, LD_qtl, gwas):
    # fit predictive model using LASSO
    h2g = her.estimate(y, "normal", A, verbose=False)

    coef, r2, logl = fit_lasso(Z_qtl, y, h2g)
    Z = gwas.beta.values / gwas.se.values

    score = np.dot(coef, Z)
    within_var = mdot([coef, LD_qtl, coef])

    return score, within_var


def main(args):
    argp = ap.ArgumentParser(description="Simulate TWAS using real genotype data",
                             formatter_class=ap.ArgumentDefaultsHelpFormatter)
    argp.add_argument("prefix", help="Prefix to PLINK-formatted data")
    argp.add_argument("--ngwas", default=100000, type=int, help="Sample size for GWAS panel")
    argp.add_argument("--nqtl", default=500, type=int, help="Sample size for eQTL panel")
    argp.add_argument("--model", choices=["10pct", "1pct", "1snp"], default="10pct", help="SNP model for generating gene expression. 10pct = 10%% of SNPs, 1pct = 1%% of SNPs, 1snp = 1 SNP")
    argp.add_argument("--eqtl-h2", default=0.1, type=float, help="The narrow-sense heritability of gene expression")
    argp.add_argument("--var-explained", default=0.01, type=float, help="Variance explained in complex trait by gene expression")
    argp.add_argument("-o", "--output", type=ap.FileType("w"), default=sys.stdout)

    args = argp.parse_args(args)

    # read in plink data
    bim, fam, G = read_plink(args.prefix, verbose=False)

    # estimate LD for population from PLINK data
    n, p = map(float, G.shape)
    p_int = int(p)
    G -= np.mean(G, axis=0)
    G /= np.std(G, axis=0)
    LD = np.dot(G.T, G) / n + np.eye(p_int) * 0.1

    # simulate trait architecture
    mapper = {"10pct":0.1 * p, "1pct":0.01 * p, "1snp":1/p}
    n_qtls = max(int(mapper[args.model]), 1) # sometimes 1snp results in 0 due to floating point arithmetic combined with int function
    c_qtls = np.random.choice(range(p_int), n_qtls)
    b_qtls = np.zeros(p_int)
    b_qtls[c_qtls] = np.random.normal(loc=0, scale=np.sqrt(args.eqtl_h2 / n_qtls), size=n_qtls)

    # simulate GWAS under assumption that expression => downstream trait
    gwas = sim_gwas(LD, args.ngwas, b_qtls, args.eqtl_h2, args.var_explained)

    # sample eQTL reference pop genotypes from MVN approx
    Z_qtl = mvn.rvs(cov=LD, size=args.nqtl)
    Z_qtl -= np.mean(Z_qtl, axis=0)
    Z_qtl /= np.std(Z_qtl, axis=0)

    # GRM and LD
    A = np.dot(Z_qtl, Z_qtl.T) / p
    LD_qtl = np.dot(Z_qtl.T, Z_qtl) / n

    # simulate gene expression
    qtl_expr = sim_trait(np.dot(Z_qtl, b_qtls), args.eqtl_h2)

    # compute the score and variance for TWAS
    score, within_var = compute_twas(Z_qtl, qtl_expr, A, LD, gwas)

    min_p_val = np.min(gwas.pval.values)
    mean_chi2 = np.mean((gwas.beta.values / gwas.se.values) ** 2)
    med_chi2 = np.median((gwas.beta.values / gwas.se.values) ** 2)

    if within_var > 0:
        z_orig_twas = score / np.sqrt(within_var)
        p_orig_twas = 2 * stats.norm.sf(np.abs(z_orig_twas))
    else:
        # on underpowered/low-h2g genes LASSO can set all weights to 0 and effectively break the variance estimate 
        z_orig_twas = 0
        p_orig_twas = 1

    # format and write output
    df = pd.DataFrame({"stat":["ngwas", "nqtl", "h2ge", "h2g", "min_gwas_p", "mean_gwas_chi2", "median_gwas_chi2", "twas_orig_t", "twas_orig_p"],
                       "values":[args.ngwas, args.nqtl, args.var_explained, args.eqtl_h2, min_p_val, mean_chi2, med_chi2, z_orig_twas, p_orig_twas]})

    df.to_csv(args.output, sep="\t", index=False)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
