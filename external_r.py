"""
This file is an example of how to define an external/custom function to fit a
predictive model of gene expression from genotype to be used by `twas_sim`. Here
we are illustrating how to call an external R script to call susieR on the simulated
data.

External modules -must- include a function named `fit` that takes as arguments:
    Z: numpy matrix of genotype
    y: numpy vector of gene expression/phenotype
    h2g: the true h2g of gene expression
    b_qtls: the true beta/effect-sizes for gene expression (i.e. eQTL)
    args: the argparse object from twas_sim; useful for pulling `args.output`
        as a prefix for temp data.

Similarly, it must return a tuple containing (coef, r2, logl):
    coef: the numpy vector for estimated eQTL weights
    r2: the predictive r2 (optional; None)
    logl: the log likelihood of the model (optional; None)

"""
import numpy as np


def fit(Z, y, h2g, b_qtls=None, args=None):
    geno_path = f"{args.output}.eqtl.genotype.txt.gz"
    pheno_path = f"{args.output}.eqtl.gexpr.txt.gz"
    coef_path = f"{args.output}.susie.coef.txt.gz"

    np.savetxt(geno_path, Z, fmt="%.5f")
    np.savetxt(pheno_path, y, fmt="%.5f")

    subprocess.run(
        f"~/miniconda3/bin/Rscript external.R {geno_path} {pheno_path} {coef_path}",
        shell=True,
        check=True,
    )

    coef = np.loadtxt(coef_path)

    return coef, None, None
