import subprocess

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
