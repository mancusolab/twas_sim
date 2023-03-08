import subprocess

import numpy as np


def fit(Z, y, h2g, b_qtls=None, args=None):

    np.savetxt(f"{args.output}.eqtl.genotype.txt.gz", Z, fmt="%.5f")
    np.savetxt(f"{args.output}.eqtl.gexpr.txt.gz", y, fmt="%.5f")

    subprocess.run(
        [
            "Rscript",
            "external.R",
            f"{args.output}.eqtl.genotype.txt.gz",
            f"{args.output}.eqtl.gexpr.txt.gz",
            f"{args.output}.susie.coef.txt.gz",
        ],
        shell=True,
    )

    coef = np.loadtxt(f"{args.output}.susie.coef.txt.gz")

    return coef, None, None
