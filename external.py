#! /usr/bin/env python
import time
import argparse as ap
import re
import sys

import limix.her as her
import numpy as np
import pandas as pd
import scipy.linalg as linalg

from numpy.linalg import multi_dot as mdot
from pandas_plink import read_plink
from scipy import stats
from scipy.stats import invgamma
from sklearn import linear_model as lm

import subprocess

def external_module(args, Z, y, h2g, b_qtls=None):
    Z_out = f"{args.output}.eqtl.genotype.txt.gz" # Z_out = "eqtl.genotype.txt.gz", example: twas_sim_loci1.eqtl.genotype.txt.gz
    np.savetxt(Z_out,Z,fmt='%.5f')
    y_out = f"{args.output}.eqtl.gexpr.txt.gz" # y_out = "eqtl.gexpr.txt.gz"
    np.savetxt(y_out,y,fmt='%.5f')

    subprocess.call("Rscript external.R str(args.IDX)", shell=True)

    with open(f"{args.output}.external_coef.txt") as coef: # with open(path, "external_coef.txt") as coef:
        coef = coef.readlines()
    return coef, None, None
