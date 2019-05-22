#! /usr/bin/env python
import argparse as ap
import functools
import os
import sys

import numpy as np


def main(args):
    argp = ap.ArgumentParser(description="Find a region constrained to have supplied number of genes.")
    argp.add_argument("loci", type=ap.FileType("r"))
    argp.add_argument("genes", type=ap.FileType("r"))
    argp.add_argument("-l", "--lgenes", type=int, default=10)
    argp.add_argument("-u", "--ugenes", type=int, default=30)
    argp.add_argument("-b", "--bound", type=float, default=5e5)
    argp.add_argument("-o", "--output", type=ap.FileType("w"), default=sys.stdout)
    argp.add_argument("--loc_output", type=ap.FileType("w"), default=sys.stdout)

    args = argp.parse_args(args)


    #1 	 10583 	 1892607
    loci = np.loadtxt(args.loci)
    n, _ = loci.shape

    # DDX11L1 1 11873 14409
    genes = np.loadtxt(args.genes, dtype=str)

    while True:
        locus = np.random.randint(0, n)
        CHR, TSS, TES = loci[locus]
        chr_flag = genes.T[1].astype(int) == CHR
        tss_flag = genes.T[2].astype(int) >= TSS
        tes_flag = genes.T[3].astype(int) <= TES
        flag = functools.reduce(np.logical_and, [chr_flag, tss_flag, tes_flag])
    
        sampled = genes[flag]

        if not (args.lgenes <= len(sampled) <= args.ugenes):
            continue

        idx = 0
        best = 0
        best_pair = None
        for jdx in range(len(sampled)):
            end = int(sampled[jdx, 3])
            start = int(sampled[idx, 2])
            if end - start < 1.5 * args.bound:
                continue
    
            if jdx - idx > best:
                best = jdx - idx 
                best_pair = (idx, jdx - 1)
    
            idx += 1
    
        if best_pair is None:
            fsampled = sampled
        else:
            fsampled = sampled[np.arange(best_pair[0], best_pair[1])]
        if len(fsampled) < args.lgenes:
            continue
        else:
            break

    for gene in fsampled:
        gene[1] = "chr" + gene[1]
        args.output.write(" ".join(gene) + os.linesep)

    start = int(fsampled[0, 2])
    end = int(fsampled[-1, 2])
    chr = fsampled[0, 1]
    out = " ".join(map(str, [chr, max(1, start - int(args.bound)), end + int(args.bound)]))
    args.loc_output.write(out + os.linesep)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
