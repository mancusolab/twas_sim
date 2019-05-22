# twas_sim
Using real genotype data, simulate a complex trait as a function of latent expression, fit eQTL weights in independent data, and perform GWAS/TWAS on complex trait.

To download the TWAS simulator first type the commands

    git clone https://github.com/mancusolab/twas_sim.git
    cd twas_sim

then,

    conda env create --file environment.yml
    conda activate twas_sim

The script `example.sh` will generate a single TWAS statistic using the simulator `sim.py`. Please be sure to update the paths in `example.sh` first. When you are done with the simulator be sure to enter the command

    conda deactivate

sim.py
    usage: sim.py [-h] [--ngwas NGWAS] [--nqtl NQTL] [--model {10pct,1pct,1snp}]
                  [--eqtl-h2 EQTL_H2] [--var-explained VAR_EXPLAINED] [-o OUTPUT]
                  prefix
    
    Simulate TWAS using real genotype data
    
    positional arguments:
      prefix                Prefix to PLINK-formatted data
    
    optional arguments:
      -h, --help            show this help message and exit
      --ngwas NGWAS         Sample size for GWAS panel (default: 100000)
      --nqtl NQTL           Sample size for eQTL panel (default: 500)
      --model {10pct,1pct,1snp}
                            SNP model for generating gene expression. 10pct = 10%
                            of SNPs, 1pct = 1% of SNPs, 1snp = 1 SNP (default:
                            10pct)
      --eqtl-h2 EQTL_H2     The narrow-sense heritability of gene expression
                            (default: 0.1)
      --var-explained VAR_EXPLAINED
                            Variance explained in complex trait by gene expression
                            (default: 0.01)
      -o OUTPUT, --output OUTPUT
