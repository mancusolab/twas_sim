# twas_sim
Using real genotype data, simulate a complex trait as a function of latent expression, fit eQTL weights in independent data, and perform GWAS/TWAS on complex trait.

To download the TWAS simulator first type the commands

    git clone https://github.com/mancusolab/twas_sim.git
    cd twas_sim

then,

    conda env create --file environment.yml
    conda activate twas_sim

The script `example.sh` will generate a single TWAS statistic using the simulator `sim.py`. Please be sure to update the paths in `example.sh` first. The script relies on PLINK-formatted genotype data. We recommend downloading [1000G](https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz) for use. When you are done with the simulator be sure to enter the command

    conda deactivate

`sim.py` is the actual simulator. Its usage is below:

    usage: sim.py [-h] [--ngwas NGWAS] [--nqtl NQTL] [--model {10pct,1pct,1snp}]
                   [--linear-model {lasso, enet, ridge}] [--eqtl-h2 EQTL_H2]
                   [--var-explained VAR_EXPLAINED] [-o OUTPUT]
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
      --linear-model {lasso, enet, ridge}
                            Linear model to predict gene expression from genotype
      --eqtl-h2 EQTL_H2     The narrow-sense heritability of gene expression
                            (default: 0.1)
      --var-explained VAR_EXPLAINED
                            Variance explained in complex trait by gene expression
                            (default: 0.01)
      -o OUTPUT, --output OUTPUT
      --seed SEED           Seed for random number generation (default: None)

The output will be a two tab-delimited reports.

The first `OUTPUT.summary.tsv` is a high-level summary that contains two columns:

| stat             | values |
| ------           | ------ |
| ngwas            | GWAS sample size |
| nqtl             | eQTL sample size  |
| nsnps            | Number of SNPs |
| h2ge             | variance explained in trait by GE |
| h2g              | Narrow-sense heritability of GE |
| mean.ldsc        | Average LD-score at the region |
| min.gwas.p       | Minimum GWAS SNP p-value |
| mean.gwas.chi2   | Mean GWAS SNP chi-sq |
| median.gwas.chi2 | Median GWAS SNP chi-sq |
| twas.z           | TWAS Z score |
| twas.p           | TWAS p-value |
| alpha            | TWAS alpha |

The second `OUTPUT.scan.tsv` is individuals statistics at each SNP. It contains the following columns:

| column              | description |
| ------              | ----------  |
| chrom               | chromosome  |
| snp                 | snp identifier |
| pos                 | bp position |
| a0                  | non-effect allele |
| a1                  | effect allele |
| maf                 | minor allele frequency |
| ld.score            | ld score (ie. sum_i r_ij^2, where r_ij is LD between snps i, j) |
| gwas.beta           | beta coefficient in GWAS |
| gwas.se             | standard error in GWAS |
| gwas.true           | true causal effect for complex trait |
| eqtl.beta           | beta coefficient in eQTL |
| eqtl.se             | standard error in eQTL |
| eqtl.model          | linear model to predict gene expression from genotype |
| eqtl.model.beta     | coefficient estimated in selected linear model |
| eqtl.true           | true causal effect for expression |
