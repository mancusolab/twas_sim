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

<<<<<<< HEAD
    usage: sim.py [-h] [--eqtl-prefix EQTL_PREFIX] [--test-prefix TEST_PREFIX]
                    [--fast-gwas-sim FAST_GWAS_SIM] [--ngwas NGWAS] [--nqtl NQTL]
                    [--model {10pct,1pct,1snp}] [--ld-ridge LD_RIDGE]
                    [--linear-model {lasso,enet,ridge}] [--eqtl-h2 EQTL_H2]
                    [--var-explained VAR_EXPLAINED] [-o OUTPUT] [--seed SEED]
                    [--sim SIM] [--locus LOCUS]
                  prefix
=======
    usage: sim.py [-h] [--eqtl-prefix EQTL_PREFIX] [--test-prefix TEST_PREFIX] [--fast-gwas-sim]
              [--ngwas NGWAS] [--nqtl NQTL] [--ncausal NCAUSAL] [--ld-ridge LD_RIDGE]
              [--linear-model {lasso,enet,ridge,trueqtl}] [--eqtl-h2 EQTL_H2]
              [--var-explained VAR_EXPLAINED] [-o OUTPUT] [--seed SEED] [--sim SIM]
              [--locus LOCUS]
              prefix
>>>>>>> dev

    Simulate TWAS using real genotype data

    positional arguments:
<<<<<<< HEAD
      prefix                Prefix to PLINK-formatted data
=======
      prefix                Prefix to PLINK-formatted data for GWAS LD information
>>>>>>> dev

    optional arguments:
      -h, --help            show this help message and exit
      --eqtl-prefix EQTL_PREFIX
<<<<<<< HEAD
                            Optional prefix to PLINK-formatted data for eQTL LD information. Otherwise use GWAS LD. (default: None)
      --test-prefix TEST_PREFIX
                            Optional prefix to PLINK-formatted data for LD information in TWAS test statistic. Otherwise use GWAS LD.
                            (default: None)
      --fast-gwas-sim FAST_GWAS_SIM
                            If 'True' (case-sensitive) then simulate GWAS summary data directly from LD (default: None)
      --ngwas NGWAS         Sample size for GWAS panel (default: 100000)
      --nqtl NQTL           Sample size for eQTL panel (default: 500)
      --model {10pct,1pct,1snp}
                            SNP model for generating gene expression. 10pct = 10% of SNPs, 1pct = 1% of SNPs, 1snp = 1 SNP (default: 10pct)
      --ld-ridge LD_RIDGE   Offset to add to LD Diagonal (default: 0.1)
      --linear-model {lasso,enet,ridge}
                            Linear model to predict gene expression from genotype. (default: lasso)
      --eqtl-h2 EQTL_H2     The narrow-sense heritability of gene expression (default: 0.1)
      --var-explained VAR_EXPLAINED
                            Variance explained in complex trait by gene expression (default: 0.01)
      -o OUTPUT, --output OUTPUT
                            Output prefix (default: None)
      --seed SEED           Seed for random number generation (default: None)
      --sim SIM             Simulation index for post-hoc analysis (default: None)
      --locus LOCUS         locus index for post-hoc analysis (default: None)

=======
                          Optional prefix to PLINK-formatted data for eQTL LD information.
                          Otherwise use GWAS LD. (default: None)
      --test-prefix TEST_PREFIX
                          Optional prefix to PLINK-formatted data for LD information in TWAS test
                          statistic. Otherwise use GWAS LD. (default: None)
      --fast-gwas-sim       If set then simulate GWAS summary data directly from LD (default: False)
      --ngwas NGWAS         Sample size for GWAS panel (default: 100000)
      --nqtl NQTL           Sample size for eQTL panel (default: 500)
      --ncausal NCAUSAL     Number of causal SNPs for gene expression/trait. Can represent explicit
                          number (e.g., 1, 10), a percentage using the 'pct' modifier (e.g.,
                          '1pct', '10pct'), or an average under a truncated Poisson model (e.g.,
                          '1avg', '10avg'). (default: 1)
      --ld-ridge LD_RIDGE   Offset to add to LD Diagonal (default: 0.1)
      --linear-model {lasso,enet,ridge,trueqtl}
                          Linear model to predict gene expression from genotype. (default: lasso)
      --eqtl-h2 EQTL_H2     The narrow-sense heritability of gene expression (default: 0.1)
      --var-explained VAR_EXPLAINED
                          Variance explained in complex trait by gene expression (default: 0.01)
      -o OUTPUT, --output OUTPUT
                          Output prefix (default: None)
      --seed SEED           Seed for random number generation (default: None)
      --sim SIM             Simulation index for post-hoc analysis (default: None)
      --locus LOCUS         locus index for post-hoc analysis (default: None)
>>>>>>> dev

The output will be a two tab-delimited reports.

The first `OUTPUT.summary.tsv` is a high-level summary that contains two columns:

| stat             | values |
| ------           | ------ |
| sim              | Simulation index |
| id               | Locus index |
| gwas.sim         | GWAS mode |
| real.time        | Real time spent on the current simulation |
| cpu.time         | CPU time spent on the current simulation |
| linear_model     | Linear model used in the current simulation |
| h2ge             | Variance explained in trait by GE |
| snp_model        | SNP model used in the current simulation |
| nsnps            | Number of SNPs |
| ngwas            | GWAS sample size |
| gwas.sim         | GWAS mode |
| nqtl             | eQTL sample size  |
| h2g              | Narrow-sense heritability of GE |
| h2g.hat          | Predicted narrow-sense heritability of GE |
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
<<<<<<< HEAD
| gwas.sim            | GWAS mode |
| ld.score            | ld score (ie. sum_i r_ij^2, where r_ij is LD between snps i, j) |
| gwas.beta           | beta coefficient in GWAS |
| gwas.se             | standard error in GWAS |
| gwas.true           | true causal effect for complex trait |
=======
| ld.score            | ld score (ie. sum_i r_ij^2, where r_ij is LD between snps i, j) |
| ld.score.causal     | ld score for causal variants |
| gwas.sim            | GWAS mode |
| gwas.true           | true causal effect for complex trait |
| gwas.beta           | beta coefficient in GWAS |
| gwas.se             | standard error in GWAS |
| eqtl.true           | true causal effect for expression |
>>>>>>> dev
| eqtl.beta           | beta coefficient in eQTL |
| eqtl.se             | standard error in eQTL |
| eqtl.model          | linear model to predict gene expression from genotype |
| eqtl.model.beta     | coefficient estimated in selected linear model |
<<<<<<< HEAD
| eqtl.true           | true causal effect for expression |
=======
>>>>>>> dev
