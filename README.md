# twas_sim
Using real genotype data, simulate a complex trait as a function of latent expression, fit eQTL weights in independent data, and perform GWAS/TWAS on complex trait.

[Installation](#Installation) | [Overview](#Overview) | [Usage](#Usage) | [Example](#Example) | [Output](#Output) | [Notes](#Notes) | [Support](#Support) | [Other Software](#Other-Software)


## Installation
To download the TWAS simulator first type the commands

    git clone https://github.com/mancusolab/twas_sim.git
    cd twas_sim

then,

    conda env create --file environment.yml
    conda activate twas_sim

## Overview
The script `example.sh` will generate a single TWAS statistic using the simulator `sim.py`. Please be sure to update the paths in `example.sh` first. The script relies on PLINK-formatted genotype data. We recommend downloading [1000G](https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz) for use. When you are done with the simulator be sure to enter the command:

```conda deactivate```

## Usage
`sim.py` is the actual simulator. Its usage is below:

    usage: sim.py [-h] [--eqtl-prefix EQTL_PREFIX] [--test-prefix TEST_PREFIX]
              [--fast-gwas-sim] [--ngwas NGWAS] [--nqtl NQTL] [--IDX IDX]
              [--ncausal NCAUSAL] [--ld-ridge LD_RIDGE]
              [--linear-model {lasso,enet,ridge,trueqtl,external}]
              [--external-module EXTERNAL_MODULE] [--eqtl-h2 EQTL_H2]
              [--h2ge H2GE] [--indep-gwas] [--h2g-gwas H2G_GWAS] [-o OUTPUT]
              [-c] [--seed SEED]
              prefix

    Simulate TWAS using real genotype data

      optional arguments:
        -h, --help            show this help message and exit
        --eqtl-prefix EQTL_PREFIX
                              Optional prefix to PLINK-formatted data for eQTL LD
                              information. Otherwise use GWAS LD. (default: None)
        --test-prefix TEST_PREFIX
                              Optional prefix to PLINK-formatted data for LD
                              information in TWAS test statistic. Otherwise use GWAS
                              LD. (default: None)
        --fast-gwas-sim       If set then simulate GWAS summary data directly from
                              LD (default: False)
        --ngwas NGWAS         Sample size for GWAS panel (default: 100000)
        --nqtl NQTL           Sample size for eQTL panel (default: 500)
        --IDX IDX             Simulation index (default: None)
        --ncausal NCAUSAL     Number of causal SNPs for gene expression/trait. Can
                              represent explicit number (e.g., 1, 10), a percentage
                              using the 'pct' modifier (e.g., '1pct', '10pct'), or
                              an average under a truncated Poisson model (e.g.,
                              '1avg', '10avg'). (default: 1)
        --ld-ridge LD_RIDGE   Offset to add to LD Diagonal (default: 0.1)
        --linear-model {lasso,enet,ridge,trueqtl,external}
                              Linear model to predict gene expression from genotype.
                              Use external to indicate an external module should be
                              loaded. (default: lasso)
        --external-module EXTERNAL_MODULE
                              Path to external Python file with custom `fit`
                              function. Only used if `--linear-module=external`.
                              E.g., if `my_module.py` contains `fit function then
                              pass in `my_module`. (default: None)
        --eqtl-h2 EQTL_H2     The SNP heritability of gene expression. (default:0.1)
        --h2ge H2GE           Phenotypic variance explained by genetic component of
                              gene expression, (default: 0.01)
        --indep-gwas          Generate GWAS effect-sizes independently from eQTLs.
                              (default: False)
        --h2g-gwas H2G_GWAS   The SNP heritability of downstream phenotype. Only
                              used when `--indep-gwas` is set. (default: 0.01)
        -o OUTPUT, --output OUTPUT
                              Output prefix (default: None)
        -c, --compress        Compress output (gzip) (default: False)
        --seed SEED           Seed for random number generation (default: None)

## Example
Key features:
* [Allele](#Allele)
* [LD](#LD)
* [Dataset](#Dataset)
* [GWAS](#GWAS)
* [TWAS](#TWAS)
* [Linear Model](#Linear-Model)
  * [Elastic Net, LASSO, and GBLUP](*Elastic Net,-LASSO,-and-GBLUP)
  * [True eQTL](*True-eQTLl)
  * [External model](*External-model)
* [Horizontal Pleiotropy](#Horizontal-Pleiotropy)

Description:
### Allele
blah
### LD
**twas_sim** supports the option to use different LD reference panels across GWAS and eQTL simulations in addition to TWAS testing. Use optional prefix ```--eqtl-prefix $Path-to-eQTL-LD-information``` and ```--test-prefix $Path-to-TWAS-LD-information``` to PLINK-formatted  eQTL and TWAS LD. Otherwise, **twas_sim** will use GWAS LD for all simulations.
### GWAS
* **Standard GWAS**
Simulates GWAS summary statistics using individual-level genotype and phenotype data.
* **Fast GWAS**
Simulates GWAS summary statistics directly using the multivariate normal distribution parameterized by LD. Use optional prefix ```--fast-gwas-sim``` for fast GWAS simulation.
### TWAS
blah
### Linear Model
Blah
* ### Elastic Net, LASSO, and GBLUP
Blah
* ### True eQTL
Blah
* ### External model
Use ```--linear-model external``` to set external predictive model and ```--external-module path-to-external-file``` to specify path to external Python file. E.g., if `my_module.py` contains `fit function then pass in `my_module`.
### Horizontal Pleiotropy
```--indep-gwas```

## Output
The output will be a two tab-delimited reports.

The first `OUTPUT.summary.tsv` is a high-level summary that contains two columns:

| stat             | values |
| ------           | ------ |
| gwas.sim         | GWAS mode |
| real.time        | Real time spent on the current simulation |
| cpu.time         | CPU time spent on the current simulation |
| linear_model     | Linear model used in the current simulation |
| h2ge             | Variance explained in trait by GE |
| snp_model        | SNP model used in the current simulation |
| nsnps            | Number of SNPs |
| ngwas            | GWAS sample size |
| nqtl             | eQTL sample size  |
| h2g              | Narrow-sense heritability of GE |
| h2g.hat          | Predicted narrow-sense heritability of GE |
| avg.ldsc         | Average LD-score at the region |
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
| ld.score.causal     | ld score for causal variants |
| gwas.sim            | GWAS mode |
| gwas.true           | true causal effect for complex trait |
| gwas.beta           | beta coefficient in GWAS |
| gwas.se             | standard error in GWAS |
| eqtl.true           | true causal effect for expression |
| eqtl.beta           | beta coefficient in eQTL |
| eqtl.se             | standard error in eQTL |
| eqtl.model          | linear model to predict gene expression from genotype |
| eqtl.model.beta     | coefficient estimated in selected linear model |

## Notes
Note

## Support
Please report any bugs or feature requests in the Issue Tracker. If users have any questions or comments, please contact Xinran Wang (xwang505@usc.edu), Zeyun Lu (zeyunlu@usc.edu), and Nicholas Mancuso (nmancuso@usc.edu).

## Other Software
Feel free to use other software developed by `Mancuso Lab <https://www.mancusolab.com/>`_:

* `MA-FOCUS <https://github.com/mancusolab/ma-focus>`_: a Bayesian fine-mapping framework using `TWAS <https://www.nature.com/articles/ng.3506>`_ statistics across multiple ancestries to identify the causal genes for complex traits.

* `SuSiE-PCA <https://github.com/mancusolab/sushie>`_: a scalable Bayesian variable selection technique for sparse principal component analysis
