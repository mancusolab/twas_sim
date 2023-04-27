# `twas_sim`

A python software leveraging real genotype data to simulate a complex trait as a function of latent expression, fit eQTL weights in independent data, and perform GWAS/TWAS on the complex trait.

`twas_sim` is described in:

> [twas_sim, a Python-based tool for simulation and power analysis of transcriptome-wide association analysis](https://doi.org/10.1093/bioinformatics/btad288)). Xinran Wang, Zeyun Lu, Arjun Bhattacharya, Bogdan Pasaniuc, Nicholas Mancuso, twas_sim, a Python-based tool for simulation and power analysis of transcriptome-wide association analysis, ***Bioinformatics***, 2023;

-------

[Installation](#Installation) | [Overview](#Overview) | [Usage](#Usage) | [Example](#Example) | [Notes](#Notes) | [Output](#Output) | [Support](#Support) | [Other Software](#Other-Software)


## Installation

To download `twas_sim`, first type the commands

    git clone https://github.com/mancusolab/twas_sim.git
    cd twas_sim

then,

    conda env create --file environment.yml
    conda activate twas_sim

## Overview

The script [example.sh](https://github.com/mancusolab/twas_sim/blob/master/example.sh) generates a single TWAS statistic using the simulator [sim.py](https://github.com/mancusolab/twas_sim/blob/master/sim.py). Please make sure to update the corresponding paths in `example.sh` first. The scripts rely on PLINK-formatted genotype data. We recommend to use data from [1000G](https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz) as reference genotypes. Make sure to enter the command when you are done with the simulator:

```conda deactivate```

#### Key features:

* **Dataset:** twas_sim first samples a genomic region uniformly at random. Then, it subsets reference genotype data to the genomic region from the previous step while removing genetic variants that are non-biallelic, MAF < 1%, HWE < 1e-5, and variant missingness > 10%. In addition, it restricts the genotype to HapMap3 variants.

* **LD reference panels:** twas_sim supports the option to use different LD reference panels across GWAS and eQTL simulations in addition to TWAS testing.

* **GWAS:** standard GWAS simulates GWAS summary statistics using individual-level genotype and phenotype data. Fast GWAS simulates GWAS summary statistics directly using the multivariate normal distribution parameterized by LD.

* **Linear model:** `twas_sim` supports predicting gene expressions using Elastic Net, LASSO, GBLUP, and true eQTL effect sizes. The dynamic import feature enables twas_sim to easily include external prediction tools.

* **TWAS:** twas_sim computes TWAS test statistics using LD, GWAS Z-scores, and estimated eQTL effect sizes.

* **Horizontal pleiotropy:** twas_sim accounts for the situation when nearby tagging genes are also tested in TWAS in addition to the causal TWAS model.

## Usage

[sim.py](https://github.com/mancusolab/twas_sim/blob/master/sim.py) is the actual simulator. Its usage is below:

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

### [example.sh](https://github.com/mancusolab/twas_sim/blob/master/example.sh)

This example script generates a single TWAS statistic using the simulator [sim.py](https://github.com/mancusolab/twas_sim/blob/master/sim.py). The simulator currently supports fitting LASSO, Elastic Net, and GBLUP prediction models to predict gene expression into GWAS. It is easily extendable with dynamic import function to include additional linear models.

<details>
<summary>show details</summary>

* First, we define GWAS sample size, eQTL sample size, eQTL model, eQTL h2g, variance explained in complex trait, and linear model:
  ```
  N=100000 # N GWAS
  NGE=500 # N EQTL
  MODEL=1 # eQTL model; see sim.py for details
  H2G=0.1 # eQTL h2g
  H2GE=0.001 # variance explained in complex trait; 0 (null) to 0.01 (huge effect) are reasonable values
  LINEAR_MODEL=enet
  ```

* Then, we call optional arguments to generate TWAS test statistics.
  * In this example, we use the first reference panel to compute GWAS LD information and the second reference panel to compute eQTL and TWAS LD information.
  ```
  python sim.py \
      $odir/twas_sim_sample1_loci${IDX} \
      --eqtl-prefix $odir/twas_sim_sample2_loci${IDX} \
      --test-prefix $odir/twas_sim_sample2_loci${IDX} \
      --ngwas $N \
      --nqtl $NGE \
      --ncausal $MODEL \
      --eqtl-h2 $H2G \
      --fast-gwas-sim \
      --IDX ${IDX}\
      --h2ge $H2GE \
      --linear-model $LINEAR_MODEL \
      --seed ${IDX} \
      --output $odir/twas_sim_loci${IDX}
  ```
</details>

### [example.external.sh](https://github.com/mancusolab/twas_sim/blob/master/example.external.sh)

This script works as a showcase of the dynamic import function mentioned above. It generates a single TWAS statistic with external python module [external_py.py](https://github.com/mancusolab/twas_sim/blob/master/external_py.py) or external R module [external_r.py](https://github.com/mancusolab/twas_sim/blob/master/external_r.py) and [external.R](https://github.com/mancusolab/twas_sim/blob/master/external.R).

<details>
<summary>show details</summary>

* First, we define GWAS sample size, eQTL sample size, eQTL model, eQTL h2g, variance explained in complex trait, and linear model:
  ```
  N=100000 # N GWAS
  NGE=500 # N EQTL
  MODEL=1 # eQTL model; see sim.py for details
  H2G=0.1 # eQTL h2g
  H2GE=0.001 # variance explained in complex trait; 0 (null) to 0.01 (huge effect) are reasonable values
  LINEAR_MODEL=external
  ```
* Then, we call optional arguments to generate TWAS test statistics.
  * In this example, we use the first reference panel to compute GWAS LD information and the second reference panel to compute eQTL and TWAS LD information.
  * twas_sim supports dynamically loading custom code (e.g., Python, R, Julia). Here, we use external R module to fit effect sizes (note: [external_r.py](https://github.com/mancusolab/twas_sim/blob/master/external_r.py) calls [external.R](https://github.com/mancusolab/twas_sim/blob/master/external.R) to call susieR on the simulated data).

  ```
  python sim.py \
      $odir/twas_sim_sample1_loci${IDX} \
      --eqtl-prefix $odir/twas_sim_sample2_loci${IDX} \
      --test-prefix $odir/twas_sim_sample2_loci${IDX} \
      --ngwas $N \
      --nqtl $NGE \
      --ncausal $MODEL \
      --eqtl-h2 $H2G \
      --fast-gwas-sim \
      --IDX ${IDX}\
      --h2ge $H2GE \
      --linear-model $LINEAR_MODEL \
      --external-module external_r \
      --seed ${IDX} \
      --output $odir/twas_sim_loci${IDX}
  ```
</details>

### [example.slurm.sh](https://github.com/mancusolab/twas_sim/blob/master/example.slurm.sh)

This is a batch script of the simulator. It generates TWAS statistic for a list of parameters specified in [slurm.params](https://github.com/mancusolab/twas_sim/blob/master/slurm.params).

<details>
<summary>show details</summary>

* First, we define a list of GWAS sample size, eQTL sample size, eQTL model, eQTL h2g, variance explained in complex trait, and linear model. The example below shows the first 4 lines of [slurm.params](https://github.com/mancusolab/twas_sim/blob/master/slurm.params):

  | ID    | N	        | NGE	     | MODEL  	| H2G	   | H2GE	    | LINEAR_MODEL |
  | ------| ------    | ------   | ------   | ------ | ------   | ------       |
  | 1	    | 50000	    | 500	     | 1    	  | 0.1	   | 0	      | enet         |
  | 2	    | 100000	  | 500	     | 1	      | 0.1	   | 0	      | enet         |
  | 3	    | 200000	  | 500	     | 1	      | 0.1	   | 0	      | enet         |
  | 4	    | 500000	  | 500	     | 1	      | 0.1	   | 0	      | enet         |
  | ...   | ...       | ...      | ...      | ...    | ...      | ...          |


* Second, we link [slurm.params](https://github.com/mancusolab/twas_sim/blob/master/slurm.params) to the shell script:
  ```
  # ID	N	NGE	MODEL	H2G	H2GE	LINEAR_MODEL
  IDX=$1
  N=$2 # N GWAS
  NGE=$3 # N EQTL
  MODEL=$4 # eQTL model; see sim.py for details
  H2G=$5 # eQTL h2g
  H2GE=$6 # h2ge in complex trait; 0 (null) to 0.01 (huge effect) are reasonable values
  LINEAR_MODEL=$7
  ```

* Then, we call optional arguments to generate TWAS test statistics for each user-defined parameter sets.
  * In this example, we use the first reference panel to compute GWAS LD information and the second reference panel to compute eQTL and TWAS LD information.
  * The first 4 lines of the [slurm.params](https://github.com/mancusolab/twas_sim/blob/master/slurm.params) generate 4 TWAS test statistics using GWAS sample size of 50K, 100K, 200K, and 500K, with all other parameters fixed.

  ```
  python sim.py \
        $odir/twas_sim_sample1_loci${IDX} \
        --eqtl-prefix $odir/twas_sim_sample2_loci${IDX} \
        --test-prefix $odir/twas_sim_sample2_loci${IDX} \
        --ngwas $N \
        --nqtl $NGE \
        --ncausal $MODEL \
        --eqtl-h2 $H2G \
        --fast-gwas-sim \
        --IDX ${IDX}\
        --h2ge $H2GE \
        --linear-model $LINEAR_MODEL \
        --seed ${IDX} \
        --output $odir/twas_sim_loci${IDX}.fast
  ```

* Here, we run ten jobs at a time and 40 jobs in total:
    
  ```
  #SBATCH --array=1-4
  ```
  ```
  start=`python -c "print( 1 + 10 * int(int($NR-1)))"`
  stop=$((start + 9))
  ```
</details>

## Notes
* **LD reference panels:** use optional prefix ```--eqtl-prefix $path-to-eQTL-LD-information``` and ```--test-prefix $path-to-TWAS-LD-information``` to point to PLINK-formatted eQTL and TWAS LD. Otherwise, twas_sim will use GWAS LD for all simulations.

* **GWAS:** twas_sim simulates standard GWAS by default ```(no optional prefix needed)```. Use optional prefix ```--fast-gwas-sim``` to simulated GWAS in the fast mode.

* **Linear Model:** use ```--linear-model enet``` to use Elastic Net model, or use ```--linear-model external``` to indicate an external model should be loaded (please see External Module below).

* **External Module:** use ```--linear-model external``` to load external predictive model and ```--external-module path-to-external-file``` to specify path to external Python file. e.g., if `my_module.py` contains `fit` function then pass in `my_module`. Please refer to script [external_py.py](https://github.com/mancusolab/twas_sim/blob/master/external_py.py) to fit external python model (OLS) and [external_r.py](https://github.com/mancusolab/twas_sim/blob/master/external_r.py) and [external.R](https://github.com/mancusolab/twas_sim/blob/master/external.R) to fit external R model ([SuSiE](https://github.com/stephenslab/susieR)).

* **Horizontal pleiotropy:** twas_sim generates GWAS effect-size using causal TWAS model by default ```(no optional prefix needed)```. Use ```--indep-gwas``` to generate GWAS effect-sizes under horizontal pleiotropy through linkage model.

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
| mean.gwas.chi2   | Mean GWAS SNP &chi;<sup>2</sup> |
| median.gwas.chi2 | Median GWAS SNP &chi;<sup>2</sup> |
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
| ld.score            | ld score (ie. sum<sub>i</sub> r<sub>ij</sub><sup>2</sup>, where r<sub>ij</sub> is LD between snps i, j) |
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

## Support

Please report any bugs or feature requests in the Issue Tracker. If users have any questions or comments, please contact Xinran Wang (xwang505@usc.edu), Zeyun Lu (zeyunlu@usc.edu), and Nicholas Mancuso (nmancuso@usc.edu).

## Other Software
Feel free to use other software developed by [Mancuso Lab](https://www.mancusolab.com/):

* [MA-FOCUS](https://github.com/mancusolab/ma-focus): a Bayesian fine-mapping framework using [TWAS](https://www.nature.com/articles/ng.3506) statistics across multiple ancestries to identify the causal genes for complex traits.

* [SuSiE-PCA](https://github.com/mancusolab/susiepca): a scalable Bayesian variable selection technique for sparse principal component analysis.

* [SuShiE](https://github.com/mancusolab/sushie): a multi-ancestry variational fine-mapping method on molecular data.
