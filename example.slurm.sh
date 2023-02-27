#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --mem=8Gb
#SBATCH --array=1-1

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
else
  NR=$SLURM_ARRAY_TASK_ID
fi

source ~/init.sh
conda activate twas_sim

# !!! change this to use SGE or the number of ind tasks per scheduler !!!
# 40 jobs in total (using slurm.params file)
# ten jobs at a time
start=`python -c "print( 1 + 10 * int(int($NR-1)))"`
stop=$((start + 9))

hapmap=HAPMAP_SNPS/
loci=ind_loci.bed
genes=glist-hg19.nodupe.autosome

# PARAMETERS
# !!! change to point to results/output directory !!!
odir=/scratch1/xwang505/TWAS/output

# !!! change to point to plink installation !!!
plink=/project/nmancuso_8/xwang505/tools/plink2

# !!! change to point to reference panel dir !!!
okg=/project/nmancuso_8/data/LDSC/1000G_EUR_Phase3_plink/

# minimum number of SNPs that need to exist for simulation
MIN_SNPS=450

# run simulation
# PARAMETERS
for IDX in `seq $start $stop`
do
  params=`sed "$((IDX+1))q;d" param_twas_sim.tsv` #change back to slurm.params
  echo "$((IDX+1)) ${params}"
  set -- junk $params
  shift

  # ID	N	NGE	MODEL	H2G	H2GE	LINEAR_MODEL
  IDX=$1
  N=$2 # N GWAS
  NGE=$3 # N EQTL
  MODEL=$4 # eQTL model; see sim.py for details
  H2G=$5 # eQTL h2g
  H2GE=$6 # h2ge in complex trait; 0 (null) to 0.01 (huge effect) are reasonable values
  LINEAR_MODEL=$7

  # get genotype
  while [[ ! -e $odir/twas_sim_loci${IDX}.bim ]]
  do
    echo "attempting ${IDX}"
    params=`shuf -n 1 $loci`
    set -- junk $params
    shift

    #1 	 10583 	 1892607
    # get the region bounds +-500kb
    chr=$1
    numchr=`echo $chr | sed 's/chr//'`
    locus_start=`python -c "print(int(max($2 - 500e3, 1)))"`
    locus_stop=`python -c "print( int(int($3) + 500e3))"`

    # grab geno from reference panel
    # keep only common hapmap3 variants
    $plink \
        --bfile $okg/1000G.EUR.QC.${numchr} \
        --chr $numchr \
        --from-bp $locus_start \
        --to-bp $locus_stop \
        --make-bed \
        --out $odir/twas_sim_loci${IDX} \
        --snps-only \
        --hwe midp 1e-5 \
        --geno 0.01 \
        --maf 0.01 \
        --allow-no-sex \
        --memory 2048 \
        --keep EUR.samples \
        --extract $hapmap/hm.$numchr.snp \
        --force-intersect

    # make sure we have at least MIN_SNPS in our data
    if [[ `wc -l $odir/twas_sim_loci${IDX}.bim | awk '{print $1}'` -lt $MIN_SNPS ]];
    then
      rm $odir/twas_sim_loci${IDX}.{bim,bed,fam}
    fi

    $plink \
    --bfile $odir/twas_sim_loci${IDX} \
    --keep EUR1.samples \
    --make-bed \
    --out $odir/twas_sim_sample1_loci${IDX}

    $plink \
    --bfile $odir/twas_sim_loci${IDX} \
    --keep EUR2.samples \
    --make-bed \
    --out $odir/twas_sim_sample2_loci${IDX}
  done

  # run simulation
  echo "running fast simulation"
  python sim.py \
      $odir/twas_sim_sample1_loci${IDX} \
      --eqtl-prefix $odir/twas_sim_sample2_loci${IDX} \
      --test-prefix $odir/twas_sim_sample2_loci${IDX} \
      --ngwas $N \
      --nqtl $NGE \
      --ncausal $MODEL \
      --eqtl-h2 $H2G \
      --fast-gwas-sim \
      --output-gexpr \
      --ld-model \
      --IDX ${IDX}\
      --h2ge $H2GE \
      --linear-model $LINEAR_MODEL \
      --seed ${IDX} \
      --output $odir/twas_sim_loci${IDX}
done

# remove temporary genotype data
rm ${odir}/twas_sim_loci${IDX}.{bed,bim,fam}
