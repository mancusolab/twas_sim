#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --mem=64Gb
#SBATCH --array=1-1920

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
else
  NR=$SLURM_ARRAY_TASK_ID
fi

eval "$(conda shell.bash hook)"
conda activate twas_sim
source ~/init.sh

start=`python -c "print( 1 + 10 *int(int($NR-1)))"`
stop=$((start + 9))

# gene complexity of sampled region
lower=5
upper=20

hapmap=/project/nmancuso_8/xwang505/twas_sim/HAPMAP_SNPS/
loci=/project/nmancuso_8/xwang505/twas_sim/ind_loci.bed
genes=/project/nmancuso_8/xwang505/twas_sim/glist-hg19.nodupe.autosome

# point to plink installation
#plink=~/bin/plink2/plink2
plink=/project/nmancuso_8/xwang505/tools/plink2

# point to reference panel dir
okg=/project/nmancuso_8/data/LDSC/1000G_EUR_Phase3_plink/

tdir=/scratch1/xwang505/tmp/

# change to point to results/output directory
odir=/scratch1/xwang505/TWAS/genotype
rdir=/scratch1/xwang505/TWAS/res

# minimum number of SNPs that need to exist for simulation
MIN_SNPS=450

# run simulation
# PARAMETERS
for IDX in `seq $start $stop`
do
  params=`sed "$((IDX+1))q;d" /project/nmancuso_8/xwang505/twas_sim/param_twas_sim.tsv`
  echo "$((IDX+1)) ${params}"
  set -- junk $params
  shift

  # SIM	ID	N	NGE	MODEL	H2G	H2GE	LINEAR_MODEL
  SIM=$1
  IDX=$2
  N=$3 # N GWAS
  NGE=$4 # N EQTL
  MODEL=$5 # eQTL model; see sim.py for details
  H2G=$6 # eQTL h2g
  H2GE=$7 # variance explained in complex trait; 0 (null) to 0.01 (huge effect) are reasonable values
  LINEAR_MODEL=$8
done

# get genotype
for IDX in `seq $start $stop`
do
  while [[ ! -e $odir/twas_sim_loci${IDX}.bim ]]
  do
    echo "attempting ${IDX}"
    python /project/nmancuso_8/xwang505/twas_sim/sample_genes.py \
    $loci \
    $genes \
    -l $lower \
    -u $upper \
    -o $tdir/gene_twas_sim_loci${IDX}.list \
    --loc_output $tdir/locus_twas_loci${IDX}.txt \
    --seed $IDX

    params=`sed "1q;d" $tdir/locus_twas_loci${IDX}.txt`
    set -- junk $params
    shift

    #1 	 10583 	 1892607
    # get the region bounds
    chr=$1
    numchr=`echo $chr | sed 's/chr//'`
    locus_start=$2
    locus_stop=$3

    # replace with your path to PLINK reference data
    $plink --bfile $okg/1000G.EUR.QC.${numchr} \
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
    --keep /project/nmancuso_8/xwang505/twas_sim/EUR.samples \
    --extract $hapmap/hm.$numchr.snp \
    --force-intersect

    if [[ `wc -l $odir/twas_sim_loci${IDX}.bim | awk '{print $1}'` -lt $MIN_SNPS ]];
    then
      rm $odir/twas_sim_loci${IDX}.{bim,bed,fam}
    fi

  done
done

for IDX in `seq $start $stop`
do
  echo "running fast simulation"
  python /project/nmancuso_8/xwang505/twas_sim/sim.py \
  ${odir}/twas_sim_loci${IDX} \
  --ngwas $N \
  --nqtl $NGE \
  --ncausal $MODEL \
  --eqtl-h2 $H2G \
  --fast-gwas-sim \
  --var-explained $H2GE \
  --linear-model $LINEAR_MODEL \
  --seed ${IDX} \
  --locus ${IDX} \
  --sim ${SIM} \
  --output ${rdir}/twas_sim${SIM}_loci${IDX}.fast

  echo "running std simulation"
  python /project/nmancuso_8/xwang505/twas_sim/sim.py \
  ${odir}/twas_sim_loci${IDX} \
  --ngwas $N \
  --nqtl $NGE \
  --ncausal $MODEL \
  --eqtl-h2 $H2G \
  --var-explained $H2GE \
  --linear-model $LINEAR_MODEL \
  --seed ${IDX} \
  --locus ${IDX} \
  --sim ${SIM} \
  --output ${rdir}/twas_sim${SIM}_loci${IDX}.std
done

#clean up!
rm ${odir}/twas_sim_loci${IDX}*
