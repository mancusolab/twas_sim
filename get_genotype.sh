#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem=4Gb
#SBATCH --array=1-200

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  IDX=$1
else
  IDX=$SLURM_ARRAY_TASK_ID
fi

source ~/init.sh

# gene complexity of sampled region
lower=5
upper=20

hapmap=/scratch1/xwang505/TWAS/twas_sim/HAPMAP_SNPS/
loci=/scratch1/xwang505/TWAS/twas_sim/ind_loci.bed
genes=/scratch1/xwang505/TWAS/twas_sim/glist-hg19.nodupe.autosome

# point to plink installation
#plink=~/bin/plink2/plink2
plink=/project/nmancuso_8/xwang505/tools/plink2

# point to reference panel dir
okg=/project/nmancuso_8/data/LDSC/1000G_EUR_Phase3_plink/

tdir=/scratch1/xwang505/tmp/

# change to point to results/output directory
odir=/project/nmancuso_8/xwang505/genotype

while [ ! -e $odir/*_loci${IDX}.bim ]
do
  echo "attempting ${IDX}"
  python /scratch1/xwang505/TWAS/twas_sim/sample_genes.py \
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
  --keep /scratch1/xwang505/TWAS/twas_sim/EUR.samples \
  --extract $hapmap/hm.$numchr.snp \
  --force-intersect
done
