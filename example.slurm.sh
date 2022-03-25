#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem=4Gb
#SBATCH --array=1-4800

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  IDX=$1
else
  IDX=$SLURM_ARRAY_TASK_ID
fi

source ~/init.sh

# gene complexity of sampled region
lower=5
upper=20

hapmap=/scratch/xwang505/twas_sim/HAPMAP_SNPS/
loci=/scratch/xwang505/twas_sim/ind_loci.bed
genes=/scratch/xwang505/twas_sim/glist-hg19.nodupe.autosome

# point to plink installation
#plink=~/bin/plink2/plink2
plink=/project/nmancuso_8/xwang505/tools/plink2

# point to reference panel dir
okg=/project/nmancuso_8/data/LDSC/1000G_EUR_Phase3_plink/

# PARAMETERS
params=`sed "$((IDX+1))q;d" /scratch/xwang505/twas_sim/param_twas_sim.tsv`
echo "$((IDX+1)) ${params}"
set -- junk $params
shift

# change to point to results/output directory
odir=/scratch/xwang505/twas_sim_tmp/

# SIM	ID	N	NGE	MODEL	H2G	H2GE	fastGWAS	LINEAR_MODEL
SIM=$1
locus=$2
N=$3 # N GWAS
NGE=$4 # N EQTL
MODEL=$5 # eQTL model; see sim.py for details
H2G=$6 # eQTL h2g
H2GE=$7 # variance explained in complex trait; 0 (null) to 0.01 (huge effect) are reasonable values
fastGWAS=$8
LINEAR_MODEL=$9


while [ ! -e $odir/twas_sim${SIM}_loci${locus}.bim ]

do
  echo "attempting ${locus}"
  python /scratch/xwang505/twas_sim/sample_genes.py \
  $loci \
  $genes \
  -l $lower \
  -u $upper \
  -o $odir/gene.twas_sim${SIM}_loci${locus}.list \
  --loc_output $odir/locus.twas_sim${SIM}_loci${locus}.txt \
  --seed $locus

  params=`sed "1q;d" $odir/locus.twas_sim${SIM}_loci${locus}.txt`
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
  --out $odir/twas_sim${SIM}_loci${locus} \
  --snps-only \
  --hwe midp 1e-5 \
  --geno 0.01 \
  --maf 0.01 \
  --allow-no-sex \
  --memory 2048 \
  --keep /scratch/xwang505/twas_sim/EUR.samples \
  --extract $hapmap/hm.$numchr.snp \
  --force-intersect

done

echo "running simulation"
python /scratch/xwang505/twas_sim/sim.py \
$odir/twas_sim${SIM}_loci${locus} \
--ngwas $N \
--nqtl $NGE \
--model $MODEL \
--eqtl-h2 $H2G \
--var-explained $H2GE \
--fast-gwas-sim $fastGWAS \
--linear-model $LINEAR_MODEL \
--seed ${locus} \
--locus ${locus} \
--sim ${SIM} \
--output /scratch/xwang505/twas_sim_results/twas_sim${SIM}_loci${locus}

rm $odir/twas_sim${SIM}_loci${locus}.bim
rm $odir/twas_sim${SIM}_loci${locus}.bed
rm $odir/twas_sim${SIM}_loci${locus}.fam
