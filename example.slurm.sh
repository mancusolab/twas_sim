#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem=8Gb
#SBATCH --array=1-4000

if [ ! $SLURM_ARRAY_TASK_ID ]; then
  NR=$1
else
  NR=$SLURM_ARRAY_TASK_ID
fi

source ~/init.sh

start=`python -c "print( 1 + 10 *int(int($NR-1)))"`
stop=$((start + 9))

for IDX in `seq $start $stop`
do
  # PARAMETERS
  params=`sed "$((IDX+1))q;d" /scratch1/xwang505/TWAS/twas_sim/param_twas_sim.tsv`
  echo "$((IDX+1)) ${params}"
  set -- junk $params
  shift

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


  OUT=/scratch1/xwang505/TWAS/res_check/twas_sim${SIM}_loci${locus}
  rm -rf $OUT*
  oloci=/project/nmancuso_8/xwang505/genotype/twas_sim_loci${locus}

  echo "running simulation"
  python /scratch1/xwang505/TWAS/twas_sim/sim.py \
  $oloci \
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
  --output $OUT

done
