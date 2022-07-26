#!/bin/bash

if [ $1 ]; then
    locus=$1
else
    locus=1
fi

eval "$(conda shell.bash hook)"
conda activate twas_sim
source ~/init.sh

# gene complexity of sampled region
lower=5
upper=20

hapmap=/project/nmancuso_8/xwang505/twas_sim/HAPMAP_SNPS/
loci=/project/nmancuso_8/xwang505/twas_sim/ind_loci.bed
genes=/project/nmancuso_8/xwang505/twas_sim/glist-hg19.nodupe.autosome

# point to plink installation
plink=/project/nmancuso_8/xwang505/tools/plink2

# point to reference panel dir
okg=/project/nmancuso_8/data/LDSC/1000G_EUR_Phase3_plink/

# PARAMETERS
# change to point to results/output directory
odir=/scratch1/xwang505/TWAS/genotype_check

# minimum number of SNPs that need to exist for simulation
MIN_SNPS=450

N=100000 # N GWAS
NGE=500 # N EQTL
MODEL=1 # eQTL model; see sim.py for details
H2G=0.1 # eQTL h2g
H2GE=0.001 # variance explained in complex trait; 0 (null) to 0.01 (huge effect) are reasonable values
LINEAR_MODEL=trueqtl

# get genotype
while [ ! -e $odir/twas_sim${locus}.bim ]
do
    echo "attempting ${locus}"
    python sample_genes.py $loci $genes -l $lower -u $upper -o $odir/gene.${locus}.list --loc_output $odir/locus.${locus}.txt

    params=`sed "1q;d" $odir/locus.${locus}.txt`
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
    --out $odir/twas_sim${locus} \
    --snps-only \
    --hwe midp 1e-5 \
    --geno 0.01 \
    --maf 0.01 \
    --allow-no-sex \
    --memory 2048 \
    --keep /project/nmancuso_8/xwang505/twas_sim/EUR.samples \
    --extract $hapmap/hm.$numchr.snp \
    --force-intersect

done

OUT=/scratch1/xwang505/TWAS/res_check/twas_sim${locus}
rm -rf $OUT*
oloci=/scratch1/xwang505/TWAS/genotype_check/twas_sim${locus}

# run simulation
echo "running fast simulation"
python /project/nmancuso_8/xwang505/twas_sim/sim.py \
$oloci \
--ngwas $N \
--nqtl $NGE \
--ncausal $MODEL \
--eqtl-h2 $H2G \
--fast-gwas-sim \
--var-explained $H2GE \
--linear-model $LINEAR_MODEL \
--output $OUT.fast

echo "running std simulation"
python /project/nmancuso_8/xwang505/twas_sim/sim.py \
$oloci \
--ngwas $N \
--nqtl $NGE \
--ncausal $MODEL \
--eqtl-h2 $H2G \
--var-explained $H2GE \
--linear-model $LINEAR_MODEL \
--output $OUT.std
