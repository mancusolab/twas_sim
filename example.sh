#!/bin/bash

if [ $1 ]; then
    locus=$1
else
    locus=1
fi

# gene complexity of sampled region
lower=5
upper=20

hapmap=~/src/twas_sim/HAPMAP_SNPS/
loci=~/src/twas_sim/ind_loci.bed
genes=~/src/twas_sim/glist-hg19.nodupe.autosome

# point to plink installation
plink=~/bin/plink2/plink2

# point to reference panel dir
okg=/project/nmancuso_8/data/LDSC/1000G_EUR_Phase3_plink/

# PARAMETERS
# change to point to results/output directory
odir=~/projects/scratch/twas_sim/

N=100000 # N GWAS
NGE=500 # N EQTL
MODEL=1 # eQTL model; see sim.py for details
H2G=0.1 # eQTL h2g
H2GE=0.001 # variance explained in complex trait; 0 (null) to 0.01 (huge effect) are reasonable values
LINEAR_MODEL=trueqtl

while [ ! -e $odir/twas_sim${SIM}_loci${locus}.bim ]
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
    --out $odir/twas_sim${SIM}_loci${locus} \
    --snps-only \
    --hwe midp 1e-5 \
    --geno 0.01 \
    --maf 0.01 \
    --allow-no-sex \
    --memory 2048 \
    --keep ~/src/twas_sim/EUR.samples \
    --extract $hapmap/hm.$numchr.snp \
    --force-intersect

done

echo "running fast simulation"
python sim.py \
$odir/twas_sim${SIM}_loci${locus} \
--ngwas $N \
--nqtl $NGE \
--ncausal $MODEL \
--eqtl-h2 $H2G \
--fast-gwas-sim \
--var-explained $H2GE \
--linear-model $LINEAR_MODEL \
--output $odir/fast

echo "running std simulation"
python sim.py \
$odir/twas_sim${SIM}_loci${locus} \
--ngwas $N \
--nqtl $NGE \
--ncausal $MODEL \
--eqtl-h2 $H2G \
--var-explained $H2GE \
--linear-model $LINEAR_MODEL \
--output $odir/std
