#!/bin/bash

if [ $1 ]; then
    locus=$1
else
    locus=1
fi

# gene complexity of sampled region
lower=5
upper=20

hapmap=HAPMAP_SNPS/
loci=ind_loci.bed
genes=glist-hg19.nodupe.autosome

# point to plink installation
plink=~/bin/plink/plink

# point to reference panel dir
okg=~/pasaniucdata/DATA/LDSC/1000G_EUR/1000G_EUR_Phase3_plink/

# PARAMETERS
# change to point to results/output directory
odir=./
N=100000 # N GWAS
NGE=500 # N EQTL
MODEL=1pct # eQTL model; see sim.py for details
H2G=0.1 # eQTL h2g
H2GE=0.001 # variance explained in complex trait; 0 (null) to 0.01 (huge effect) are reasonable values


while [ ! -e $odir/${locus}.bim ]
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
    $plink --bfile $okg/1000G.EUR.QC.${numchr} --chr $numchr --from-bp $locus_start --to-bp $locus_stop --make-bed \
        --out $odir/${locus} --snps-only --hwe midp 1e-5 --geno 0.01 --maf 0.01 --allow-no-sex \
        --memory 2048 --keep EUR.samples --extract $hapmap/hm.$numchr.snp #--silent
done

echo "running simulation"
python sim.py \
    $odir/${locus} \
    --ngwas $N \
    --nqtl $NGE \
    --model $MODEL \
    --eqtl-h2 $H2G \
    --var-explained $H2GE \
    --output ${locus}
