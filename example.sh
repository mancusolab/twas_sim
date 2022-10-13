#!/bin/bash

if [ $1 ]; then
    IDX=$1
else
    IDX=1
fi

hapmap=HAPMAP_SNPS/
loci=ind_loci.bed
genes=glist-hg19.nodupe.autosome

# PARAMETERS
# !!! change to point to results/output directory !!!
odir=/scratch1/xwang505/TWAS/

# !!! change to point to plink installation !!!
plink=/project/nmancuso_8/xwang505/tools/plink2

# !!! change to point to reference panel dir !!!
okg=/project/nmancuso_8/data/LDSC/1000G_EUR_Phase3_plink/

# !!! change to values for your study !!!
N=100000 # N GWAS
NGE=500 # N EQTL
MODEL=1 # eQTL model; see sim.py for details
H2G=0.1 # eQTL h2g
H2GE=0.001 # variance explained in complex trait; 0 (null) to 0.01 (huge effect) are reasonable values
LINEAR_MODEL=trueqtl

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
    --keep EUR.samples \
    --extract $hapmap/hm.$numchr.snp \
    --force-intersect

    # make sure we have at least MIN_SNPS in our data
    if [[ `wc -l $odir/twas_sim_loci${IDX}.bim | awk '{print $1}'` -lt $MIN_SNPS ]];
    then
      rm $odir/twas_sim_loci${IDX}.{bim,bed,fam}
    fi
done

# geno was pulled, run simulation
echo "running fast simulation"
python sim.py \
    ${odir}/twas_sim_loci${IDX} \
    --ngwas $N \
    --nqtl $NGE \
    --ncausal $MODEL \
    --eqtl-h2 $H2G \
    --fast-gwas-sim \
    --var-explained $H2GE \
    --linear-model $LINEAR_MODEL \
    --output ${odir}/twas_sim_loci${IDX}

# remove temporary genotype data
rm ${odir}/twas_sim_loci${IDX}.{bed,bim,fam}
