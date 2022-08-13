#!/bin/bash

if [ $1 ]; then
    IDX=$1
else
    IDX=1
fi

# gene complexity of sampled region
lower=5
upper=20

hapmap=/scratch1/xwang505/TWAS/twas_sim/HAPMAP_SNPS/
loci=/scratch1/xwang505/TWAS/twas_sim/ind_loci.bed
genes=/scratch1/xwang505/TWAS/twas_sim/glist-hg19.nodupe.autosome

# point to plink installation
plink=/project/nmancuso_8/xwang505/tools/plink2

# point to reference panel dir
okg=/project/nmancuso_8/data/LDSC/1000G_EUR_Phase3_plink/

# PARAMETERS
# change to point to results/output directory
odir=/scratch1/xwang505/TWAS/genotype_check
rdir=/scratch1/xwang505/TWAS/res_check

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
    python sample_genes.py $loci $genes -l $lower -u $upper -o $odir/gene.${IDX}.list --loc_output $odir/locus.${IDX}.txt

    params=`sed "1q;d" $odir/locus.${IDX}.txt`
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
    --out $odir/twas_sim${SIM}_loci${IDX} \
    --snps-only \
    --hwe midp 1e-5 \
    --geno 0.01 \
    --maf 0.01 \
    --allow-no-sex \
    --memory 2048 \
    --keep /scratch1/xwang505/TWAS/twas_sim/EUR.samples \
    --extract $hapmap/hm.$numchr.snp \
    --force-intersect

    if [[ `wc -l $odir/twas_sim_loci${IDX}.bim | awk '{print $1}'` -lt $MIN_SNPS ]];
    then
      rm $odir/twas_sim_loci${IDX}.{bim,bed,fam}
    fi
done

# run simulation
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
--output ${rdir}/twas_sim${SIM}_loci${IDX}.std

#clean up!
rm ${odir}/twas_sim_loci${IDX}*
