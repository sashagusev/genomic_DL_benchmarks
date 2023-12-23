#!/bin/bash

# Assume plink is in the command line
###module load plink/1.90
# Chromosome to run on
CHR=$1
DIR="../data"

# this is where results will go:
rm -f  chr${CHR}.bonf.tab

# iterate over each gene
cat $DIR/gene_list.hg19.bed | awk -v c=$CHR '$1 == "chr"c' | while read line; do

# identify the 200kb locus
p0=`echo $line | awk '{ print $2 - 100e3 }'`
p1=`echo $line | awk '{ print $2 + 100e3 }'`
g=`echo $line | tr ':' '\t' | awk '{ if($NF == "+" || $NF == "-" ) print $(NF-1); else print $NF }'`

echo $p0 $p1 $g

# Generate the gene expression phenotype

paste \
<(cat $DIR/geuvadis_peer_normalized.transpose.tab  | head -n1 | tr '\t' '\n' | awk '{ print 0,$1 }' | tr ' ' '\t') \
<(cat $DIR/geuvadis_peer_normalized.transpose.tab  | awk -vg=$g 'toupper($1) == g' | cut -f 2- | tr '\t' '\n') > $CHR.$g.gexp.pheno

# Generate the phenotypes for the predictors (from which we'll get the weights)

paste \
<(cat $DIR/predixcan_197KB_sample_preds.transpose.tab  | head -n1 | tr '\t' '\n' | awk '{ print 0,$1 }' | tr ' ' '\t') \
<(cat $DIR/predixcan_197KB_sample_preds.transpose.tab  | awk -vg=$g 'toupper($1) == g' | cut -f 2- | tr '\t' '\n') > $CHR.$g.predixcan.pheno

paste \
<(cat $DIR/enformer_sample_preds.transpose.tab  | head -n1 | tr '\t' '\n' | awk '{ print 0,$1 }' | tr ' ' '\t') \
<(cat $DIR/enformer_sample_preds.transpose.tab  | awk -vg=$g 'toupper($1) == g' | cut -f 2- | tr '\t' '\n') > $CHR.$g.enformer.pheno

# For each phenotype, run a QTL analysis on this locus
for p in gexp predixcan enformer; do
plink --bfile geuv_chr${CHR} --chr $CHR --from-bp $p0 --to-bp $p1 --pheno $CHR.$g.$p.pheno --assoc --out $CHR.$g.$p
rm -f $CHR.$g.$p.pheno 
done

# Now compute the weighted Bonferroni correction
# pull out the promoter SNP
snp=`grep -i -w $g $DIR/gene_list.hg19.1kg_phase3_maf01.closest.snp | awk '{ print $1 }'`
if [ -f $CHR.$g.gexp.qassoc ]; then
R --slave --args $CHR.$g $snp < compute_bonferroni.R >> chr${CHR}.bonf.tab
fi

done
