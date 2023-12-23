# Data download and preparation

All relevant data has already been generated but the code is here for reproducibility. 

Download the large files from here: https://www.dropbox.com/scl/fo/y6m00ojoadwaaoo9kjeaa/h?rlkey=jv5zb7wq6ujwtsqzq7oldn67j&dl=0

*Download meta data from Huang et al repository*

```
wget https://raw.githubusercontent.com/ni-lab/personalized-expression-benchmark/main/predictions/geuvadis_peer_normalized.csv
wget https://raw.githubusercontent.com/ni-lab/personalized-expression-benchmark/main/predictions/predixcan_197KB_sample_preds.csv
wget https://raw.githubusercontent.com/ni-lab/personalized-expression-benchmark/main/data/gene_list.csv
```

*Download predictions from Huang et al.*

```
for t in basenji2 enformer expecto xpresso; do
wget https://raw.githubusercontent.com/ni-lab/personalized-expression-benchmark/main/predictions/${t}_sample_preds.csv
done
```

*Clean up and transpose the data to make it easier to load for eQTL analysis*

```
R --slave < convert.R
```

*Identify promoter variants and extract them*

Assumes `bedtools` and `plink` (1.90) are in the command line:

```
# Convert the reference gene list to a bed file
cat gene_list.csv  | tr ',' '\t' | awk '{ print "chr"$2,$3,$3+1,$1":"$4 }' | tr ' ' '\t' | sort -k1,1 -k2,2n > gene_list.hg19.bed

# Use bedtools to identify the nearest SNP to each gene
cat geuv*.bim  | awk '{ print "chr"$1,$4,$4+1,$2 }' | tr ' ' '\t' | sort -k1,1 -k2,2n | bedtools closest -b - -a gene_list.hg19.bed -d > gene_list.hg19.1kg_phase3_maf01.closest.bed

# Subset to the relevant columns
cat gene_list.hg19.1kg_phase3_maf01.closest.bed | awk '{ print $8,$4 }' | tr ':' '\t' | tr ' ' '\t' > gene_list.hg19.1kg_phase3_maf01.closest.snp

# Make an extract file for plink
cat gene_list.hg19.1kg_phase3_maf01.closest.bed | awk '{ print $8 }' | sort | uniq > gene_list.hg19.1kg_phase3_maf01.closest.bed.extract

# Use plink to extract these positions (assume we have one merged genotype file)
plink --bfile phase3_maf01 --keep geuvadis_peer_normalized.keep --extract gene_list.hg19.1kg_phase3_maf01.closest.bed.extract --recode A --out geuv_promoter_snp.raw --allow-extra-chr

# Make a snp - gene mapping
cat gene_list.hg19.1kg_phase3_maf01.closest.bed | awk '{ print $8,$4 }' | tr ':' '\t' | tr ' ' '\t' > geuv_promoter_snp_mapping.tab
```

*Generate a csv file with the promoter SNP for each gene*

This will use the above files to make a new `promoter_snp.csv` that contains the promoter SNP for each gene and is structured like the other prediction CSVs
 
```
R --slave < promoter_snp.csv
```
