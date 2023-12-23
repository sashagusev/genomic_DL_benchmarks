# load in the true expression
true_gexp = (read.csv("geuvadis_peer_normalized.csv",as.is=T,head=T,row.names=1))

# --- generate a matrix using promoter SNPs 
nearest_snp = read.table("geuv_promoter_snp.raw",as.is=T,head=T)
snp_names = gsub("_[A-Z]+","",colnames(nearest_snp))
m = match( rownames(true_gexp) , nearest_snp$IID )
nearest_snp = nearest_snp[m,]
gene_map = read.table("geuv_promoter_snp_mapping.tab",as.is=T)

m = match(colnames(true_gexp),tolower(gene_map[,2]))
snp_map1 = gene_map[ m , 1 ]
m = match(colnames(true_gexp),tolower(gene_map[,3]))
snp_map2 = gene_map[ m , 1 ]
snp_map2[ is.na(snp_map2) ] = snp_map1[ is.na(snp_map2) ]

m = match(snp_map2,snp_names)
keep = !is.na(m)
snp_map2 = snp_map2[keep]
true_gexp = true_gexp[ , keep ]
m = match(snp_map2,snp_names)
snps = nearest_snp[ , m ]

rownames(snps) = rownames(true_gexp)
colnames(snps) = colnames(true_gexp)

write.csv(snps,file="promoter_snp.csv")
# ---
