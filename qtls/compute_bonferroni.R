args = (commandArgs(trailingOnly=TRUE))

inp = args[1]
snp = args[2]

# Load in the eQTL results
eqtl = read.table(paste(inp,".gexp.qassoc",sep=''),as.is=T,head=T)

# Standard bonferroni correction
p.bonf = min(eqtl$P*nrow(eqtl))

# Load in the predixcan results
twas = read.table(paste(inp,".predixcan.qassoc",sep=''),as.is=T,head=T)
min_twas = which.min(twas$P)

# Weighted bonferroni correction
p.twas = min( eqtl$P[-min_twas]*(nrow(eqtl)-1)*2 , eqtl$P[min_twas]*2 )

# Load in the enformer results
enf = read.table(paste(inp,".enformer.qassoc",sep=''),as.is=T,head=T)
min_enf = which.min(enf$P)

# Weighted bonferroni correction
p.enformer = min( eqtl$P[-min_enf]*(nrow(eqtl)-1)*2 , eqtl$P[min_enf]*2 )

# Load in the promoter SNP and use it for weighted Bonferroni correction
min_prom = which(eqtl$SNP == snp)
if ( length(min_prom) != 0 ) {
p.promoter = min( eqtl$P[-min_prom]*(nrow(eqtl)-1)*2 , eqtl$P[min_prom]*2 )
} else {
p.promoter = NA
}

# Print the results
cat( inp , p.bonf , p.promoter , p.enformer , p.twas , '\n' , sep='\t' )
