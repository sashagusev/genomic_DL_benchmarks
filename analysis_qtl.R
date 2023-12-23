library("qvalue")

# read in the gene-wise results
tab = read.table("qtls/all_genes_bonferroni.tab",as.is=T)
colnames(tab) = c("Gene","Flat","Promoter","Enformer","Predixcan")
methods = colnames(tab)[2:ncol(tab)]

result_pi1 = rep(NA,length(methods))
names(result_pi1) = methods
result_hits = result_pi1

# compute the pi1 and number of Bonferroni significant hits
for ( i in 1:length(methods) ) {
  p = tab[,methods[i]]
  p[p>1] = 1
  p = p[!is.na(p)]
  qv = qvalue(p)
  result_pi1[i] = 1 - qv$pi0
  result_hits[i] = sum(p<0.05/length(p))
}

# Visualize

pdf("plot_eqtl.pdf",height=4,width=8)
par(mfrow=c(1,2))
methods = c("Flat","Promoter","Enformer","Predixcan")
barplot(result_pi1,border=NA,col="royalblue",las=2,ylim=c(0,1),ylab="FDR pi1 (non-null genes)")
barplot(result_hits,names=methods,border=NA,col="royalblue",las=2,ylab="# Bonferroni significant genes")
dev.off()