# load in the true expression
true_gexp = (read.csv("data/geuvadis_peer_normalized.csv",as.is=T,head=T,row.names=1))

# --- other predictors
method_files = c("data/promoter_snp.csv",
				 "data/enformer_sample_preds.csv",
                 "data/xpresso_sample_preds.csv",
                 "data/predixcan_197KB_sample_preds.csv",
                 "data/expecto_sample_preds.csv",
                 "data/basenji2_sample_preds.csv")
method_names = c("Promoter","Enformer","Xpresso","Predixcan","Expecto","Basenji2")
# ---


# --- harmonize the data
cat("Harmonizing data\n")

predictions = list()
keep = rep(TRUE,ncol(true_gexp))
for ( i in 1:length(method_files) ) {
  predictions[[i]] = read.csv(method_files[i],as.is=T,head=T,row.names=1)
  # reorder the rows
  predictions[[i]] = predictions[[i]][ match(rownames(true_gexp),rownames(predictions[[i]])) , ]
  keep = keep & !is.na(match(colnames(true_gexp),colnames(predictions[[i]])))
  cat( method_names[i] , sum(keep) , '\n' )
}

# reorder so all the columns match
true_gexp = true_gexp[,keep]
for ( i in 1:length(method_files) ) {
  m = match(colnames(true_gexp),colnames(predictions[[i]]))
  predictions[[i]] = predictions[[i]][ , m ]
}

# --- finally compute the prediction rsq
cat("Computing gene Rsqs\n")

rsq_methods = matrix(NA,nrow=ncol(true_gexp),ncol=length(method_names))
colnames(rsq_methods) = method_names
rsq_joint = matrix(NA,nrow=ncol(true_gexp),ncol=1)
colnames(rsq_joint) = "Predixcan+\nEnformer"
pval_joint = rsq_joint

for ( i in 1:ncol(true_gexp) ) {
  for ( m in 1:length(method_names)) {
    reg = summary(lm( as.numeric(true_gexp[,i]) ~ as.numeric(predictions[[m]][,i]) ))
    rsq_methods[i,m] = reg$adj.r.sq    
  }
  
  # joint prediction using predixcan + enformer
  reg = summary(lm( as.numeric(true_gexp[,i]) ~ as.numeric(predictions[[ which(method_names == "Enformer") ]][,i]) + as.numeric(predictions[[ which(method_names == "Predixcan") ]][,i]) ))
  rsq_joint[i] = reg$adj.r.sq
  pval_joint[i] = reg$coef[2,4]

  if ( i %% 100 == 0 ) {
    cat (i ,'\n' )
  }
}

# Add the joint result
rsq_methods = cbind(rsq_methods,rsq_joint)

# --- Plot
cat("Plotting figures\n")

# Compute the means
mean_rsq = apply(rsq_methods,2,mean)
se_rsq = apply(rsq_methods,2,sd)/sqrt(nrow(rsq_methods))
# Reorder by magnitude
ord = order(mean_rsq)
mean_rsq = mean_rsq[ord]
se_rsq = se_rsq[ord]

pdf("rsq_barplot.pdf")
xval = barplot(mean_rsq,border=NA,col="royalblue",las=2,ylim=c(0,0.1),ylab="Mean Adjusted R2")
arrows(xval,mean_rsq-1.96*se_rsq,xval,mean_rsq+1.96*se_rsq,len=0)
abline(h=mean_rsq[ "Promoter" ],lty=2)
legend("topleft",legend=c("Promoter SNP"),bty="n",lty=2)
dev.off()

pdf("rsq_scatterplot.pdf")
plot(rsq_methods[,"Predixcan"],rsq_methods[,"Predixcan+\nEnformer"],cex=0.5,pch=19,col="#00000010",las=1,xlab="Adjusted R2 Predixcan",ylab="Adjusted R2 Predixcan + Enformer")
top = pval_joint < 0.05/length(pval_joint)
points( rsq_methods[top,"Predixcan"],rsq_methods[top,"Predixcan+\nEnformer"],col="royalblue",cex=0.5)
dev.off()
