gexp = read.csv("geuvadis_peer_normalized.csv",as.is=T,head=T,row.names=1)
write.table( t(gexp) , quote=F , sep='\t' , row.names=T , col.names=T , file="geuvadis_peer_normalized.transpose.tab" )

gexp = read.csv("predixcan_197KB_sample_preds.csv",as.is=T,head=T,row.names=1)
write.table( t(gexp) , quote=F , sep='\t' , row.names=T , col.names=T , file="predixcan_197KB_sample_preds.transpose.tab" )

gexp = read.csv("enformer_sample_preds.csv",as.is=T,head=T,row.names=1)
write.table( t(gexp) , quote=F , sep='\t' , row.names=T , col.names=T , file="enformer_sample_preds.tab" )
