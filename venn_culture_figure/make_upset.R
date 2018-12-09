library(UpSetR)
data=read.table("upsetr_inputs.txt",header=TRUE,sep='\t',row.names=1)
upset(data,order.by="freq")
