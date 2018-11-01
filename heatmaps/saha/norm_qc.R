#following suggestions here: https://support.bioconductor.org/p/62954/
rm(list=ls())
library(DESeq2)
library(limma)
library(preprocessCore)
#load the counts 
data=read.table("hydrogel.soft.stiff.saha.txt",header=TRUE,sep='\t')
rownames(data)=paste(data$Chrom,data$Start,data$End,sep="_")
data$Chrom=NULL
data$Start=NULL
data$End=NULL

#1. quantile normalization
print(mean(data))
data_qnorm=normalize.quantiles(as.matrix(data))
print(mean(data))

#2. rlogTransform 
data_qnorm_rlog=rlog(round(data_qnorm))


#3. Remove batch effect 
batches=data.frame(read.table("batches_hydrogel.soft.stiff.saha.txt",header=TRUE,sep='\t'))
batches$Stiffness=factor(batches$Stiffness)
batches$Saha=factor(batches$Saha)
batches$Replicate=NULL
batches$Sample=NULL 
mod1=model.matrix(~Stiffness+Saha,data=batches)
#we just keep the terms for stiffness and saha 
#normalized=removeBatchEffect(data_qnorm_rlog,design=mod1)
normalized=removeBatchEffect()

#4. save table of normalized values
row.names(normalized)=row.names(data)
colnames(normalized)=colnames(data)
write.table(normalized,file="custom_normalized.txt",quote=FALSE,col.names=TRUE,row.names=TRUE,sep='\t')