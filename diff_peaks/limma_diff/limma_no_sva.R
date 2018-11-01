rm(list=ls())
library(ggplot2)
library(sva)
library(data.table)
library(limma)

data=data.frame(read.table('hydrogel.soft.stiff.saha.txt',header=TRUE,sep='\t'))
rownames(data)=paste(data$Chrom,data$Start,data$End,sep="_")
data$Chrom=NULL
data$Start=NULL 
data$End=NULL
data=voom(counts=data,normalize.method = "quantile")

batches=data.frame(read.table("batches_hydrogel.soft.stiff.saha.txt",header=TRUE,sep='\t'))
rownames(batches)=batches$Replicate
batches$Stiffness=factor(batches$Stiffness)
batches$Saha=factor(batches$Saha)
batches$Sample=factor(batches$Sample)

filtered=data$E
mod2=model.matrix(~0+Sample,data=batches)
fit <- lmFit(filtered,mod2)
cont.matrix=makeContrasts(soft_vs_stiff="Sample100Pa - Sample2000Pa",
                          saha_soft_vs_stiff="Samplesaha_100Pa - Samplesaha_2000Pa",
                          soft_saha_vs_nosaha="Samplesaha_100Pa - Sample100Pa",
                          stiff_saha_vs_nosaha="Samplesaha_2000Pa - Sample2000Pa",
                          sahasoft_vs_stiff="Samplesaha_100Pa - Sample2000Pa",
                          sahastiff_vs_soft="Samplesaha_2000Pa - Sample100Pa",
                          levels=mod2)
fit2=contrasts.fit(fit,cont.matrix)
e=eBayes(fit2)
comparisons=colnames(cont.matrix)
for(i in seq(1,length(comparisons)))
{
  tab<-topTable(e, number=nrow(e),coef=i,p.value = 0.05)
  if(nrow(tab)>0){
  names(tab)[1]=comparisons[i]
  write.table(tab,file=paste("diff_",comparisons[i],".tsv",sep=""),quote=FALSE,sep='\t',row.names = TRUE,col.names = TRUE)
  png(paste("volcano_peaks",comparisons[i],".png",sep=""))
  volcanoplot(e,coef=i,highlight =0,names=rownames(tab),main=comparisons[i])
  dev.off() 
  }
  else{
    write.table(tab,file=paste("diff_",comparisons[i],".tsv",sep=""),quote=FALSE,sep='\t',row.names = FALSE,col.names = FALSE)
  }
}
