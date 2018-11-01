rm(list=ls())
library(ggplot2)
library(sva)
library(data.table)
library(limma)

data=data.frame(read.table('2D.3D.soft.tcps.mg.txt',header=TRUE,sep='\t'))
rownames(data)=paste(data$Chrom,data$Start,data$End,sep="_")
data$Chrom=NULL
data$Start=NULL 
data$End=NULL
data=voom(counts=data,normalize.method = "quantile")

batches=data.frame(read.table("batches_mamgland.txt",header=TRUE,sep='\t'))
rownames(batches)=batches$Replicate
batches$Sample=factor(batches$Sample)

filtered=data$E
mod2=model.matrix(~0+Sample,data=batches)
fit <- lmFit(filtered,mod2)
cont.matrix=makeContrasts(mam_vs_3dsoft="SampleMamgland - Sample100Pa3D",
                          mam_vs_2dsoft="SampleMamgland - Sample2DSoft",
                          mam_vs_2dtcps="SampleMamgland - Sample2DTCPS",
                          soft3d_vs_soft2d="Sample2DSoft - Sample100Pa3D",
                          soft3d_vs_tcps2d="Sample100Pa3D - Sample2DTCPS",
                          soft2d_vs_tcps2d="Sample2DSoft - Sample2DTCPS",
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
