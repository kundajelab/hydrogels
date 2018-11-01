rm(list=ls())
library(ggplot2)
library(sva)
library(DESeq2)
library("BiocParallel")
parallelFlag=TRUE
register(MulticoreParam(20))

data=data.frame(read.table('background_subtracted.txt',header=TRUE,sep='\t'))
rownames(data)=paste(data$Chrom,data$Start,data$End,sep="_")
data$Chrom=NULL
data$Start=NULL 
data$End=NULL
data[data<0]=0
data=data[rowSums(data)>0,]

#load the normalization factors
normfactors=data.frame(read.table("norm_factors.txt",header=TRUE,sep='\t'))
normfactors=normfactors[1:nrow(data),]
rownames(normfactors)=rownames(data)

batches=data.frame(read.table("batches_hydrogel.soft.v.stiff.txt",header=TRUE,sep='\t'))
rownames(batches)=batches$Replicate
batches$Stiffness=factor(batches$Stiffness)
batches$Saha=factor(batches$Saha)
batches$Sample=factor(batches$Sample)

#svaseq

mod0=model.matrix(~1,data=batches)
mod1=model.matrix(~as.factor(Sample),data=batches)
sva.obj=svaseq(as.matrix(data),mod1,mod0)
sur_var=data.frame(sva.obj$sv)
names(sur_var)=c("sv1","sv2")
batches=cbind(batches,sur_var)

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = batches,
                              design = ~Sample + sv2)
			      
#provide the custom normalization factors
sizeFactors(dds)=c(0.0783,0.0545,0.0826,0.049,0.0422,0.0466)			      
normalizationFactors(dds)=as.matrix(normfactors)

dds <- DESeq(dds,parallel = parallelFlag)

#get the results for various contrasts 
res=results(dds,contrast=c("Sample","100Pa","2000Pa"),parallel=TRUE)
res$logPadj=-10*log10(res$padj)
res=as.data.frame(res)
resfilt=na.omit(res)
numsig=sum(resfilt$padj < 0.05)
print(numsig) 
png("deseq_truerep_sva_nosahaSoft_vs_nosahaStiff_customNorm.png")
ggplot(data=res,
       aes(x=res$log2FoldChange,
           y=res$logPadj))+
  geom_point()+
  xlab("Log2(FC)")+
  ylab("-10*log10(P-adj)")+
  ggtitle(paste("100Pa, NO saha vs 2000Pa, NO saha\nTrue Reps -- single comparison, with SVAseq,",numsig,"peaks w/ padj<0.05"))
dev.off()
write.table(resfilt,file="true_rep_sva_0.05_diff_pairwise_customNorm.tsv",quote=FALSE,sep='\t',row.names=TRUE,col.names=TRUE)
