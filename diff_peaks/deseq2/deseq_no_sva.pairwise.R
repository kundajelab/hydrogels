rm(list=ls())
library(ggplot2)
library(sva)
library(data.table)
library(DESeq2)
library("BiocParallel")
parallelFlag=TRUE
register(MulticoreParam(4))

#load the read counts
data=data.frame(read.table('hydrogel.atac.soft.v.stiff.counts.txt',header=TRUE,sep='\t'))
rownames(data)=paste(data$Chrom,data$Start,data$End,sep="_")
data$Chrom=NULL
data$Start=NULL 
data$End=NULL


batches=data.frame(read.table("batches_hydrogel.soft.v.stiff.txt",header=TRUE,sep='\t'))
rownames(batches)=batches$Replicate
batches$Stiffness=factor(batches$Stiffness)
batches$Saha=factor(batches$Saha)
batches$Sample=factor(batches$Sample)

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = batches,
                              design = ~ Sample)

dds <- DESeq(dds,parallel = parallelFlag)

#get the results for various contrasts 
res=results(dds,contrast=c("Sample","100Pa","2000Pa"),parallel=TRUE)
res$logPadj=-10*log10(res$padj)
res=as.data.frame(res)
resfilt=na.omit(res)
numsig=sum(resfilt$padj<0.05)
ggplot(data=res,
       aes(x=res$log2FoldChange,
           y=res$logPadj))+
  geom_point()+
  xlab("Log2(FC)")+
  ylab("-10*log10(P-adj)")+
  ggtitle(paste("100Pa, NO saha vs 2000Pa, NO saha\nTrue Reps -- single comparison, no SVAseq,",numsig,"peaks w/ padj<0.05"))

