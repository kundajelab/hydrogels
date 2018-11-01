rm(list=ls())
library(ggplot2)
library(sva)
library(data.table)
library(DESeq2)
library("BiocParallel")
parallelFlag=TRUE
register(MulticoreParam(4))

data=data.frame(read.table('hydrogel.soft.stiff.saha.ppr.txt',header=TRUE,sep='\t'))
rownames(data)=paste(data$Chrom,data$Start,data$End,sep="_")
data$Chrom=NULL
data$Start=NULL 
data$End=NULL
#data=voom(counts=data,normalize.method = "quantile")

batches=data.frame(read.table("batches_hydrogel.soft.stiff.saha.ppr.txt",header=TRUE,sep='\t'))
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
ggplot(data=res,
       aes(x=res$log2FoldChange,
           y=res$logPadj))+
  geom_point()+
  xlab("Log2(FC)")+
  ylab("-10*log10(P-adj)")+
  ggtitle("100Pa, NO saha vs 2000Pa, NO saha\nPseudo Reps, no SVAseq, 99 peaks w/ padj<0.05")

