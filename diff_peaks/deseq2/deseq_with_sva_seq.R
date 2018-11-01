rm(list=ls())
library(ggplot2)
library(sva)
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

#svaseq
data=data[rowSums(data)>0,]
mod0=model.matrix(~1,data=batches)
mod1=model.matrix(~0+Stiffness+Saha,data=batches)
sva.obj=svaseq(as.matrix(data),mod1,mod0)
sur_var=data.frame(sva.obj$sv)
names(sur_var)=c("sv1")
batches=cbind(batches,sur_var)

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = batches,
                              design = ~ Sample + sv1)

dds <- DESeq(dds,parallel = parallelFlag)

#get the results for various contrasts 
res=results(dds,contrast=c("Sample","100Pa","2000Pa"),parallel=TRUE)
res$logPadj=-10*log10(res$padj)
res=as.data.frame(res)
resfilt=na.omit(res)
numsig=sum(resfil$pad < 0.05)
png("deseq_pseudorep_sva_nosahaSoft_vs_nosahaStiff.png")
ggplot(data=res,
       aes(x=res$log2FoldChange,
           y=res$logPadj))+
  geom_point()+
  xlab("Log2(FC)")+
  ylab("-10*log10(P-adj)")+
  ggtitle(paster("100Pa, NO saha vs 2000Pa, NO saha\nPseudo Reps, with SVAseq,",numsig,"peaks w/ padj<0.05"))
dev.off()

