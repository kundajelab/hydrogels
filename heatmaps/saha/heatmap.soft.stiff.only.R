rm(list=ls())
library(gplots)
library(preprocessCore)
library(DESeq2)
require(gtools)
require(RColorBrewer)
library(cba)
library(limma)
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

inputf="hydrogel.soft.stiff.saha.limma.diffpeaks.txt"
#inputf="hydrogel.soft.stiff.saha.minus.sv.diffpeaks.txt"
#inputf="hydrogel.soft.stiff.saha.sample.diffpeaks.txt"
#inputf="hydrogel.soft.stiff.saha.diffpeaks.txt"
data=read.table(inputf,header=TRUE,sep='\t',stringsAsFactors = FALSE,row.names = 1)
data=as.matrix(data)
data=data[,1:6]
batches=data.frame(read.table("batches_with_svs.txt",header=TRUE,sep='\t'))
batches=batches[1:6,]
batches$Saha=NULL
batches$Stiffness=factor(batches$Stiffness)
batches$Sample=factor(batches$Sample)
mod1=model.matrix(~Sample,data=batches)
covariates=batches[c("sv1","sv2","sv3")]
data=removeBatchEffect(data,covariates = covariates,design=mod1)
#We split the fold change matrix into 1% quantiles 
quantile.range <- quantile(data, probs = seq(0, 1, 0.01))
#we scale the breaks in the heatmap color palette according to the quantiles. 
palette.breaks <- seq(quantile.range["5%"], quantile.range["95%"], 0.1)
indices=seq(1,1726,10)
palette.breaks=palette.breaks[indices]
colnames(data)=c("100Pa_Rep1","100Pa_Rep2","100Pa_Rep3","2000Pa_Rep1","2000Pa_Rep2","2000Pa_Rep3")

heatmap.2(data,
          scale     = "none",
          col       = rev(colorRampPalette(brewer.pal(10, "RdBu"))(length(palette.breaks) - 1)),
          distfun   = function(x) dist(x,method="euclidean"),
          hclustfun = function(x) hclust(x, method="ward.D"),
          Rowv=TRUE,
          Colv=FALSE,
          trace="none",
          density.info="none",
          cexCol = 0.9,
          margin=c(10,10),
          key.xlab = "Normalized Signal",
          breaks = palette.breaks,
          labRow="")

# heatmap.2(data,
#           scale     = "none",
#           col       = rev(colorRampPalette(brewer.pal(10, "RdBu"))(length(palette.breaks) - 1)),
#           distfun   = function(x) dist(x,method="euclidean"),
#           hclustfun = function(x) hclust(x, method="ward.D"),
#           Rowv=TRUE,
#           Colv=FALSE,
#           trace="none",
#           cexCol = 0.9,
#           margins=c(15,5),
#           breaks = palette.breaks,
#           labRow="")

#compute row z-scores 

zscores=t(apply(data,1,scale))
zdist=dist(zscores,method='euclidean')
zdist_clustered=hclust(zdist,method="ward.D")
ordered=order.optimal(zdist,zdist_clustered$merge)
sorted=zscores[ordered$order,]
colnames(sorted)=colnames(data)
rownames(sorted)=rownames(data)
quantile.range <- quantile(sorted, probs = seq(0, 1, 0.01))
#we scale the breaks in the heatmap color palette according to the quantiles. 
palette.breaks <- seq(quantile.range["5%"], quantile.range["95%"], 0.1)

colnames(sorted)=c("100Pa_Rep1","100Pa_Rep2","100Pa_Rep3","2000Pa_Rep1","2000Pa_Rep2","2000Pa_Rep3")
heatmap.2(sorted,
          scale     = "none",
          col       = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),
          Rowv=FALSE,
          Colv=FALSE,
          trace="none",
          cexCol = 0.9,
          margins=c(15,5),
          labRow="")
