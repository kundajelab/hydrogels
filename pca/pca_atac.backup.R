rm(list=ls())
library(ggplot2)
library(sva)
library(data.table)
library(limma)
c25 <- c("dodgerblue2","#E31A1C", # red
         "green4",
         "#6A3D9A", # purple
         "#FF7F00", # orange
         "black","gold1",
         "skyblue2","#FB9A99", # lt pink
         "palegreen2",
         "#CAB2D6", # lt purple
         "#FDBF6F", # lt orange
         "gray70", "khaki2",
         "maroon","orchid1","deeppink1","blue1","steelblue4",
         "darkturquoise","green1","yellow4","yellow3",
         "darkorange4","brown")

data=data.frame(read.table('../age.atac.counts.txt.filtered',header=TRUE,sep='\t'))
rownames(data)=paste(data$Chrom,data$Start,data$End,sep="_")
data$Chrom=NULL
data$Start=NULL 
data$End=NULL
batches=data.frame(read.table("../batches.atac.txt.filtered",header=TRUE,sep='\t'))
rownames(batches)=batches$Replicate
batches$Day=factor(batches$Day)
batches$Age=factor(batches$Age)
batches$Sample=factor(batches$Sample)
batches$Batch=factor(batches$Batch)

mod0=model.matrix(~1,data=batches)
mod1=model.matrix(~Day+Age,data=batches)
v=voom(counts=data,design=mod1,normalize.method = "quantile")
sva.obj=sva(v$E,mod1,mod0,method="irw")
sur_var=data.frame(sva.obj$sv)

summary(lm(sva.obj$sv ~ batches$Age+batches$Day))
full.design.sv=cbind(mod1,sur_var)
v=voom(counts=data,design=full.design.sv)
fit <- lmFit(v)
sv_contribs=coefficients(fit)[,7:16] %*% t(fit$design[,7:16])
filtered=v$E-sv_contribs
write.table(filtered,"filtered_cpm_peaks.txt",quote=FALSE,sep='\t')

spearman_cor=cor(filtered,method="spearman")
write.csv(spearman_cor,"atac.spearman_cor.csv")
heatmap(spearman_cor,Rowv=NA,Colv=NA)
pearson_cor=cor(filtered,method="pearson")
heatmap(pearson_cor,Rowv=NA,Colv=NA)
write.csv(pearson_cor,"atac.pearson_cor.csv")
data.pca=prcomp(t(filtered),scale=FALSE,center=FALSE)
var_explained=as.character(round(100*data.pca$sdev^2/sum(data.pca$sdev^2),2))

png('atac.varexplained.png',height=10,width=10,units='in',res=300)
barplot(100*data.pca$sdev^2/sum(data.pca$sdev^2),las=2,ylab="% Variance Explained",xlab="Principal Component",ylim=c(0,100))
text(1:12,100*data.pca$sdev^2/sum(data.pca$sdev^2),labels=var_explained,pos=3,cex=1)
dev.off() 

day_labels=batches[rownames(data.pca$x),]$Day
age_labels=batches[rownames(data.pca$x),]$Age

#p1 vs p2 


png(filename = "age.atac.pca.1vs2.png",height=10, width=10, units='in', res=300)
plot(data.pca$x[,c(1,2)],col=c25[day_labels],pch=16)
text(data.pca$x[,c(1,2)],labels=rownames(data.pca$x),pos=3, cex=0.5)
title("PC1 vs PC2")
dev.off() 
png(filename = "age.atac.pca.2vs3.png",height=10, width=10, units='in', res=300)
plot(data.pca$x[,c(2,3)],col=c25[day_labels],pch=16)
text(data.pca$x[,c(2,3)],labels=rownames(data.pca$x),pos=3, cex=0.5)
title("PC2 vs PC3")
dev.off() 
png(filename = "age.atac.pca.1vs3.png",height=10, width=10, units='in', res=300)
plot(data.pca$x[,c(1,3)],col=c25[day_labels],pch=16)
text(data.pca$x[,c(1,3)],labels=rownames(data.pca$x),pos=3, cex=0.5)
title("PC1 vs PC3")
dev.off() 


#GET DIFFERENTIAL GENES ACROSS OLD VS YOUNG FOR SAME TIMEPOINT & BETWEEN SAME AGE FOR ADJACENT TIMEPOINTS

mod2=model.matrix(~0+Sample,data=batches)
full.design.sv=cbind(mod2,sur_var)
v=voom(counts=data,design=full.design.sv,normalize="quantile",plot=T)
fit <- lmFit(v)
cont.matrix=makeContrasts(d0_Y_vs_d0_O="Sampled0_Y - Sampled0_Old",
                          d1_Y_vs_d1_O="Sampled1_Y - Sampled1_Old",
                          d3_Y_vs_d3_O="Sampled3_Y - Sampled3_Old",
                          d5_Y_vs_d5_O="Sampled5_Y - Sampled5_Old",
                          d7_Y_vs_d7_O="Sampled7_Y - Sampled7_Old",
                          d1_Y_vs_d0_Y="Sampled1_Y - Sampled0_Y",
                          d3_Y_vs_d1_Y="Sampled3_Y - Sampled1_Y",
                          d5_Y_vs_d3_Y="Sampled5_Y - Sampled3_Y",
                          d7_Y_vs_d5_Y="Sampled7_Y - Sampled5_Y",
                          d1_O_vs_d0_O="Sampled1_Old - Sampled0_Old",
                          d3_O_vs_d1_O="Sampled3_Old - Sampled1_Old",
                          d5_O_vs_d3_O="Sampled5_Old - Sampled3_Old",
                          d7_O_vs_d5_O="Sampled7_Old - Sampled5_Old",
                          levels=full.design.sv)
fit2=contrasts.fit(fit,cont.matrix)
e=eBayes(fit2)
comparisons=colnames(cont.matrix)
for(i in seq(1,length(comparisons)))
{
  tab<-topTable(e, number=nrow(e),coef=i,lfc=1,p.value  = 0.01)
  names(tab)[1]=comparisons[i]
  tab$Gene=rownames(tab)
  write.table(tab,file=paste("differential_peaks_",comparisons[i],".tsv",sep=""),quote=TRUE,sep='\t',row.names = FALSE)
  png(paste("volcano_peaks",comparisons[i],".png",sep=""))
  volcanoplot(e,coef=i,highlight =0,names=rownames(tab),main=comparisons[i])
  dev.off() 
}
