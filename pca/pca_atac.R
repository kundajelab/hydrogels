rm(list=ls())
library(ggplot2)
library(sva)
library(limma)

#read all data
alldata=read.table('hydrogel.soft.stiff.saha.txt',header=TRUE,sep='\t')
rownames(alldata)=paste(alldata$Chrom,alldata$Start,alldata$End,sep='_')
alldata$Chrom=NULL
alldata$Start=NULL 
alldata$End=NULL

#use limma voom to normalize the read counts 
alldata=voom(alldata,normalize.method = "quantile")

data.pca=prcomp(t(alldata$E),scale=FALSE,center=FALSE)
var_explained=round(100*data.pca$sdev^2/sum(data.pca$sdev^2),2)
var_explained_df=data.frame(seq(1,length(var_explained)),var_explained)
names(var_explained_df)=c("PC","VarianceExplained")
p1=ggplot(data=var_explained_df,
       aes(x=var_explained_df$PC,
           y=var_explained_df$VarianceExplained,
           label=as.character(var_explained_df$VarianceExplained)))+
  geom_bar(stat="identity")+
  geom_text(size=3,color='red')+
  xlab("PC")+
  ylab("% Variance Explained")

pca_df=data.frame(data.pca$x)
p2=ggplot(data=pca_df,aes(x=pca_df$PC1,y=pca_df$PC2,label=rownames(pca_df)))+
  geom_point() +
  geom_text(color="#000000",alpha=0.5,size=3)+
  xlab("PC1")+
  ylab("PC2")+
  theme_bw(20)

p3=ggplot(data=pca_df,aes(x=pca_df$PC2,y=pca_df$PC3,label=rownames(pca_df)))+
  geom_point() +
  geom_text(color="#000000",alpha=0.5,size=3)+
  xlab("PC2")+
  ylab("PC3")+
  theme_bw(20)

p4=ggplot(data=pca_df,aes(x=pca_df$PC1,y=pca_df$PC3,label=rownames(pca_df)))+
  geom_point() +
  geom_text(color="#000000",alpha=0.5,size=3)+
  xlab("PC1")+
  ylab("PC3")+
  theme_bw(20)

source("~/helpers.R")
multiplot(p1,p2,p3,p4,cols=2)

