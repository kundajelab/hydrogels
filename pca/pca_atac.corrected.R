rm(list=ls())
library(ggplot2)
library(sva)
library(limma)

#read normalized data from deseq2
data=read.table('hydrogel.soft.stiff.saha.limma.txt',header=TRUE,sep='\t',row.names=1)
#remove batch effects 
batches=data.frame(read.table("batches_with_svs.txt",header=TRUE,sep='\t'))
batches$Stiffness=factor(batches$Stiffness)
batches$Saha=factor(batches$Saha)
mod1=model.matrix(~Sample,data=batches)
covariates=batches[c("sv1","sv2","sv3")]
data=removeBatchEffect(data,covariates = covariates,design=mod1)



data.pca=prcomp(t(data),scale=TRUE,center=FALSE)
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
pca_df$Saha=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)
pca_df$Stiffness=c("Soft","Soft","Soft","Stiff","Stiff","Stiff","Soft","Soft","Soft","Stiff","Stiff","Stiff")
pca_df$Saha=factor(pca_df$Saha)
p2=ggplot(data=pca_df,aes(x=pca_df$PC1,y=pca_df$PC2,color=pca_df$Stiffness,shape=pca_df$Saha))+
  geom_point(size=10) +
  xlab("PC1: 48.13% Var. Explained")+
  ylab("PC2: 19.86% Var. Explained")+
  theme_bw(20)+
  scale_color_manual(values=c("#fc8d62","#66c2a5"),name="Soft v. Stiff")+
  scale_shape_manual(values=c("o","x"),name="SAHA vs No SAHA")

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

