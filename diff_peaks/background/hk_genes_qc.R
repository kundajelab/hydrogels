rm(list=ls())
library(ggplot2)
library(limma)
data=read.table("hydrogel.soft.stiff.saha.HK.txt",header=FALSE,sep='\t')
data=data[,4:15]
normd=voom(counts=data,normalize.method = "quantile")$E
normd=as.data.frame(normd)

norm_soft=rowMeans(normd[,1:3])
norm_stiff=rowMeans(normd[,4:6])
norm_saha_soft=rowMeans(normd[,7:9])
norm_saha_stiff=rowMeans(normd[,10:12])

norm_collapsed=data.frame(norm_soft,norm_stiff,norm_saha_soft,norm_saha_stiff)
names(norm_collapsed)=c("SoftNoSAHA","StiffNoSAHA",'SoftSAHA','StiffSAHA')
data=norm_collapsed
p1= ggplot(data=data)+
    geom_point(aes(x=data$SoftNoSAHA,y=data$StiffNoSAHA,alpha=0.1))+
    ggtitle("(Soft No SAHA) vs (Stiff No SAHA)")+
    xlab("Soft, No SAHA")+
    ylab("Stiff, No SAHA")

p2= ggplot(data=data)+
  geom_point(aes(x=data$SoftSAHA,y=data$StiffSAHA,alpha=0.1))+
  ggtitle("(Soft SAHA) vs (Stiff SAHA)")+
  xlab("Soft, SAHA")+
  ylab("Stiff, SAHA")


p3= ggplot(data=data)+
  geom_point(aes(x=data$SoftNoSAHA,y=data$SoftSAHA,alpha=0.1))+
  ggtitle("(Soft No SAHA) vs (Soft SAHA)")+
  xlab("Soft, No SAHA")+
  ylab("Soft, SAHA")


p4= ggplot(data=data)+
  geom_point(aes(x=data$StiffNoSAHA,y=data$StiffSAHA,alpha=0.1))+
  ggtitle("(Stiff No SAHA) vs (Stiff SAHA)")+
  xlab("Stiff, No SAHA")+
  ylab("Stiff, SAHA")

source('~/helpers.R')
multiplot(p1,p2,p3,p4,cols=2)


