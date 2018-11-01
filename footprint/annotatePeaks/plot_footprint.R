rm(list=ls())
library(ggplot2)
library(reshape2)

soft=read.table("soft.txt",header=TRUE,sep='\t')[,1:2]
stiff=read.table("stiff.txt",header=TRUE,sep='\t')[,2]
diff=read.table("diffpeaks.txt",header=TRUE,sep='\t')[,2]
data=data.frame(soft,stiff,diff)
names(data)=c("Position","Soft","Stiff","Diff")

ggplot(data=data)+
  geom_line(aes(x=data$Position,y=data$Soft,color="Soft"))+
  geom_line(aes(x=data$Position,y=data$Stiff,color="Stiff"))+
  geom_line(aes(x=data$Position,y=data$Diff,color="Differential Peaks"))+
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'))+
  xlab("Distance from Peak")+
  ylab("Sites per base per peak")+
  theme_bw(20)+
  ggtitle("Sp1 Binding Motif")


