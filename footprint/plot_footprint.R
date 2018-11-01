rm(list=ls())
library(ggplot2)
data=read.table("summary.txt",header=TRUE,sep='\t')
data$Position=data$Position-100
data$Soft=data$soft.5p+data$soft.3p
data$Stiff=data$stiff.5p+data$stiff.3p
data$Stiff=data$Stiff*1.47
ggplot(data=data)+
  geom_line(aes(x=data$Position,y=data$Soft,color="Soft"))+
  geom_line(aes(x=data$Position,y=data$Stiff,color="Stiff"))+
  scale_color_manual(values=c('#000000','#e41a1c'))+
  xlab("Distance from Motif Center")+
  ylab("Average Cuts")+
  theme_bw(20)+
  ggtitle("Sp1")+
  ylim(c(0,2.0))
