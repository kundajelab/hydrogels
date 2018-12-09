rm(list=ls())
library(ggplot2)
data=read.table("summary.txt",header=TRUE,sep='\t')
data$Position=data$Position-100
data$Soft=data$soft.neg.cuts.bed+data$soft.pos.cuts.bed
data$Stiff=data$stiff.neg.cuts.bed+data$stiff.pos.cuts.bed
data$Stiff=data$Stiff*1.47
ggplot(data=data)+
  geom_line(aes(x=data$Position,y=data$Soft,color="Soft"))+
  geom_line(aes(x=data$Position,y=data$Stiff,color="Stiff"))+
  scale_color_manual(values=c('#000000','#e41a1c'))+
  xlab("Distance from Motif Center")+
  ylab("Average Cuts")+
  theme_bw()+
  ggtitle("Sp1")+
  ylim(c(0,2.0))
subset_data=subset(data,select=c("Position","Soft","Stiff"))
write.table(subset_data,file="sp1.footprint.tsv",row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')