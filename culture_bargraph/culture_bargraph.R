rm(list=ls())
library(ggplot2)
library(reshape2)
data=read.table("bargraph_inputs.csv",header=TRUE,sep='\t')
data$Term=factor(data$Term,levels=data$Term)
names(data)=c("Term","2D Soft + In Vivo","2D TCPS + In Vivo","3D Soft + In Vivo")
data=melt(data)
p1=ggplot(data=data,
          aes(x=data$Term,
              y=data$value,
              group=data$variable,
              fill=data$variable))+
  geom_bar(stat='identity',position='dodge')+
  xlab("")+
  ylab("-log10(Enrichment)")+
  scale_fill_manual(name="",values=c('#d95f02','#7570b3','#1b9e77'))+
  coord_flip()+
  theme_bw()

svg("culture_bargraphs.svg")
print(p1)
dev.off() 
print(p1)