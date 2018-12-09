rm(list=ls())
data=read.table("motif_rank.tsv",header=TRUE,sep='\t')
data$logPvalue=-log10(data$Pvalue)
data$Rank=factor(data$Rank,levels=c("8","7","6","5","4","3","2","1"))
library(ggplot2)
ggplot(data=data,aes(x=data$Rank,
                     y=data$logPvalue))+
  geom_bar(fill='red',stat='identity')+
  xlab("Motif Rank")+
  ylab("motif enrichment\n-log10(p value)")+
  coord_flip()+
  theme_bw(20)


  