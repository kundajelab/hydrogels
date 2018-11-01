rm(list=ls())
library(ggplot2)
library(sva)
library(DESeq2)
library("BiocParallel")
parallelFlag=TRUE
register(MulticoreParam(20))

data=data.frame(read.table('background_subtracted.txt',header=TRUE,sep='\t'))
rownames(data)=paste(data$Chrom,data$Start,data$End,sep="_")
data$Chrom=NULL
data$Start=NULL 
data$End=NULL
data[data<0]=0
data=data[rowSums(data)>0,]

#load the normalization factors
normfactors=data.frame(read.table("norm_factors.txt",header=TRUE,sep='\t'))
normfactors=normfactors[1:nrow(data),]
rownames(normfactors)=rownames(data)

batches=data.frame(read.table("batches_hydrogel.soft.stiff.saha.txt",header=TRUE,sep='\t'))
rownames(batches)=batches$Replicate
batches$Stiffness=factor(batches$Stiffness)
batches$Saha=factor(batches$Saha)
batches$Sample=factor(batches$Sample)

#svaseq
mod0=model.matrix(~1,data=batches)
mod1=model.matrix(~Stiffness+Saha,data=batches)
sva.obj=svaseq(as.matrix(data),mod1,mod0)
sur_var=data.frame(sva.obj$sv)
names(sur_var)=c("sv1","sv2","sv3")
batches=cbind(batches,sur_var)
#drop the columns we don't want for analysis
batches$Stiffness=NULL
batches$Saha=NULL
batches$Replicate=NULL

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = batches,
                              design = ~Sample + sv1 + sv2 + sv3 )
			      
#provide the custom normalization factors
sizeFactors(dds)=c(0.0412995158,0.0273144907,0.0438329811,0.0240900163,0.0208288154,0.023301384,0.0157312261,0.020082931,0.0292282123,0.0275870297,0.0163430175,0.0208949729)
normalizationFactors(dds)=as.matrix(normfactors)

dds <- DESeq(dds,parallel = parallelFlag)

#extract the normalized counts for purposes of PCA
normcounts=counts(dds,normalized=TRUE)

#approach 1: remove SV's 
sv_contribs=coefficients(dds)[,5:7] %*% t(attr(dds,"modelMatrix")[,5:7])
filtered=normcounts - sv_contribs
write.table(filtered,"hydrogel.soft.stiff.saha.minus.sv.txt",quote=FALSE,sep='\t',col.names=TRUE,row.names=TRUE)

#approach 2: extract sample contribs
sample_contribs=coefficients(dds)[,1:4] %*% t(attr(dds,"modelMatrix")[,1:4])
write.table(sample_contribs,"hydrogel.soft.stiff.saha.sample.txt",quote=FALSE,sep='\t',col.names=TRUE,row.names=TRUE)

#approach 3: apply limma batch effect removal
library(limma)
tokeep=model.matrix(~Sample,data=batches)
normalized=removeBatchEffect(normcounts,design=tokeep)
write.table(normalized,"hydrogel.soft.stiff.saha.limma.txt",quote=FALSE,sep='\t',col.names=TRUE,row.names=TRUE)



#get the results for various contrasts
sample1=c("100Pa","saha_100Pa","saha_100Pa","saha_2000Pa","saha_100Pa","saha_2000Pa")
sample2=c("2000Pa","saha_2000Pa","2000Pa","100Pa","100Pa","2000Pa")
numcomparisons=length(sample1)
for(i in seq(1,numcomparisons)){
	res=results(dds,contrast=c("Sample",sample1[i],sample2[i]),parallel=TRUE)
	res$logPadj=-10*log10(res$padj)
	res=as.data.frame(res)
	res=na.omit(res)
	numsig=sum(res$padj <= 0.05)
	sigsubset=res[res$padj<=0.05,]
	outtable=paste("diffpeaks",sample1[i],sample2[i],"tsv",sep='.')
	write.table(sigsubset,file=outtable,quote=FALSE,sep='\t',row.names=TRUE,col.names=TRUE)
	print(paste(sample1[i],sample2[i],numsig))
	outpng=paste("volcano",sample1[i],sample2[i],"svg",sep='.')
	outlabel=paste(sample1[i],'vs.',sample2[i],"Diff. peaks:",numsig)
	res$color=res$padj<=0.05
	svg(outpng,width=5,height=5,pointsize=12)
	ggplot(data=res,
       	aes(x=res$log2FoldChange,
		 y=res$logPadj,
		 color=res$color))+
		 geom_point()+	   
		 xlab("Log2(FC)")+
		 ylab("-10*log10(P-adj)")+
  		 ggtitle(outlabel)+
		 scale_color_manual(values=c("#000000","#ff0000"),name="P.adj<0.05")
		 dev.off()
}