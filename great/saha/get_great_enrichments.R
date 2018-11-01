rm(list=ls())
tasks=c("diffpeaks.100Pa.2000Pa.bed",
	"diffpeaks.saha_100Pa.2000Pa.bed",
	"diffpeaks.saha_2000Pa.100Pa.bed",
	"diffpeaks.saha_100Pa.100Pa.bed",
	"diffpeaks.saha_100Pa.saha_2000Pa.bed",
	"diffpeaks.saha_2000Pa.2000Pa.bed")
background=read.table("background.bed",header=FALSE,sep='\t')
names(background)=c("chr","start","end")
library(rGREAT,quietly=TRUE)
set.seed(123)
for(taskname in tasks){
	     bed=read.table(taskname,header=FALSE,sep='\t')
	     names(bed)=c("chr","start","end")
	     job=submitGreatJob(bed,species="hg19",bg=background,request_interval = 30)
	     #browser()
	     #necessary to avoid the job query timing out. Also from rstudio use ctl+shift+enter to execute with echo. 
	     tb = getEnrichmentTables(job)
	     names(tb)
	     go_tables = getEnrichmentTables(job, category = c("GO"))
	     go_mol=go_tables[1]$`GO Molecular Function`
	     go_proc=go_tables[2]$`GO Biological Process`
	     go_cell=go_tables[3]$`GO Cellular Component`

	     #Apply filters on 

	     go_mol=go_mol[go_mol$Hyper_Adjp_BH < 0.2,]
	     #go_mol=go_mol[go_mol$Binom_Adjp_BH < 0.2,]

	     go_proc=go_proc[go_proc$Hyper_Adjp_BH<0.2,]
	     #go_proc=go_proc[go_proc$Binom_Adjp_BH<0.2,]

	     go_cell=go_cell[go_cell$Hyper_Adjp_BH<0.2,]
	     #go_cell=go_cell[go_cell$Binom_Adjp_BH<0.2,]

pathway_tables=getEnrichmentTables(job,category=c("Pathway Data"))
panther=pathway_tables[1]$`PANTHER Pathway`
biocyc=pathway_tables[2]$`BioCyc Pathway`
msigdb=pathway_tables[3]$`MSigDB Pathway`

panther=panther[panther$Hyper_Adjp_BH<0.2,]
#panther=panther[panther$Binom_Adjp_BH<0.2,]

biocyc=biocyc[biocyc$Hyper_Adjp_BH<0.2,]
#biocyc=biocyc[biocyc$Binom_Adjp_BH<0.2,]

msigdb=msigdb[msigdb$Hyper_Adjp_BH<0.2,]
#msigdb=msigdb[msigdb$Binom_Adjp_BH<0.2,]

#combine all results 
sig_results=rbind(panther,biocyc,msigdb,go_mol,go_proc,go_cell)
write.csv(sig_results,file=paste(taskname,'great',sep='.'),quote=TRUE,sep='\t')
}
