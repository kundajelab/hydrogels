import pandas as pd
data=pd.read_table("start_codons_of_HK.txt",header=None,sep='\t')
outf=open("5kb_upstream_start_codons_of_HK.txt",'w')
chrom=data[0]
start_codon=data[3]
shift=5000
bed_start=start_codon-shift
bed_end=start_codon
for i in range(len(chrom)):
    outf.write(chrom[i]+'\t'+str(bed_start[i])+'\t'+str(bed_end[i])+'\n')
    
