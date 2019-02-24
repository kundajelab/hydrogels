data=open("finalized.soft.stiff.diff.peaks.pc1.soft.stiff.highest.2.5.bed",'r').read().strip().split('\n')
outf=open("finalized.soft.stiff.diff.peaks.pc1.soft.stiff.highest.2.5.with.names.bed",'w')
for line in data:
    tokens=line.split('\t')
    outf.write(line+'\t'+'_'.join(tokens)+'\n')
    
