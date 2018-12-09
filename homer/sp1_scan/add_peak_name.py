data=open("hydrogel.soft.stiff.saha.bed",'r').read().strip().split('\n')
outf=open('tmp.bed','w')
for line in data:
    outf.write(line+'\t'+'_'.join(line.split('\t'))+'\n')
    
