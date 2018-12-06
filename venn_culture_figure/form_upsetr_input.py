x2dsoft=open("2dsoft.idr.bed",'r').read().strip().split('\n')
x2dtcps=open("2dtcps.idr.bed",'r').read().strip().split('\n')
x3dsoft=open('3dsoft.idr.bed','r').read().strip().split('\n')
mamgland=open('mamgland.idr.bed','r').read().strip().split('\n')
outf=open('upsetr_inputs.txt','w')
outf.write('Peak\t2DSoft\t2DTCPS\t3DSoft\tInVivo\n')
peakset=set([]) 
x2dsoft_dict=dict()
x2dtcps_dict=dict()
x3dsoft_dict=dict()
mamgland_dict=dict()
for line in x2dsoft:
    line=line.replace('\t','_')
    peakset.add(line) 
    x2dsoft_dict[line]=1
for line in x2dtcps:
    line=line.replace('\t','_')
    peakset.add(line) 
    x2dtcps_dict[line]=1
for line in x3dsoft:
    line=line.replace('\t','_')
    peakset.add(line) 
    x3dsoft_dict[line]=1
for line in mamgland:
    line=line.replace('\t','_')
    peakset.add(line) 
    mamgland_dict[line]=1

peakset=list(peakset)
for entry in peakset:
    outf.write(entry)
    if entry in x2dsoft_dict:
        outf.write('\t1')
    else:
        outf.write('\t0')
    if entry in x2dtcps_dict:
        outf.write('\t1')
    else:
        outf.write('\t0')
    if entry in x3dsoft_dict:
        outf.write('\t1')
    else:
        outf.write('\t0')
    if entry in mamgland_dict:
        outf.write('\t1')
    else:
        outf.write('\t0')
    outf.write('\n')
    




