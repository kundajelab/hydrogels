for f in mam.2d.soft mam.2d.tcps mam.3d
do
    annotatePeaks.pl $f.bed hg19 -go go.$f 
done
