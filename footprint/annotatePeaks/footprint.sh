annotatePeaks.pl 100Pa_3D_ppr.naive_overlap.filt.narrowPeak hg19 -m sp1.motif -size 1000 -hist 5 > soft.txt
annotatePeaks.pl 2000Pa_3D_ppr.naive_overlap.filt.narrowPeak hg19 -m sp1.motif -size 1000 -hist 5 > stiff.txt
annotatePeaks.pl diffpeaks.100Pa.2000Pa.bed hg19 -m sp1.motif -size 1000 -hist 5 > diffpeaks.txt


