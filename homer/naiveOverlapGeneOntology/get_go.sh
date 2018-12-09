for sample in 100Pa_3D_ppr.naive_overlap.filt.narrowPeak.bed 2000Pa_3D_ppr.naive_overlap.filt.narrowPeak.bed
do    
    annotatePeaks.pl $sample hg19 -genomeOntology genomeOntology.$sample
done



