#We want to calculate the motif enrichment in the subset of peaks that did not revert to soft-like phenotypes after treatment with SAHA.
bedtools intersect -v -a finalized.soft.stiff.diff.peaks.pc1.soft.stiff.highest.2.5.bed -b reverted.bed > nonreverted.bed
#Pad to 350 bp 
python ~/anna_utils/seq_utils/pad.py --input_bed nonreverted.bed --desired_length 350 --output_bed nonreverted.padded.bed --chromsizes /mnt/data/annotations/by_release/hg19.GRCh37/hg19.chrom.sizes

#Extract the fasta sequence corresponding to the bed file
fastaFromBed -fi /mnt/data/annotations/by_release/hg19.GRCh37/hg19.genome.fa -bed nonreverted.padded.bed -fo nonreverted.padded.fa
