#Pad to 350 bp 
python ~/anna_utils/seq_utils/pad.py --input_bed stiff.diffpeaks.bed --desired_length 350 --output_bed stiff.diffpeaks.padded.bed --chromsizes /mnt/data/annotations/by_release/hg19.GRCh37/hg19.chrom.sizes

#Extract the fasta sequence corresponding to the bed file
fastaFromBed -fi /mnt/data/annotations/by_release/hg19.GRCh37/hg19.genome.fa -bed stiff.diffpeaks.padded.bed -fo stiff.diffpeaks.padded.fa
