#Pad to 350 bp 
python ~/anna_utils/seq_utils/pad.py --input_bed soft.diffpeaks.bed --desired_length 350 --output_bed soft.diffpeaks.padded.bed --chromsizes /mnt/data/annotations/by_release/hg19.GRCh37/hg19.chrom.sizes

#Extract the fasta sequence corresponding to the bed file
fastaFromBed -fi /mnt/data/annotations/by_release/hg19.GRCh37/hg19.genome.fa -bed soft.diffpeaks.padded.bed -fo soft.diffpeaks.padded.fa
