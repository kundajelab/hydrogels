#Peaks that go down 
#Pad to 350 bp 
python ~/anna_utils/seq_utils/pad.py --input_bed diff.down.with.saha.bed --desired_length 350 --output_bed stiff.diffpeaks.down.padded.bed --chromsizes /mnt/data/annotations/by_release/hg19.GRCh37/hg19.chrom.sizes

#Extract the fasta sequence corresponding to the bed file
fastaFromBed -fi /mnt/data/annotations/by_release/hg19.GRCh37/hg19.genome.fa -bed stiff.diffpeaks.down.padded.bed -fo stiff.diffpeaks.down.padded.fa


#Peaks that go up 
#Pad to 350 bp 
python ~/anna_utils/seq_utils/pad.py --input_bed diff.up.with.saha.bed --desired_length 350 --output_bed stiff.diffpeaks.up.padded.bed --chromsizes /mnt/data/annotations/by_release/hg19.GRCh37/hg19.chrom.sizes

#Extract the fasta sequence corresponding to the bed file
fastaFromBed -fi /mnt/data/annotations/by_release/hg19.GRCh37/hg19.genome.fa -bed stiff.diffpeaks.up.padded.bed -fo stiff.diffpeaks.up.padded.fa
