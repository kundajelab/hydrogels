for f in pc1.soft.stiff.highest.2.5.padded pc1.soft.stiff.highest.2.padded pc1.soft.stiff.highest.3.padded pc1.soft.stiff.highest.4.padded
do
    fastaFromBed -fi /mnt/data/annotations/by_release/hg19.GRCh37/hg19.genome.fa -bed $f.bed -fo $f.fa
    echo $f
done

