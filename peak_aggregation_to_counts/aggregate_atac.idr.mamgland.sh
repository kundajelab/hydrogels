#get a merged peak file
rm idr.optimal_set.mamgland.bed
rm idr.optimal_set.mamgland.sorted.bed
rm idr.optimal_set.mamgland.sorted.merged.bed

for peak_file in `cut -f2 mamgland.atac.peaks.idr.txt`
do
    zcat $peak_file >> idr.optimal_set.mamgland.bed
done
bedtools sort -i idr.optimal_set.mamgland.bed > idr.optimal_set.mamgland.sorted.bed
bedtools merge -i idr.optimal_set.mamgland.sorted.bed > idr.optimal_set.mamgland.sorted.merged.bed
echo "generated merged peak file!"

#get the counts for each sample within each peak
export numfiles=`cat atac.tagalign.mamgland.txt| wc -l`
echo $numfiles
for i in $(seq 1 $numfiles)
do
    cur_sample_name=`head -n $i atac.tagalign.mamgland.txt | tail -n1 | cut -f1`
    echo $cur_sample_name > counts.$cur_sample_name.txt
    cur_tagalign=`head -n $i atac.tagalign.mamgland.txt | tail -n1 | cut -f2`
    echo $cur_sample_name
    bedtools coverage -counts -a idr.optimal_set.mamgland.sorted.merged.bed -b $cur_tagalign | cut -f4 >> counts.$cur_sample_name.txt
done
paste counts.*.txt > mamgland.idr.atac.counts.txt
#clean up temporary files
rm counts.*.txt
echo "Chrom\\tStart\\tEnd" >> idr.optimal_set.mamgland.sorted.merged.header.bed
cat idr.optimal_set.mamgland.sorted.merged.bed >> idr.optimal_set.mamgland.sorted.merged.header.bed
paste idr.optimal_set.mamgland.sorted.merged.header.bed mamgland.idr.atac.counts.txt > tmp
mv tmp mamgland.idr.atac.counts.txt
