#get a merged peak file
rm idr.optimal_set.bed
rm idr.optimal_set.sorted.bed
rm idr.optimal_set.sorted.merged.bed

for peak_file in `cut -f2 atac.peaks.idr.txt`
do
    zcat $peak_file >> idr.optimal_set.bed
done
bedtools sort -i idr.optimal_set.bed > idr.optimal_set.sorted.bed
bedtools merge -i idr.optimal_set.sorted.bed > idr.optimal_set.sorted.merged.bed
echo "generated merged peak file!"

#get the counts for each sample within each peak
export numfiles=`cat atac.tagalign.soft.stiff.saha.txt| wc -l`
echo $numfiles
for i in $(seq 1 $numfiles)
do
    cur_sample_name=`head -n $i atac.tagalign.soft.stiff.saha.txt | tail -n1 | cut -f1`
    echo $cur_sample_name > counts.$cur_sample_name.txt
    cur_tagalign=`head -n $i atac.tagalign.soft.stiff.saha.txt | tail -n1 | cut -f2`
    echo $cur_sample_name
    bedtools coverage -counts -a idr.optimal_set.sorted.merged.bed -b $cur_tagalign | cut -f4 >> counts.$cur_sample_name.txt
done
paste counts.*.txt > hydrogel.idr.atac.counts.txt
#clean up temporary files
rm counts.*.txt
echo "Chrom\\tStart\\tEnd" >> idr.optimal_set.sorted.merged.header.bed
cat idr.optimal_set.sorted.merged.bed >> idr.optimal_set.sorted.merged.header.bed
paste idr.optimal_set.sorted.merged.header.bed hydrogel.idr.atac.counts.txt > tmp
mv tmp hydrogel.idr.atac.counts.txt
