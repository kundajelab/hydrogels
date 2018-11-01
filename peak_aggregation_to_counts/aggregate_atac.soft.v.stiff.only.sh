#get a merged peak file
#rm naive_overlap.optimal_set.bed
#rm naive_overlap.optimal_set.sorted.bed
#rm naive_overlap.optimal_set.sorted.merged.bed

for peak_file in `cut -f2 atac.peaks.naiveo.soft.v.stiff.txt`
do
    zcat $peak_file >> naive_overlap.optimal_set.soft.v.stiff.bed
done
bedtools sort -i naive_overlap.optimal_set.soft.v.stiff.bed > naive_overlap.optimal_set.soft.v.stiff.sorted.bed
bedtools merge -i naive_overlap.optimal_set.sorted.bed > naive_overlap.optimal_set.soft.v.stiff.sorted.merged.bed

echo "generated merged peak file!"

#get the counts for each sample within each peak
export numfiles=`cat atac.tagalign.soft.v.stiff.txt| wc -l`
echo $numfiles
for i in $(seq 1 $numfiles)
do
    cur_sample_name=`head -n $i atac.tagalign.soft.v.stiff.txt | tail -n1 | cut -f1`
    echo $cur_sample_name > counts.$cur_sample_name.txt
    cur_tagalign=`head -n $i atac.tagalign.soft.v.stiff.txt | tail -n1 | cut -f2`
    echo $cur_sample_name
    bedtools coverage -counts -a naive_overlap.optimal_set.soft.v.stiff.sorted.merged.bed -b $cur_tagalign | cut -f4 >> counts.$cur_sample_name.txt
done
paste counts.*.txt > hydrogel.atac.soft.v.stiff.counts.txt
#clean up temporary files
rm counts.*.txt
paste naive_overlap.optimal_set.soft.v.stiff.sorted.merged.bed hydrogel.atac.soft.v.stiff.counts.txt > tmp
mv tmp hydrogel.atac.soft.v.stiff.counts.txt
