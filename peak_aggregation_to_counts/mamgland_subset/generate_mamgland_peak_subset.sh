for peak_file in `cut -f2 atac.peaks.naiveo.txt`
do
    zcat $peak_file >> naive_overlap.optimal_set.mamgland.bed
done
bedtools sort -i naive_overlap.optimal_set.mamgland.bed > naive_overlap.optimal_set.sorted.mamgland.bed
bedtools merge -i naive_overlap.optimal_set.sorted.mamgland.bed > naive_overlap.optimal_set.sorted.merged.mamgland.bed

echo "generated merged peak file!"

