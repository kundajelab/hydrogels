cut -f1,2,3 ../../diff_peaks/hydrogel.soft.stiff.saha.txt > background.bed

for sample in diffpeaks.100Pa.2000Pa diffpeaks.saha_100Pa.100Pa diffpeaks.saha_100Pa.2000Pa diffpeaks.saha_100Pa.saha_2000Pa diffpeaks.saha_2000Pa.100Pa diffpeaks.saha_2000Pa.2000Pa
do
    cut -f1 ../../diff_peaks/deseq2_custom_norm/$sample.tsv | sed --expression='s/\_/\t/g' | sort | grep -v "baseMean" | bedtools sort -i stdin > $sample.bed
    echo $sample
done
