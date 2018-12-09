for f in idr.fa reverted.fa	 
do
    meme-chip -db /mnt/data/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme -oc meme.$f $f -bfile background.markov
done
