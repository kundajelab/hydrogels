for f in pc1.soft.stiff.highest.2.5.padded.fa pc1.soft.stiff.highest.2.padded.fa pc1.soft.stiff.highest.3.padded.fa pc1.soft.stiff.highest.4.padded.fa
do
    meme-chip -db /mnt/data/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme -o meme.$f $f
done
