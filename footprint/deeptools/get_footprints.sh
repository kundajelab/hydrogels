#run compute matrix to collect the data needed for plotting
#computeMatrix scale-regions --samplesLabel "Soft" "Stiff" \
#	      -p 30 \
#	      -S /srv/scratch/annashch/hydrogels/100Pa_3D/out/signal/macs2/pooled_rep/100Pa_1_TAAGGCGA_1_pf.PE2SE.nodup.tn5_pooled.pf.pval.signal.bigwig /srv/scratch/annashch/hydrogels/2000Pa_3D/out/signal/macs2/pooled_rep/2000Pa_1_TCCTGAGC_1_pf.nodup.tn5_pooled.pf.pval.signal.bigwig \
#	      -R sp1.intersection.bed \
#	      -o sp1.2.matrix.mat.gz \
#	      --binSize=2 \
#	      --beforeRegionStartLength 150 \
#	      --afterRegionStartLength 150 \
#	      --regionBodyLength 200

#plotHeatmap -m sp1.2.matrix.mat.gz -out sp1.2.svg

#
computeMatrix reference-point --samplesLabel "Soft" "Stiff" \
	      -p 30 \
	      -S /srv/scratch/annashch/hydrogels/100Pa_3D/out/signal/macs2/pooled_rep/100Pa_1_TAAGGCGA_1_pf.PE2SE.nodup.tn5_pooled.pf.pval.signal.bigwig /srv/scratch/annashch/hydrogels/2000Pa_3D/out/signal/macs2/pooled_rep/2000Pa_1_TCCTGAGC_1_pf.nodup.tn5_pooled.pf.pval.signal.bigwig \
	      -R sp1.intersection.pos.bed \
	      --referencePoint center \
	      -o sp1.matrix.pos.mat.gz \
	      --upstream 500 \
	      --downstream 500 \
	      --binSize=1

plotHeatmap -m sp1.matrix.pos.mat.gz -out sp1.pos.svg

computeMatrix reference-point --samplesLabel "Soft" "Stiff" \
	      -p 30 \
	      -S /srv/scratch/annashch/hydrogels/100Pa_3D/out/signal/macs2/pooled_rep/100Pa_1_TAAGGCGA_1_pf.PE2SE.nodup.tn5_pooled.pf.pval.signal.bigwig /srv/scratch/annashch/hydrogels/2000Pa_3D/out/signal/macs2/pooled_rep/2000Pa_1_TCCTGAGC_1_pf.nodup.tn5_pooled.pf.pval.signal.bigwig \
	      -R sp1.intersection.neg.bed \
	      --referencePoint center \
	      -o sp1.matrix.neg.mat.gz \
	      --upstream 500 \
	      --downstream 500 \
	      --binSize=1

plotHeatmap -m sp1.matrix.neg.mat.gz -out sp1.neg.svg


