#!/bin/bash
#for sample in pc1.soft.stiff.highest.2.5.padded.bed pc1.soft.stiff.highest.2.padded.bed pc1.soft.stiff.highest.3.padded.bed pc1.soft.stiff.highest.4.padded.bed sp1.bed combo.bed
#do
#    findMotifsGenome.pl $sample hg19 homer.$sample -bg background.bed &
#    findMotifsGenome.pl $sample hg19 homer.nobackground.$sample &
#done


#for sample in  pc1.soft.stiff.highest.2.padded.bed
#do
    #findMotifsGenome.pl $sample hg19 homer.100.nobackground.$sample -size 100
#    findMotifsGenome.pl $sample hg19 homer.given.nobackground.$sample -size given
#done


for sample in  combo.bed
do
    findMotifsGenome.pl $sample hg19 homer.$sample -bg background.bed  &
    findMotifsGenome.pl $sample hg19 homer.nobackground.$sample &
    findMotifsGenome.pl $sample hg19 homer.100.nobackground.$sample -size 100 &
    
done
