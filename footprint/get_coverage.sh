# stiff positive 
bedtools coverage -a sp1.intersection.pos.bed -b 2000Pa_3D.pos.5p.bed -d  > stiff.pos.cuts.bed &

#stiff negative
bedtools coverage -a sp1.intersection.neg.bed -b 2000Pa_3D.neg.3p.bed -d > stiff.neg.cuts.bed &


# soft positive 
bedtools coverage -a sp1.intersection.pos.bed -b 100Pa_3D.pos.5p.bed -d > soft.pos.cuts.bed &

#soft negative
bedtools coverage -a sp1.intersection.neg.bed -b 100Pa_3D.neg.3p.bed -d > soft.neg.cuts.bed &
