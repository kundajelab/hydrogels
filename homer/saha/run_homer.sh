#!/bin/bash
for sample in diffpeaks.100Pa.2000Pa diffpeaks.saha_100Pa.100Pa diffpeaks.saha_100Pa.2000Pa diffpeaks.saha_100Pa.saha_2000Pa diffpeaks.saha_2000Pa.100Pa diffpeaks.saha_2000Pa.2000Pa
do
    findMotifsGenome.pl $sample.bed hg19 homer.$sample -bg background.bed & 
done

