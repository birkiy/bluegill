#!/bin/bash


wget https://www.encodeproject.org/files/ENCFF331JNU/@@download/ENCFF331JNU.bigWig # A549 GR - 0h DEX - bigwig
wget https://www.encodeproject.org/files/ENCFF809IAK/@@download/ENCFF809IAK.bigWig # A549 GR - 4h DEX - bigwig 

wget https://www.encodeproject.org/files/ENCFF995VAB/@@download/ENCFF995VAB.bed.gz # narrowPeaks - 4h DEX - bigwig

gunzip ENCFF995VAB.bed.gz


cut -f1-3 ENCFF995VAB.bed \
| awk '{print $0"\tGR"}' \
> GR.hg38.bed


cat data/GR.hg38.bed <(shuf -n 10000 random.hg38.bed ) \
| sort -k1,1 -k2,2n \
| grep -v "_" \
> data/ALL.hg38.bed

