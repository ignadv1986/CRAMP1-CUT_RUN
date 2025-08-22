computeMatrix reference-point \
  --referencePoint TSS \
  -b 1000 -a 1000 \
  -R histones_TSS.bed H1_TSS.bed \
  -S CRAMP1_WT_merged.bw CRAMP1_SANT_merged.bw CRAMP1_D1_merged.bw \
  --samplesLabel WT SANT D1 \
  --skipZeros \
  --missingDataAsZero \
  -o TSS_matrix.gz

plotProfile -m TSS_matrix.gz -out TSS_profile.png --perGroup --samplesLabel WT SANT D1 --plotTitle "Signal around H1 TSS ±1kb"

