bedtools intersect -v \
  -a CRAMP1_WT_merged.bed \
  -b CRAMP1_D1_merged.bed \
  > WT_only_vs_D1.bed

bedtools intersect -v \
  -a CRAMP1_WT_merged.bed \
  -b CRAMP1_SANT_merged.bed \
  > WT_only_vs_SANT.bed

#And find shared peaks

bedtools intersect -u \
  -a CRAMP1_WT_merged.bed \
  -b CRAMP1_D1_merged.bed \
  > WT_and_D1_shared.bed

bedtools intersect -u \
  -a CRAMP1_WT_merged.bed \
  -b CRAMP1_SANT_merged.bed \
  > WT_and_SANT_shared.bed

