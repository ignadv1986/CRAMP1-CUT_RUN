#!/bin/bash
set -e

# Define the TSS set (change this filename for each script)
TSS_SET="H1_TSS.bed"
TSS_NAME=$(basename "$TSS_SET" .bed)

# Samples chosen (replicates only)
samples=(
  "CRAMP1_WT_R1_peaks_peaks.narrowPeak"
  "CRAMP1_WT_R3_peaks_peaks.narrowPeak"
  "CRAMP1_SANT_R1_peaks_peaks.narrowPeak"
  "CRAMP1_SANT_R2_peaks_peaks.narrowPeak"
  "CRAMP1_D1_R2_peaks_peaks.narrowPeak"
  "CRAMP1_D1_R3_peaks_peaks.narrowPeak"
)

# Output file
output="peak_score_summary_${TSS_NAME}.tsv"
echo -e "Sample\tTSS_Set\tNum_Overlaps\tMean_Score\tMedian_Score\tTotal_Score" > "$output"

echo "Processing TSS set: $TSS_NAME"

for peakfile in "${samples[@]}"; do
  sample=$(basename "$peakfile" .narrowPeak)
  echo "  Sample: $sample"

  # Intersect peaks with current TSS set
  bedtools intersect -wa -a "$peakfile" -b "$TSS_SET" > "${sample}_${TSS_NAME}_overlaps.bed"

  scores=$(cut -f5 "${sample}_${TSS_NAME}_overlaps.bed")
  num_overlaps=$(echo "$scores" | wc -l)

  if [ "$num_overlaps" -eq 0 ]; then
    mean_score="NA"
    median_score="NA"
    total_score="0"
  else
    mean_score=$(echo "$scores" | awk '{sum+=$1; count++} END {if(count>0) print sum/count; else print "NA"}')
    median_score=$(echo "$scores" | sort -n | awk '{
      a[NR]=$1
    } END {
      if (NR%2==1) {print a[(NR+1)/2]} else {print (a[NR/2]+a[NR/2+1])/2}
    }')
    total_score=$(echo "$scores" | awk '{sum+=$1} END {print sum}')
  fi

  echo -e "${sample}\t${TSS_NAME}\t${num_overlaps}\t${mean_score}\t${median_score}\t${total_score}" >> "$output"
done

echo "🎉 Done! Summary saved to $output"

