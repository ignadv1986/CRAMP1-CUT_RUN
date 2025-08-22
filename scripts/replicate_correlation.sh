#!/bin/bash

run_correlation() {
  local group_name=$1
  shift
  local bw_files=("$@")

  echo "📊 Running replicate correlation for $group_name..."

  local npz_file="${group_name}_replicates.npz"
  local corr_plot="${group_name}_replicates_correlation.png"

  multiBigwigSummary bins -b "${bw_files[@]}" -o "$npz_file"

  plotCorrelation -in "$npz_file" -c spearman -p heatmap -o "$corr_plot" \
                  --plotTitle "${group_name} replicate correlation" \
                  --whatToPlot heatmap --skipZeros

  echo "✅ Done for $group_name. Output: $corr_plot"
}

# Run replicate correlation for each group with your actual files
run_correlation "WT" CRAMP1_WT_R1.clean.bw CRAMP1_WT_R2.clean.bw CRAMP1_WT_R3.clean.bw
run_correlation "SANT" CRAMP1_SANT_R1.clean.bw CRAMP1_SANT_R2.clean.bw CRAMP1_SANT_R3.clean.bw
run_correlation "D1" CRAMP1_D1_R1.clean.bw CRAMP1_D1_R2.clean.bw CRAMP1_D1_R3.clean.bw

