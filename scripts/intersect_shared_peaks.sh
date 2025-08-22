#!/bin/bash

# Make sure you're in the correct directory
INPUT_DIR="/Volumes/dhb572/Sequencing_data/CRAMP1_CUTRUN/alignment"
cd "$INPUT_DIR" || exit 1

# List of your sample groups
samples=("CRAMP1_WT" "CRAMP1_D1" "CRAMP1_SANT")

for sample in "${samples[@]}"; do
  echo "Processing $sample"

  # Define replicate file names
  rep1="${sample}_R1_vs_GFP_R1.stringent.bed"
  rep2="${sample}_R2_vs_GFP_R2.stringent.bed"
  rep3="${sample}_R3_vs_GFP_R3.stringent.bed"

  # Check that all replicate files exist
  if [[ -f "$rep1" && -f "$rep2" && -f "$rep3" ]]; then
    # Intersect all three replicates to find shared peaks
    bedtools intersect -a "$rep1" -b "$rep2" | bedtools intersect -a - -b "$rep3" > "${sample}_shared_peaks.bed"
    echo "✅ Created ${sample}_shared_peaks.bed"
  else
    echo "⚠️ Missing one or more replicates for $sample"
  fi
done

echo "Finished intersecting all samples."

