#!/bin/bash

# Versions and replicates
versions=("WT" "D1" "SANT")
replicates=("R1" "R2" "R3")

# Output directory
mkdir -p macs2_output

# Loop over all combinations
for v in "${versions[@]}"; do
  for r in "${replicates[@]}"; do

    # Define file names
    treatment="CRAMP1_${v}_${r}.filtered.cleaned.bam"
    control="GFP_${r}.target.dedup.sorted.bam"
    outname="CRAMP1_${v}_${r}_peaks"

    # Check files exist
    if [[ -f "$treatment" && -f "$control" ]]; then
      echo "Calling peaks for $treatment with control $control"

      # Call MACS2
      macs2 callpeak -t "$treatment" \
                     -c "$control" \
                     -f BAM \
                     -g hs \
                     -n "$outname" \
                     --outdir macs2_output \
                     --nomodel \
                     --shift -100 \
                     --extsize 200 \
                     --keep-dup all \
                     -q 0.05

    else
      echo "❌ Missing file(s) for $v $r: $treatment or $control not found!"
    fi

  done
done

