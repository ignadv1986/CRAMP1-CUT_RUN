#!/bin/bash

# Define the three conditions
conditions=("WT" "D1" "SANT")

# Loop through each condition
for condition in "${conditions[@]}"; do
    echo "Processing condition: $condition"

    # Find all relevant narrowPeak files
    replicate_files=$(ls CRAMP1_${condition}_R*_peaks_peaks.narrowPeak)

    # Combine all replicates into one temporary file
    cat $replicate_files > CRAMP1_${condition}_all.tmp

    # Sort and merge
    sort -k1,1 -k2,2n CRAMP1_${condition}_all.tmp > CRAMP1_${condition}_all_sorted.bed
    bedtools merge -i CRAMP1_${condition}_all_sorted.bed > CRAMP1_${condition}_merged.bed

    # Clean up
    rm CRAMP1_${condition}_all.tmp

    echo "Merged peak file for $condition saved as: CRAMP1_${condition}_merged.bed"
    echo "----"
done

