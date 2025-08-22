#!/bin/bash

set -euo pipefail

# Function to generate consensus peaks for a given sample type
generate_consensus() {
    sample=$1
    echo "Processing ${sample}..."

    # Get all replicates for this sample
    files=($(ls CRAMP1_${sample}_R*_peaks_peaks.narrowPeak))

    if [ ${#files[@]} -ne 3 ]; then
        echo "Expected 3 replicates for ${sample}, found ${#files[@]}."
        exit 1
    fi

    # Extract pairwise intersections
    bedtools intersect -a "${files[0]}" -b "${files[1]}" > "${sample}_1_2.bed"
    bedtools intersect -a "${files[0]}" -b "${files[2]}" > "${sample}_1_3.bed"
    bedtools intersect -a "${files[1]}" -b "${files[2]}" > "${sample}_2_3.bed"

    # Combine and merge intersections to get consensus
    cat "${sample}_1_2.bed" "${sample}_1_3.bed" "${sample}_2_3.bed" \
        | cut -f1-3 \
        | sort -k1,1 -k2,2n \
        | bedtools merge > "CRAMP1_${sample}_consensus.bed"

    # Clean up temporary files
    rm "${sample}_1_2.bed" "${sample}_1_3.bed" "${sample}_2_3.bed"

    echo "Finished CRAMP1_${sample}_consensus.bed"
}

# Run for each condition
generate_consensus WT
generate_consensus SANT
generate_consensus D1

# Merge all consensus peak sets into a final set
cat CRAMP1_WT_consensus.bed CRAMP1_SANT_consensus.bed CRAMP1_D1_consensus.bed \
    | sort -k1,1 -k2,2n \
    | bedtools merge > CRAMP1_all_conditions_consensus.bed

echo "✅ All done! Final consensus: CRAMP1_all_conditions_consensus.bed"

