#!/bin/bash

# Parameters
GENOME_SIZE=2913022398  # hg38 effective genome size; change if needed
BIN_SIZE=10             # bin size for bigWig resolution
NORM=CPM                # or RPGC if using deepTools with scaling factors
THREADS=4               # adjust to your system

# Loop over all replicate BAMs
for bam in CRAMP1_*_R*.filtered.cleaned.bam; do
    base=$(basename "$bam" .filtered.cleaned.bam)
    echo "Generating bigWig for $base"
    bamCoverage \
        -b "$bam" \
        -o "${base}.bw" \
        --binSize $BIN_SIZE \
        --normalizeUsing $NORM \
        --effectiveGenomeSize $GENOME_SIZE \
        --extendReads \
        --numberOfProcessors $THREADS
done

# Also process merged BAMs
for bam in CRAMP1_*_merged.bam; do
    base=$(basename "$bam" .bam)
    echo "Generating bigWig for $base"
    bamCoverage \
        -b "$bam" \
        -o "${base}.bw" \
        --binSize $BIN_SIZE \
        --normalizeUsing $NORM \
        --effectiveGenomeSize $GENOME_SIZE \
        --extendReads \
        --numberOfProcessors $THREADS
done

