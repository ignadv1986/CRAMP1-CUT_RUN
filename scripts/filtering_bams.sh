#!/bin/bash

set -e  # Exit on any error

# Conditions and replicates
conds=("WT" "SANT")
reps=("R1" "R2" "R3")

# Loop through all samples
for cond in "${conds[@]}"; do
  for rep in "${reps[@]}"; do
    input="CRAMP1_${cond}_${rep}.target.dedup.sorted.bam"
    filtered_tmp="CRAMP1_${cond}_${rep}.filtered.tmp.bam"
    header="CRAMP1_${cond}_${rep}.header.sam"
    sq_lines="CRAMP1_${cond}_${rep}.sq.sam"
    meta_lines="CRAMP1_${cond}_${rep}.meta.sam"
    new_header="CRAMP1_${cond}_${rep}.new_header.sam"
    final_bam="CRAMP1_${cond}_${rep}.filtered.cleaned.bam"

    echo "Processing $input..."

    # Step 1: Filter reads from autosomes and sex chromosomes
    samtools view -bh -o "$filtered_tmp" "$input" chr{1..22} chrX chrY

    # Step 2: Extract full header
    samtools view -H "$input" > "$header"

    # Step 3: Keep only @SQ lines for desired chromosomes
    grep -E "^@SQ\s+SN:(chr([1-9]|1[0-9]|2[0-2]|X|Y))\b" "$header" > "$sq_lines"
    grep -v "^@SQ" "$header" > "$meta_lines"
    cat "$meta_lines" "$sq_lines" > "$new_header"

    # Step 4: Reheader the BAM
    samtools reheader "$new_header" "$filtered_tmp" > "$final_bam"

    # Step 5: Index
    samtools index "$final_bam"

    # Optional: Clean up temporary files
    rm "$filtered_tmp" "$header" "$sq_lines" "$meta_lines" "$new_header"

    echo "✔️ Finished $final_bam"
  done
done

echo "✅ All BAMs processed!"

