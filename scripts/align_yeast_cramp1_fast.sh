#!/bin/bash

# Directory containing gzipped FASTQ files
RAW_DIR="raw_data"

# Yeast Bowtie2 index (no extensions)
YEAST_INDEX="yeast_index"

# Output file for mapped read counts
echo -e "Sample\tYeast_Reads" > yeast_spikein_counts.tsv

# List of sample names (excluding _R1/_R2)
samples=(
  CUTRUN_CRAMP1_GBN_rep1
  CUTRUN_CRAMP1_GBN_rep2
  CUTRUN_CRAMP1_GBN_rep3
  CUTRUN_CRAMP1_SANT_rep1
  CUTRUN_CRAMP1_SANT_rep2
  CUTRUN_CRAMP1_SANT_rep3
  CUTRUN_CRAMP1_WT_rep1
  CUTRUN_CRAMP1_WT_rep2
  CUTRUN_CRAMP1_WT_rep3
)

for sample in "${samples[@]}"; do
  fq1="${RAW_DIR}/${sample}_R1.fastq.gz"
  fq2="${RAW_DIR}/${sample}_R2.fastq.gz"
  sam="${sample}_yeast.sam"
  bam="${sample}_yeast_sorted.bam"

  echo "🔄 Processing $sample..."

  # Fast alignment to yeast genome
  bowtie2 -x "$YEAST_INDEX" -1 "$fq1" -2 "$fq2" -S "$sam" --fast -p 4

  # Convert SAM to sorted BAM
  samtools view -bS "$sam" | samtools sort -o "$bam"

  # Index the BAM
  samtools index "$bam"

  # Count mapped reads (exclude unmapped)
  count=$(samtools view -c -F 4 "$bam")

  # Record in output table
  echo -e "${sample}\t${count}" >> yeast_spikein_counts.tsv

  # Clean up intermediate SAM
  rm "$sam"

  echo "✅ Done: $sample → $count yeast reads"
done

echo "🎉 All alignments complete. Summary saved to: yeast_spikein_counts.tsv"

