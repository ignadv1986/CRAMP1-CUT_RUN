# CRAMP1-CUT_RUN
CUT&RUN analysis of CRAMP1 protein in human cells

## Project summary

In this project, we analyzed the DNA-binding ability of wild-type CRAMP1 and two distinct mutants (∆SANT and ∆D1) using CUT&RUN. The experiment was performed in CRAMP1 knockout human osteosarcoma U2OS cells reconstituted with inducible expression of GFP-tagged CRAMP1 WT, ∆SANT, or ∆D1. Cells expressing an empty GFP vector served as a control.

The experimental procedures are detailed in [this publication](https://www.sciencedirect.com/science/article/pii/S1097276525003090?via%3Dihub)

Here, we present the bioinformatics workflow used to process and analyze the resulting sequencing data, highlighting key quality control steps, peak calling, and downstream interpretation.

**Note:** Raw sequencing data processing (adapter trimming with fastp, sequence mapping with bowtie2 and removal of duplicated sequences with Picard MarkDuplicates) was performed by the sequencing facility using standard pipelines. This portfolio focuses on the subsequent steps starting with BAM files.

---
## Background

Through CRISPR/Cas9 screens, we identified the previously uncharacterized protein CRAMP1 as a regulator of sensitivity to Topoisomerase 2 inhibitors, a class of anticancer drugs. Subsequent experiments revealed that CRAMP1 localizes to nuclear condensates known as histone locus bodies (HLBs), sites of histone gene transcription. Surprisingly, CRAMP1 knockout did not affect core histone levels but caused approximately a 50% reduction in linker histone H1 levels. The GBD1 domain of CRAMP1 (here termed D1 for simplicity) regulates both localization to HLBs and H1 expression, while the SANT domain is dispensable for localization but essential for maintaining correct H1 levels. To investigate whether CRAMP1 binds histone gene promoters and to determine the roles of these domains, we performed CUT&RUN on GFP-CRAMP1 expressing cell lines.

---
## Goals

- Perform quality control and peak analysis of CUT&RUN data.
- Visualize and interpret peaks.
- Determine CRAMP1 ability to bind the promoters of histone H1 genes.
- Assess differences between WT and mutant proteins. 
---
## Tools

- **bash, conda** - workflow/environment control.
- **bowtie2** - sequence alignment.
- **samtools** - sample quality check.
- **bamPEFragmentSize** - fragment size distribution.
- **bamCoverage** - conversion to BigWig.
- **MACS2** - peak calling.
---
## Workflow
1. **Map reads to spike-in DNA** Mapping of .bam files to Saccharomyces cerevisae DNA using bowtie2. Calculate scaling factors (max yeast reads)/(sample yeast reads)
2. **Sample quality evaluation** Evaluate the quality of the samples after mapping to human genome using samtools in a conda environment.
3. **Fragment size determination** Using bamtools to confirm the right sized of the sequenced fragments.
4. 
5. 
