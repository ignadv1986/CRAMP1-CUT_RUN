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
- **samTools** - sample quality check.
- **bamPEFragmentSize** - fragment size distribution.
- **bedtools** - .bedgraph files generation.
- **deepTools** (bigwigCompare, multiBigwigSummary, plotCorrelation, computeMatrix, plotProfile, multiBigWigSummary)
- **USCS-bigWIgtobedGraph** - conversion of .bedgraph files into .bw.
- **MACS2** - peak calling.
- **SeqMonk** - track visualization.
---
## Workflow
1. **Sample Quality Assessment**
- Aligned reads were quality-checked using samtools.
- Fragment size distributions were confirmed with bamtools to ensure expected CUT&RUN fragment profiles.
2. **Spike-In Normalization**
- Sequencing reads were mapped to *Saccharomyces cerevisiae* spike-in DNA using bowtie2.
- A scaling factor was calculated based on yeast read counts to normalize for technical variation across samples.
3. **Filtering Reads**
- Only reads mapping to canonical chromosomes (1–22, X, and Y) were retained for further analysis.
4. **Visualization and Replicate Assessment**
- Initial visualization of mapped reads was done in SeqMonk.
5. **Generation and Normalization of Coverage Tracks**
- BedGraph files were generated from BAMs with spike-in scaling applied.
- BedGraph files were converted to BigWig format.
- BigWig files were normalized to negative control samples using deepTools’ bigwigCompare to highlight specific signal over background.
6. **Replicate Assessment**
- Replicate concordance was evaluated using deepTools (multiBigwigSummary and plotCorrelation) on normalized BigWig files.
7. **Peak Calling and Intersection Analysis**
- Peaks were called on both merged replicates and individual samples using MACS2.
- Unique and shared peaks between conditions were identified using bedtools intersect.
8. **Binding Profile Analysis at Regions of Interest**
- Binding at transcription start sites (TSS) was quantified using deepTools’ computeMatrix with reference to TSS.
- Profiles were plotted with plotProfile to compare binding patterns across different gene subsets.
