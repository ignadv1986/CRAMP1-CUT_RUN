library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
peak_file <- "WT_unique_peaks.bed"
peak_gr <- readPeakFile(peak_file)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peak_annot <- annotatePeak(peak_gr, TxDb=txdb, tssRegion=c(-3000, 3000), annoDb="org.Hs.eg.db")
plotAnnoPie(peak_annot)
