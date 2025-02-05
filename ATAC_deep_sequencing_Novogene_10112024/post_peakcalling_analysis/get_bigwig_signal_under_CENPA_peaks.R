library(MotifDb)
library(motifmatchr)
library(Biostrings)
library(clusterProfiler)
library(pheatmap)
library(readr)
library(ggsci)
library(ggrepel)
library(chromPlot)
library(rtracklayer)
library(BiocParallel)

# Set up parallel backend with 20 cores
bp <- MulticoreParam(workers = 20)

# Read in CENPA ChIP-seq peaks
CENPA <- lapply(dir("/lustre/fs4/risc_lab/store/jyeung/for_Hide/Sep2024/Xenbase/ChIPseq/peaks", full.names = TRUE), read.table)
CENPA_peaks <- list()
for (i in seq_along(CENPA)) {
  CENPA_peaks[[i]] <- GRanges(seqnames = CENPA[[i]]$V1, ranges = IRanges(start = CENPA[[i]]$V2, end = CENPA[[i]]$V3), ID = CENPA[[i]]$V4)
}
# Merge replicates into a single peak set and remove scaffold chromosomes
CENPA_peaks <- reduce(unlist(GRangesList(CENPA_peaks)))
CENPA_peaks <- CENPA_peaks[grepl("Chr", seqnames(CENPA_peaks))]

# Split each peak into its own GRanges object
CENPA_peaks_list <- split(CENPA_peaks, seq_along(CENPA_peaks))

# Read in signal from CENPA ChIP-seq bigwigs
CENPA_bw <- lapply(dir("/lustre/fs4/risc_lab/store/jyeung/for_Hide/Sep2024/Xenbase/ChIPseq/bigwigs", full.names = TRUE, pattern = "CENPA"), import.bw)

# Define the bigwig GRanges to use
gr <- CENPA_bw[[1]]

# Parallelized function to find overlaps for each GRanges in CENPA_peaks_list
overlap_results <- bplapply(CENPA_peaks_list, function(gr_list) {
  hits <- findOverlaps(gr_list, gr)
  if (length(hits) > 0) {
    gr[subjectHits(hits)]  # Extract overlapping rows from BigWig GRanges
  } else {
    GRanges()  # Return an empty GRanges if no overlaps are found
  }
}, BPPARAM = bp)

save(overlap_results, file = "/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/workspaces/CENPA_peaks_rep1_signal_bw.RData")

# Parallelized summation of scores
summed_score <- bplapply(overlap_results, function(overlap_results) {
  sum(overlap_results$score)
}, BPPARAM = bp)

save(summed_score, file = "/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/workspaces/CENPA_peaks_rep1_summed_score.RData")

# Add scores to CENPA peaks and save
CENPA_peaks_rep1 <- CENPA_peaks
CENPA_peaks_rep1$score <- unlist(summed_score)
save(CENPA_peaks_rep1, file = "/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/workspaces/CENPA_peaks_rep1_GR_score.RData")
