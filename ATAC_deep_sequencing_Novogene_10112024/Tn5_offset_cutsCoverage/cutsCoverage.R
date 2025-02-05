library(readr)
library(GenomicAlignments)
library(Rsamtools)
library(MotifDb)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(soGGi)
library(rtracklayer)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)


# name bams
load("/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/workspaces/workspace1.RData")

scalefactor <- 1/dds_tiles$sizeFactor
process_bams <- function(bams, fileName) {
  # Initialize lists
  read1 <- list()
  read2 <- list()
  Firsts <- list()
  First_Pos_toCut <- list()
  First_Neg_toCut <- list()
  Seconds <- list()
  Second_Pos_toCut <- list()
  Second_Neg_toCut <- list()
  test_toCut <- list()
  
  # Extract read1 and read2 from bams
  for (i in 1:length(bams)) {
    read1[[i]] <- first(bams[[i]])
    read2[[i]] <- second(bams[[i]])
  }
  
  # Process reads and calculate cut sites
  for (i in 1:length(read1)) {
    # Read 1
    Firsts[[i]] <- resize(granges(read1[[i]]), fix = "start", 1)
    # Account for Tn5 shift
    First_Pos_toCut[[i]] <- shift(granges(Firsts[[i]][strand(read1[[i]]) == "+"]), 4)
    First_Neg_toCut[[i]] <- shift(granges(Firsts[[i]][strand(read1[[i]]) == "-"]), -5)
    
    # Read 2
    Seconds[[i]] <- resize(granges(read2[[i]]), fix = "start", 1)
    # Account for Tn5 shift
    Second_Pos_toCut[[i]] <- shift(granges(Seconds[[i]][strand(read2[[i]]) == "+"]), 4)
    Second_Neg_toCut[[i]] <- shift(granges(Seconds[[i]][strand(read2[[i]]) == "-"]), -5)
    
    # Combine all cut sites
    test_toCut[[i]] <- c(First_Pos_toCut[[i]], First_Neg_toCut[[i]], Second_Pos_toCut[[i]], Second_Neg_toCut[[i]])
  }
  
  # Convert cut sites to RleList
  cutsCoverage <- lapply(test_toCut, coverage)
  names(cutsCoverage) <- names(bams)
  
  # Save the cutsCoverage data
  save(cutsCoverage, file=paste0("/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Tn5_offset_cutsCoverage/", fileName, "_cutsCoverage.RData"))
  
}

export_normalized_bigwig <- function(cutsCoverage, scalefactor,output_dir) {
  # Normalize cutsCoverage using scalefactor
  normcutsCoverage <- list()
  for (i in 1:length(cutsCoverage)) {
    normcutsCoverage[[i]] <- cutsCoverage[[i]]*scalefactor[i]
  }
  
  # Export BigWig files
  for (i in 1:length(cutsCoverage)) {
    export.bw(cutsCoverage[[i]], con = paste0(output_dir, "/", names(cutsCoverage)[i], "_cutsCoverage.bw"))
    export.bw(normcutsCoverage[[i]], con = paste0(output_dir, "/", names(cutsCoverage)[i], "_scaleFactornorm_cutsCoverage.bw"))
  }
}

Interphase_bams <- lapply(bamdir[1:3], readGAlignmentPairs)
save(Interphase_bams, file="/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Tn5_offset_cutsCoverage/mergedbams_Interphase_reps_GAlignments.RData")
process_bams(Interphase_bams, fileName="Interphase_reps")

Mitotic_bams <- lapply(bamdir[4:6], readGAlignmentPairs)
save(Mitotic_bams, file="/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Tn5_offset_cutsCoverage/mergedbams_Mitotic_reps_GAlignments.RData")
process_bams(Mitotic_bams, fileName="Mitotic_reps")

#need to load back in cutsCoverage workspace first
load("/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Tn5_offset_cutsCoverage/Interphase_reps_cutsCoverage.RData")
names(cutsCoverage) <- samplenames[1:3]
export_normalized_bigwig(cutsCoverage, output_dir ="/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Tn5_offset_cutsCoverage", scalefactor = scalefactor[1:3])

# need to load back in cutsCoverage workspace first
load("/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Tn5_offset_cutsCoverage/Mitotic_reps_cutsCoverage.RData")
names(cutsCoverage) <- samplenames[4:6]
export_normalized_bigwig(cutsCoverage, output_dir ="/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/ATAC_deep_sequencing_Novogene_10112024/Tn5_offset_cutsCoverage", scalefactor = scalefactor[4:6])

