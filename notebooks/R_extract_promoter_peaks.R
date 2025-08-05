# ─────────────────────────────────────────────────────────────────────────────
# Script: extract_promoter_peaks.R
#
# Purpose:
#   1. Iterate over CTCF timepoint annotation files (t1–t4) and extract promoter peaks.
#   2. Extract promoter peaks from the GPS2 silencing annotation file.
#   3. Define a reusable function to extract promoter peaks from any annotation file.
#   4. Write out promoter-only BED files (seqnames, start, end) for downstream analysis.
# ─────────────────────────────────────────────────────────────────────────────

# ========== Cell 1: Setup & CTCF Timepoints ==========
# Directory containing annotation .txt files for CTCF timepoints
annotation_dir1 <- "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/annotation"
setwd(annotation_dir1)

# Timepoint identifiers corresponding to CTCF replicates
timepoints <- c("t1", "t2", "t3", "t4")

for (tp in timepoints) {
  # Build input/output filenames
  infile  <- paste0("CTCF_", tp, "_bf_annotated.txt")
  outfile <- paste0("CTCF_", tp, "_promoter_only.bed")
  
  # Read the annotation table
  anno <- read.delim(infile, header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
  
  # Subset rows where the annotation column contains “Promoter”
  promoter_peaks <- subset(anno, grepl("Promoter", annotation))
  
  # Select only the BED columns (seqnames, start, end)
  bed <- promoter_peaks[, c("seqnames", "start", "end")]
  
  # Write out as a BED file without headers or row names
  write.table(bed,
              file      = outfile,
              sep       = "\t",
              quote     = FALSE,
              row.names = FALSE,
              col.names = FALSE)
}

# ========== Cell 2: GPS2 Silencing Promoter Extraction ==========
# Single-file extraction for GPS2 silencing annotation
annotated_file <- "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Silencing/GPS2/results/annotation/sictl_annotated.txt"
output_bed     <- "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Silencing/GPS2/results/annotation/sictl_promoter_peaks.bed"

anno <- read.delim(annotated_file,
                   header = TRUE,
                   sep    = "\t",
                   stringsAsFactors = FALSE)

promoter_peaks <- subset(anno, grepl("Promoter", annotation))
bed             <- promoter_peaks[, c("seqnames", "start", "end")]

write.table(bed,
            file      = output_bed,
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# ========== Cell 3: Reusable Function Definition ==========
# Define a function to extract promoter peaks from any annotation file
extract_promoter_peaks <- function(infile, outfile) {
  anno <- read.delim(infile,
                     header = TRUE,
                     sep    = "\t",
                     stringsAsFactors = FALSE)
  
  promoter_peaks <- subset(anno, grepl("Promoter", annotation))
  bed             <- promoter_peaks[, c("seqnames", "start", "end")]
  
  write.table(bed,
              file      = outfile,
              sep       = "\t",
              quote     = FALSE,
              row.names = FALSE,
              col.names = FALSE)
}

# Example usage:
extract_promoter_peaks(
  infile  = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/results/annotation/ATF4_d6_basal_mm39_annotated.txt",
  outfile = "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/results/annotation/ATF4_d6_basal_promoter.bed"
)
