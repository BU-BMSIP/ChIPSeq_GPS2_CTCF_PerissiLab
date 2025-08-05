# ───────────────────────────────────────────────────────────────
# Script: generate_common_peaks.R
#
# Purpose:
#   1. Load an annotated peak table.
#   2. Read a list of “common” gene symbols (GPS2 day6 ∩ CTCF).
#   3. Filter the annotation to promoter peaks for those common genes.
#   4. Output a BED file (seqnames, start, end) for visualization.
# ───────────────────────────────────────────────────────────────


# 1) Load required package
library(dplyr)

# 2) Define file paths (modify here to reuse with different files)
annot_file         <- "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Adipocyte_differentiation/GPS2/results/annotation/gps2_day6_annotated.txt"
common_genes_file  <- "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Adipocyte_differentiation/GPS2/results/annotation/GPS2d6_CTCF_common_genes.txt"
output_bed         <- "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Adipocyte_differentiation/GPS2/results/annotation/gps2_day6_filtered_by_CTCF_GPS2_d6_common_genes_promoter_only.bed"

# 3) Read inputs
anno_df       <- read.delim(annot_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
common_genes  <- readLines(common_genes_file)

# 4) Filter to promoter annotations only
promoter_df   <- anno_df %>% 
  filter(grepl("Promoter", annotation))

# 5) Further filter to rows whose SYMBOL is in the common gene list
filtered_df   <- promoter_df %>% 
  filter(SYMBOL %in% common_genes)

# 6) Select only the BED columns (chr, start, end) and write to file
bed_df        <- filtered_df %>% 
  select(seqnames, start, end)

write.table(
  bed_df,
  file      = output_bed,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

cat("BED file written to:", output_bed, "\n")
