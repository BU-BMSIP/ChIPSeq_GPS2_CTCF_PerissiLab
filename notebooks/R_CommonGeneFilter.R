# ──────────────────────────────────────────────────────
# Script: CommonGeneFilter.R
# Purpose: 
#   1. Read an annotated peak table.
#   2. Read a list of “common” gene symbols (GPS2 d6 & CTCF).
#   3. Filter the annotation to only those common genes.
#   4. Save both a filtered annotation table and a small BED file 
#      (chr, start, end) for visualization.
# ──────────────────────────────────────────────────────

# 1) Read the full annotation table from a tab-delimited file
anno_file <- "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Adipocyte_differentiation/GPS2/results/annotation/gps2_day0_annotated.txt"
anno_df   <- read.delim(anno_file, stringsAsFactors = FALSE)

# 2) Read the list of common gene symbols into a character vector
gene_file    <- "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Adipocyte_differentiation/GPS2/results/annotation/GPS2d6_CTCF_common_genes.txt"
common_genes <- read.table(gene_file, stringsAsFactors = FALSE)[, 1]

# 3) Filter the annotation data frame to rows whose SYMBOL is in common_genes
filtered_df <- subset(anno_df, SYMBOL %in% common_genes)

# 4) Save the filtered annotation as a tab-delimited file (with headers)
output_file <- "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Adipocyte_differentiation/GPS2/results/annotation/gps2_day0_filtered_by_CTCF_GPS2_d6_common_genes.txt"
write.table(
  filtered_df,
  file      = output_file,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

# 5) Extract only the BED‐relevant columns (chr, start, end) and save as .bed
bed_df <- filtered_df[, c("seqnames", "start", "end")]
write.table(
  bed_df,
  file      = sub(".txt", ".bed", output_file),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
