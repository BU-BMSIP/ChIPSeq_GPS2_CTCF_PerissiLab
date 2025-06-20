# Load required package
library(dplyr)

# Define file paths
annot_file <- "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/GPS2_CTCFpeaks_merged_annotated.txt"
common_genes_file <- "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/annotation/GPS2_CTCF_common_genes.txt"
output_bed <- "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/annotation/GPS2_CTCF_common_gene_peaks_centered_±3kb.bed"

# Read in the peak annotation table
anno <- read.table(annot_file, sep = "\t", header = TRUE)

# Read in the list of GPS2-CTCF common target genes
common_genes <- readLines(common_genes_file)

# Filter peaks that are associated with the common genes
filtered <- subset(anno, SYMBOL %in% common_genes)

# Calculate center of each peak and define ±3kb regions
center <- round((filtered$start + filtered$end) / 2)
start_new <- pmax(center - 3000, 0)
end_new <- center + 3000

# Construct a data frame
bed_df <- data.frame(
  chr = filtered$seqnames,
  start = start_new,
  end = end_new,
  name = filtered$SYMBOL,
  score = ".",
  strand = filtered$strand
)

# Write to BED file (tab-separated, no column names or row numbers)
write.table(bed_df, file = output_bed,
            sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

cat("BED file centered on GPS2+CTCF peak regions written to:", output_bed, "\n")
