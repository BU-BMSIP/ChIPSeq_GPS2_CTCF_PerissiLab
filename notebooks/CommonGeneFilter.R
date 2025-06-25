# Read annotated peaks
anno_file <- "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/results/annotation/ATF4_d6_basal_mm39_annotated.txt"
anno_df <- read.delim(anno_file, stringsAsFactors = FALSE)

# Read common gene list
gene_file <- "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/results/annotation/GPS2_ATF4_CTCF_common_genes.txt"
common_genes <- read.table(gene_file, stringsAsFactors = FALSE)[, 1]

# Filter
filtered_df <- subset(anno_df, SYMBOL %in% common_genes)

# Save filtered result
output_file <- "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/results/annotation/ATF4_filtered_by_ATF4_CTCF_GPS2_common_genes.txt"
write.table(filtered_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
#The BED format is used for visualization (only chr/start/end is reserved)
bed_df <- filtered_df[, c("seqnames", "start", "end")]
write.table(bed_df, file = sub(".txt", ".bed", output_file), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
