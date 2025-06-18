annotation_dir1 <- "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/annotation"
setwd(annotation_dir1)

timepoints <- c("t1", "t2", "t3", "t4")

for (tp in timepoints) {
  infile <- paste0("CTCF_", tp, "_bf_annotated.txt")
  outfile <- paste0("CTCF_", tp, "_promoter_only.bed")
  
  anno <- read.delim(infile, header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
  
  promoter_peaks <- subset(anno, grepl("Promoter", anno$annotation))
  
  bed <- promoter_peaks[, c("seqnames", "start", "end")]

  write.table(bed, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

annotated_file <- "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Silencing/GPS2/results/annotation/sictl_annotated.txt"
output_bed <- "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Silencing/GPS2/results/annotation/sictl_promoter_peaks.bed"
anno <- read.delim(annotated_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)


promoter_peaks <- subset(anno, grepl("Promoter", anno$annotation))


bed <- promoter_peaks[, c("seqnames", "start", "end")]

write.table(bed, file = output_bed,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
