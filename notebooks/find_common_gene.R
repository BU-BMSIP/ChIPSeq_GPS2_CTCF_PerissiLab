# load library
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(org.Mm.eg.db)

txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene

# ========== GPS2 (sictl) promoter peak annotate ==========
gps2_peaks <- readPeakFile("/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Silencing/GPS2/results/annotation/sictl_promoter_peaks.bed")
gps2_anno <- annotatePeak(gps2_peaks, TxDb=txdb, tssRegion=c(-2000, 2000), annoDb="org.Mm.eg.db")
gps2_df <- as.data.frame(gps2_anno)
gps2_genes <- unique(na.omit(gps2_df$SYMBOL))
write.table(gps2_genes, file="/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/annotation/GPS2_gene_symbols.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


# ========== annotate ctcf all time point ==========
ctcf_files <- c("/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/annotation/CTCF_t1_promoter_only.bed", 
                "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/annotation/CTCF_t2_promoter_only.bed", 
                "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/annotation/CTCF_t3_promoter_only.bed", 
                "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/annotation/CTCF_t4_promoter_only.bed")

ctcf_all_genes <- c()

for (file in ctcf_files) {
  peak <- readPeakFile(file)
  anno <- annotatePeak(peak, TxDb=txdb, tssRegion=c(-2000, 2000), annoDb="org.Mm.eg.db")
  df <- as.data.frame(anno)
  ctcf_all_genes <- c(ctcf_all_genes, na.omit(df$SYMBOL))
}

ctcf_genes <- unique(ctcf_all_genes)
write.table(ctcf_genes, file="/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/annotation/CTCF_all_gene_symbols.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


# ========== find common genes ==========
common_genes <- intersect(gps2_genes, ctcf_genes)
write.table(common_genes, file="/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/annotation/GPS2_CTCF_common_genes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
cat("Number of common genes：", length(common_genes), "\n")


# ========== annotate ATF4 promoter peaks ==========
output_dir <- "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/results/annotation"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.table(atf4_genes, file=file.path(output_dir, "ATF4_gene_symbols.txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)

atf4_peaks <- readPeakFile("/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/results/annotation/ATF4_d6_basal_promoter.bed")
atf4_anno <- annotatePeak(atf4_peaks, TxDb=txdb, tssRegion=c(-2000, 2000), annoDb="org.Mm.eg.db")
atf4_df <- as.data.frame(atf4_anno)
atf4_genes <- unique(na.omit(atf4_df$SYMBOL))
write.table(atf4_genes, file="/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/results/annotation/ATF4_gene_symbols.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

# ========== find common genes between ATF4 and CTCF ==========
common_genes_atf4_ctcf <- intersect(atf4_genes, ctcf_genes)
write.table(common_genes_atf4_ctcf, file="/projectnb/perissilab/Xinyu/GPS2_CHIPseq/ATF4_3T3L1/results/annotation/ATF4_CTCF_common_genes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

cat("Number of common genes between ATF4 and CTCF：", length(common_genes_atf4_ctcf), "\n")

