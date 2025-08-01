# Load required package
library(dplyr)

# Define file paths
annot_file <- "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Adipocyte_differentiation/GPS2/results/annotation/gps2_day6_annotated.txt"
common_genes_file <- "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Adipocyte_differentiation/GPS2/results/annotation/GPS2d6_CTCF_common_genes.txt"
output_bed <- "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/Adipocyte_differentiation/GPS2/results/annotation/gps2_day6_filtered_by_CTCF_GPS2_d6_common_genes_promoter_only.bed"

# 如果你想以后快速换文件，只需改上面三行

# 读取文件
anno <- read.delim(annot_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
common_genes <- readLines(common_genes_file)

# 先过滤 promoter
promoter_only <- anno %>% filter(grepl("Promoter", annotation))

# 再按基因名过滤
filtered <- promoter_only %>% filter(SYMBOL %in% common_genes)

# 输出
bed_df <- filtered %>% select(seqnames, start, end)
write.table(bed_df, output_bed, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

cat("生成完成：", output_bed, "\n")
