library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

# setting your output dir
output_dir <- "/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/Enrichment"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# read gene symbol list
genes <- read.table("/projectnb/perissilab/Xinyu/GPS2_CHIPseq/CTCF_3T3L1/results/annotation/GPS2_CTCF_common_genes.txt", stringsAsFactors = FALSE)[,1]

# transfer into Entrez ID
gene_entrez <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
entrez_ids <- gene_entrez$ENTREZID

#GO Biological Process
ego <- enrichGO(gene = entrez_ids,
                OrgDb = org.Mm.eg.db,
                ont = "BP",     # Biological Process
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2,
                readable = TRUE)

# barplot
barplot(ego, showCategory = 20, title = "GO BP Enrichment")

#  KEGG Enrichment Analysis
ekegg <- enrichKEGG(gene = entrez_ids,
                    organism = 'mmu', 
                    pvalueCutoff = 0.05)

# dotplot
dotplot(ekegg, showCategory = 20, title = "KEGG Pathway Enrichment")

write.csv(as.data.frame(ego), file = file.path(output_dir, "GO_BP_enrichment.csv"), row.names = FALSE)
write.csv(as.data.frame(ekegg), file = file.path(output_dir, "KEGG_enrichment.csv"), row.names = FALSE)

png(file.path(output_dir, "KEGG_dotplot.png"), width=10, height=6)
dotplot(ekegg, showCategory = 20, title = "KEGG Pathway Enrichment")
dev.off()
png(file.path(output_dir, "GO_BP_barplot.png"), width=1000, height=600)
barplot(ego, showCategory = 20, title = "GO BP Enrichment")
dev.off()

