library(Seurat)
library(DESeq2)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)

sce <- readRDS("cd45p_integrated_annotated.rds")

# Macrophage / monocyte populations
celltype_of_interest <- "Monocytes"
sub <- subset(sce, celltype == celltype_of_interest)

# Build pseudobulk
counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
donor <- sub$sample
condition <- sub$condition

# aggregate counts by donor
pb_counts <- sapply(unique(donor), function(d) {
  cells <- which(donor == d)
  Matrix::rowSums(counts[, cells, drop = FALSE])
})

# metadata for DESeq2
pb_meta <- data.frame(
  donor = unique(donor),
  condition = condition[match(unique(donor), donor)]
)
rownames(pb_meta) <- pb_meta$donor

# Run DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = pb_counts,
  colData = pb_meta,
  design = ~ condition
)

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "cirrhotic", "healthy"))
res <- as.data.frame(res)
res$gene <- rownames(res)

# Pathway / mechanism analysis
up_genes <- res %>%
  filter(padj < 0.05 & log2FoldChange > 0.5) %>%
  pull(gene)

entrez <- bitr(up_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

ego <- enrichGO(
  gene = entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

head(ego)
dotplot(ego, showCategory = 20)

# Generate GO-informed gene target ranking

ego_df <- as.data.frame(ego)

# Create a gene table
gene_table <- ego_df %>%
  select(Description, p.adjust, geneID) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  rename(gene = geneID)

# Create a gene rank based on pathways each gene appears
gene_rank <- gene_table %>%
  group_by(gene) %>%
  summarise(
    n_pathways = n(),
    best_p = min(p.adjust)
  ) %>%
  arrange(best_p)

head(gene_rank, 20)



