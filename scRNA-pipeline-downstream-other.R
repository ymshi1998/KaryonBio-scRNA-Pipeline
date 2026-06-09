library(Seurat)
library(DESeq2)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)

sce <- readRDS("cd45n_integrated_annotated.rds")

# 3 Categories:
# - Stromal (HSC / fibroblast / myofibroblast-like): Fibroblasts, Smooth_muscle_cells, HSC_CD34+, HSC_-G-CSF
# - Myeloid (Macrophage / Monocyte): Macrophage, Monocyte
# - Endothelial: Endothelial_cells

# Separate all the cell types to those 3 categores and others
sce$celltype_big <- "Other"

sce$celltype_big[sce$celltype %in% c(
  "Fibroblasts", "Smooth_muscle_cells", "HSC_CD34+", "HSC_-G-CSF"
)] <- "Stromal"

sce$celltype_big[sce$celltype %in% c(
  "Macrophage", "Monocyte"
)] <- "Myeloid"

sce$celltype_big[sce$celltype %in% c(
  "Endothelial_cells"
)] <- "Endothelial"

# Function to run pseudobulk + DESeq2
run_pseudobulk_DE <- function(sce, celltype_name) {
  
  message("Running pseudobulk for: ", celltype_name)
  
  sub <- subset(sce, celltype_big == celltype_name)
  
  counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
  donor <- sub$donor
  condition <- sub$condition
  
  # aggregate counts by donor
  pb_counts <- sapply(unique(donor), function(d) {
    cells <- which(donor == d)
    Matrix::rowSums(counts[, cells, drop = FALSE])
  })
  
  # metadata
  pb_meta <- data.frame(
    donor = unique(donor),
    condition = condition[match(unique(donor), donor)]
  )
  rownames(pb_meta) <- pb_meta$donor
  
  # DESeq2
  dds <- DESeqDataSetFromMatrix(
    countData = pb_counts,
    colData = pb_meta,
    design = ~ condition
  )
  
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", "cirrhotic", "healthy"))
  res <- as.data.frame(res)
  res$gene <- rownames(res)
  
  return(res)
}

# Run pseudobulk + DESeq2 to the 3 categories
res_stromal <- run_pseudobulk_DE(sce, "Stromal")
res_myeloid <- run_pseudobulk_DE(sce, "Myeloid")
res_endo    <- run_pseudobulk_DE(sce, "Endothelial")

# Run GO for the combination of 3 categories
geneList <- list(
  Stromal = res_stromal %>% filter(padj < 0.05 & log2FoldChange > 0.5) %>% pull(gene),
  Myeloid = res_myeloid %>% filter(padj < 0.05 & log2FoldChange > 0.5) %>% pull(gene),
  Endothelial = res_endo %>% filter(padj < 0.05 & log2FoldChange > 0.5) %>% pull(gene)
)

# Convert Gene ID to Entrez ID
convert_to_entrez <- function(genes) {
  bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
}

geneList_entrez <- list(
  Stromal = convert_to_entrez(geneList$Stromal),
  Myeloid = convert_to_entrez(geneList$Myeloid),
  Endothelial = convert_to_entrez(geneList$Endothelial)
)

cc <- compareCluster(
  geneCluster = geneList_entrez,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

dotplot(cc, showCategory = 20)

# Generate GO-informed gene target ranking

ego_df <- as.data.frame(cc)

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



