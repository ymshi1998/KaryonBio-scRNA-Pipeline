library(Seurat)
library(SingleR)
library(celldex)
library(glmGamPoi)

# Preprocessing Function
process_sample <- function(path, sample_id, condition, donor){
  # Read 10x Original Data
  data <- Read10X(path)
  
  # Create Seurat Object
  sce <- CreateSeuratObject(counts = data, project = sample_id, min.cells = 3, min.features = 200)
  
  # QC
  # Calculate mitochondria percentage
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  
  # Filter out low quality cells
  sce <- subset(sce, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
  
  sce$sample <- sample_id
  sce$condition <- condition
  sce$donor <- donor
  
  # Normalization
  sce <- SCTransform(sce, vars.to.regress = "percent.mt", verbose = FALSE)
  return(sce)
}

# Load all cd45- samples
samples_cd45n <- list(
  healthy1_cd45n_A = process_sample("human_healthy1_cd45-A/", "healthy1-cd45n-A", "healthy", "healthy1"),
  healthy1_cd45n_B = process_sample("human_healthy1_cd45-B/", "healthy1-cd45n-B", "healthy", "healthy1"),
  healthy2_cd45n   = process_sample("human_healthy2_cd45-/", "healthy2-cd45n", "healthy", "healthy2"),
  healthy3_cd45n_A = process_sample("human_healthy3_cd45-A/", "healthy3-cd45n-A", "healthy", "healthy3"),
  healthy3_cd45n_B = process_sample("human_healthy3_cd45-B/", "healthy3-cd45n-B", "healthy", "healthy3"),
  healthy4_cd45n   = process_sample("human_healthy4_cd45-/", "healthy4-cd45n", "healthy", "healthy4"),
  cirrhotic1_cd45n_A    = process_sample("human_cirrhotic1_cd45-A/", "cirrhotic1-cd45n-A", "cirrhotic", "cirrhotic1"),
  cirrhotic1_cd45n_B    = process_sample("human_cirrhotic1_cd45-B/", "cirrhotic1-cd45n-B", "cirrhotic", "cirrhotic1"),
  cirrhotic2_cd45n      = process_sample("human_cirrhotic2_cd45-/", "cirrhotic2-cd45n", "cirrhotic", "cirrhotic2"),
  cirrhotic3_cd45n      = process_sample("human_cirrhotic3_cd45-/", "cirrhotic3-cd45n", "cirrhotic", "cirrhotic3")
)

# Multi-Samples Integration
features <- SelectIntegrationFeatures(samples_cd45n)
samples_cd45n <- PrepSCTIntegration(samples_cd45n, anchor.features = features)

anchors <- FindIntegrationAnchors(
  object.list = samples_cd45n,
  normalization.method = "SCT",
  anchor.features = features
)

sce_cd45n <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT"
)

VlnPlot(
  sce_cd45n,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample",
  pt.size = 0.1,
  ncol = 3
)


# Switch to integrated assay
DefaultAssay(sce_cd45n) <- "integrated"

# PCA + UMAP
sce_cd45n <- RunPCA(sce_cd45n, verbose = FALSE)
ElbowPlot(sce_cd45n)
sce_cd45n <- RunUMAP(sce_cd45n, dims = 1:30)

# Find Clusters
sce_cd45n <- FindNeighbors(sce_cd45n, dims = 1:30)
sce_cd45n <- FindClusters(sce_cd45n, resolution = 0.5)

# Find Marker genes
markers <- FindAllMarkers(sce_cd45n, only.pos = TRUE)
head(markers)

backup_sce <- sce_cd45n

# Cell Type Annotation
sce_cd45n <- JoinLayers(sce_cd45n, assay = "RNA")

DefaultAssay(sce_cd45n) <- "RNA"
sce_cd45n <- NormalizeData(sce_cd45n)
expr <- GetAssayData(sce_cd45n, assay = "RNA", layer = "data")

# Immune for CD45-
ref <- HumanPrimaryCellAtlasData()
pred <- SingleR(
  test = expr,
  ref = ref,
  labels = ref$label.main
)

sce_cd45n$celltype <- pred$labels

# Visualization
DimPlot(sce_cd45n, group.by = "condition")
DimPlot(sce_cd45n, group.by = "celltype", split.by = "condition")

# Save to file
saveRDS(sce_cd45n, file = "cd45n_integrated_annotated.rds")

