library(Seurat)
library(SingleR)
library(celldex)
library(glmGamPoi)

# Preprocessing Function
process_sample <- function(path, sample_id, condition){
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
  
  # Normalization
  sce <- SCTransform(sce, vars.to.regress = "percent.mt", verbose = FALSE)
  return(sce)
}

# Load all cd45+ samples
samples <- list(
  healthy1_cd45p = process_sample("human_healthy1_cd45+/", "healthy1-cd45p", "healthy"),
  healthy2_cd45p = process_sample("human_healthy2_cd45+/", "healthy2-cd45p", "healthy"),
  healthy3_cd45p = process_sample("human_healthy3_cd45+/", "healthy3-cd45p", "healthy"),
  healthy4_cd45p = process_sample("human_healthy4_cd45+/", "healthy4-cd45p", "healthy"),
  healthy5_cd45p = process_sample("human_healthy5_cd45+/", "healthy5-cd45p", "healthy"),
  cirrhotic1_cd45p = process_sample("human_cirrhotic1_cd45+/", "cirrhotic1-cd45p", "cirrhotic"),
  cirrhotic2_cd45p = process_sample("human_cirrhotic2_cd45+/", "cirrhotic2-cd45p", "cirrhotic"),
  cirrhotic3_cd45p = process_sample("human_cirrhotic3_cd45+/", "cirrhotic3-cd45p", "cirrhotic"),
  cirrhotic4_cd45p = process_sample("human_cirrhotic4_cd45+/", "cirrhotic4-cd45p", "cirrhotic"),
  cirrhotic5_cd45p = process_sample("human_cirrhotic5_cd45+/", "cirrhotic5-cd45p", "cirrhotic")
)

# Multi-Samples Integration
features <- SelectIntegrationFeatures(samples)
samples <- PrepSCTIntegration(samples, anchor.features = features)

anchors <- FindIntegrationAnchors(
  object.list = samples,
  normalization.method = "SCT",
  anchor.features = features
)

sce <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT"
)

VlnPlot(
  sce,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample",
  pt.size = 0.1,
  ncol = 3
)


# Switch to integrated assay
DefaultAssay(sce) <- "integrated"

# PCA + UMAP
sce <- RunPCA(sce, verbose = FALSE)
ElbowPlot(sce)
sce <- RunUMAP(sce, dims = 1:30)

# Find Clusters
sce <- FindNeighbors(sce, dims = 1:30)
sce <- FindClusters(sce, resolution = 0.5)

# Find Marker genes
markers <- FindAllMarkers(sce, only.pos = TRUE)
head(markers)

backup_sce <- sce

# Cell Type Annotation
sce <- JoinLayers(sce, assay = "RNA")

DefaultAssay(sce) <- "RNA"
sce <- NormalizeData(sce)
expr <- GetAssayData(sce, assay = "RNA", layer = "data")

# Immune for CD45+
ref <- MonacoImmuneData()
pred <- SingleR(
  test = expr,
  ref = ref,
  labels = ref$label.main
)

sce$celltype <- pred$labels

# Visualization
DimPlot(sce, group.by = "condition")
DimPlot(sce, group.by = "celltype", split.by = "condition")

# Save to file
saveRDS(sce, file = "cd45p_integrated_annotated.rds")

