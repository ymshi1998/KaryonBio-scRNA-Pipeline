# KaryonBio-scRNA-Pipeline
This is a Mini-Pipeline for Discovery and Prioritization of Cell-Type-Specific Biomarkers in Human Liver Fibrosis.

### Introduction
The primary dataset of this pipeline is GSE136103 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136103), the primary article is 'Resolving the fibrotic niche of human liver cirrhosis at single cell level' (https://pubmed.ncbi.nlm.nih.gov/31597160/). This is an end-to-end analysis workflow that
- identifies fibrosis-associated cell populations
- discovers cell-type-specific fibrosis biomarkers 
- prioritizes a short list of translational biomarker or therapeutic target candidates.


#### Files
**'scRNA-pipeline-human-immune.R'**, focuses on CD45+ human liver cell scRNA analysis.

**'scRNA-pipeline-human-parenchymal.R'**, focuses on CD45- human liver cell scRNA analysis.


#### Instructions to Run the Script
<mark>CRITICAL</mark>: Please download and unzip the dataset file from the provided link. After that, you will see a list of healthy/cirrhotic human liver cell files, healthy/cirrhotic human blood files, and healthy/cirrhotic rodent files. Each sample sequencing generates 3 corresponding file: barcodes.tsv.gz, genes.tsv.gz, and matrix.mtx.gz.

For each cell type and sequencing batch, for example, `GSM4041151_healthy1_cd45-A_barcodes.tsv.gz`, `GSM4041151_healthy1_cd45-A_genes.tsv.gz`, and `GSM4041151_healthy1_cd45-A_matrix.mtx.gz`, create a folder parallel to the running R script named `human_healthy1_cd45-A`. Paste the three files inside, and rename them as `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz` in order to enable Seurat to load and create Seurat object. 

After you do this step, you will have a list of folders, each containing `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz`. You will now be able to run all the R scripts.
