
######## 03_Seurat_QC_Integration.R #######

library(Seurat)
library(dplyr)
library(future)

# Enable parallel processing
plan(multisession, workers = 10)
options(future.globals.maxSize = 20 * 1024 ^ 3) # Set max RAM to 20GB


# 1. Spatial Transcriptomics (ST-seq) Pipeline
# 1.1 Load Raw RDS files (Assuming objects were created and stored in a list)
# Replace these with your actual file paths
st_files <- c("S01_raw.rds", "S02_raw.rds", "S11_raw.rds", "S12_raw.rds")
st_list <- lapply(st_files, readRDS)
names(st_list) <- c("S01", "S02", "S11", "S12")

# 1.2 Quality Control (Based on Methods)
st_list <- lapply(X = st_list, FUN = function(x) {
 # Calculate mitochondrial percentage (adjust pattern to your species, e.g., "^mt-")
 x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^SgalMit") 
 
 # Filter bins: nFeature 100-3200, nCount 100-7500, mt < 10%
 x <- subset(x, subset = nFeature_Spatial >= 100 & nFeature_Spatial <= 3200 & 
              nCount_Spatial >= 100 & nCount_Spatial <= 7500 & 
              percent.mt < 10)
 
 # SCTransform normalization
 x <- SCTransform(x, assay = "Spatial", verbose = FALSE, conserve.memory = TRUE)
 return(x)
})

# 1.3 RPCA Integration for ST-seq
# Select 2,000 integration features as stated in methods
st_features <- SelectIntegrationFeatures(object.list = st_list, nfeatures = 2000)
st_list <- PrepSCTIntegration(object.list = st_list, anchor.features = st_features)

# Note: RPCA requires PCA to be run on each individual object first
st_list <- lapply(X = st_list, FUN = function(x) {
 x <- RunPCA(x, features = st_features, verbose = FALSE)
})

# Find integration anchors using RPCA
st_anchors <- FindIntegrationAnchors(object.list = st_list, 
                                     normalization.method = "SCT",
                                     anchor.features = st_features, 
                                     reduction = "rpca", 
                                     k.anchor = 7, 
                                     k.score = 50)

# Integrate data
st_integrated <- IntegrateData(anchorset = st_anchors, 
                               normalization.method = "SCT", 
                               new.assay.name = "integrated")

# 1.4 Dimensionality Reduction & Clustering
DefaultAssay(st_integrated) <- "integrated"
st_integrated <- RunPCA(st_integrated, npcs = 30, verbose = FALSE)

# Parameters from Methods: top 15 PCs, resolution 0.4
st_integrated <- RunUMAP(st_integrated, reduction = "pca", dims = 1:15)
st_integrated <- FindNeighbors(st_integrated, reduction = "pca", dims = 1:15)
st_integrated <- FindClusters(st_integrated, resolution = 0.4, random.seed = 42)

# 1.5 Sub-clustering of CI cells (Assuming CI cells are cluster '3', adjust accordingly)
# Method: FindSubCluster with resolution 0.2
ci_cluster_id <- "3" 
st_integrated <- FindSubCluster(st_integrated, 
                                cluster = ci_cluster_id, 
                                graph.name = "integrated_snn", 
                                subcluster.name = "CI_subclusters",
                                resolution = 0.2, 
                                algorithm = 1)

saveRDS(st_integrated, "ST_integrated_final.rds")
message("ST-seq integration and clustering complete.")


# 2. Single-Nuclei RNA-seq (SN-seq) Pipeline
# 2.1 Load Raw SN RDS files
sn_files <- c("SN_Rep1_raw.rds", "SN_Rep2_raw.rds")
sn_list <- lapply(sn_files, readRDS)

# 2.2 Quality Control for SN-seq
sn_list <- lapply(X = sn_list, FUN = function(x) {
 x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^SgalMit") 
 
 # Filter nuclei: nFeature 200-4000, nCount 200-6000, mt < 10%
 x <- subset(x, subset = nFeature_RNA >= 200 & nFeature_RNA <= 4000 & 
              nCount_RNA >= 200 & nCount_RNA <= 6000 & 
              percent.mt < 10)
 
 x <- SCTransform(x, assay = "RNA", verbose = FALSE, conserve.memory = TRUE)
 return(x)
})

# 2.3 RPCA Integration for SN-seq
sn_features <- SelectIntegrationFeatures(object.list = sn_list, nfeatures = 2000)
sn_list <- PrepSCTIntegration(object.list = sn_list, anchor.features = sn_features)

sn_list <- lapply(X = sn_list, FUN = function(x) {
 x <- RunPCA(x, features = sn_features, verbose = FALSE)
})

sn_anchors <- FindIntegrationAnchors(object.list = sn_list, 
                                     normalization.method = "SCT",
                                     anchor.features = sn_features, 
                                     reduction = "rpca", 
                                     k.anchor = 7, 
                                     k.score = 50)

sn_integrated <- IntegrateData(anchorset = sn_anchors, 
                               normalization.method = "SCT", 
                               new.assay.name = "integrated")

# 2.4 Dimensionality Reduction & Clustering
DefaultAssay(sn_integrated) <- "integrated"
sn_integrated <- RunPCA(sn_integrated, npcs = 50, verbose = FALSE)

# Parameters from Methods: top 30 PCs, resolution 0.5
sn_integrated <- RunUMAP(sn_integrated, reduction = "pca", dims = 1:30)
sn_integrated <- FindNeighbors(sn_integrated, reduction = "pca", dims = 1:30)
sn_integrated <- FindClusters(sn_integrated, resolution = 0.5, random.seed = 42)

saveRDS(sn_integrated, "SN_integrated_final.rds")
message("SN-seq integration and clustering complete.")