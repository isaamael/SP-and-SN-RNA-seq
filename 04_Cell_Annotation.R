

######## 04_Cell_Annotation.R #######

library(Seurat)
library(COSG)
library(UCell)
library(MetaNeighbor)
library(SingleCellExperiment)
library(BiocParallel)

st_obj <- readRDS("ST_integrated_final.rds")
sn_obj <- readRDS("SN_integrated_final.rds")

# Identify cluster-specific markers using COSG
# Based on methods: top 100 specifically expressed genes
DefaultAssay(st_obj) <- "integrated"
Idents(st_obj) <- "celltype"

st_markers <- cosg(st_obj, 
                   groups = 'all', 
                   assay = 'integrated',
                   slot = 'data', 
                   mu = 10, 
                   remove_lowly_expressed = TRUE,
                   expressed_pct = 0.1, 
                   n_genes_user = 100)

# Format COSG output to dataframe
marker_genes <- character()
marker_clusters <- character()
for (cluster_name in names(st_markers$names)) {
 genes <- st_markers$names[[cluster_name]]
 marker_genes <- c(marker_genes, genes)
 marker_clusters <- c(marker_clusters, rep(cluster_name, length(genes)))
}
st_marker_df <- data.frame(Gene_ID = marker_genes, Cluster = marker_clusters)
write.table(st_marker_df, "ST_COSG_markers.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Extract marker lists for UCell scoring
marker_list <- lapply(names(st_markers$names), function(col) st_markers$names[[col]])
names(marker_list) <- names(st_markers$names)

# SN-seq cell type matching using UCell
# Parameters set according to methods: maxRank = 1000, w_neg = 0.5
bpparam <- SnowParam(workers = 8)
sn_obj <- AddModuleScore_UCell(sn_obj, 
                               features = marker_list, 
                               maxRank = 1000, 
                               assay = "integrated", 
                               slot = "data", 
                               w_neg = 0.5, 
                               BPPARAM = bpparam)

# SN-seq cell type matching using MetaNeighbor
# Methods: ST-seq data served as the pre-trained reference model
common_genes <- intersect(rownames(st_obj), rownames(sn_obj))

st_sce <- as.SingleCellExperiment(st_obj, assay = "integrated")
sn_sce <- as.SingleCellExperiment(sn_obj, assay = "integrated")

pretrained_model <- trainModel(var_genes = common_genes,
                               dat = st_sce,
                               study_id = st_sce$orig.ident,
                               cell_type = st_sce$celltype)

aurocs <- MetaNeighborUS(trained_model = pretrained_model,
                         var_genes = common_genes,
                         dat = sn_sce,
                         study_id = sn_sce$orig.ident,
                         cell_type = sn_sce$seurat_clusters,
                         one_vs_best = TRUE,
                         fast_version = TRUE)

write.table(aurocs, "MetaNeighbor_AUROC_ST_to_SN.txt", sep="\t", quote=FALSE)