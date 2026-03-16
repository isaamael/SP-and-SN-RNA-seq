
######## 05_Spatial_Modeling.R #######
library(Seurat)
library(scDist)
library(miloR)
library(miloDE)
library(SingleCellExperiment)
library(dplyr)
library(BiocParallel)

st_obj <- readRDS("ST_integrated_final.rds")
bpparam <- SnowParam(workers = 8)

# Calculate transcriptional distances between states using scDist
# Methods: top 10 PCs, batch included as random effects
scD_input <- list(
 cont = as.data.frame(st_obj@assays$integrated@data),
 metadata = as.data.frame(st_obj@meta.data)
)

scdist_res <- scDist(normalized_counts = scD_input$cont,
                     meta.data = scD_input$metadata,
                     d = 10, 
                     fixed.effects = "treat",
                     random.effects = "orig.ident",
                     clusters = "celltype")

write.table(scdist_res$results, "scDist_results.txt", sep="\t", quote=FALSE)

# Differential abundance (DA) testing using miloR
# Methods: top 30 PCs, sampling fraction 0.1, k=20
st_sce <- as.SingleCellExperiment(st_obj, assay = "integrated")
st_milo <- Milo(st_sce)

st_milo <- buildGraph(st_milo, k = 20, d = 30, transposed = TRUE, reduced.dim = "PCA")
st_milo <- makeNhoods(st_milo, prop = 0.1, k = 20, d = 30, refined = TRUE, reduced_dims = "PCA")
st_milo <- countCells(st_milo, meta.data = as.data.frame(colData(st_milo)), samples="orig.ident")

exp_design <- data.frame(
 sample = unique(st_obj$orig.ident),
 treatment = ifelse(grepl("S0", unique(st_obj$orig.ident)), "CK", "LTSS")
)
rownames(exp_design) <- exp_design$sample

st_milo <- calcNhoodDistance(st_milo, d = 30, reduced.dim = "PCA")
da_results <- testNhoods(st_milo, design = ~ treatment, design.df = exp_design)

# Filter out neighborhoods with dominant cell type < 80%
da_results <- annotateNhoods(st_milo, da_results, coldata_col = "celltype")
da_results <- da_results %>% filter(celltype_fraction >= 0.8)

write.table(da_results, "miloR_DA_results.txt", sep="\t", quote=FALSE)

# Differential expression (DE) testing using MiloDE
# Methods: k=10, d=30, p=0.2
st_milo_de <- buildGraph(st_milo, k = 10, d = 30, transposed = TRUE, reduced.dim = "PCA")
st_milo_de <- makeNhoods(st_milo_de, prop = 0.2, k = 10, d = 30, refined = TRUE, reduced_dims = "PCA")

# Calculate AUC to assess treatment perturbation
stat_auc <- calc_AUC_per_neighbourhood(st_milo_de, 
                                       sample_id = "orig.ident", 
                                       condition_id = "treat", 
                                       BPPARAM = bpparam)

# Filter neighborhoods (AUC < 0.5 removed as per methods)
valid_nhoods <- stat_auc$Nhood[stat_auc$auc >= 0.5 & !is.na(stat_auc$auc)]

# Perform DE testing on perturbed neighborhoods
de_stat <- de_test_neighbourhoods(st_milo_de,
                                  sample_id = "orig.ident",
                                  design = ~ treat,
                                  covariates = "treat",
                                  subset_nhoods = valid_nhoods,
                                  output_type = "SCE",
                                  BPPARAM = bpparam)

# Extract and merge final MiloDE results
stat_de_magnitude <- rank_neighbourhoods_by_DE_magnitude(de_stat)
nhood_stat_ct <- data.frame(Nhood = 1:ncol(nhoods(st_milo_de)), Nhood_center = colnames(nhoods(st_milo_de)))
nhood_stat_ct <- annotateNhoods(st_milo_de, nhood_stat_ct, coldata_col = "celltype")

milo_de_summary <- merge(stat_auc, stat_de_magnitude, by=c("Nhood", "Nhood_center"))
milo_de_summary <- merge(milo_de_summary, nhood_stat_ct, by=c("Nhood", "Nhood_center"))

# Apply logFC > 0.25 and FDR < 0.1 thresholds (can be applied downstream or saved directly)
write.table(milo_de_summary, "MiloDE_summary_results.txt", sep="\t", quote=FALSE)


