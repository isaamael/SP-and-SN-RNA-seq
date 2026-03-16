
########### 06_CARD_Deconvolution.R #########
library(Seurat)
library(CARD)
library(dplyr)

st_obj <- readRDS("ST_integrated_final.rds")
sn_obj <- readRDS("SN_integrated_final.rds")

# Prepare spatial data inputs
st_counts <- st_obj@assays$Spatial@counts
st_coords <- st_obj@images$sample1@coordinates[, c("row", "col")]
colnames(st_coords) <- c("x", "y")

# Prepare single-nuclei reference inputs
sn_counts <- sn_obj@assays$RNA@counts
sn_meta <- sn_obj@meta.data[, c("orig.ident", "celltype")]
colnames(sn_meta) <- c("sampleInfo", "cellType")
sn_meta$cellType <- factor(sn_meta$cellType)

# Create CARD object and run deconvolution
# Methods: Genes with > 50 counts and detected in at least 5 spots
card_obj <- createCARDObject(
 sc_count = sn_counts,
 sc_meta = sn_meta,
 spatial_count = st_counts,
 spatial_location = st_coords,
 ct.varname = "cellType",
 ct.select = unique(sn_meta$cellType),
 sample.varname = "sampleInfo",
 minCountGene = 50,
 minCountSpot = 5
) 

card_obj <- CARD_deconvolution(CARD_object = card_obj)

# Add CARD predictions as metadata to the Seurat object
st_obj <- AddMetaData(st_obj, metadata = as.data.frame(card_obj@Proportion_CARD))
write.table(card_obj@Proportion_CARD, "CARD_Predicted_Proportions.txt", sep="\t", quote=FALSE)

# Model goodness-of-fit validation (RMSE calculation)
basis_matrix <- card_obj@algorithm_matrix$B
prop_matrix <- card_obj@Proportion_CARD
common_genes <- intersect(rownames(basis_matrix), rownames(card_obj@spatial_countMat))

basis_matrix <- basis_matrix[common_genes, ]
obs_matrix <- card_obj@spatial_countMat[common_genes, rownames(prop_matrix)]
reconstructed_matrix <- basis_matrix %*% t(prop_matrix)

# Calculate RMSE for each spot to rule out structural artifacts
residuals <- obs_matrix - reconstructed_matrix
rmse_per_spot <- sqrt(colMeans(residuals^2))
rmse_df <- data.frame(Spot = names(rmse_per_spot), RMSE = rmse_per_spot)
write.table(rmse_df, "CARD_RMSE_Validation.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Identify "juxtaposed cells" based on Euclidean distance
# Methods: 10 non-CI cells with the minimum Euclidean distance to a CI
ci_spots <- rownames(st_obj@meta.data[st_obj@meta.data$celltype == "CI", ])
dist_matrix <- as.matrix(dist(st_coords))

# Extract the 10 nearest neighbors for each CI spot
n_neighbors <- 10
juxtaposed_list <- lapply(ci_spots, function(spot) {
 distances <- dist_matrix[spot, ]
 sorted_indices <- order(distances)
 # Exclude self and select top n
 closest_spots <- names(distances[sorted_indices])[2:(n_neighbors + 1)]
 return(data.frame(CI_Spot = spot, Juxtaposed_Spot = closest_spots))
})

juxtaposed_df <- do.call(rbind, juxtaposed_list)
juxtaposed_df <- juxtaposed_df %>%
 left_join(st_obj@meta.data[, c("celltype")], by = c("Juxtaposed_Spot" = "row.names"))

write.table(juxtaposed_df, "Spatial_Juxtaposed_Cells.txt", sep="\t", quote=FALSE, row.names=FALSE)
