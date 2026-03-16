
########### 07_scWGCNA_Analysis.R #########

library(Seurat)
library(scWGCNA)
library(dplyr)
library(reshape2)

# Load MiloDE results
milo_results <- read.table("MiloDE_summary_results.txt", header = TRUE, sep = "\t")

# Filter highly perturbed neighborhoods and extract significant DEGs
# Methods: NHs with dominant cell type > 80%, AUC > 0.7, DEGs in >= 2 NHs
nh_valid <- milo_results %>% filter(celltype_fraction >= 0.8 & auc > 0.7)
de_genes_freq <- milo_results %>% group_by(gene) %>% summarise(n_nhoods = sum(pval_corrected < 0.1, na.rm = TRUE))
sig_genes <- de_genes_freq$gene[de_genes_freq$n_nhoods >= 2]

milo_filtered <- milo_results %>% filter(gene %in% sig_genes & Nhood %in% nh_valid$Nhood)

# Construct gene-by-neighborhood matrix
# Methods: Entries represent log2FC, non-DEGs set to zero
milo_filtered$logFC[milo_filtered$pval_corrected >= 0.1] <- 0
gene_nh_matrix <- dcast(data = milo_filtered, formula = gene ~ Nhood, value.var = "logFC")
rownames(gene_nh_matrix) <- gene_nh_matrix$gene
gene_nh_matrix <- gene_nh_matrix[, -1]

# Create pseudo-Seurat object for scWGCNA
wgcna_obj <- CreateSeuratObject(counts = gene_nh_matrix)
DefaultAssay(wgcna_obj) <- "RNA"
wgcna_obj <- FindVariableFeatures(wgcna_obj)
wgcna_obj[["RNA"]]@scale.data <- as.matrix(wgcna_obj[["RNA"]]@data)

# Run scWGCNA framework
# Methods: Soft threshold = 10, min module size = 15, cut height = 0.25
modules_wgcna <- run.scWGCNA(p.cells = wgcna_obj, 
                             s.cells = wgcna_obj, 
                             is.pseudocell = FALSE, 
                             features = rownames(wgcna_obj),
                             power = 10,
                             minModuleSize = 15,
                             mergeCutHeight = 0.25,
                             merging = TRUE)

# Compile and save module statistics
cluster_list <- lapply(seq_along(modules_wgcna$module.genes), function(i) {
 data.frame(module = i, gene = modules_wgcna$module.genes[[i]])
})
module_df <- do.call(rbind, cluster_list)

write.table(module_df, "scWGCNA_Gene_Modules.txt", sep="\t", quote=FALSE, row.names=FALSE)
saveRDS(modules_wgcna, "scWGCNA_Network_Object.rds")