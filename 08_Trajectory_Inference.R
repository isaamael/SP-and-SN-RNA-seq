
########### 08_Trajectory_Inference.R #########

library(Seurat)
library(slingshot)
library(monocle)
library(SCP)
library(dplyr)
library(BiocParallel)

st_obj <- readRDS("ST_integrated_final.rds")
bpparam <- SnowParam(workers = 8)

# Subset dataset to major mesophyll cells and CIs
traj_obj <- subset(st_obj, idents = c("PMC", "MC1", "MC2", "CMC_0", "CMC_1"))

# Global lineage trajectory inference using Slingshot
# Methods: subCI1 (CMC_1) designated as starting cell type, interpret in reverse
slingshot_res <- RunSlingshot(srt = traj_obj, 
                              group.by = "sub_cluster", 
                              reduction = "UMAP", 
                              start = "CMC_1", 
                              reverse = TRUE)

# Fit Generalized Additive Models (GAMs) to identify pseudotime-dependent genes
# Methods: Top 2,000 highly variable genes modeled against pseudotime
slingshot_res <- FindVariableFeatures(slingshot_res, nfeatures = 2000)
hvg_genes <- VariableFeatures(slingshot_res)

slingshot_res <- RunDynamicFeatures(srt = slingshot_res, 
                                    assay = "Spatial", 
                                    slot = "data", 
                                    lineages = c("Lineage1", "Lineage2", "Lineage3"), 
                                    features = hvg_genes, 
                                    BPPARAM = bpparam)

# Extract and save dynamic features
dynamic_genes <- c()
for (lin in c("Lineage1", "Lineage2", "Lineage3")) {
 dynamic_genes <- unique(c(dynamic_genes, slingshot_res@tools[[paste0("DynamicFeatures_", lin)]]$DynamicFeatures$features))
}
write.table(data.frame(Gene = dynamic_genes), "Slingshot_GAM_Dynamic_Genes.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Monocle 2 pseudotime analysis between CI sub-clusters
# Methods: Isolate subCI0 and subCI1, filter 1500 highly variable genes
sub_ci <- subset(st_obj, idents = c("CMC_0", "CMC_1"))
sub_ci <- FindVariableFeatures(sub_ci, nfeatures = 1500)
monocle_genes <- VariableFeatures(sub_ci)

expr_matrix <- as.matrix(GetAssayData(sub_ci, slot = "data"))[monocle_genes, ]
cell_meta <- new('AnnotatedDataFrame', data = sub_ci@meta.data)
gene_meta <- new('AnnotatedDataFrame', data = data.frame(gene_short_name = rownames(expr_matrix), row.names = rownames(expr_matrix)))

cds <- newCellDataSet(expr_matrix,
                      phenoData = cell_meta,
                      featureData = gene_meta,
                      expressionFamily = uninormal())

# Dimensionality reduction and cell ordering via DDRTree
cds <- estimateSizeFactors(cds)
cds <- setOrderingFilter(cds, monocle_genes)
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree', norm_method = "non")
cds <- orderCells(cds)

# Model gene expression against pseudotime
# Methods: differentialGeneTest used to identify top 100 pseudotime-dependent genes
diff_test_res <- differentialGeneTest(cds[monocle_genes, ], 
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)", 
                                      cores = 8)

diff_test_res <- diff_test_res %>% arrange(qval)
top_100_genes <- head(diff_test_res, 100)

write.table(top_100_genes, "Monocle2_Top100_Pseudotime_Genes.txt", sep="\t", quote=FALSE, row.names=FALSE)
saveRDS(cds, "Monocle2_CI_Subclusters.rds")