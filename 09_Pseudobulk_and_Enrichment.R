
########### 09_Pseudobulk_and_Enrichment.R #########
library(Seurat)
library(DESeq2)
library(clusterProfiler)
library(dplyr)

# 1. Pseudo-bulk Differential Expression Analysis
st_obj <- readRDS("ST_integrated_final.rds")

# Aggregate expression to sample level
pseudo_counts <- AggregateExpression(st_obj, 
                                     assays = "Spatial", 
                                     slot = "counts", 
                                     return.seurat = FALSE, 
                                     group.by = "orig.ident")$Spatial

colnames(pseudo_counts) <- c("S01_CK", "S02_CK", "S11_LTSS", "S12_LTSS")

# Setup metadata for DESeq2
colData <- data.frame(
 row.names = colnames(pseudo_counts),
 condition = factor(c("CK", "CK", "LTSS", "LTSS"), levels = c("CK", "LTSS"))
)

# Run DESeq2 pipeline
dds <- DESeqDataSetFromMatrix(countData = pseudo_counts, colData = colData, design = ~ condition)
dds <- DESeq(dds, fitType = 'local', minReplicatesForReplace = 7, parallel = FALSE) 

res <- results(dds)
res_df <- as.data.frame(res)

# Filter DEGs (Methods: Fold Change >= 2, adjusted P-value < 0.05)
# FIXED: Changed from 'pvalue' to 'padj' to match manuscript
res_deg <- res_df[which(abs(res_df$log2FoldChange) >= 1 & res_df$padj < 0.05), ]
res_deg <- res_deg[order(res_deg$padj, decreasing = FALSE), ]

write.table(res_deg, file = "Pseudobulk_DEGs.txt", quote = FALSE, sep = "\t")


# 2. Functional Enrichment Analysis (GO & KEGG)

# Load custom databases (Paths should be relative in GitHub repository)
Kterm2gene <- read.table("term2gene_kegg.db", col.names = c("KO", "GID"))
Kterm2name <- read.table("erm2name_kegg.db", sep = "\t", col.names = c("KO", "Pathway"))

Gterm2gene <- read.table("term2gene_go.db", col.names = c("GO", "GID"))
Gterm2name <- read.table("term2name_go.db", sep = "\t", col.names = c("GO", "term", "ONTOLOGY"))
Gterm2nameBP <- Gterm2name[Gterm2name$ONTOLOGY == "BP", ]

# Define target gene list (e.g., DEGs from Pseudobulk or MiloDE)
target_genes <- rownames(res_deg)

# KEGG Enrichment
# Methods: adjusted significance threshold of 0.1
kegg_res <- enricher(gene = target_genes,
                     TERM2GENE = Kterm2gene,
                     TERM2NAME = Kterm2name,
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.1, 
                     qvalueCutoff = 0.1)

write.table(as.data.frame(kegg_res), file = "Enrichment_KEGG_Results.txt", sep="\t", quote=FALSE, row.names=FALSE)

# GO-BP Enrichment
go_res <- enricher(gene = target_genes,
                   TERM2GENE = Gterm2gene,
                   TERM2NAME = Gterm2nameBP,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.1, 
                   qvalueCutoff = 0.1)

write.table(as.data.frame(go_res), file = "Enrichment_GO_BP_Results.txt", sep="\t", quote=FALSE, row.names=FALSE)


