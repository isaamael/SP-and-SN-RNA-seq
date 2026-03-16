######## 02_Bulk_DESeq2.R #######
library(DESeq2)

# 1. Load data
# Note: Ensure the working directory is correctly set or use relative paths
counts_file <- "All_gene_counts.list"
count_df <- read.table(counts_file, header = TRUE, stringsAsFactors = FALSE)

# Remove 'geneLength' column if present and format the matrix
if("geneLength" %in% colnames(count_df)) {
 count_df <- count_df[, !names(count_df) %in% "geneLength"]
}

# Set gene IDs as rownames and convert to numeric matrix
countData <- as.matrix(count_df[, -1])
rownames(countData) <- count_df$ID

# 2. Prepare sample metadata (colData)
# Assuming the order in matrix is MCK_1, MCK_2, MST_1, MST_2, PCK_1...
sample_conditions <- factor(c(rep("MCK", 2), rep("MST", 2), rep("PCK", 2), rep("PST", 2)))
colData <- data.frame(
 condition = sample_conditions,
 row.names = colnames(countData)
)

# 3. Create DESeq2 object and run analysis
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData, 
                              design = ~ condition)

dds <- DESeq(dds)

# 4. Extract results for MCK vs MST
# To compare other groups, adjust the 'contrast' argument
res <- results(dds, contrast = c("condition", "MST", "MCK"))

# 5. Filter Differentially Expressed Genes (DEGs)
# Thresholds defined in Methods: Fold Change >= 2 (log2FC >= 1) and padj < 0.05
fc_threshold <- 2  
fdr_threshold <- 0.05      

res_filtered <- res[which(abs(res$log2FoldChange) >= log2(fc_threshold) & res$padj < fdr_threshold), ]
res_filtered <- res_filtered[order(res_filtered$padj), ] # Order by significance

DEG_df <- as.data.frame(res_filtered)
DEG_df <- na.omit(DEG_df)

# 6. Save results
write.table(DEG_df, file = "MCK_vs_MST_DEGs.txt", sep = "\t", quote = FALSE, row.names = TRUE)

message("DESeq2 analysis complete. Number of DEGs: ", nrow(DEG_df))

