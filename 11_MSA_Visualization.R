
########### 11_MSA_Visualization.R #########

library(Biostrings)
library(ggmsa)
library(ggplot2)

# 1. Load trimmed alignment sequences
seqs <- readAAStringSet("mulit_cax1_pep.aln.tri")

# 2. Rename sequence IDs to species names
# Define species prefix mapping
species_map <- c(
 "Brara." = "Brassica rapa",
 "AT2G|AtCAX" = "Arabidopsis thaliana",
 "Medtr" = "Medicago truncatula",
 "Glyma" = "Glycine max",
 "VIT_" = "Vitis vinifera",
 "Cucsa." = "Cucumis sativus",
 "Nta" = "Nicotiana tabacum",
 "SlCAX" = "Solanum lycopersicum",
 "PGSC0003DMG" = "Solanum tuberosum",
 "Bevul." = "Beta vulgaris",
 "AH" = "Amaranthus hypochondriacus",
 "Spov3_chr" = "Spinacia oleracea",
 "Bradi" = "Brachypodium distachyon",
 "Sobic." = "Sorghum bicolor",
 "Zm00001d" = "Zea mays",
 "LOC_Os" = "Oryza sativa"
)

# Apply renaming logic
new_names <- names(seqs)
for (pattern in names(species_map)) {
 new_names[grepl(pattern, new_names)] <- species_map[[pattern]]
}
names(seqs) <- new_names

# 3. Sort sequences based on predefined evolutionary order
evolutionary_order <- c(
 "Arabidopsis thaliana", "Brassica rapa", "Medicago truncatula", "Glycine max",
 "Vitis vinifera", "Cucumis sativus", "Nicotiana tabacum", "Solanum lycopersicum",
 "Solanum tuberosum", "Beta vulgaris", "Amaranthus hypochondriacus", "Spinacia oleracea",
 "Brachypodium distachyon", "Sorghum bicolor", "Zea mays", "Oryza sativa"
)

# Keep only matched sequences and order them
matched_idx <- match(evolutionary_order, names(seqs))
matched_idx <- matched_idx[!is.na(matched_idx)]
seqs_sorted <- seqs[matched_idx]

# 4. Visualize the conserved N-terminal autoinhibitory domain using ggmsa
# Coordinates (start=20, end=40) target the specific S-cluster region
p_msa <- ggmsa(seqs_sorted, 
               start = 20, end = 40, 
               color = "Chemistry_AA",
               font = "mono",                 
               char_width = 0.7,              
               seq_name = TRUE,
               border = "grey25",
               consensus_views = TRUE,
               use_dot = FALSE,               
               disagreement = FALSE,           
               ignore_gaps = TRUE) + 
 labs(title = "CAX1 orthologs N-terminal MSA") +
 theme_msa() + 
 theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
       axis.text.x = element_text(size = 8))

# Save publication-ready plot
ggsave("Figure_CAX1_MSA_Aligned.pdf", plot = p_msa, width = 8, height = 4, dpi = 600, bg = "white", device = "pdf", useDingbats = FALSE)




