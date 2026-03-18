---

### README.md


This repository contains the custom scripts and analysis pipelines used in the manuscript "SlCAX3 drives the formation of crystal idioblasts for tomato ion compartmentalization under salt stress".



**Overview**  
The codebase covers the complete downstream analysis of spatial transcriptomics (ST-seq), single-nuclei RNA-seq (SN-seq), and bulk RNA-seq data, including:
- Data quality control, SCTransform normalization, and RPCA integration.
- Spatial modeling (MiloDE, miloR, scDist).
- Spatial deconvolution (CARD) and cell-cell interaction (CellChat).
- Trajectory inference (Slingshot, Monocle2) and gene co-expression networks (scWGCNA).




**Data Availability**  
The raw sequencing data and processed matrices (rds and count matrices) required to run these scripts have been deposited in the NCBI Gene Expression Omnibus (GEO) under accession number xxxxxx.




**System Requirements and Dependencies**  
All downstream analyses were performed in `R` (version 4.x.x). The major R packages used include:
- `Seurat` (v4.3.0) - *Note: Seurat v4 is required for this pipeline.*
- `CARD` (v1.1)
- `CellChat` (v1.6.1)
- `scWGCNA` / `hdWGCNA`
- `miloR` & `MiloDE`
- `monocle` (v2.x) & `slingshot`
- `DESeq2`, `clusterProfiler`, `scDist`, `UCell`, `MetaNeighbor`




**Usage**  
The scripts are numbered sequentially (`01` to `11`) according to the logical flow of the analysis described in the **Materials and Methods** section of the manuscript. 
To reproduce the findings, users should download the processed data from GEO, adjust the file paths in the scripts accordingly, and run them sequentially.




**License**  
This project is licensed under the MIT License.



**Contact**
For any questions regarding the code or data, please contact at isaac2607388620@gmail.com.
