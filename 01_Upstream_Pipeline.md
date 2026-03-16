# Upstream Data Processing Pipeline


## 1. Single-Nuclei RNA-seq (SN-seq) Data Processing

```bash
# Example command for SN-seq processing
cellranger count --id=SN_Rep1 \
                 --transcriptome=/path/to/ref_tomato_sl4 \
                 --fastqs=/path/to/raw_fastq/ \
                 --sample=SN_Rep1 \
                 --localcores=16 \
                 --localmem=64
				 
## 2.  Bulk RNA-seq Data Processing				 
# Step 1: Quality control and adapter trimming
fastp -i raw_R1.fq.gz -I raw_R2.fq.gz \
      -o clean_R1.fq.gz -O clean_R2.fq.gz \
      -w 8 -q 20 -u 40

# Step 2: Sequence alignment to the tomato reference genome (Sl4)
hisat2 -p 16 -x /path/to/hisat2_index/tomato_sl4 \
       -1 clean_R1.fq.gz -2 clean_R2.fq.gz \
       -S aligned_reads.sam

# Step 3: Convert SAM to sorted BAM
samtools view -bS aligned_reads.sam | samtools sort -@ 8 -o sorted_reads.bam
samtools index sorted_reads.bam

# Step 4: Gene expression quantification
featureCounts -T 16 -p -t exon -g gene_id \
              -a /path/to/annotation/tomato_sl4.gtf \
              -o gene_counts.txt sorted_reads.bam