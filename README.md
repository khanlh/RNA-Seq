# RNA-seq Nextflow Pipeline with Nextflow 

## 1. Basic Biological Background  

**RNA-seq (RNA sequencing)** is a next-generation sequencing (NGS) technology used to measure gene expression levels across the transcriptome.  

The general workflow:  
1. Extract RNA from biological samples  
2. Convert RNA to cDNA and prepare sequencing libraries  
3. Sequence cDNA to produce **FASTQ files** (reads + quality scores)  
4. Process data to quantify **gene expression** and detect potential variants  

**Applications:**  
- Compare gene expression between conditions (e.g., disease vs. control)  
- Identify novel genes or splicing isoforms  
- Detect sequence variants from transcriptome data  

---

## 2. Pipeline Steps  

### a. Quality Control & Trimming – Fastp  
- Removes adapters, low-quality bases, or Ns.  
- Uses **sliding window filtering** and quality distribution.  
- Generates QC reports (HTML, JSON).  

### b. Alignment – HISAT2  
- Aligns reads to a reference genome.  
- Based on **Burrows–Wheeler Transform (BWT)** and **FM-index**.  
- Optimized for RNA-seq with **spliced alignment** (exon–exon junctions).  

### c. Read Counting – featureCounts  
- Assigns aligned reads to genomic features (exons/genes).  
- Algorithm: read → mapped position → overlap with exon regions → increment count.  
- Produces a **counts matrix** (genes × samples).  

### d. Differential Expression (Downstream)  
- Use R packages (**DESeq2**, **edgeR**) on counts matrix.  
- **Log2 Fold Change (log2FC):**  

$$
\text{log2FC} = \log_2 \left( \frac{\text{Expression}_{\text{Condition A}}}{\text{Expression}_{\text{Condition B}}} \right)
$$

  **Interpretation:**  
  - log2FC = 1 → 2× higher in A vs. B  
  - log2FC = -1 → 2× lower in A vs. B  
  - log2FC ≈ 0 → no difference  

### e. Variant Calling (Optional)  
- With **samtools/bcftools**, SNPs and INDELs can be called.  
- RNA-seq variant calling has biases (coverage, allele-specific expression).  

---

## 3. System Requirements  

- [Nextflow](https://www.nextflow.io/) ≥ 22.10.0  
- Tools available in the PATH:  
  - fastp  
  - hisat2  
  - subread (featureCounts)  
  - samtools, bcftools  

Quick install with Conda:  
```bash
conda create -n rnaseq \
  fastp hisat2 subread samtools bcftools -c bioconda -c conda-forge
conda activate rnaseq
