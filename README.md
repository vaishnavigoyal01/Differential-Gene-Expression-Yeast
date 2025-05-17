🧬 NGS-Data-Analysis

This repository contains scripts, results, and documentation for a comprehensive Next-Generation Sequencing (NGS) Data Analysis pipeline, focusing on Differential Gene Expression (DGE) using command-line tools and R-based packages.

📌 Objective To perform a step-by-step analysis of RNA-Seq data to identify differentially expressed genes between experimental conditions, using tools such as FastQC, HISAT2, featureCounts, and DESeq2.

🗂️ Repository Structure NGS-Data-Analysis/ ├── fastqc-output/ # Quality check reports from FastQC and MultiQC ├── feature-count/ # Read count files from featureCounts ├── mapped/ # Aligned reads (SAM/BAM files) ├── reads/ # Input reads (FASTQ files) ├── reference/ # Reference genome files and index ├── summary/ # Final summary results (plots, DEGs, tables) ├── trim/ # Trimmed reads from Trimmomatic ├── README.md # Project documentation

🔧 Tools Used

| Tool | Description |
| FastQC | Quality control of raw reads | | MultiQC | Quality control results summary report | | Trimmomatic | Trimming adapters and low-quality bases | | HISAT2 | Read alignment to reference genome | | Samtools | File conversion and indexing | | featureCounts | Counting mapped reads to genes | | R (DESeq2) | Differential expression analysis |

🧪 Workflow Overview

Raw Reads Quality Check

Tool: FastQC
Output: HTML reports (fastqc-output/)
Trimming Low-Quality Reads

Tool: Trimmomatic
Output: Cleaned FASTQ files (trim/)
Read Mapping

Tool: HISAT2
Reference: Stored in reference/
Output: SAM/BAM files (mapped/)
Counting Reads

Tool: featureCounts
Output: Count matrix files (feature-count/)
Differential Expression Analysis

Tool: DESeq2 in R
Output: DEGs, visualizations (summary/)
📚 References FastQC Documentation HISAT2 Manual featureCounts Manual DESeq2 Bioconductor
