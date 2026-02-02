⚠️ **WORK IN PROGRESS**

#  Part III – Variant Visualization

## Table of Contents

- [Introduction](#introduction)
- [Integrative Genomics Viewer (IGV)](#integrative-genomics-viewer-(igv))
- [Bash shell scripting for NGS-DNA Somatic analysis](#v-bash-shell-scripting-for-ngs-dna-somatic-analysis)
- [Quality control (QC)](#quality-control-qc--01_qcsh)
- [Trimming + QC](#trimming--qc--02_trimsh)
- [Alignment and BAM preprocessing](#alignment-and-bam-preprocessing--03_align__bam_preprocesssh)
    - [Folder structure: From QC – Trimming/Filtering – Alignment + BAM preprocessing](#folder-structure-from-qc---trimmingfiltering---alignment--bam-preprocessing)



## Introduction

The last part of this basic tutorial on DNA-NGS analysis of a small FASTQ dataset consists of visualizing the reads mapped to the reference genome and confirm the called variants shown in [Part II – Somatic analysis](README_somatic_analysis_Part2-3.md) using the dataset `SRR30536566`. A software widely used for this purpose is called **IGV**.

In this section then, I will show you:

- what IGV is

- how to install it

- how to load the data

- and how IGV helps assessing whether variants are true from artifacts

>**Note**: This tutorial will not show complex gene alterations like copy number variations (CNV), translocations or other types of structural variant (SV).

**Documentation**

- <https://igv.org/>: Check "**Citing IGV**". Recommendation reference: <http://cancerres.aacrjournals.org/content/77/21/e31.long>
- <https://igv.org/doc/desktop/>: Installation
- **Best practices for variant calling in clinical sequencing **. DOI: <https://doi.org/10.1186/s13073-020-00791-w>
- Dr. Sarah Dean, PhD - <https://www.youtube.com/watch?v=VyNpu3ubAGY>: Excellent tutorial!


## Integrative Genomics Viewer (IGV)

IGV is a visualization tool, allowing the observation of mapped reads to the reference genome from NGS datasets and the understanding of different types of variant calls from genomes more intuitively. In this sense, it's possible to explore different type of mutations, including single nucleotide variants (SNV) such as SNP and Indels, SV, but also RNA interference screens, gene expression, methylation and genomic annotations ([Robison, et al. 2011](https://www.nature.com/articles/nbt.1754)).


### IGV: Installation

Go to <https://igv.org/doc/desktop/#DownloadPage/> 

1. Click on the type of computer you have to downloading IGV. For those working with Mac, pay attention to the Chip your computer has (Apple or Intel Chip). Because IGV works with Java, I would recommend to download the package of IGV that **includes** Java.

2. Once downloaded, unpack it. Then move the icon to Applications (for those having Mac).

3. Open IGV - double click

### Loading the data

Ideally, in order to have a comprehensive view of the analysed DNA-NGS dataset, there should be loaded 5 different types of files. The files used for the visualization are those derived from the analysis of dataset `SRR30536566`, which was done in  👉 [Part II – Somatic analysis](README_somatic_analysis_Part2-3.md). The files are:

- **.bam**: `SRR30536566.sorted.markdup.md.bam`



- **.bam.bai**: `SRR30536566.sorted.markdup.md.bam.bai`



- **.vcf.gz**: `SRR30536566.postfiltered.vcf.gz`



- **.vcf.gz.tbi**: `SRR30536566.postfiltered.vcf.gz.tbi`


- **.gtf.gz**: `gencode.v38.annotation.gtf.gz`


- **.bed** or **sorted.bed**: `crc_panel_7genes_sorted.hg38.bed`


| Format File | File from `SRR30536566` | Purpose in IGV | Why It's Needed |
|-------------|-------------------------|----------------|-----------------|
| **BAM** | `SRR30536566.sorted.markdup.md.bam` | Shows **aligned sequencing reads** mapped to the reference genome with coverage depth and read details | Core visualization file - shows actual read alignments, base qualities, and how reads support variant calls |
| **BAM Index** | `SRR30536566.sorted.markdup.md.bam.bai` | **Index** for the BAM file allowing fast navigation to any genomic position | Without this index, IGV would be extremely slow or unable to load large BAM files |
| **VCF** (compressed) | `SRR30536566.postfiltered.vcf.gz` | Displays **variant calls** as colored bars showing variant positions, alleles, and quality metrics | Allows visual confirmation of variant calls against aligned reads |
| **VCF Index** | `SRR30536566.postfiltered.vcf.gz.tbi` | **Index** for the VCF file enabling quick navigation to variant positions | Essential for random access in compressed VCF files |
| **GTF** (compressed) | `gencode.v38.annotation.gtf.gz` | Shows **gene annotations** - gene boundaries, exons, transcripts as separate tracks | Provides biological context (is variant in a gene? in an exon? which transcript?) |
| **BED** | `crc_panel_7genes_sorted.hg38.bed` | Highlights **targeted sequencing regions** as a colored track | Shows which genomic regions were actually sequenced in your panel |



