⚠️ **WORK IN PROGRESS**

#  Part III – Variant Visualization

## Table of Contents

- [Introduction](#introduction)
- [Integrative Genomics Viewer (IGV)](#integrative-genomics-viewer-igv)
- [Bash shell scripting for NGS-DNA Somatic analysis](#v-bash-shell-scripting-for-ngs-dna-somatic-analysis)
- [Quality control (QC)](#quality-control-qc--01_qcsh)
- [Trimming + QC](#trimming--qc--02_trimsh)
- [Alignment and BAM preprocessing](#alignment-and-bam-preprocessing--03_align__bam_preprocesssh)
    - [Folder structure: From QC – Trimming/Filtering – Alignment + BAM preprocessing](#folder-structure-from-qc---trimmingfiltering---alignment--bam-preprocessing)



## Introduction

The last part of this basic tutorial on DNA-NGS analysis of a small FASTQ dataset consists of visualizing the reads mapped to the reference genome and confirming the called variants shown in [Part II – Somatic analysis](README_somatic_analysis_Part2-3.md) using the dataset `SRR30536566`. A software widely used for this purpose is called **IGV**.

In this section then, I will show you:

- what IGV is

- how to install it

- how to load the data

- and how IGV helps assess whether variants are real or artifacts

> **Note**: This tutorial will not cover complex genomic alterations such as copy number variations (CNVs), translocations, or other types of structural variants (SVs).

**Documentation**

- <https://igv.org/>: Check "**Citing IGV**". Recommendation reference: <http://cancerres.aacrjournals.org/content/77/21/e31.long>
- <https://igv.org/doc/desktop/>: Installation
- **Best practices for variant calling in clinical sequencing**. DOI: <https://doi.org/10.1186/s13073-020-00791-w>
- Dr. Sarah Dean, PhD - <https://www.youtube.com/watch?v=VyNpu3ubAGY>: Excellent tutorial!


## Integrative Genomics Viewer (IGV)

IGV is an open-source tool for visualization of NGS data, allowing the observation of mapped reads to the reference genome from NGS datasets and the understanding of different types of variant calls from genomes more intuitively. In this sense, it's possible to explore different type of mutations, including single nucleotide variants (SNVs) and small insertions/deletions (indels), SV, but also RNA interference screens, gene expression, methylation and genomic annotations ([Robison, et al. 2011](https://www.nature.com/articles/nbt.1754)).


### IGV: Installation

Go to <https://igv.org/doc/desktop/#DownloadPage/> 

1. Click on the type of computer you have to download IGV. For those working with Mac, pay attention to the chip architecture your computer has (Apple or Intel Chip). Because IGV works with Java, I would recommend to download the package of IGV that **includes** Java.

2. Once downloaded, unpack it. Then move the icon to Applications (for those having Mac).

3. Open IGV - double click

### Type of files

Ideally, in order to have a comprehensive view of the analyzed DNA-NGS dataset, there should be loaded 5 different types of files. The files used for the visualization are those derived from the analysis of dataset `SRR30536566`, which was done in  👉 [Part II – Somatic analysis](README_somatic_analysis_Part2-3.md). **Table 1** shows a full list of files necessary for IGV visualization.

**Table 1A**: Files to upload into IGV.

| Format File | File from `SRR30536566` | Purpose in IGV | Why It's Needed |
|-------------|-------------------------|----------------|-----------------|
| **Reference FASTA** (optional but recommended) | `Homo_sapiens_assembly38.fasta` | Provides the **reference genome sequence** against which the reads were aligned | IGV must use the same reference genome build (and contig naming) that was used for alignment. Loading your exact FASTA ensures coordinate compatibility and correct reference base display.<br>⚠️ **Problem**: It's not easy to load into IGV. |
| **BAM** | `SRR30536566.sorted.markdup.md.bam` | Shows **aligned sequencing reads** mapped to the reference genome with coverage depth and read details | Core visualization file - shows actual read alignments, base quality scores, chromosome coordinates and how reads support variant calls. BAM + Index together are the minimum files for IGV visualization. |
| **BAM Index** | `SRR30536566.sorted.markdup.md.bam.bai` | **Index** for the BAM file allowing fast navigation/access to any genomic position | Without this index, IGV would be extremely slow or unable to load large BAM files. |
| **VCF** (compressed) | `SRR30536566.postfiltered.vcf.gz` | Displays **variant calls** as colored bars showing variant positions, alleles, and quality metrics | Allows visual confirmation of variant calls against aligned reads. |
| **VCF Index** | `SRR30536566.postfiltered.vcf.gz.tbi` | **Index** for the VCF file enabling quick navigation to variant positions | Essential for random access in compressed VCF files. |
| **GTF** (compressed) | `gencode.v38.annotation.gtf.gz` | Shows **gene annotations** - gene boundaries (start/stop genomic positions), codons & exons, gene strand (pos/neg), transcripts as separate tracks.<br> Some **.gtf** files provide information about the name and location of the specific primers used for targeted amplicon generation | Provides biological context (is variant in a gene? in an exon? which transcript?) |
| **BED**<br>(browser extensible data) | `crc_panel_7genes_sorted.hg38.bed` | Highlights **targeted sequencing regions** as a colored track | Shows which genomic regions were actually sequenced in your panel.<br>Also shows start/stop regions of targeted genes. Usually, this file is provided by the producer of the gene panel to show the gene regions targeted by the primers they created. |


**Table 1B**: Germline Database File (**Optional**) to upload into IGV.

| Format File | File from Resources | Purpose in IGV | Important Notes |
|-------------|---------------------|----------------|-----------------|
| **VCF** (germline database) | `af-only-gnomad.hg38.vcf.gz`<br>(from `/somatic_resources/`) | Shows **population allele frequencies** at known germline variant positions | **⚠️ WARNING**: This shows common germline variants in the general population, NOT a patient-matched normal. It helps identify variants that are likely germline (common in population) vs. potentially somatic (rare). |

>**Note**: VCF, GTF and BED files are **supporting files**, that means, they will help to understand better the visualization.

### Using Population Databases as Reference

Since no matched normal sample is available for `SRR30536566`, you can use population databases as a reference:

1. **gnomAD VCF** (`af-only-gnomad.hg38.vcf.gz`):
   - Shows **population allele frequencies** at millions of positions
   - **Purpose**: Identify variants that are common in the general population (likely germline)
   - **Limitation**: Cannot prove a variant is somatic, only suggests it might be germline if common

2. **What to look for**:
   - If your variant overlaps a gnomAD position with high allele frequency (>1%), it's likely a germline variant
   - If no gnomAD entry exists or frequency is very low (<0.1%), it *could* be somatic
   - **Important**: Rare germline variants exist, so absence from gnomAD doesn't guarantee somatic status

3. **Clinical context**:
   - In clinical practice, **this is insufficient** - you need patient-matched normal

4. For research/tutorial purposes, it provides useful context

### Opening and Loading files into IGV

When opening IGV, the window will be almost completely blank as shown in **Figure 1 (left)**

1. Click on **IGV** icon
2. Load "**hg38**" human reference genome
  - **Genomes** > **Download Hosted Genome** → Pop-up window called **Hosted Genomes**
  - Type "Human"
  - Select: **Human (hg38 1kg/GATK)**
3. Load **BAM**, **BED** and **VCF** files
  - **File** > **Load from File...**
4. Type the name of one of the seven genes of dataset `SRR30536566`. For example "**NRAS**".

**Reminder**: There are seven genes from **Targeted sequencing panels**. See 👉 [Part II – Somatic analysis – Variant calling with Mutect2](README_somatic_analysis_Part2-3.md#variant-calling-with-mutect2--04_mutect2sh)

```bash
Final CRC 7-gene BED file contents:
Chr Gene_Start  Gene_End  Gene_Name
chr1	114704468	114716771	NRAS
chr10	87863624	87971930	PTEN
chr12	25205245	25250936	KRAS
chr15	66386836	66491544	MAP2K1
chr19	49635291	49640143	RRAS
chr3	179148113	179240093	PIK3CA
chr7	140719326	140924929	BRAF
```

The IGV will look more or less like in **Figure 1 (right)**. 


**Figure 1**: First glance at IGV (left panel) and after loading BAM and supporting files. Each file has its own track in IGV.

![Figure 1: First glance at IGV (left panel and after loading IGV with the files)](images/IGV_starting_1.png)


> **Note**: The **GTF** file was not loaded into IGV because it requires more than 8 GB of available RAM.
  
>**Key message**: Ideally, **IGV must use the same reference genome that was used for alignment**. Unfortunately, IGV doesn't easily handle .fasta reference genomes as the hosted genomes. There's a way to do so but it implies the use of **igvtools**; however, this is out of the scope of this tutorial.


### Assessing the validity of called variants

### 1. NRAS

From step number 4:

5. Click on "Go" and IGV will display **NRAS**. In **Figure 2**, it's possible to see that **NRAS**:

- is located on **chr1**

- has 7 exons represented as blue thick lines (thin lines are introns) according to the **hg38**

- has a read coverage for mainly exons that extends a bit to the beginning and end of introns

- has a SNV located in the **exon 3** as the **VCF** panel indicates a single base (two blue vertical lines)

- is full covered by the **BED** file

- is downstream the **CSDE1** gene

**Figure 2**: **NRAS** gene structure overview at glance.

![Figure 2](images/IGV_NRAS_zoom_out_1_labeled.png)

To visualize in more detail the SNP:

6. Zoom in by increasing sliding bar to **+** and center at exon 3, where the SNV is.

7. Click on the two blue lines from the **VCF** track.

It will pop up two windows:

- **Genotype Information**: Sample-level, from FORMAT column of **VCF**.

   - Header
   
   ```bash
   Sample: DMBEL-EIDR-071
   Genotype: G/T
   Type: HET
   ```
   ✔ Heterozygous NRAS Q61 mutation
   
   ✔ Matches expectation
   
   - Key fields explained
   
   **AD**: `648,115`
   Allele Depth
   
   - 648 reads supporting reference (G)
   - 115 reads supporting alternate (T)
   
   **Allele Frequency: 0.154** → Calculated as: `115 / (648 + 115) ≈ 0.154`
   
   **Depth**: 763 → Total reads at this position (good coverage).
   
   **Strand Bias**: `328,320,61,54` ✔ Alt allele well balanced → trustworthy variant
   
   ```bash
   Ref Forward: 328
   Ref Reverse: 320
   Alt Forward: 61
   Alt Reverse: 54
   ```
   
   **Quality**: `-1`
   This is normal for somatic callers (Mutect2)
   

- **Variant attributes**: Site-level (INFO fields)

   - Header

   ```bash
   Reference: G*
   Alternate: T
   Type: SNP
   ```

   - **TLOD**: `323.24`  ✔ extremely strong
   Tumor Log Odds → **This is the most important somatic quality metric**
   It measures confidence that variant ≠ sequencing noise
     - 10 is decent
     - 100 is excellent
   
   - **AS_SB_TABLE**: `[328, 320|61, 54]`  ✔ No strand bias
   Same strand-bias numbers as above:
   
   ```bash
   Ref: 328 F / 320 R
   Alt: 61 F / 54 R
   ```

   - **MMQ**: `[60, 60]` ✔ Maximum MQ → excellent alignment
   Median Mapping Quality
     - Ref reads: 60
     - Alt reads: 60

  - **MBQ**: `[41, 41]` ✔ Both alleles high quality
  Median Base Quality
  
  - **GERMQ**: `93`
  Germline quality score
    - High value means **unlikely germline**
    - Supports somatic interpretation
    
  - **POPAF**: `5.6` Estimated population allele frequency (% × 10⁻⁴ scale depending on caller)
  Low value → not common in population
  
  - **AS_FilterStatus**: SITE
  Variant passed all site-level filters.

8. Click on the amino acid **Q** on track "Sequence" at the same position of the SNV

It will pop up one window:

```markdown
name: NRAS
location: chr1:114704469-114716771 (-)
id: NM_002524.5
--------------
Exon number: 3
Amino acid coding number: 61
chr1:114713800-114713978
https://www.ncbi.nlm.nih.gov/gene/?term=NM_002524.5
```

**Figure 3** shows the pop up windows with the information about variant attribute, genotype and amino acid of NRAS's SNV.

![Figure 3](images/IGV_NRAS_zoom_in_2_labeled.png)

### IGV-based variant interpretation
***Visual inspection in IGV (reference genome: hg38 1kg/GATK) confirmed a heterozygous somatic single-nucleotide variant in the NRAS gene (NM_002524.5), located in exon 3 at genomic position chr1:114,713,909 (G>T). This variant affects codon 61, resulting in an amino-acid substitution from glutamine (Q) to lysine (K) (p.Gln61Lys, Q61K).***
***The variant is supported by high read depth (~760×), a variant allele frequency of ~15%, balanced forward and reverse strand representation, high mapping and base qualities, and a strong somatic log-odds score (TLOD = 323), with no evidence of strand bias or sequencing artifact. These features are consistent with a high-confidence somatic NRAS Q61K mutation.***


**Visit** 👉 [Part II – Somatic analysis – Final clinical report table](README_somatic_analysis_Part2-3.md#final-clinical-report-table) and compared the information about NRAS (G>T | Q61) variant with the one provided by IGV.








   
