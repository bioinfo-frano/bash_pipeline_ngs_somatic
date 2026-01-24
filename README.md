# Bash_pipeline_NGS

Welcome to my **DNA-NGS tutorial** 👋  

This repository shows you how to analyze **DNA sequencing datasets**  
using **Bash script pipelines**, starting from raw FASTQ data  
and running everything on a **home workstation or laptop**.

You can follow the tutorial step by step, or jump directly to the analysis part.

---

## Tutorial structure

### 1️⃣ Part I – Preparation & setup  
Learn how to prepare a clean and reproducible environment:
- Folder structure
- Reference genome setup and integrity checks
- SRA data selection and download
- Conda environments and tool installation

➡️ **Start here:**  
👉 [Part I – Preparation & setup](README_setup.md)

---

### 2️⃣ Part II – Somatic analysis (Bash pipelines)  
Perform a **somatic DNA-NGS analysis** following GATK best practices:
- FASTQ processing and QC
- Alignment and BAM processing
- Somatic variant calling with **Mutect2**
- Variant filtering and annotation

➡️ **Go to analysis:**  
👉 [Part II – Somatic analysis](README_somatic_analysis.md)

---

## 🔮 Future extensions

This repository is designed to grow. Planned additions include:

- **Part III – Somatic - IGV analysis**
  - Learn how to use IGV using the outputted sorted .BAM, .BAI and .VCF files from 👉 [Part II – Somatic analysis](README_somatic_analysis.md)
  - Check potential artifacts and confirm annotated variants

- **Part IV – Germline analysis**
  - Additional datasets
  - Pipeline optimizations and best practices

---

## 🧬 Target audience

This tutorial is intended for:
- Bioinformatics students
- Life scientists learning NGS analysis
- Researchers who want a **transparent, Bash-only workflow**

---

## 📌 Notes

- The pipeline is optimized for **educational clarity**, not HPC clusters
- All steps are reproducible and runnable on a local machine
- Real public datasets from **SRA** are used

---

Happy sequencing analysis! 🚀
