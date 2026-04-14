# Dietary Fiber-Driven Microbiome-Lipid Axis Maintains Brain Homeostasis through Promoting Adult Neurogenesis


**Repository for the code and data used in the manuscript revision**

## 📄 Manuscript Abstract

Brain homeostasis is fundamental to cognitive function, yet its disruption by neuroinflammation, neuronal injury, or degenerative processes leads to cognitive impairment across diverse neurological conditions. Here we demonstrate that a structurally defined, highly soluble, short-chain-enriched inulin (HSI) effectively maintains brain homeostasis through a microbiome-dependent mechanism. Using the 5xFAD mouse model of amyloid pathology, we show that HSI rescues recognition memory, reduces plaque burdens, and dampens microglial and astrocytic activation. Single-nucleus RNA sequencing revealed a shift toward neurogenesis-linked transcriptional programs in neuronal populations, accompanied by suppression of stress and inflammatory pathways in microglia. Multi-compartment microbiome-metabolome profiling (gut, circulation, and brain) converged on lysophosphatidylethanolamine (LPE) 18:2 as the only gut-derived lipid consistently elevated after HSI supplementation. Notably, systemic administration of LPE 18:2 at the onset of pathology improved cognition, lowered amyloid deposition and microglial cytokine responses, and promoted adult hippocampal neurogenesis toward mature neurons. Mechanistically, LPE 18:2 binds to the orphan GPCR GPR37 to enhance Wnt/β-catenin signaling through facilitating LRP6 maturation, thereby driving neural stem cell proliferation and neuronal differentiation. These findings delineate a diet-driven fiber-microbiome-lipid-receptor axis that maintains brain homeostasis and nominate LPE 18:2 as a tractable, diet-linked mediator with broad therapeutic potential for conditions involving cognitive impairment and neuroinflammation.


---

## 📁 Repository Structure

```
├── scripts/                     # Scripts used for data processing and analysis
│   ├── scRNAseq_analysis.R      # Main Seurat pipeline for single-cell analysis
│   └── DE_analysis.R            # Differential expression analysis code
│
├── README.md                    # This file
└── requirements.txt             # R package dependencies (see below)
```

---

## 🔧 Setup and Requirements

### For Single Cell RNA Sequencing Data Analysis:

We used [R (≥4.2)](https://www.r-project.org/) with the following R packages:
* `Seurat` ≥ 5.0
* `tidyverse`
* `ggplot2`
* `flowCore`
* `ComplexHeatmap`
* `Gprofiler2`
* `clusterProfiler`
* `enrichplot`
* `hdwgcna`



---

## 🧬 Data Availability

All datasets (processed and/or raw) will be made available upon acceptance via a public repository (e.g., GEO, FlowRepository, Zenodo).

Please refer to the `data/README.md` file in each subfolder for specific descriptions and usage notes.

---



## 🖋️ Citation

If you use this code, please cite:

**\[Author list TBD]**
"Dietary Fiber-Driven Microbiome-Lipid Axis Maintains Brain Homeostasis through Promoting Adult Neurogenesis"


---

## 📬 Contact

For questions regarding the code or data, please contact:

**Yiqiao Wang** – \[[yiqiaowang@jnu.edu.cn](mailto:yiqiaowang@jnu.edu.cn)]

**Yuxi Guo** – \[[yuxi.guo@ki.se](mailto:yuxi.guo@ki.se)]



