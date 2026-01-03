# Deciphering Asymmetric Regulatory Logic in 17q-Amplified Breast Cancer

![Status](https://img.shields.io/badge/Status-Research_Prototype-blue)
![Languages](https://img.shields.io/badge/Languages-R_%7C_Python_%7C_Perl_%7C_Shell_%7C_Cypher-orange)
![Graph DB](https://img.shields.io/badge/Database-Neo4j-green)

## 🧬 Project Overview

This repository contains the complete analysis pipeline for the research project: **"Deciphering Asymmetric Regulatory Logic in 17q Amplified Breast Cancer via Boolean Modeling and BRCA1 Stratification."**

### Scientific Objective
Chromosomal amplification at 17q is a driver of breast cancer, yet linear correlation methods fail to capture the asymmetric ("if-then") regulatory logic governing this system. This project utilizes a **Boolean Implication Network** approach to map these directed interactions across transcriptomic (RNA-seq) and epigenomic (DNA Methylation) layers.

By stratifying patients based on **17q Copy Number Status** (GAIN vs. DIS) and **BRCA1 Expression** (Up vs. Down), this workflow uncovers hidden topological architectures—such as "lock-in" regulatory loops and network bottlenecks—that predict therapeutic vulnerabilities.

---

## 📂 Repository Structure & Workflow

The code is organized by analytical stage, processing data from raw TCGA downloads to graph database visualization.

### 1. Data Acquisition
**Directory:** `tcga-brca-integration/`
* `tcga_brca_integration.r`: Uses `TCGAbiolinks` to download and integrate RNA-seq, DNA Methylation, and CNV data for the BRCA cohort, filtering for complete multi-omics profiles ($n=429$ tumor samples).

### 2. Differential Analysis
**Directory:** `differential-analysis/`
Performs `limma`-based differential expression (DEG) and methylation (DMP) analysis to select significant features for network construction.
* **Chromosome 17q Global Contrasts:** `deg_chr17q_script.r` and `dmp_chr17q_script.r` analyze GAIN-Chr17q vs. DIS-Chr17q vs. Control.
* **BRCA1 Stratification:** `deg_brca1_variant_script.r` and `dmp_brca1_variant_script.r` analyze specific *upBRCA1* vs. *downBRCA1* subsets.

### 3. Feature Annotation & Preparation
**Directories:** `gene-probe-anno/` & `pre-stepminer/`
* `gene_annotation_refactor.r`: Maps Ensembl and Probe IDs to gene symbols and chromosomal locations.
* `pre-stepminer/`: Refactors and formats expression/methylation matrices for specific cohorts (GAIN, DIS, upBRCA1, downBRCA1) to match the input requirements of the StepMiner algorithm.

### 4. Boolean Network Inference (Core Pipeline)
**Directory:** `stepminer–booleannet-workflow/`
The core computational engine for discretizing continuous data and inferring logic.
* `StepMiner_algorithm.ipynb` & `stepminer-1.1.jar`: Discretizes continuous omics data into Boolean states (Low, Intermediate, High).
* `booleannet_pipeline.sh`: Shell script orchestrating the `BooleanNet` algorithm to identify asymmetric implications (e.g., *High → High*, *Low → Low*).
* `extract_exp.pl` / `extract_met.pl`: Perl scripts to parse and format the raw BooleanNet output.

### 5. Statistical Validation
**Directory:** `permutation-test/`
Implements rigorous permutation testing ($K=50$ iterations) to control False Discovery Rate (FDR).
* `permutation_test_script_exp.sh` / `_met.sh`: Randomizes patient samples while preserving data distribution to generate null models.
* `extract_permutation.pl`: Extracts statistics to calculate empirical FDR and determine optimal thresholds.

### 6. Multi-Omics Integration & Network Analysis
**Directories:** `interlayer-construction/` & `network_analysis/`
* **Interlayer Construction:** Scripts like `brca1_interlayer_integration.r` map DNA methylation nodes to gene expression nodes (Same-Gene Interlayer Mapping) to model cis-regulatory effects.
* **Network Analysis:** Scripts like `brca1_network_ud_contrast.r` and `gain_dis_network_gd_contrast.r` analyze topology (centrality, degree distribution) across the different biological contrasts.

### 7. Graph Database Management
**Directory:** `neo4j-cypher/`
* `core_neo4j_query_script.cypher`: Cypher queries to load nodes/edges into Neo4j and perform graph-based queries (e.g., identifying Ego-Networks centered on BRCA1).
* `pre-prosessing-neo4j.txt`: Guidelines for formatting CSVs for Neo4j import.

---

## 💻 Prerequisites & Dependencies

To run the full pipeline, the following tools are required:

* **R (v4.5.2+):** `TCGAbiolinks`, `limma`, `dplyr`, `igraph`.
* **Python (v3.12) & Jupyter:** Required for running the StepMiner notebook.
* **Java Runtime Environment (JRE):** Required for `stepminer-1.1.jar`.
* **Perl:** Required for output parsing scripts (`.pl`).
* **Neo4j:** For graph database management and visualization.

---

## 📊 Key Findings

* **Asymmetric Dominance:** The regulatory landscape of 17q-amplified tumors is dominated by asymmetric subset implications (*High $\rightarrow$ High*), reflecting a "lock-in" of oncogenic states.
* **BRCA1 Topology:**
    * *upBRCA1* networks preserve regulatory heterogeneity.
    * *downBRCA1* networks show program consolidation, serving as a topological proxy for Homologous Recombination Deficiency (HRD).
* **Epigenetic Decoupling:** We identified a paradox where hypermethylation of 17q genes (e.g., *MAPT*, *BRIP1*) correlates with transcriptional upregulation, driven by gene dosage effects.

---

## 📜 Citation

...

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
