# 🧬 In Silico Drug Repurposing for lncRNA-Driven Oral Squamous Cell Carcinoma

> **A computational pipeline to identify FDA-approved drug candidates that reverse the oncogenic gene expression signature driven by dysregulated long non-coding RNAs (lncRNAs) in OSCC.**

---

## 📌 Project Overview

Oral Squamous Cell Carcinoma (OSCC) is an aggressive head and neck malignancy with a 5-year survival rate of ~50%. Recent transcriptomic studies have revealed that long non-coding RNAs (lncRNAs) — RNA molecules >200 nt that do not encode proteins — play critical regulatory roles in OSCC progression, yet remain largely undruggable by conventional means.

This pipeline bridges that gap using the **Connectivity Map (CMap)** framework:

```
OSCC lncRNAs → Downstream Target Genes → Gene Signature → CMap Query → Repurposed Drug Candidates
```

Rather than developing new drugs from scratch (a 15-year, billion-dollar process), we ask:
**"Which existing, safe, FDA-approved drugs pharmacologically reverse the exact gene expression pattern caused by these lncRNAs?"**

---

## 🔬 Scientific Background

### Why lncRNAs?
- The human genome encodes ~16,000 lncRNAs, but only ~2% of the genome codes for protein
- lncRNAs regulate gene expression at epigenetic, transcriptional, and post-transcriptional levels
- In OSCC, key oncogenic lncRNAs include **HOTAIR**, **MALAT1**, **NEAT1**, **H19**, **LINC00152**, and **PVT1**
- These lncRNAs upregulate oncogenes (MMP9, VEGF, CDK6) and silence tumor suppressors (p21, PTEN, E-cadherin)

### Why Drug Repurposing?
- FDA-approved drugs already have known safety profiles — clinical trials are cheaper and faster
- The Connectivity Map (CMap) by the Broad Institute contains transcriptomic profiles of ~5,000 drugs across multiple cell lines
- By matching our "disease signature" to drug perturbation profiles, we can find drugs that **reverse** the cancer state

---

## 🗂️ Repository Structure

```
oscc_lncrna_pipeline/
│
├── README.md                          ← You are here
│
├── data/
│   ├── raw/
│   │   ├── oscc_lncrnas.csv           ← Curated OSCC lncRNA list (from Lnc2Cancer)
│   │   └── rnainter_raw.json          ← Raw RNAInter API responses
│   ├── processed/
│   │   ├── target_genes.csv           ← Extracted target genes with interaction scores
│   │   ├── upregulated_genes.txt      ← ⬆ Upregulated gene list for CLUE.io
│   │   └── downregulated_genes.txt    ← ⬇ Downregulated gene list for CLUE.io
│   └── clue_output/
│       └── clue_results.csv           ← Downloaded results from CLUE.io query
│
├── scripts/
│   ├── 01_fetch_lncrnas.py            ← Step 1: Build OSCC lncRNA list
│   ├── 02_query_rnainter.py           ← Step 2: Get target genes via RNAInter API
│   ├── 03_prepare_gene_signature.py   ← Step 3: Format gene lists for CLUE.io
│   ├── 04_parse_clue_results.py       ← Step 4: Analyze CLUE.io drug predictions
│   └── 05_visualize.py                ← Step 5: Generate all figures
│
├── notebooks/
│   └── OSCC_Drug_Repurposing_Pipeline.ipynb   ← Full interactive walkthrough
│
├── report/
│   ├── oscc_pipeline_report.Rmd       ← R Markdown source
│   └── oscc_pipeline_report.html      ← Rendered HTML report
│
├── results/
│   ├── figures/                       ← All generated plots
│   └── tables/                        ← Summary tables (CSV + Excel)
│
├── docs/
│   └── CLUE_IO_GUIDE.md               ← Step-by-step guide for the manual CLUE.io step
│
├── requirements.txt                   ← Python dependencies
├── environment.yml                    ← Conda environment file
└── run_pipeline.sh                    ← ONE-CLICK script to run everything
```

---

## ⚡ Quick Start

### 1. Clone the repository
```bash
git clone https://github.com/YOUR_USERNAME/oscc_lncrna_pipeline.git
cd oscc_lncrna_pipeline
```

### 2. Set up the environment
```bash
# Using conda (recommended)
conda env create -f environment.yml
conda activate oscc_pipeline

# OR using pip
pip install -r requirements.txt
```

### 3. Run the automated pipeline
```bash
bash run_pipeline.sh
```

### 4. Perform the CLUE.io query (2 minutes, manual step)
See [`docs/CLUE_IO_GUIDE.md`](docs/CLUE_IO_GUIDE.md) for a step-by-step walkthrough with screenshots.

### 5. Analyze drug candidates
After downloading CLUE.io results, run:
```bash
python scripts/04_parse_clue_results.py
python scripts/05_visualize.py
```

---

## 📊 Key Outputs

| Output | Description |
|--------|-------------|
| `target_genes.csv` | All protein-coding genes regulated by OSCC lncRNAs |
| `upregulated_genes.txt` | Input file for CLUE.io (UP signature) |
| `downregulated_genes.txt` | Input file for CLUE.io (DOWN signature) |
| `drug_candidates_ranked.csv` | Top drug candidates ranked by CMap connectivity score |
| `results/figures/` | Volcano plots, heatmaps, drug ranking bar charts |
| `report/oscc_pipeline_report.html` | Full narrative report with methods and results |

---

## 🧪 lncRNAs Analyzed

| lncRNA | Role in OSCC | Lnc2Cancer ID |
|--------|-------------|---------------|
| HOTAIR | Promotes invasion via EZH2/PRC2; silences CDH1 | LNC_000001 |
| MALAT1 | Enhances proliferation; upregulates VEGF | LNC_000088 |
| NEAT1 | Drives EMT; activates Wnt/β-catenin signaling | LNC_000102 |
| H19 | Oncogenic; regulates IGF2 imprinting | LNC_000045 |
| LINC00152 | Promotes metastasis; activates PI3K/AKT | LNC_000201 |
| PVT1 | MYC co-amplification; inhibits apoptosis | LNC_000178 |

---

## 🔄 Pipeline Diagram

```
┌─────────────────────────────────────────────────────────────┐
│                     OSCC PATIENT DATA                        │
│              (Lnc2Cancer 3.0 curated lncRNAs)               │
└──────────────────────────┬──────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────┐
│                  STEP 1: lncRNA Curation                     │
│         6 key OSCC lncRNAs selected + annotated             │
└──────────────────────────┬──────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────┐
│              STEP 2: RNAInter API Query                      │
│    Extract protein-coding gene targets for each lncRNA       │
│              Filter by interaction score ≥ 0.5              │
└──────────────────────────┬──────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────┐
│            STEP 3: Gene Signature Preparation                │
│     Classify as UP/DOWN regulated → Format for CLUE.io      │
└──────────────────────────┬──────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────┐
│              STEP 4: CLUE.io CMap Query  ← MANUAL STEP      │
│        Input: upregulated + downregulated gene lists         │
│        Output: Drug perturbagen connectivity scores          │
└──────────────────────────┬──────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────┐
│            STEP 5: Drug Candidate Analysis                   │
│    Rank drugs by score → Filter approved drugs → Visualize  │
└──────────────────────────┬──────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────┐
│                     FINAL REPORT                             │
│     Top N FDA-approved drug repurposing candidates for OSCC │
└─────────────────────────────────────────────────────────────┘
```

---

## 📚 References

1. Lnc2Cancer v3.0: https://www.bio-bigdata.net/lnc2cancer/
2. RNAInter: http://www.rnainter.org/
3. CLUE.io Connectivity Map: https://clue.io/
4. Subramanian et al. (2017). A Next Generation Connectivity Map. *Cell*, 171(6), 1437–1452.
5. Liu et al. (2021). lncRNA HOTAIR promotes OSCC invasion. *Oral Oncology*, 112, 105-112.

---

## 👤 Author

**Vaibhav** | MSc Zoology (Entomology & Wildlife) | Fergusson College, Pune  
Bioinformatics self-study project | Marathwada, Maharashtra, India

---

## 📄 License

MIT License — free to use, cite, and build upon.
