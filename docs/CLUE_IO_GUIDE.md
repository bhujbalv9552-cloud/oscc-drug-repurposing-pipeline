# 📋 CLUE.io Step-by-Step Guide
## The Manual Query Step (takes ~5 minutes)

This guide walks you through the only step that requires a web browser.
Everything else in the pipeline runs automatically.

---

## Why CLUE.io?

CLUE.io hosts the **Connectivity Map (CMap)** — a database built by the Broad Institute (MIT/Harvard) containing the gene expression profiles of ~5,000 drugs and genetic perturbations across multiple cancer cell lines.

By comparing your "disease signature" (what genes OSCC turns ON and OFF) to the drug profiles, CMap tells you which drug would REVERSE your disease state — like finding the key that fits your lock.

This is **real translational medicine** — the same method used in published Nature Medicine papers.

---

## Before You Start

Make sure Step 3 of the pipeline has run. You need these two files:
- `data/processed/upregulated_genes.txt`
- `data/processed/downregulated_genes.txt`

Open both files in Notepad (Windows) or TextEdit (Mac) and keep them handy.

---

## Step-by-Step Instructions

### Step 1: Create a Free Account
1. Open your browser and go to: **https://clue.io**
2. Click **"Sign Up"** in the top right
3. Register with your email address (free, no credit card)
4. Verify your email and log in

---

### Step 2: Go to the Query Page
1. After logging in, click **"Query"** in the top navigation bar
   OR go directly to: **https://clue.io/query**
2. You will see the "Query the CMap" interface

---

### Step 3: Set Up Your Query

**Query Name:** Give it a name like `OSCC_lncRNA_repurposing`

**Query Type:** Select **"Gene Expression"**

**Cell Line:** 
- Select **"A375"** (melanoma line) — this has the best CMap coverage
- Alternatively: **"MCF7"** (breast) or **"PC3"** (prostate)
- NOTE: There is no oral cancer line in CMap, but A375 is the standard choice for cross-cancer repurposing studies. This is scientifically valid and you must mention this in your methods section.

---

### Step 4: Paste Your Gene Lists

**UP Genes (Upregulated in disease):**
1. Open `data/processed/upregulated_genes.txt`
2. Select all text (Ctrl+A) and copy (Ctrl+C)
3. Paste into the **"UP"** text box on CLUE.io

**DOWN Genes (Downregulated in disease):**
1. Open `data/processed/downregulated_genes.txt`
2. Select all text (Ctrl+A) and copy (Ctrl+C)
3. Paste into the **"DOWN"** text box on CLUE.io

---

### Step 5: Submit the Query
1. Review your inputs — both UP and DOWN boxes should have genes
2. Click the **"Submit Query"** button
3. You will receive an email notification when results are ready
4. **Wait time:** Usually 10–45 minutes depending on server load

---

### Step 6: Download Your Results
1. When notified, return to https://clue.io/query
2. Click on your completed query under **"My Queries"**
3. You will see a ranked list of drug perturbagens
4. Click **"Download"** → select **CSV format**
5. Rename the downloaded file to: `clue_results.csv`
6. Move it to: `data/clue_output/clue_results.csv`

---

### Step 7: Continue the Pipeline
After saving the file, return to your terminal and:

```bash
python scripts/04_parse_clue_results.py
python scripts/05_visualize.py
```

OR simply run: `bash run_pipeline.sh` (it will pause and wait for you at this step)

---

## How to Interpret the Results

| Score Range | Meaning |
|-------------|---------|
| −90 to −100 | **Excellent** reversal candidate — drug strongly opposes OSCC signature |
| −75 to −90  | **Strong** candidate — investigate further |
| −50 to −75  | **Moderate** — worth noting |
| −50 to 0    | Weak — unlikely to be useful |
| 0 to +100   | **AVOID** — drug would WORSEN the cancer state |

Focus your analysis on **FDA-approved drugs** with scores below −75.
These are the most clinically actionable findings.

---

## What to Write in Your Methods Section

> "Gene expression connectivity analysis was performed using the Broad Institute's Connectivity Map (CMap) CLUE.io platform (Subramanian et al., 2017). The disease gene signature, comprising upregulated and downregulated genes derived from lncRNA–target interactions in OSCC, was queried against the Level 5 CMap compound perturbagen dataset using the A375 cell line. Drugs with a normalized connectivity score below −75 were considered candidate repurposing agents."

---

## Troubleshooting

**"No results found for your gene list"**  
→ Some genes may not be in the CMap reference. Reduce list to top 50 genes.

**"Query failed"**  
→ CMap uses HUGO gene symbols. Make sure genes like "TP53" not "p53"

**Results seem random**  
→ Check that UP and DOWN lists are not swapped

---

## Reference to Cite

Subramanian, A., et al. (2017). A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. *Cell*, 171(6), 1437–1452.e17. https://doi.org/10.1016/j.cell.2017.10.049
