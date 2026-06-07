#!/bin/bash
# ═══════════════════════════════════════════════════════════════════════════════
#  OSCC lncRNA Drug Repurposing Pipeline — One-Click Runner
#  Run this from the project root: bash run_pipeline.sh
# ═══════════════════════════════════════════════════════════════════════════════

echo ""
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║   OSCC lncRNA → Drug Repurposing Pipeline                   ║"
echo "║   Author: Vaibhav | MSc Zoology | Fergusson College, Pune   ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo ""

# Check Python
if ! command -v python &> /dev/null; then
    echo "❌ Python not found. Please install Python 3.8+ and try again."
    exit 1
fi

# Install dependencies if needed
echo "📦 Checking dependencies..."
pip install -q pandas matplotlib requests openpyxl

echo ""
echo "▶ STEP 1/5 — Curating OSCC lncRNA list..."
python scripts/01_fetch_lncrnas.py
echo ""

echo "▶ STEP 2/5 — Querying target genes (RNAInter)..."
python scripts/02_query_rnainter.py
echo ""

echo "▶ STEP 3/5 — Preparing gene signature for CLUE.io..."
python scripts/03_prepare_gene_signature.py
echo ""

echo "════════════════════════════════════════════════════════════"
echo "  ⏸ MANUAL STEP REQUIRED — CLUE.io Query"
echo "  See: docs/CLUE_IO_GUIDE.md for full instructions"
echo ""
echo "  1. Visit https://clue.io/query"
echo "  2. Paste contents of data/processed/upregulated_genes.txt"
echo "  3. Paste contents of data/processed/downregulated_genes.txt"
echo "  4. Submit and download results"
echo "  5. Save to: data/clue_output/clue_results.csv"
echo "════════════════════════════════════════════════════════════"
echo ""
read -p "  Press ENTER when you have downloaded the CLUE.io results... "
echo ""

echo "▶ STEP 4/5 — Parsing drug candidates..."
python scripts/04_parse_clue_results.py
echo ""

echo "▶ STEP 5/5 — Generating figures..."
python scripts/05_visualize.py
echo ""

echo "╔══════════════════════════════════════════════════════════════╗"
echo "║   ✅ PIPELINE COMPLETE                                       ║"
echo "║                                                              ║"
echo "║   Results:  results/tables/drug_candidates_ranked.csv        ║"
echo "║   Figures:  results/figures/                                 ║"
echo "║   Report:   report/oscc_pipeline_report.html (run R Markdown)║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo ""
