"""
Script 03: Gene Signature Preparation for CLUE.io
==================================================
PURPOSE:
    Takes the target gene interaction table from Step 2 and
    prepares it in the exact format required by CLUE.io's
    Connectivity Map query interface.

    CLUE.io needs:
      1. A list of UPREGULATED genes   (the "disease UP signature")
      2. A list of DOWNREGULATED genes (the "disease DOWN signature")

    These represent the gene expression state OF THE CANCER CELL
    that we want to reverse with a drug.

    Logic:
      - If a gene is NET UPREGULATED in OSCC → goes into UP list
      - If a gene is NET DOWNREGULATED in OSCC → goes into DOWN list
      - CLUE.io will find drugs whose perturbation OPPOSES this pattern

HOW TO RUN:
    python scripts/03_prepare_gene_signature.py

OUTPUT:
    data/processed/upregulated_genes.txt    ← Paste into CLUE.io
    data/processed/downregulated_genes.txt  ← Paste into CLUE.io
    data/processed/gene_signature_summary.csv
"""

import pandas as pd
import os

os.makedirs("data/processed", exist_ok=True)

print("=" * 60)
print("  STEP 3: Gene Signature Preparation")
print("=" * 60)

# ── LOAD TARGET GENES ─────────────────────────────────────────────────────────
df = pd.read_csv("data/processed/target_genes.csv")

# ── RESOLVE CONFLICTS (same gene appearing in both directions) ────────────────
# If a gene is regulated by multiple lncRNAs in different directions,
# use the majority vote and the mean score

gene_summary = (
    df.groupby(["target_gene", "net_effect_in_OSCC"])
    .agg(
        n_lncrnas=("lncrna_name", "nunique"),
        mean_score=("interaction_score", "mean"),
        supporting_lncrnas=("lncrna_name", lambda x: ", ".join(sorted(set(x))))
    )
    .reset_index()
)

# For genes with conflicting directions, keep the one with higher score sum
gene_direction = (
    gene_summary
    .sort_values("mean_score", ascending=False)
    .drop_duplicates("target_gene", keep="first")
    .copy()
)

# ── SEPARATE UP vs DOWN ───────────────────────────────────────────────────────
up_genes   = gene_direction[gene_direction.net_effect_in_OSCC == "up"].sort_values("mean_score", ascending=False)
down_genes = gene_direction[gene_direction.net_effect_in_OSCC == "down"].sort_values("mean_score", ascending=False)

# ── CLUE.io RECOMMENDS: 50–150 genes per list for optimal results ─────────────
# We'll keep high-confidence genes (score ≥ 0.75)
up_filtered   = up_genes[up_genes.mean_score >= 0.75]
down_filtered = down_genes[down_genes.mean_score >= 0.75]

# ── SAVE PLAIN TEXT GENE LISTS (exactly what you paste into CLUE.io) ─────────
with open("data/processed/upregulated_genes.txt", "w") as f:
    for gene in up_filtered["target_gene"]:
        f.write(gene + "\n")

with open("data/processed/downregulated_genes.txt", "w") as f:
    for gene in down_filtered["target_gene"]:
        f.write(gene + "\n")

# ── SAVE SUMMARY TABLE ────────────────────────────────────────────────────────
gene_direction.to_csv("data/processed/gene_signature_summary.csv", index=False)

# ── PRINT SUMMARY ─────────────────────────────────────────────────────────────
print(f"\n  Gene Signature Summary:")
print(f"  ┌─────────────────────────────────────────────┐")
print(f"  │  Total unique target genes  : {gene_direction['target_gene'].nunique():3d}           │")
print(f"  │  Upregulated in OSCC        : {len(up_filtered):3d}           │")
print(f"  │  Downregulated in OSCC      : {len(down_filtered):3d}           │")
print(f"  └─────────────────────────────────────────────┘")

print(f"\n  ⬆ TOP UPREGULATED GENES (OSCC disease signature):")
for _, row in up_filtered.head(10).iterrows():
    print(f"      {row['target_gene']:10s}  score={row['mean_score']:.2f}  driven by: {row['supporting_lncrnas']}")

print(f"\n  ⬇ TOP DOWNREGULATED GENES (tumor suppressor losses):")
for _, row in down_filtered.head(10).iterrows():
    print(f"      {row['target_gene']:10s}  score={row['mean_score']:.2f}  driven by: {row['supporting_lncrnas']}")

print(f"""
  ✅ Files ready for CLUE.io:
      data/processed/upregulated_genes.txt
      data/processed/downregulated_genes.txt

  ────────────────────────────────────────────────
  📋 NEXT STEP — CLUE.io QUERY (2 minutes):
  ────────────────────────────────────────────────
  1. Go to: https://clue.io/query
  2. Create a free account (required)
  3. Select: "Gene Expression" query type
  4. Paste contents of upregulated_genes.txt into UP box
  5. Paste contents of downregulated_genes.txt into DOWN box
  6. Select cell line: "A375" (melanoma) or "PC3" (prostate)
     — these have the best CMap coverage
  7. Click "Submit Query"
  8. Wait ~10–30 minutes for results
  9. Download results as CSV
  10. Save to: data/clue_output/clue_results.csv
  11. Then run: python scripts/04_parse_clue_results.py
  ────────────────────────────────────────────────
  See docs/CLUE_IO_GUIDE.md for screenshots!
""")
print("=" * 60)
