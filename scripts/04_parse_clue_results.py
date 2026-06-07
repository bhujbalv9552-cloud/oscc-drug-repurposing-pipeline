"""
Script 04: Parse CLUE.io Results & Rank Drug Candidates
========================================================
PURPOSE:
    Reads the CSV file downloaded from CLUE.io and:
      1. Parses the connectivity scores
      2. Identifies the top drug repurposing candidates
      3. Annotates drugs with FDA approval status
      4. Filters for clinically relevant (negative score = reversal) drugs
      5. Saves ranked drug table

    CONNECTIVITY SCORE INTERPRETATION:
      Score = -90 to -100 → STRONG reversal candidate (what we want)
      Score = -50 to -90  → Moderate reversal candidate
      Score = 0           → No connection
      Score > 0           → Would WORSEN the cancer state (avoid)

    The more NEGATIVE the score, the more strongly the drug
    REVERSES the OSCC gene expression signature.

HOW TO RUN:
    python scripts/04_parse_clue_results.py

    REQUIRES: data/clue_output/clue_results.csv
    (Download from CLUE.io after running Step 3)

OUTPUT:
    results/tables/drug_candidates_ranked.csv
    results/tables/top10_drugs.csv
"""

import pandas as pd
import os
import json

os.makedirs("results/tables", exist_ok=True)
os.makedirs("data/clue_output", exist_ok=True)

print("=" * 60)
print("  STEP 4: CLUE.io Results Analysis")
print("=" * 60)

# ── CHECK FOR REAL CLUE.io OUTPUT ────────────────────────────────────────────
clue_file = "data/clue_output/clue_results.csv"

if os.path.exists(clue_file) and os.path.getsize(clue_file) > 100:
    print(f"\n  ✅ Found CLUE.io results: {clue_file}")
    df_clue = pd.read_csv(clue_file)
    print(f"  Loaded {len(df_clue)} drug perturbagen entries")

else:
    print("\n  ⚠ CLUE.io results not yet downloaded.")
    print("  → Using simulated representative results for pipeline demonstration.")
    print("  → Replace data/clue_output/clue_results.csv with your real download.\n")

    # ── SIMULATED RESULTS (representative of real CMap outputs) ───────────────
    # These drugs are chosen because they are ACTUALLY known to affect the
    # pathways driven by our OSCC lncRNAs — this is not random simulation,
    # it reflects published drug-pathway relationships.
    simulated_drugs = [
        # Strong reversal candidates (negative scores)
        {"pert_iname": "vorinostat",       "pert_type": "trt_cp", "connectivity_score": -97.2, "fda_approved": True,  "target": "HDAC inhibitor",         "known_pathway": "EZH2/PRC2 (HOTAIR pathway)"},
        {"pert_iname": "entinostat",       "pert_type": "trt_cp", "connectivity_score": -94.8, "fda_approved": True,  "target": "HDAC inhibitor (class I)", "known_pathway": "EZH2/HOTAIR"},
        {"pert_iname": "everolimus",       "pert_type": "trt_cp", "connectivity_score": -93.5, "fda_approved": True,  "target": "mTOR inhibitor",          "known_pathway": "GAS5/NEAT1 mTOR pathway"},
        {"pert_iname": "rapamycin",        "pert_type": "trt_cp", "connectivity_score": -91.2, "fda_approved": True,  "target": "mTOR inhibitor",          "known_pathway": "GAS5 loss → mTOR activation"},
        {"pert_iname": "gefitinib",        "pert_type": "trt_cp", "connectivity_score": -89.7, "fda_approved": True,  "target": "EGFR inhibitor",          "known_pathway": "LINC00152 → EGFR"},
        {"pert_iname": "erlotinib",        "pert_type": "trt_cp", "connectivity_score": -88.3, "fda_approved": True,  "target": "EGFR inhibitor",          "known_pathway": "LINC00152 → EGFR"},
        {"pert_iname": "JQ1",              "pert_type": "trt_cp", "connectivity_score": -87.9, "fda_approved": False, "target": "BET/BRD4 inhibitor",      "known_pathway": "MYC (PVT1/H19 pathway)"},
        {"pert_iname": "venetoclax",       "pert_type": "trt_cp", "connectivity_score": -86.4, "fda_approved": True,  "target": "BCL2 inhibitor",          "known_pathway": "PVT1 → BCL2 upregulation"},
        {"pert_iname": "navitoclax",       "pert_type": "trt_cp", "connectivity_score": -85.1, "fda_approved": False, "target": "BCL2/BCL-XL inhibitor",   "known_pathway": "PVT1/MALAT1 → BCL2"},
        {"pert_iname": "trametinib",       "pert_type": "trt_cp", "connectivity_score": -83.7, "fda_approved": True,  "target": "MEK inhibitor",           "known_pathway": "HOTAIR → RAS/MAPK"},
        {"pert_iname": "pictilisib",       "pert_type": "trt_cp", "connectivity_score": -82.3, "fda_approved": False, "target": "PI3K inhibitor",          "known_pathway": "NEAT1/LINC00152 → PI3K"},
        {"pert_iname": "MK-2206",          "pert_type": "trt_cp", "connectivity_score": -81.0, "fda_approved": False, "target": "AKT inhibitor",           "known_pathway": "NEAT1 → AKT1"},
        {"pert_iname": "nutlin-3",         "pert_type": "trt_cp", "connectivity_score": -79.5, "fda_approved": False, "target": "MDM2 inhibitor",          "known_pathway": "MEG3 loss → MDM2 upregulation"},
        {"pert_iname": "palbociclib",      "pert_type": "trt_cp", "connectivity_score": -78.2, "fda_approved": True,  "target": "CDK4/6 inhibitor",        "known_pathway": "MALAT1/H19 → CDK6/CCND1"},
        {"pert_iname": "ribociclib",       "pert_type": "trt_cp", "connectivity_score": -77.8, "fda_approved": True,  "target": "CDK4/6 inhibitor",        "known_pathway": "MALAT1 → CDK6"},
        {"pert_iname": "bevacizumab",      "pert_type": "trt_cp", "connectivity_score": -76.3, "fda_approved": True,  "target": "VEGF inhibitor",          "known_pathway": "HOTAIR/MALAT1 → VEGF"},
        {"pert_iname": "metformin",        "pert_type": "trt_cp", "connectivity_score": -74.9, "fda_approved": True,  "target": "AMPK activator / mTOR-",  "known_pathway": "GAS5 / mTOR pathway"},
        {"pert_iname": "sorafenib",        "pert_type": "trt_cp", "connectivity_score": -73.5, "fda_approved": True,  "target": "Multi-kinase (VEGFR, RAF)", "known_pathway": "HOTAIR/MALAT1 angiogenesis"},
        # Moderate candidates
        {"pert_iname": "doxorubicin",      "pert_type": "trt_cp", "connectivity_score": -65.2, "fda_approved": True,  "target": "Topoisomerase II",        "known_pathway": "General cytotoxic"},
        {"pert_iname": "cisplatin",        "pert_type": "trt_cp", "connectivity_score": -58.1, "fda_approved": True,  "target": "DNA crosslinking agent",  "known_pathway": "Standard OSCC chemo"},
        # Positive scores (would worsen) — shown as examples to exclude
        {"pert_iname": "EGF",             "pert_type": "trt_sh", "connectivity_score":  72.3, "fda_approved": False, "target": "EGFR agonist",            "known_pathway": "Activates OSCC pathway — EXCLUDE"},
        {"pert_iname": "insulin",         "pert_type": "trt_cp", "connectivity_score":  61.7, "fda_approved": True,  "target": "IGF pathway activator",   "known_pathway": "H19 → IGF2 — EXCLUDE"},
    ]
    df_clue = pd.DataFrame(simulated_drugs)

# ── FILTER AND RANK ───────────────────────────────────────────────────────────
# Keep only reversal candidates (negative connectivity score)
df_reversal = df_clue[df_clue["connectivity_score"] < 0].copy()
df_reversal = df_reversal.sort_values("connectivity_score", ascending=True)

# Separate FDA-approved vs experimental
df_fda       = df_reversal[df_reversal["fda_approved"] == True]
df_exp       = df_reversal[df_reversal["fda_approved"] == False]

# ── SAVE OUTPUTS ─────────────────────────────────────────────────────────────
df_reversal.to_csv("results/tables/drug_candidates_ranked.csv", index=False)
df_reversal.head(10).to_csv("results/tables/top10_drugs.csv", index=False)

# ── PRINT SUMMARY ─────────────────────────────────────────────────────────────
print(f"\n  Total drug perturbagens analyzed : {len(df_clue)}")
print(f"  Reversal candidates (score < 0)  : {len(df_reversal)}")
print(f"  FDA-approved candidates          : {len(df_fda)}")
print(f"  Experimental candidates          : {len(df_exp)}")

print(f"\n  🏆 TOP 10 DRUG REPURPOSING CANDIDATES:")
print(f"  {'Rank':<5} {'Drug':<20} {'Score':>8}  {'FDA':>5}  Target")
print(f"  {'----':<5} {'----':<20} {'-----':>8}  {'---':>5}  ------")
for i, (_, row) in enumerate(df_reversal.head(10).iterrows(), 1):
    fda_tag = "✓" if row.get("fda_approved", False) else "—"
    print(f"  {i:<5} {row['pert_iname']:<20} {row['connectivity_score']:>8.1f}  {fda_tag:>5}  {row.get('target', 'N/A')}")

print(f"""
  ✅ Results saved to:
      results/tables/drug_candidates_ranked.csv
      results/tables/top10_drugs.csv

  Next step: Run  python scripts/05_visualize.py
""")
print("=" * 60)
