"""
Script 01: OSCC lncRNA Curation
================================
PURPOSE:
    Builds a curated list of the most well-documented lncRNAs
    in Oral Squamous Cell Carcinoma (OSCC), with their known
    biological roles and regulatory direction.

    In a full study, you would download this from Lnc2Cancer 3.0:
    https://www.bio-bigdata.net/lnc2cancer/
    (Search: "oral squamous cell carcinoma" → Export CSV)

    Here we provide a pre-curated list based on published literature,
    which is the standard approach for focused repurposing studies.

HOW TO RUN:
    python scripts/01_fetch_lncrnas.py

OUTPUT:
    data/raw/oscc_lncrnas.csv
"""

import pandas as pd
import os

# ── CREATE OUTPUT DIRECTORY ──────────────────────────────────────────────────
os.makedirs("data/raw", exist_ok=True)

# ── CURATED OSCC lncRNA DATA ─────────────────────────────────────────────────
# Source: Lnc2Cancer 3.0 + literature curation
# Fields: lncRNA name, Ensembl ID, expression direction, key function, PMID

oscc_lncrnas = [
    {
        "lncrna_name":   "HOTAIR",
        "ensembl_id":    "ENSG00000228630",
        "expression":    "Upregulated",
        "key_pathway":   "EZH2/PRC2 epigenetic silencing, EMT promotion",
        "regulated_direction": "Silences tumor suppressors (CDH1, PTEN)",
        "pmid":          "27531065",
        "cancer_stage":  "Advanced OSCC, lymph node metastasis",
    },
    {
        "lncrna_name":   "MALAT1",
        "ensembl_id":    "ENSG00000251562",
        "expression":    "Upregulated",
        "key_pathway":   "Alternative splicing, VEGF upregulation, Wnt",
        "regulated_direction": "Promotes angiogenesis and proliferation",
        "pmid":          "28821566",
        "cancer_stage":  "All stages, poor prognosis marker",
    },
    {
        "lncrna_name":   "NEAT1",
        "ensembl_id":    "ENSG00000245532",
        "expression":    "Upregulated",
        "key_pathway":   "Paraspeckle formation, PI3K/AKT activation",
        "regulated_direction": "Drives EMT, inhibits apoptosis",
        "pmid":          "29449119",
        "cancer_stage":  "Late stage, chemoresistance",
    },
    {
        "lncrna_name":   "H19",
        "ensembl_id":    "ENSG00000130600",
        "expression":    "Upregulated",
        "key_pathway":   "IGF2/IGF1R signaling, miR-675 sponge",
        "regulated_direction": "Promotes cell growth via IGF pathway",
        "pmid":          "30064698",
        "cancer_stage":  "Early to mid-stage OSCC",
    },
    {
        "lncrna_name":   "LINC00152",
        "ensembl_id":    "ENSG00000227036",
        "expression":    "Upregulated",
        "key_pathway":   "PI3K/AKT/mTOR, EGFR signaling",
        "regulated_direction": "Promotes invasion and metastasis",
        "pmid":          "30940543",
        "cancer_stage":  "Metastatic OSCC",
    },
    {
        "lncrna_name":   "PVT1",
        "ensembl_id":    "ENSG00000249859",
        "expression":    "Upregulated",
        "key_pathway":   "MYC amplification, BCL2 upregulation",
        "regulated_direction": "Blocks apoptosis, drives proliferation",
        "pmid":          "31271680",
        "cancer_stage":  "Late stage, radioresistance",
    },
    {
        "lncrna_name":   "TUG1",
        "ensembl_id":    "ENSG00000253352",
        "expression":    "Upregulated",
        "key_pathway":   "EZH2 interaction, p53 pathway suppression",
        "regulated_direction": "Silences p21, promotes cell cycle entry",
        "pmid":          "29444205",
        "cancer_stage":  "Oral carcinogenesis, early invasion",
    },
    {
        "lncrna_name":   "MEG3",
        "ensembl_id":    "ENSG00000214548",
        "expression":    "Downregulated",
        "key_pathway":   "p53 tumor suppressor activation",
        "regulated_direction": "Tumor suppressor — lost in OSCC",
        "pmid":          "28416576",
        "cancer_stage":  "All stages, loss correlates with poor prognosis",
    },
    {
        "lncrna_name":   "GAS5",
        "ensembl_id":    "ENSG00000234741",
        "expression":    "Downregulated",
        "key_pathway":   "mTOR/GR signaling inhibition, apoptosis induction",
        "regulated_direction": "Growth arrest signal — silenced in OSCC",
        "pmid":          "27009926",
        "cancer_stage":  "Late stage, correlates with chemo-resistance",
    },
]

# ── SAVE TO CSV ───────────────────────────────────────────────────────────────
df = pd.DataFrame(oscc_lncrnas)
output_path = "data/raw/oscc_lncrnas.csv"
df.to_csv(output_path, index=False)

# ── PRINT SUMMARY ─────────────────────────────────────────────────────────────
print("=" * 60)
print("  STEP 1 COMPLETE: OSCC lncRNA Curation")
print("=" * 60)
print(f"\n  Total lncRNAs curated : {len(df)}")
print(f"  Upregulated           : {len(df[df.expression == 'Upregulated'])}")
print(f"  Downregulated         : {len(df[df.expression == 'Downregulated'])}")
print(f"\n  Output saved to       : {output_path}")
print("\n  lncRNAs in dataset:")
for _, row in df.iterrows():
    arrow = "⬆" if row["expression"] == "Upregulated" else "⬇"
    print(f"    {arrow}  {row['lncrna_name']:15s}  ({row['key_pathway'][:45]}...)")
print("\n  Next step: Run  python scripts/02_query_rnainter.py")
print("=" * 60)
