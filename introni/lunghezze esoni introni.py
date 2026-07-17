#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from pathlib import Path

# ============================================================
# CONFIG
# ============================================================

BASE_DIR = Path(r"C:\Users\Chiara Panico\OneDrive - Politecnico di Torino\Desktop\verifica codici")

EXON_INTRON_FILE = BASE_DIR / "1_Coordinate_esoni_introni.csv"
REGULATORY_FILE  = BASE_DIR / "coordinate_REGULATORY_complete.csv"

OUT_FILE = BASE_DIR / "feature_lengths_for_rho_table.csv"

GENES = ["ASH1L", "NSD1", "NSD2", "NSD3", "SETD2"]

# ============================================================
# HELPERS
# ============================================================

def clean_columns(df):
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    return df


def interval_length(df):
    """
    Coordinate genomiche 1-based inclusive:
    length = feature_end - feature_start + 1
    """
    return pd.to_numeric(df["feature_end"]) - pd.to_numeric(df["feature_start"]) + 1


def sum_by_subset(df, subset_name):
    sub = df[df["subset"] == subset_name].copy()
    if sub.empty:
        return {}
    sub["length_bp"] = interval_length(sub)
    return sub.groupby("gene_symbol")["length_bp"].sum().to_dict()


def sum_by_label(df, label_name):
    sub = df[df["feature_label"] == label_name].copy()
    if sub.empty:
        return {}
    sub["length_bp"] = interval_length(sub)
    return sub.groupby("gene_symbol")["length_bp"].sum().to_dict()


def gene_lengths(df):
    """
    Gene length = gene_end - gene_start + 1
    """
    tmp = df[["gene_symbol", "gene_start", "gene_end"]].drop_duplicates().copy()
    tmp["Gene_length"] = (
        pd.to_numeric(tmp["gene_end"]) -
        pd.to_numeric(tmp["gene_start"]) +
        1
    )
    return tmp.groupby("gene_symbol")["Gene_length"].first().to_dict()


# ============================================================
# MAIN
# ============================================================

df_exint = clean_columns(pd.read_csv(EXON_INTRON_FILE))
df_reg   = clean_columns(pd.read_csv(REGULATORY_FILE))

rows = []

# Dizionari lunghezze
gene_len = gene_lengths(df_exint)

exon_mane    = sum_by_subset(df_exint, "EXON_MANE")
intron_mane  = sum_by_subset(df_exint, "INTRON_MANE")
exon_union   = sum_by_subset(df_exint, "EXON_UNION")
intron_union = sum_by_subset(df_exint, "INTRON_UNION")

promoter       = sum_by_label(df_reg, "Promoter")
enhancer       = sum_by_label(df_reg, "Enhancer")
ctcf           = sum_by_label(df_reg, "CTCF")
emar_only      = sum_by_label(df_reg, "EMAR_ONLY")
emar_and_other = sum_by_label(df_reg, "EMAR_AND_OTHER")
other_genic    = sum_by_label(df_reg, "OTHER_GENIC")

for gene in GENES:
    rows.append({
        "Gene": gene,
        "Gene_length": gene_len.get(gene, 0),

        "Exon_MANE_length": exon_mane.get(gene, 0),
        "Exon_Union_length": exon_union.get(gene, 0),
        "Intron_MANE_length": intron_mane.get(gene, 0),
        "Intron_Union_length": intron_union.get(gene, 0),

        "Promoter_length": promoter.get(gene, 0),
        "Enhancer_length": enhancer.get(gene, 0),
        "CTCF_length": ctcf.get(gene, 0),
        "EMAR_only_length": emar_only.get(gene, 0),
        "EMAR_and_Other_length": emar_and_other.get(gene, 0),
        "Other_genic_length": other_genic.get(gene, 0),
    })

out = pd.DataFrame(rows)
out.to_csv(OUT_FILE, index=False)

print(out.to_string(index=False))
print(f"\nFile salvato in:\n{OUT_FILE}")