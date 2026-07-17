#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calcola la copertura per-gene in bp e in percentuale per le principali categorie
usando DUE file di input:
- uno per regulatory
- uno per esoni/introni

Scrive SOLO un file riassuntivo con una riga per gene.

Note:
- usa la lunghezza del gene puro: gene_end - gene_start + 1
- per ogni categoria calcola la copertura UNICA in bp (unione degli intervalli),
  quindi senza doppio conteggio di overlap interni alla stessa categoria
- categorie diverse possono comunque sovrapporsi tra loro, quindi le percentuali
  non devono sommare a 100
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Dict, Iterable, List, Tuple
import pandas as pd
import numpy as np


# =========================
# CONFIG
# =========================
REGULATORY_CSV = Path(r"C:/Users/user/OneDrive - Politecnico di Torino/Desktop/verifica codici/ANALISI REGULATORY FEATURES/coordinate_REGULATORY_complete.csv")
GENIC_CSV = Path(r"C:/Users/user/OneDrive - Politecnico di Torino/Desktop/verifica codici/1_Coordinate_esoni_introni.csv")
OUTPUT_SUMMARY_CSV = Path(r"C:/Users/user/OneDrive - Politecnico di Torino/Desktop/verifica codici/gene_coverage_summary_two_inputs.csv")
ROUND_PCT = 4

REQUIRED_COLS = [
    "gene_symbol", "gene_start", "gene_end",
    "feature_start", "feature_end",
    "subset", "feature_type", "feature_label", "logic_name",
]


# =========================
# HELPERS
# =========================
def norm_text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip().str.lower()


def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Unisce intervalli genomici 1-based inclusivi."""
    cleaned = []
    for start, end in intervals:
        if pd.isna(start) or pd.isna(end):
            continue
        start = int(start)
        end = int(end)
        if end < start:
            start, end = end, start
        cleaned.append((start, end))

    if not cleaned:
        return []

    cleaned.sort(key=lambda x: (x[0], x[1]))
    merged = [cleaned[0]]

    for start, end in cleaned[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + 1:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))

    return merged


def covered_bp(intervals: Iterable[Tuple[int, int]]) -> int:
    merged = merge_intervals(intervals)
    return int(sum(end - start + 1 for start, end in merged))


def load_and_merge_inputs(regulatory_csv: Path, genic_csv: Path) -> pd.DataFrame:
    if not regulatory_csv.exists():
        raise FileNotFoundError(f"File non trovato: {regulatory_csv}")
    if not genic_csv.exists():
        raise FileNotFoundError(f"File non trovato: {genic_csv}")

    df_reg = pd.read_csv(regulatory_csv)
    df_gen = pd.read_csv(genic_csv)

    for col in REQUIRED_COLS:
        if col not in df_reg.columns:
            df_reg[col] = pd.NA
        if col not in df_gen.columns:
            df_gen[col] = pd.NA

    df = pd.concat([
        df_reg[REQUIRED_COLS].copy(),
        df_gen[REQUIRED_COLS].copy(),
    ], ignore_index=True)

    df = df.drop_duplicates().reset_index(drop=True)
    return df


# =========================
# MAIN
# =========================
def main() -> None:
    df = load_and_merge_inputs(REGULATORY_CSV, GENIC_CSV)

    subset = norm_text(df["subset"])
    feature_type = norm_text(df["feature_type"])
    feature_label = norm_text(df["feature_label"])
    logic_name = norm_text(df["logic_name"])

    category_rules: Dict[str, Callable[[pd.DataFrame], pd.Series]] = {
        "enhancer": lambda x: feature_label.eq("enhancer"),
        "promoter": lambda x: feature_label.eq("promoter"),
        "ctcf": lambda x: feature_label.eq("ctcf"),
        "emar_only": lambda x: feature_label.eq("emar_only") | logic_name.eq("emar_only"),
        "emar_and_other": lambda x: feature_label.eq("emar_and_other") | logic_name.eq("emar_and_other"),
        "other_genic": lambda x: (
            subset.eq("other_genic") |
            feature_type.eq("other_genic") |
            feature_label.eq("other_genic") |
            logic_name.eq("other_genic")
        ),
        "regulatory_all": lambda x: (
            subset.eq("regulatory") |
            feature_type.eq("regulatory") |
            logic_name.eq("regulatory")
        ),
        "exon_mane": lambda x: subset.eq("exon_mane") | logic_name.eq("exon_mane"),
        "intron_mane": lambda x: subset.eq("intron_mane") | logic_name.eq("intron_mane"),
        "exon_union": lambda x: subset.eq("exon_union") | logic_name.eq("exon_union"),
        "intron_union": lambda x: subset.eq("intron_union") | logic_name.eq("intron_union"),
    }

    summary_rows = []

    for gene, gdf in df.groupby("gene_symbol", dropna=False):
        gene_start = int(gdf["gene_start"].min())
        gene_end = int(gdf["gene_end"].max())
        gene_length_bp = gene_end - gene_start + 1

        row = {
            "gene_symbol": gene,
            "gene_start": gene_start,
            "gene_end": gene_end,
            "gene_length_bp": gene_length_bp,
        }

        for category_name, rule in category_rules.items():
            mask = rule(gdf)
            cat_df = gdf.loc[mask, ["feature_start", "feature_end"]].copy()
            bp = covered_bp(cat_df.itertuples(index=False, name=None))
            pct = (bp / gene_length_bp * 100.0) if gene_length_bp > 0 else np.nan
            row[f"{category_name}_bp"] = bp
            row[f"{category_name}_pct"] = round(pct, ROUND_PCT)

        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows).sort_values("gene_symbol").reset_index(drop=True)
    summary_df.to_csv(OUTPUT_SUMMARY_CSV, index=False)

    print(f"Creato: {OUTPUT_SUMMARY_CSV}")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
