#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Classifica le righe EMAR in:
- emar_only
- emar_and_other

partendo da file per gene con colonne tipo:
feature_type ; ID ; feature_start ; feature_end ; info

Formato atteso:
- separatore: ;
- colonne minime:
    feature_type
    feature_start
    feature_end

Output:
- un nuovo file per ogni gene, con feature_type aggiornato
- colonna aggiuntiva: emar_overlap_with
"""

import os
import glob
from typing import Optional, Set

import pandas as pd


# ===================== CONFIG =====================

INPUT_DIR = r"C:\Users\user\OneDrive - Politecnico di Torino\Desktop\verifica codici\coordinate"
OUTPUT_DIR = r"C:\Users\user\OneDrive - Politecnico di Torino\Desktop\verifica codici\coordinate_per_gene_emar_classified"

FILE_PATTERN = "*.csv"
SEP = ";"

# categorie che rendono una EMAR = emar_and_other
OTHER_FEATURES = ("promoter", "enhancer", "ctcf")


# ===================== HELPERS =====================

def normalize_text(x) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip().lower()


def clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    df.columns = [str(c).strip() for c in df.columns]
    return df


def is_emar_type(feature_type: str) -> bool:
    """
    Considera come EMAR tutte le righe che contengono emar:
    - emar
    - emar_only
    - emar_and_other
    """
    ft = normalize_text(feature_type)
    return "emar" in ft


def get_other_category(feature_type: str) -> Optional[str]:
    ft = normalize_text(feature_type)
    for cat in OTHER_FEATURES:
        if cat in ft:
            return cat
    return None


def overlaps(start1: int, end1: int, start2: int, end2: int) -> bool:
    """
    Overlap tra intervalli chiusi [start, end].
    Basta anche 1 bp in comune.
    """
    return max(start1, start2) <= min(end1, end2)


def classify_single_file(input_csv: str, output_csv: str) -> None:
    df = pd.read_csv(input_csv, sep=SEP, engine="python")
    df = clean_columns(df)

    required_cols = ["feature_type", "feature_start", "feature_end"]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Nel file {input_csv} manca la colonna obbligatoria '{col}'")

    # Ripulisce spazi nei valori testuali
    df["feature_type"] = df["feature_type"].astype(str).str.strip()

    # Converte coordinate
    df["feature_start"] = pd.to_numeric(df["feature_start"], errors="coerce")
    df["feature_end"] = pd.to_numeric(df["feature_end"], errors="coerce")

    bad_mask = df["feature_start"].isna() | df["feature_end"].isna()
    if bad_mask.any():
        n_bad = int(bad_mask.sum())
        print(f"[WARN] {os.path.basename(input_csv)}: {n_bad} righe con coordinate non valide verranno ignorate")
        df = df.loc[~bad_mask].copy()

    df["feature_start"] = df["feature_start"].astype(int)
    df["feature_end"] = df["feature_end"].astype(int)

    # Aggiunge colonna informativa
    df["emar_overlap_with"] = ""

    # Estrae righe promoter/enhancer/ctcf
    other_rows = []
    for idx, row in df.iterrows():
        cat = get_other_category(row["feature_type"])
        if cat is not None:
            other_rows.append({
                "cat": cat,
                "start": int(row["feature_start"]),
                "end": int(row["feature_end"]),
            })

    total_emar = 0
    total_only = 0
    total_and_other = 0

    for idx, row in df.iterrows():
        ftype = row["feature_type"]

        if not is_emar_type(ftype):
            continue

        total_emar += 1
        emar_start = int(row["feature_start"])
        emar_end = int(row["feature_end"])

        overlap_cats: Set[str] = set()

        for other in other_rows:
            if overlaps(emar_start, emar_end, other["start"], other["end"]):
                overlap_cats.add(other["cat"])

        if overlap_cats:
            df.at[idx, "feature_type"] = "emar_and_other"
            df.at[idx, "emar_overlap_with"] = ";".join(sorted(overlap_cats))
            total_and_other += 1
        else:
            df.at[idx, "feature_type"] = "emar_only"
            df.at[idx, "emar_overlap_with"] = ""
            total_only += 1

    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    df.to_csv(output_csv, sep=SEP, index=False)

    print(f"\nOK: {os.path.basename(input_csv)}")
    print(f"  salvato -> {output_csv}")
    print(f"  EMAR totali: {total_emar}")
    print(f"  emar_only: {total_only}")
    print(f"  emar_and_other: {total_and_other}")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    files = sorted(glob.glob(os.path.join(INPUT_DIR, FILE_PATTERN)))
    if not files:
        raise FileNotFoundError(f"Nessun file trovato in {INPUT_DIR} con pattern {FILE_PATTERN}")

    print(f"Trovati {len(files)} file")

    for input_csv in files:
        base = os.path.basename(input_csv)
        stem, ext = os.path.splitext(base)
        output_csv = os.path.join(OUTPUT_DIR, f"{stem}_emar_classified{ext}")
        classify_single_file(input_csv, output_csv)

    print("\n✅ Classificazione completata su tutti i file.")


if __name__ == "__main__":
    main()