#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import pandas as pd

# =========================
# CONFIG
# =========================

COORD_FILE = r"C:\Users\Chiara Panico\OneDrive - Politecnico di Torino\Desktop\verifica codici\Coordinate_COMPLETE_5geni.csv"

INTRON_FILES = {
    "ASH1L": r"C:\Users\Chiara Panico\OneDrive - Politecnico di Torino\Desktop\verifica codici\analisi introni\introni\introni_ASH1L.csv",
    "NSD1":  r"C:\Users\Chiara Panico\OneDrive - Politecnico di Torino\Desktop\verifica codici\analisi introni\introni\introni_NSD1.csv",
    "NSD2":  r"C:\Users\Chiara Panico\OneDrive - Politecnico di Torino\Desktop\verifica codici\analisi introni\introni\introni_NSD2.csv",
    "NSD3":  r"C:\Users\Chiara Panico\OneDrive - Politecnico di Torino\Desktop\verifica codici\analisi introni\introni\intron_NSD3.csv",
    "SETD2": r"C:\Users\Chiara Panico\OneDrive - Politecnico di Torino\Desktop\verifica codici\analisi introni\introni\SETD2.csv",
}

OUT_FILE = r"C:\Users\Chiara Panico\OneDrive - Politecnico di Torino\Desktop\verifica codici\Coordinate_INTRONI_5geni.csv"

FLANK_BP = 1000

# Se True, blocca il programma se la lunghezza originale non coincide con quella
# ricavata dalle coordinate normalizzate
STRICT_LENGTH_CHECK = True

# =========================
# HELPERS
# =========================

def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [
        str(c).strip().replace("\ufeff", "").strip()
        for c in df.columns
    ]
    return df


def find_column(df: pd.DataFrame, candidates):
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    return None


def parse_intron_label(label: str) -> str:
    """
    Normalizza etichette tipo:
    - Intron1-2
    - Intron 1-2
    - intron 1-2
    Restituisce ad es. 'Intron 1-2'
    """
    s = str(label).strip()
    s = re.sub(r"\s+", " ", s)
    s = re.sub(r"(?i)^intron\s*", "Intron ", s)
    s = re.sub(r"(?i)^introne\s*", "Intron ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s


def intron_feature_id(gene: str, label: str) -> str:
    """
    Es: ASH1L_INTRON_1_2
    """
    m = re.search(r"(\d+)\s*-\s*(\d+)", label)
    if m:
        return f"{gene}_INTRON_{m.group(1)}_{m.group(2)}"
    safe = re.sub(r"[^A-Za-z0-9]+", "_", label.upper()).strip("_")
    return f"{gene}_{safe}"


# =========================
# MAIN
# =========================

def main():
    coord = pd.read_csv(COORD_FILE)
    coord = normalize_columns(coord)

    required_coord_cols = [
        "gene_symbol", "gene_id", "gene_chr", "gene_start", "gene_end", "gene_strand"
    ]
    for c in required_coord_cols:
        if c not in coord.columns:
            raise ValueError(f"Manca la colonna '{c}' in {COORD_FILE}")

    # una riga riassuntiva per gene
    gene_info = (
        coord.groupby("gene_symbol", as_index=False)
        .agg({
            "gene_id": "first",
            "gene_chr": "first",
            "gene_start": "min",
            "gene_end": "max",
            "gene_strand": "first"
        })
    )

    all_rows = []

    for gene, path in INTRON_FILES.items():
        if not os.path.exists(path):
            print(f"[WARN] file non trovato: {path}")
            continue

        df = pd.read_csv(path, encoding="utf-8-sig")
        df = normalize_columns(df)

        label_col = find_column(df, ["Intron", "introne", "Exon / Intron", "Exon/Intron"])
        start_col = find_column(df, ["Start", "start"])
        end_col = find_column(df, ["End", "end"])
        length_col = find_column(df, ["Length", "length", "Lunghezza", "lunghezza"])

        if label_col is None or start_col is None or end_col is None or length_col is None:
            raise ValueError(
                f"Colonne non riconosciute in {path}. Colonne trovate: {list(df.columns)}"
            )

        ginfo = gene_info[gene_info["gene_symbol"] == gene]
        if ginfo.empty:
            raise ValueError(f"Gene {gene} non trovato in {COORD_FILE}")

        ginfo = ginfo.iloc[0]
        gene_id = ginfo["gene_id"]
        gene_chr = ginfo["gene_chr"]
        gene_start = int(ginfo["gene_start"])
        gene_end = int(ginfo["gene_end"])
        gene_strand = int(ginfo["gene_strand"])

        ext_start = gene_start - FLANK_BP

        for _, row in df.iterrows():
            raw_label = row[label_col]
            raw_start = pd.to_numeric(row[start_col], errors="coerce")
            raw_end = pd.to_numeric(row[end_col], errors="coerce")
            raw_length = pd.to_numeric(row[length_col], errors="coerce")

            if pd.isna(raw_start) or pd.isna(raw_end) or pd.isna(raw_length):
                continue

            # Normalizzazione coordinate: sempre ordine genomico crescente
            feature_start = int(min(raw_start, raw_end))
            feature_end = int(max(raw_start, raw_end))
            feature_length = int(raw_length)

            # Lunghezza attesa dalle coordinate genomiche inclusive
            expected_length = feature_end - feature_start + 1

            label = parse_intron_label(raw_label)
            feature_id = intron_feature_id(gene, label)

            if feature_length != expected_length:
                msg = (
                    f"[ERR] Lunghezza non coerente per {gene} - {label}\n"
                    f"      file: {os.path.basename(path)}\n"
                    f"      start originale = {raw_start}\n"
                    f"      end originale   = {raw_end}\n"
                    f"      start normaliz. = {feature_start}\n"
                    f"      end normaliz.   = {feature_end}\n"
                    f"      length file     = {feature_length}\n"
                    f"      length attesa   = {expected_length}"
                )
                if STRICT_LENGTH_CHECK:
                    raise ValueError(msg)
                else:
                    print(msg)

            local_start0 = feature_start - ext_start
            local_end0_exclusive = feature_end - ext_start

            all_rows.append({
                "subset": "INTRON",
                "gene_symbol": gene,
                "gene_id": gene_id,
                "gene_chr": gene_chr,
                "gene_start": gene_start,
                "gene_end": gene_end,
                "gene_strand": gene_strand,
                "feature_id": feature_id,
                "feature_type": "intron",
                "feature_label": label,
                "feature_chr": gene_chr,
                "feature_start": feature_start,
                "feature_end": feature_end,
                "feature_strand": gene_strand,
                "feature_length": feature_length,
                "local_start0": local_start0,
                "local_end0_exclusive": local_end0_exclusive,
                "logic_name": "intron",
                "mane_transcript_id": "",
                "info": ""
            })

    out = pd.DataFrame(all_rows)

    # stesso ordine colonne del file coordinate completo + feature_length
    target_cols = [
        "subset", "gene_symbol", "gene_id", "gene_chr", "gene_start", "gene_end",
        "gene_strand", "feature_id", "feature_type", "feature_label", "feature_chr",
        "feature_start", "feature_end", "feature_strand", "feature_length",
        "local_start0", "local_end0_exclusive", "logic_name", "mane_transcript_id", "info"
    ]
    out = out[target_cols].sort_values(["gene_symbol", "feature_start", "feature_end"])

    out.to_csv(OUT_FILE, index=False)
    print(f"Creato: {OUT_FILE}")
    print(out.head(10).to_string(index=False))


if __name__ == "__main__":
    main()