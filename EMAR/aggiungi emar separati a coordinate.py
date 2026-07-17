#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Integra le coordinate EMAR_ONLY, EMAR_AND_OTHER e OTHER_GENIC dai file per-gene
nel file coordinate principale.

Assunzione:
- i file per-gene contengono già feature_type con valori tipo:
    * emar_only
    * emar_and_other
    * lncRNA / miRNA / pseudogene / ecc. per OTHER_GENIC
- qui non si ricalcola nulla: si fa solo import + adattamento formato

Logica duplicati:
- usa la tripla (feature_start, feature_end, feature_label_lower)

Convenzione coordinate:
- local_start0         = feature_start_1b - ext_start
- local_end0_exclusive = feature_end_1b   - ext_start

Strand:
- EMAR_ONLY / EMAR_AND_OTHER: usa gene_strand
- OTHER_GENIC: feature_strand lasciato vuoto; info salvato as-is
"""

import os
import pandas as pd

# ===================== CONFIGURAZIONE =====================

COORD_FILE  = r"C:\Users\user\OneDrive - Politecnico di Torino\Desktop\verifica codici\Coordinate_COMPLETE_5geni.csv"
FLANK_BP    = 1000
OUTPUT_FILE = r"C:\Users\user\OneDrive - Politecnico di Torino\Desktop\verifica codici\pippo.csv"

GENE_FILES = {
    "ASH1L": r"C:\Users\user\OneDrive - Politecnico di Torino\Desktop\verifica codici\coordinate_per_gene_emar_classified\ASH1L_coordinate emar_emar_classified.csv",
    "NSD1":  r"C:\Users\user\OneDrive - Politecnico di Torino\Desktop\verifica codici\coordinate_per_gene_emar_classified\NSD1_emar_classified.csv",
    "NSD2":  r"C:\Users\user\OneDrive - Politecnico di Torino\Desktop\verifica codici\coordinate_per_gene_emar_classified\NSD2_emar_classified.csv",
    "NSD3":  r"C:\Users\user\OneDrive - Politecnico di Torino\Desktop\verifica codici\coordinate_per_gene_emar_classified\nsd3_emar_classified.csv",
    "SETD2": r"C:\Users\user\OneDrive - Politecnico di Torino\Desktop\verifica codici\coordinate_per_gene_emar_classified\setd2_emar_classified.csv",
}

EMAR_ONLY_KEYWORDS = ["emar_only", "emar only"]
EMAR_AND_OTHER_KEYWORDS = ["emar_and_other", "emar and other", "emar+other", "emar + other"]

OTHER_GENIC_KEYWORDS = [
    "lncrna", "lnc_rna", "snrna", "mirna", "microrna",
    "pseudogene", "pseudogeno", "scrna", "scarna",
    "processed pseudo", "novel transcript"
]

# ==========================================================


def load_main_coords(path):
    df = pd.read_csv(path, sep=",", skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    df = df.apply(lambda col: col.str.strip() if col.dtype == "object" else col)

    if "info" not in df.columns:
        df["info"] = ""

    return df


def remove_old_emar_rows(df_main):
    """
    Rimuove eventuali righe EMAR/EMAR_ONLY/EMAR_AND_OTHER già presenti,
    così il file finale viene ricostruito pulito.
    """
    flabel = df_main.get("feature_label", pd.Series("", index=df_main.index)).astype(str).str.strip().str.lower()
    logic = df_main.get("logic_name", pd.Series("", index=df_main.index)).astype(str).str.strip().str.lower()

    mask_remove = (
        flabel.isin(["emar", "emar_only", "emar_and_other"]) |
        logic.isin(["emar", "emar_only", "emar_and_other"])
    )

    removed = int(mask_remove.sum())
    if removed > 0:
        print(f"Rimosse {removed} righe EMAR/EMAR_ONLY/EMAR_AND_OTHER già presenti dal file principale.")

    return df_main.loc[~mask_remove].copy(), removed


def get_gene_info(df_main, gene):
    rows = df_main[df_main["gene_symbol"] == gene]
    if rows.empty:
        raise ValueError(f"Gene '{gene}' non trovato nel file coordinate principale.")

    r = rows.iloc[0]
    return {
        "gene_symbol": gene,
        "gene_id": r["gene_id"],
        "gene_chr": r["gene_chr"],
        "gene_start": int(r["gene_start"]),
        "gene_end": int(r["gene_end"]),
        "gene_strand": int(r["gene_strand"]),
    }


def load_gene_file(path):
    df = pd.read_csv(path, sep=";", skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    df = df.apply(lambda col: col.str.strip() if col.dtype == "object" else col)

    rename_map = {}
    for c in df.columns:
        cl = c.lower()
        if cl in ("gene", "gene_symbol"):
            rename_map[c] = "gene_symbol"
        elif cl == "feature_type":
            rename_map[c] = "feature_type"
        elif cl == "feature_start":
            rename_map[c] = "feature_start"
        elif cl == "feature_end":
            rename_map[c] = "feature_end"
        elif cl == "info":
            rename_map[c] = "info"

    df = df.rename(columns=rename_map)

    if "info" not in df.columns:
        df["info"] = ""

    return df


def classify_row(feature_type):
    ft = str(feature_type).strip().lower()

    if any(kw in ft for kw in EMAR_ONLY_KEYWORDS):
        return "EMAR_ONLY"

    if any(kw in ft for kw in EMAR_AND_OTHER_KEYWORDS):
        return "EMAR_AND_OTHER"

    if any(kw in ft for kw in OTHER_GENIC_KEYWORDS):
        return "OTHER_GENIC"

    return None


def calc_local_coords(fs, fe, gene_start, flank=1000):
    ext_start = gene_start - flank
    return fs - ext_start, fe - ext_start


def build_existing_key_set(df_main, gene):
    """
    Tripla (start, end, label_lower) per controllo duplicati.
    """
    sub = df_main[df_main["gene_symbol"] == gene].copy()
    sub["feature_start"] = pd.to_numeric(sub["feature_start"], errors="coerce")
    sub["feature_end"] = pd.to_numeric(sub["feature_end"], errors="coerce")
    sub = sub.dropna(subset=["feature_start", "feature_end"])

    return set(zip(
        sub["feature_start"].astype(int),
        sub["feature_end"].astype(int),
        sub["feature_label"].astype(str).str.strip().str.lower(),
    ))


def build_new_row(label, fs, fe, local_s, local_e, ginfo, info_raw=""):
    if label == "EMAR_ONLY":
        subset = "REGULATORY"
        feature_type = "regulatory"
        feature_label = "EMAR_ONLY"
        logic_name = "emar_only"
        feature_strand = ginfo["gene_strand"]

    elif label == "EMAR_AND_OTHER":
        subset = "REGULATORY"
        feature_type = "regulatory"
        feature_label = "EMAR_AND_OTHER"
        logic_name = "emar_and_other"
        feature_strand = ginfo["gene_strand"]

    elif label == "OTHER_GENIC":
        subset = "OTHER_GENIC"
        feature_type = "other_genic"
        feature_label = "OTHER_GENIC"
        logic_name = "other_genic"
        feature_strand = ""

    else:
        raise ValueError(f"Label non gestita: {label}")

    return {
        "subset": subset,
        "gene_symbol": ginfo["gene_symbol"],
        "gene_id": ginfo["gene_id"],
        "gene_chr": ginfo["gene_chr"],
        "gene_start": ginfo["gene_start"],
        "gene_end": ginfo["gene_end"],
        "gene_strand": ginfo["gene_strand"],
        "feature_id": "",
        "feature_type": feature_type,
        "feature_label": feature_label,
        "feature_chr": ginfo["gene_chr"],
        "feature_start": fs,
        "feature_end": fe,
        "feature_strand": feature_strand,
        "local_start0": local_s,
        "local_end0_exclusive": local_e,
        "logic_name": logic_name,
        "mane_transcript_id": "",
        "info": info_raw,
    }


def process_gene(gene, gene_file_path, df_main):
    ginfo = get_gene_info(df_main, gene)
    existing_keys = build_existing_key_set(df_main, gene)

    ext_start = ginfo["gene_start"] - FLANK_BP
    ext_end = ginfo["gene_end"] + FLANK_BP
    seq_len = ext_end - ext_start

    df_gene = load_gene_file(gene_file_path)

    new_rows = []
    stats = {
        "emar_only_found": 0,
        "emar_only_added": 0,
        "emar_only_dup": 0,
        "emar_and_other_found": 0,
        "emar_and_other_added": 0,
        "emar_and_other_dup": 0,
        "other_found": 0,
        "other_added": 0,
        "other_dup": 0,
        "invalid": 0,
    }

    for _, row in df_gene.iterrows():
        label = classify_row(row["feature_type"])
        if label is None:
            continue

        if label == "EMAR_ONLY":
            stats["emar_only_found"] += 1
        elif label == "EMAR_AND_OTHER":
            stats["emar_and_other_found"] += 1
        elif label == "OTHER_GENIC":
            stats["other_found"] += 1

        try:
            fs = int(float(str(row["feature_start"]).replace(" ", "")))
            fe = int(float(str(row["feature_end"]).replace(" ", "")))
        except (ValueError, KeyError):
            print(f"   [WARN] Coordinate non valide per '{row['feature_type']}' — salto.")
            stats["invalid"] += 1
            continue

        if fs >= fe:
            print(f"   [WARN] start >= end ({fs}>={fe}) per '{row['feature_type']}' — salto.")
            stats["invalid"] += 1
            continue

        key = (fs, fe, label.lower())
        if key in existing_keys:
            if label == "EMAR_ONLY":
                stats["emar_only_dup"] += 1
            elif label == "EMAR_AND_OTHER":
                stats["emar_and_other_dup"] += 1
            else:
                stats["other_dup"] += 1
            continue

        local_s, local_e = calc_local_coords(fs, fe, ginfo["gene_start"], FLANK_BP)

        if local_e <= 0 or local_s >= seq_len:
            print(f"   [WARN] '{row['feature_type']}' ({fs}-{fe}) fuori dalla sequenza — salto.")
            stats["invalid"] += 1
            continue

        local_s = max(0, local_s)
        local_e = min(seq_len, local_e)

        info_raw = str(row.get("info", "")).strip()
        new_rows.append(build_new_row(label, fs, fe, local_s, local_e, ginfo, info_raw))
        existing_keys.add(key)

        if label == "EMAR_ONLY":
            stats["emar_only_added"] += 1
        elif label == "EMAR_AND_OTHER":
            stats["emar_and_other_added"] += 1
        else:
            stats["other_added"] += 1

    return new_rows, stats


def main():
    print("Carico file coordinate principale...")
    df_main = load_main_coords(COORD_FILE)
    print(f"   Righe iniziali: {len(df_main)}")

    df_main, removed_rows = remove_old_emar_rows(df_main)

    all_new_rows = []
    all_stats = {}

    for gene, gene_file_path in GENE_FILES.items():
        print(f"\nGene: {gene}")
        if not os.path.exists(gene_file_path):
            print(f"   [WARN] File non trovato: {gene_file_path} — salto.")
            continue

        try:
            new_rows, stats = process_gene(gene, gene_file_path, df_main)
        except ValueError as e:
            print(f"   [ERR] {e} — salto.")
            continue

        all_new_rows.extend(new_rows)
        all_stats[gene] = stats

        print(
            f"   EMAR_ONLY      -> trovati={stats['emar_only_found']} | "
            f"aggiunti={stats['emar_only_added']} | duplicati={stats['emar_only_dup']}"
        )
        print(
            f"   EMAR_AND_OTHER -> trovati={stats['emar_and_other_found']} | "
            f"aggiunti={stats['emar_and_other_added']} | duplicati={stats['emar_and_other_dup']}"
        )
        print(
            f"   OTHER_GENIC    -> trovati={stats['other_found']} | "
            f"aggiunti={stats['other_added']} | duplicati={stats['other_dup']}"
        )
        if stats["invalid"]:
            print(f"   Invalide/fuori sequenza: {stats['invalid']}")

    if all_new_rows:
        df_new = pd.DataFrame(all_new_rows)
        df_out = pd.concat([df_main, df_new], ignore_index=True)
    else:
        df_out = df_main.copy()
        print("\n[INFO] Nessuna nuova riga da aggiungere.")

    df_out.to_csv(OUTPUT_FILE, index=False)

    emar_only_tot = sum(s["emar_only_added"] for s in all_stats.values())
    emar_and_other_tot = sum(s["emar_and_other_added"] for s in all_stats.values())
    other_tot = sum(s["other_added"] for s in all_stats.values())

    print(f"\nFile salvato: {OUTPUT_FILE}")
    print(f"   Righe EMAR vecchie rimosse : {removed_rows}")
    print(f"   EMAR_ONLY aggiunti         : {emar_only_tot}")
    print(f"   EMAR_AND_OTHER aggiunti    : {emar_and_other_tot}")
    print(f"   OTHER_GENIC aggiunti       : {other_tot}")
    print(f"   Righe totali finali        : {len(df_out)}")


if __name__ == "__main__":
    main()