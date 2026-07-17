#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Pipeline unico (stage 2 + 3) per analisi su tanti geni, con:
- gestione strand del gene (GeneStrand = +1 / -1 salvata nei file)
- creazione delle mappe (mask) in entrambe le direzioni:
    * orientation="genomic": Position 0..N-1 segue coordinate genomiche crescenti
    * orientation="reversed": stessa mappa ma invertita (Position 0..N-1 parte dall'estremo opposto)

Requisiti:
- Coordinate.csv (output dello script 1) deve contenere almeno:
    gene_symbol, gene_start, gene_end, gene_strand,
    feature_start, feature_end, feature_type, feature_label
- se presenti, vengono usati anche:
    subset, logic_name

Output:
- mask/genomic/... + mask/reversed/... (con le stesse strutture di file attese dallo sliding)
- sliding_{gene}_{feature}_win10k_step{step}_{region}_{orientation}.csv
"""

import os
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# ===================== CONFIGURAZIONE =====================

FEATURE_TABLE = r"C:\Users\user\OneDrive - Politecnico di Torino\Desktop\verifica codici\ANALISI REGULATORY FEATURES\coordinate_REGULATORY_complete.csv"

# oppure lascia [] e fai leggere da genes.txt (uno per riga)
GENES: List[str] = []  # se vuoto, legge da GENE_FILE_FALLBACK

GENE_FILE_FALLBACK = r"C:\Users\user\OneDrive - Politecnico di Torino\Desktop\verifica codici\0_genes.txt"

FLANK_BP = 1000

# Directory base (il codice creerà due sottocartelle: genomic/ e reversed/)
BASE_DIR = r"C:\Users\user\OneDrive - Politecnico di Torino\Desktop\verifica codici\2_masks_sliding_regulatory"

# Nomi file per le mappe per-base
FILE_PATTERNS: Dict[str, str] = {
    "EXON_MANE":      "mappa_EXON_MANE_{gene}.csv",
    "INTRON_MANE":    "mappa_INTRON_MANE_{gene}.csv",
    "EXON_UNION":     "mappa_EXON_UNION_{gene}.csv",
    "INTRON_UNION":   "mappa_INTRON_UNION_{gene}.csv",
    "REG_ALL":        "mappa_REGULATORY_{gene}.csv",
    "PROMOTER":       os.path.join("regulatory_labels", "mappa_promoter_{gene}.csv"),
    "ENHANCER":       os.path.join("regulatory_labels", "mappa_enhancer_{gene}.csv"),
    "CTCF":           os.path.join("regulatory_labels", "mappa_ctcf_{gene}.csv"),
    "EMAR_ONLY":      os.path.join("emar_labels", "mappa_emar_only_{gene}.csv"),
    "EMAR_AND_OTHER": os.path.join("emar_labels", "mappa_emar_and_other_{gene}.csv"),
    "OTHER_GENIC":    os.path.join("other_genic_labels", "mappa_other_genic_{gene}.csv"),
}

# Sliding window
WINDOW_SIZE = 10_000
STEPS = [10]

REGIONS = {
    "full": {"start_pos": 0, "trim_end": 0, "inner_trim": 0},
    "winTrim500": {"start_pos": 0, "trim_end": 0, "inner_trim": 500},
}

OUT_DIR = "3_slidingresults_CON_Introni_EMAR"


# ===================== HELPERS =====================

def _load_genes() -> List[str]:
    if GENES:
        return GENES
    if not os.path.exists(GENE_FILE_FALLBACK):
        raise FileNotFoundError(f"GENES è vuota e non trovo {GENE_FILE_FALLBACK}")
    with open(GENE_FILE_FALLBACK, "r", encoding="utf-8") as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def assegna_categorie(row) -> List[str]:
    """
    Decide in quali mask accendere le basi per una data riga.

    Restituisce una lista di categorie tra:
      ["EXON_MANE", "INTRON_MANE", "EXON_UNION", "INTRON_UNION",
       "REG_ALL", "PROMOTER", "ENHANCER", "CTCF",
       "EMAR_ONLY", "EMAR_AND_OTHER", "OTHER_GENIC"]

    Regole:
    - per esoni/introni usa soprattutto subset e logic_name, se presenti
    - mantiene anche fallback su feature_type / feature_label
    - EMAR_ONLY e EMAR_AND_OTHER vengono assegnati da feature_label / logic_name
    - OTHER_GENIC viene assegnato da subset / feature_type / feature_label / logic_name
    """
    cats: List[str] = []

    subset = str(row.get("subset", "")).lower().strip()
    ftype = str(row.get("feature_type", "")).lower().strip()
    flabel = str(row.get("feature_label", "")).lower().strip()
    logic = str(row.get("logic_name", "")).lower().strip()

    # ---- ESONI ----
    if subset == "exon_mane" or logic == "exon_mane":
        cats.append("EXON_MANE")

    if subset == "exon_union" or logic == "exon_union":
        cats.append("EXON_UNION")

    # fallback vecchio
    if "exon" in ftype:
        cats.append("EXON_UNION")
        if "mane" in logic or "mane" in flabel or subset == "exon_mane":
            cats.append("EXON_MANE")

    # ---- INTRONI ----
    if subset == "intron_mane" or logic == "intron_mane":
        cats.append("INTRON_MANE")

    if subset == "intron_union" or logic == "intron_union":
        cats.append("INTRON_UNION")

    if "intron" in ftype:
        if "mane" in logic or "mane" in flabel or subset == "intron_mane":
            cats.append("INTRON_MANE")
        if "union" in logic or "union" in flabel or subset == "intron_union":
            cats.append("INTRON_UNION")

    # ---- REGOLATORY classiche ----
    if subset == "regulatory" or "regulat" in ftype or any(kw in flabel for kw in ["promoter", "enhancer", "ctcf"]):
        cats.append("REG_ALL")

    if "promoter" in ftype or "promoter" in flabel:
        cats.append("PROMOTER")
    if "enhancer" in ftype or "enhancer" in flabel:
        cats.append("ENHANCER")
    if "ctcf" in ftype or "ctcf" in flabel:
        cats.append("CTCF")

    # ---- EMAR esplicite dal file coordinate ----
    if (
        "emar_only" in ftype or
        "emar_only" in flabel or
        logic == "emar_only"
    ):
        cats.append("EMAR_ONLY")
        cats.append("REG_ALL")

    if (
        "emar_and_other" in ftype or
        "emar_and_other" in flabel or
        logic == "emar_and_other"
    ):
        cats.append("EMAR_AND_OTHER")
        cats.append("REG_ALL")

    # ---- OTHER_GENIC ----
    if (
        subset == "other_genic" or
        "other_genic" in ftype or
        "other_genic" in flabel or
        logic == "other_genic"
    ):
        cats.append("OTHER_GENIC")

    return sorted(set(cats))
    """
    Decide in quali mask accendere le basi per una data riga.

    Restituisce una lista di categorie tra:
      ["EXON_MANE", "INTRON_MANE", "EXON_UNION", "INTRON_UNION",
       "REG_ALL", "PROMOTER", "ENHANCER", "CTCF",
       "EMAR_ONLY", "EMAR_AND_OTHER"]

    Regole:
    - per esoni/introni usa soprattutto subset e logic_name, se presenti
    - mantiene anche fallback su feature_type / feature_label
    - EMAR_ONLY e EMAR_AND_OTHER vengono assegnati solo da nome/label
    """
    cats: List[str] = []

    subset = str(row.get("subset", "")).lower().strip()
    ftype = str(row.get("feature_type", "")).lower().strip()
    flabel = str(row.get("feature_label", "")).lower().strip()
    logic = str(row.get("logic_name", "")).lower().strip()

    # ---- ESONI / INTRONI da subset esplicito ----
    if subset == "exon_mane":
        cats.append("EXON_MANE")
    elif subset == "intron_mane":
        cats.append("INTRON_MANE")
    elif subset == "exon_union":
        cats.append("EXON_UNION")
    elif subset == "intron_union":
        cats.append("INTRON_UNION")

    # ---- Fallback: logic_name ----
    if logic == "exon_mane":
        cats.append("EXON_MANE")
    elif logic == "intron_mane":
        cats.append("INTRON_MANE")
    elif logic == "exon_union":
        cats.append("EXON_UNION")
    elif logic == "intron_union":
        cats.append("INTRON_UNION")

    # ---- Fallback ulteriore su tipo/label ----
    # utile se il file non ha subset ma contiene tracce nel nome
    if "exon" in ftype:
        if "mane" in logic or "mane" in flabel or subset == "exon_mane":
            cats.append("EXON_MANE")
        if "union" in logic or "union" in flabel or subset == "exon_union":
            cats.append("EXON_UNION")

    if "intron" in ftype:
        if "mane" in logic or "mane" in flabel or subset == "intron_mane":
            cats.append("INTRON_MANE")
        if "union" in logic or "union" in flabel or subset == "intron_union":
            cats.append("INTRON_UNION")

    # ---- REGOLATORY classiche ----
    if subset == "regulatory" or "regulat" in ftype or any(kw in flabel for kw in ["promoter", "enhancer", "ctcf"]):
        cats.append("REG_ALL")

    if "promoter" in ftype or "promoter" in flabel:
        cats.append("PROMOTER")
    if "enhancer" in ftype or "enhancer" in flabel:
        cats.append("ENHANCER")
    if "ctcf" in ftype or "ctcf" in flabel:
        cats.append("CTCF")

    # ---- EMAR esplicite dal file coordinate ----
    if "emar_only" in ftype or "emar_only" in flabel:
        cats.append("EMAR_ONLY")
        cats.append("REG_ALL")

    if "emar_and_other" in ftype or "emar_and_other" in flabel:
        cats.append("EMAR_AND_OTHER")
        cats.append("REG_ALL")

    return sorted(set(cats))


def autodetect_value_col(df: pd.DataFrame) -> str:
    candidates = [c for c in df.columns if c.lower().startswith("is")]
    if candidates:
        return candidates[0]
    for c in df.columns:
        if c != "Position" and pd.api.types.is_numeric_dtype(df[c]):
            return c
    raise ValueError("Impossibile auto-rilevare la colonna dei valori (Is* o numerica).")


def load_mask(path: str, value_col: Optional[str]) -> Tuple[np.ndarray, np.ndarray, str, Optional[int]]:
    """
    Ritorna pos, vals, col_usata, gene_strand (se presente)
    """
    df = pd.read_csv(path)
    if "Position" not in df.columns:
        raise ValueError(f"Nel file {path} manca la colonna 'Position'.")

    col = value_col or autodetect_value_col(df)
    if col not in df.columns:
        raise ValueError(f"Colonna '{col}' non trovata in {path}. Colonne: {list(df.columns)}")

    pos = df["Position"].to_numpy()
    vals = pd.to_numeric(df[col], errors="coerce").fillna(0).to_numpy()

    gstrand = None
    if "GeneStrand" in df.columns:
        try:
            gstrand = int(pd.to_numeric(df["GeneStrand"], errors="coerce").dropna().iloc[0])
        except Exception:
            gstrand = None

    return pos, vals, col, gstrand


def sliding_stats_single_mask(
    vals: np.ndarray,
    window_size: int,
    step: int,
    start_pos: int,
    end_excl: int,
    inner_trim: int = 0,
) -> pd.DataFrame:
    n = vals.size
    start_pos = max(start_pos, 0)
    end_excl = min(end_excl, n)

    empty_cols = [
        "WindowStart", "WindowEnd",
        "FeatureCount", "BlockCount",
        "MappingSum", "Fraction_Feature",
    ]

    if start_pos >= end_excl:
        return pd.DataFrame(columns=empty_cols)

    region_len = end_excl - start_pos
    if region_len < window_size:
        return pd.DataFrame(columns=empty_cols)

    if inner_trim * 2 >= window_size:
        raise ValueError(f"inner_trim troppo grande: {inner_trim}*2 >= {window_size}")

    effective_window_size = window_size - 2 * inner_trim

    b = (vals > 0)
    ps_bool = np.cumsum(np.insert(b.astype(np.int32), 0, 0))
    ps_vals = np.cumsum(np.insert(vals.astype(np.int64), 0, 0))

    run_start = np.zeros(n, dtype=bool)
    if n > 0 and b[0]:
        run_start[0] = True
    if n > 1:
        run_start[1:] = b[1:] & (~b[:-1])
    ps_run = np.cumsum(np.insert(run_start.astype(np.int32), 0, 0))

    starts_outer = np.arange(start_pos, end_excl - window_size + 1, step)
    ends_outer_excl = starts_outer + window_size
    if starts_outer.size == 0:
        return pd.DataFrame(columns=empty_cols)

    inner_starts = starts_outer + inner_trim
    inner_ends_excl = ends_outer_excl - inner_trim

    def sum_range(ps, s, e):
        return ps[e] - ps[s]

    feature_counts = sum_range(ps_bool, inner_starts, inner_ends_excl)
    mapping_sum = sum_range(ps_vals, inner_starts, inner_ends_excl)
    block_counts = sum_range(ps_run, inner_starts, inner_ends_excl)
    frac = feature_counts.astype(float) / float(effective_window_size)

    return pd.DataFrame({
        "WindowStart": inner_starts,
        "WindowEnd": inner_ends_excl - 1,
        "FeatureCount": feature_counts.astype(int),
        "BlockCount": block_counts.astype(int),
        "MappingSum": mapping_sum.astype(int),
        "Fraction_Feature": frac,
    })


# ===================== STAGE 2: MASK =====================

def build_masks_both_orientations(df: pd.DataFrame, genes: List[str]) -> Dict[str, str]:
    """
    Crea le mappe per-base per ogni gene e categoria in:
      BASE_DIR/genomic/...
      BASE_DIR/reversed/...
    Restituisce un dict orient->base_dir effettiva.
    """
    needed_cols = [
        "gene_symbol", "gene_start", "gene_end", "gene_strand",
        "feature_start", "feature_end",
        "feature_type", "feature_label"
    ]
    for c in needed_cols:
        if c not in df.columns:
            raise ValueError(f"Manca la colonna '{c}' nel file {FEATURE_TABLE}")

    orient_dirs = {
        "genomic": os.path.join(BASE_DIR, "genomic"),
        "reversed": os.path.join(BASE_DIR, "reversed"),
    }
    for d in orient_dirs.values():
        os.makedirs(d, exist_ok=True)

    for gene in genes:
        df_gene = df[df["gene_symbol"] == gene].copy()
        if df_gene.empty:
            print(f"[WARN] Nessuna riga per gene {gene}, salto.")
            continue

        gene_start = int(df_gene["gene_start"].min())
        gene_end = int(df_gene["gene_end"].max())

        try:
            gene_strand = int(pd.to_numeric(df_gene["gene_strand"], errors="coerce").dropna().iloc[0])
        except Exception:
            gene_strand = 0

        ext_start = gene_start - FLANK_BP
        ext_end = gene_end + FLANK_BP
        length = ext_end - ext_start

        if length <= 0:
            print(f"[ERR] Gene {gene}: length <= 0 (start={ext_start}, end={ext_end}), salto.")
            continue

        print(f"\n🧬 Gene {gene}: strand={gene_strand:+d} | ext_start={ext_start}, ext_end={ext_end}, N={length}")

        positions = np.arange(length, dtype=np.int64)
        masks = {key: np.full(length, -1, dtype=np.int8) for key in FILE_PATTERNS.keys()}

        for _, row in df_gene.iterrows():
            cats = assegna_categorie(row)
            if not cats:
                continue

            # Coordinate input: 1-based inclusive
            # Conversione a 0-based half-open:
            #   start_0b = start_1b - 1
            #   end_0b   = end_1b
            f_start = int(row["feature_start"]) - 1
            f_end = int(row["feature_end"])

            idx_start = max(0, f_start - ext_start)
            idx_end = min(length, f_end - ext_start)

            if idx_start >= idx_end:
                continue

            for cat in cats:
                masks[cat][idx_start:idx_end] = 1

        for orientation, base_out in orient_dirs.items():
            for cat, pattern in FILE_PATTERNS.items():
                arr = masks[cat]

                if orientation == "reversed":
                    arr_out = arr[::-1]
                    pos_out = positions[::-1]
                else:
                    arr_out = arr
                    pos_out = positions

                rel_path = pattern.format(gene=gene)
                out_path = os.path.join(base_out, rel_path)
                os.makedirs(os.path.dirname(out_path), exist_ok=True)

                df_out = pd.DataFrame({
                    "Position": pos_out,
                    "IsFeature": arr_out.astype(int),
                    "GeneStrand": gene_strand,
                    "Orientation": orientation,
                })
                df_out.to_csv(out_path, index=False)

        print(f"  ✓ salvate mask in: {orient_dirs['genomic']} e {orient_dirs['reversed']}")

    print("\n✅ Creazione mask per-base (entrambe le direzioni) completata.")
    return orient_dirs


# ===================== STAGE 3: SLIDING =====================

def run_sliding_for_orientations(genes: List[str], orient_dirs: Dict[str, str]) -> None:
    os.makedirs(OUT_DIR, exist_ok=True)

    feature_config = FILE_PATTERNS

    for gene in genes:
        print(f"\n🧬 Sliding → Gene: {gene}")
        for feat_name, pattern in feature_config.items():
            for orientation, base_dir in orient_dirs.items():
                rel_path = pattern.format(gene=gene)
                path = os.path.join(base_dir, rel_path)

                if not os.path.exists(path):
                    print(f"  [WARN] {feat_name} ({orientation}): file mancante → {path} (salto)")
                    continue

                try:
                    _, vals, used_col, gstrand = load_mask(path, None)
                except Exception as e:
                    print(f"  [ERR] {feat_name} ({orientation}): errore nel file {path} → {e}")
                    continue

                N = vals.size
                strand_txt = f"{gstrand:+d}" if isinstance(gstrand, int) else "NA"
                print(f"  ▶ {feat_name} | orientation={orientation} | strand={strand_txt} | N={N} | col='{used_col}'")

                for region_name, cfg in REGIONS.items():
                    start = cfg["start_pos"]
                    trim = cfg["trim_end"]
                    inner_trim = cfg["inner_trim"]
                    end_excl = N - trim if trim > 0 else N

                    if start >= end_excl:
                        continue

                    for step in STEPS:
                        df_out = sliding_stats_single_mask(
                            vals=vals,
                            window_size=WINDOW_SIZE,
                            step=step,
                            start_pos=start,
                            end_excl=end_excl,
                            inner_trim=inner_trim,
                        )
                        if df_out.empty:
                            continue

                        out_name = f"sliding_{gene}_{feat_name}_win10k_step{step}_{region_name}_{orientation}.csv"
                        out_path = os.path.join(OUT_DIR, out_name)
                        df_out.to_csv(out_path, index=False)

    print("\n✅ Sliding completato (entrambe le direzioni).")


# ===================== MAIN =====================

def main():
    genes = _load_genes()
    df = pd.read_csv(FEATURE_TABLE)

    orient_dirs = build_masks_both_orientations(df, genes)
    run_sliding_for_orientations(genes, orient_dirs)


if __name__ == "__main__":
    main()