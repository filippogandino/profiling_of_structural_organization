#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calcola correlazioni Pearson tra MappingSum (sliding windows) e H_value,
gestendo PIÙ BIN di H organizzati in sottocartelle.

Aggiunge:
- p-value
- intervalli di confidenza del coefficiente di Pearson
- correzioni per test multipli

Struttura consigliata:
H_BASE_DIR/
  H_50_56/
    GENE_STEP.csv
  H_57_60/
    GENE_STEP.csv

SLIDING_DIR/
  sliding_GENE_FEATURE_win10k_step10_full_genomic.csv
  sliding_GENE_FEATURE_win10k_step10_full_reversed.csv
  ...

Output:
- correlations_*.csv con colonne:
    H_bin, Gene, Feature, Region, Orientation, Step,
    Pearson_r, P_value, CI_low, CI_high, N_windows, Sliding_file,
    P_Bonferroni, P_Holm, P_FDR_BH, Reject_Bonferroni, Reject_Holm, Reject_FDR_BH

NOTE IMPORTANTI
---------------
1) Allineamento H vs MappingSum
   Se H e sliding hanno lunghezze diverse, vengono tagliate entrambe
   alla lunghezza minima, come facevi manualmente in Excel.

2) Pearson
   Usiamo scipy.stats.pearsonr perché restituisce sia r che il p-value,
   e permette di calcolare l'intervallo di confidenza.

3) Correzioni multiple
   Applichiamo:
   - Bonferroni
   - Holm
   - FDR Benjamini-Hochberg
"""

import os
import glob
import re
from typing import Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests


# =====================================
# CONFIGURAZIONE
# =====================================

H_BASE_DIR = r"C:\Users\user\OneDrive - Politecnico di Torino\Desktop\verifica codici\5_KL_5geniverifica"

H_BINS = {
    "0.50-0.56": "H_50_56",
    "0.57-0.60": "H_57_60",
}

SLIDING_DIR = r"C:\Users\user\OneDrive - Politecnico di Torino\Desktop\verifica codici\3_slidingresults_consensus_pc80"
OUTPUT_FILE = r"6_correlations_pc80.csv"

# Livello di confidenza per CI di Pearson
CI_LEVEL = 0.95

# Alpha per i test multipli
ALPHA = 0.05

# Metodi di correzione da applicare
APPLY_BONFERRONI = True
APPLY_HOLM = True
APPLY_FDR_BH = True

# Se True, salta correlazioni con meno di 3 punti
REQUIRE_AT_LEAST_3_POINTS = True

# Pattern file sliding con orientation finale
pattern_sliding = re.compile(
    r"sliding_(?P<gene>[^_]+)_(?P<feature>.+)_win10k_step(?P<step>\d+)_(?P<region>[^_]+)_(?P<orientation>genomic|reversed)\.csv$"
)

# Pattern file H
pattern_h = re.compile(r"(?P<gene>[^_]+)_(?P<step>\d+)\.csv$")


# =====================================
# FUNZIONI DI SUPPORTO
# =====================================

def load_hmetric_file(path: str) -> pd.Series:
    """
    Legge un file H e restituisce una Series numerica chiamata 'H_value'.

    Usa esclusivamente la colonna contenente '<H<'.
    """
    df = pd.read_csv(path)

    h_col = None
    for c in df.columns:
        if "<H<" in str(c):
            h_col = c
            break

    if h_col is None:
        raise ValueError(
            "Colonna con pattern '<H<' non trovata nel file H. "
            f"Colonne disponibili: {list(df.columns)}"
        )

    s = pd.to_numeric(df[h_col], errors="coerce").dropna().reset_index(drop=True)
    return s.rename("H_value")


def load_all_hmetrics_for_bin(h_dir: str):
    """
    Carica tutti i file H di una singola cartella bin.
    Ritorna dict: (gene, step) -> Series(H_value)
    """
    hmetrics = {}
    for path in glob.glob(os.path.join(h_dir, "*.csv")):
        fname = os.path.basename(path)
        m = pattern_h.match(fname)
        if not m:
            print(f"[INFO] File H ignorato (nome non matcha pattern): {fname}")
            continue

        gene = m.group("gene")
        step = int(m.group("step"))

        try:
            s = load_hmetric_file(path)
        except Exception as e:
            print(f"[ERR] Errore nel parsing H file {fname}: {e}")
            continue

        hmetrics[(gene, step)] = s

    return hmetrics


def load_sliding_mappingsum(path: str) -> Optional[pd.Series]:
    """
    Legge un file sliding e restituisce la colonna 'MappingSum' come Series numerica pulita.
    """
    try:
        df = pd.read_csv(path)
    except Exception as e:
        print(f"[ERR] Errore nel leggere {os.path.basename(path)}: {e}")
        return None

    if "MappingSum" not in df.columns:
        print(f"[WARN] In {os.path.basename(path)} manca la colonna 'MappingSum'. Salto.")
        return None

    s_map = pd.to_numeric(df["MappingSum"], errors="coerce").dropna().reset_index(drop=True)
    return s_map


def align_for_corr(s_h: pd.Series, s_map: pd.Series) -> Tuple[pd.Series, pd.Series, int]:
    """
    Allinea le due serie per la correlazione tagliando alla lunghezza minima.

    Il disallineamento atteso è fisiologico: l'ultima finestra KL viene scartata
    se la sequenza residua è < 10kbp, quindi len(s_h) può essere leggermente
    inferiore a len(s_map). Il troncamento avviene sempre dalla coda (tail).

    Returns: (s_h_aligned, s_map_aligned, n_usati)
    """
    n_h   = len(s_h)
    n_map = len(s_map)
    n     = min(n_h, n_map)

    discarded_h   = max(0, n_h   - n)
    discarded_map = max(0, n_map - n)

    return s_h.iloc[:n].reset_index(drop=True), s_map.iloc[:n].reset_index(drop=True), n


def alignment_diagnostics(gene: str, step: int, feature: str, orientation: str,
                           n_h: int, n_map: int) -> None:
    """
    Stampa un report diagnostico sull'allineamento tra serie H e sliding.
    Da chiamare ogni volta che le due lunghezze differiscono.
    """
    n_used     = min(n_h, n_map)
    discarded  = abs(n_h - n_map)
    longer     = "H-metric" if n_h > n_map else "sliding/mask"
    pct        = discarded / max(n_h, n_map) * 100

    print(
        f"  [ALIGN] gene={gene} | feat={feature} | orient={orientation} | step={step}\n"
        f"          H-metric={n_h} finestre | sliding={n_map} finestre\n"
        f"          → usate={n_used} | scartate dalla coda={discarded} ({pct:.2f}%)"
        f" [dalla serie più lunga: {longer}]"
    )


def compute_pearson_with_stats(
    s_h: pd.Series,
    s_map: pd.Series,
    ci_level: float = 0.95
) -> dict:
    """
    Calcola Pearson r, p-value e CI.
    Gestisce casi problematici restituendo NaN dove necessario.
    """
    out = {
        "Pearson_r": np.nan,
        "P_value": np.nan,
        "CI_low": np.nan,
        "CI_high": np.nan,
    }

    if len(s_h) != len(s_map):
        raise ValueError("Le due serie devono avere la stessa lunghezza dopo l'allineamento.")

    n = len(s_h)
    if n < 2:
        return out

    x = pd.to_numeric(s_map, errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(s_h, errors="coerce").to_numpy(dtype=float)

    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]

    if len(x) < 2:
        return out

    if REQUIRE_AT_LEAST_3_POINTS and len(x) < 3:
        return out

    # Serie costanti -> Pearson non informativo
    if np.all(x == x[0]) or np.all(y == y[0]):
        return out

    try:
        res = pearsonr(x, y)
        out["Pearson_r"] = float(res.statistic)
        out["P_value"] = float(res.pvalue)

        try:
            ci = res.confidence_interval(confidence_level=ci_level)
            out["CI_low"] = float(ci.low)
            out["CI_high"] = float(ci.high)
        except Exception:
            pass

    except Exception:
        pass

    return out


def apply_multiple_testing_corrections(df: pd.DataFrame, alpha: float = 0.05) -> pd.DataFrame:
    """
    Aggiunge correzioni multiple ai p-value.
    """
    out = df.copy()

    if "P_value" not in out.columns:
        return out

    valid_mask = out["P_value"].notna()
    pvals = out.loc[valid_mask, "P_value"].to_numpy(dtype=float)

    if pvals.size == 0:
        return out

    if APPLY_BONFERRONI:
        reject, p_corr, _, _ = multipletests(pvals, alpha=alpha, method="bonferroni")
        out.loc[valid_mask, "P_Bonferroni"] = p_corr
        out.loc[valid_mask, "Reject_Bonferroni"] = reject

    if APPLY_HOLM:
        reject, p_corr, _, _ = multipletests(pvals, alpha=alpha, method="holm")
        out.loc[valid_mask, "P_Holm"] = p_corr
        out.loc[valid_mask, "Reject_Holm"] = reject

    if APPLY_FDR_BH:
        reject, p_corr, _, _ = multipletests(pvals, alpha=alpha, method="fdr_bh")
        out.loc[valid_mask, "P_FDR_BH"] = p_corr
        out.loc[valid_mask, "Reject_FDR_BH"] = reject

    return out


# =====================================
# MAIN
# =====================================

def main():
    # ===============================
    # Carica tutti gli H per ogni bin
    # ===============================
    all_h_by_bin = {}
    print("➡ Carico file H-metrics per ogni bin...")

    for bin_label, subdir in H_BINS.items():
        h_dir = os.path.join(H_BASE_DIR, subdir)
        if not os.path.isdir(h_dir):
            print(f"❌ Cartella bin non trovata: {h_dir}")
            raise SystemExit(1)

        hmetrics = load_all_hmetrics_for_bin(h_dir)
        if not hmetrics:
            print(f"❌ Nessun H-metric caricato per bin {bin_label} in {h_dir}")
            raise SystemExit(1)

        all_h_by_bin[bin_label] = hmetrics
        print(f"  ✓ {bin_label}: caricati {len(hmetrics)} file H (gene,step)")

    # ===============================
    # Scorri sliding files
    # ===============================
    print("\n➡ Carico sliding files e calcolo correlazioni...")

    results = []
    already_reported_alignment = set()

    for path in glob.glob(os.path.join(SLIDING_DIR, "sliding_*.csv")):
        fname = os.path.basename(path)
        m = pattern_sliding.match(fname)
        if not m:
            print(f"[INFO] Sliding file ignorato (nome non matcha pattern): {fname}")
            continue

        gene = m.group("gene")
        feature = m.group("feature")
        step = int(m.group("step"))
        region = m.group("region")
        orientation = m.group("orientation")

        s_map = load_sliding_mappingsum(path)
        if s_map is None:
            continue

        for bin_label, hmetrics in all_h_by_bin.items():
            key = (gene, step)
            if key not in hmetrics:
                continue

            s_h = hmetrics[key]

            if len(s_h) != len(s_map):
                alignment_diagnostics(
                    gene=gene, step=step, feature=feature,
                    orientation=orientation,
                    n_h=len(s_h), n_map=len(s_map)
                )
                if (gene, step) not in already_reported_alignment:
                    already_reported_alignment.add((gene, step))

            s_h_aligned, s_map_aligned, n = align_for_corr(s_h, s_map)

            stats = compute_pearson_with_stats(
                s_h=s_h_aligned,
                s_map=s_map_aligned,
                ci_level=CI_LEVEL,
            )

            if pd.isna(stats["Pearson_r"]):
                continue

            results.append({
                "H_bin": bin_label,
                "Gene": gene,
                "Feature": feature,
                "Region": region,
                "Orientation": orientation,
                "Step": step,
                "Pearson_r": stats["Pearson_r"],
                "P_value": stats["P_value"],
                "CI_low": stats["CI_low"],
                "CI_high": stats["CI_high"],
                "N_windows": n,
                "Sliding_file": fname,
            })

    # ===============================
    # Salva risultati
    # ===============================
    if not results:
        print("\n❌ Nessuna correlazione calcolata.")
        raise SystemExit(2)

    out = pd.DataFrame(results).sort_values(
        ["H_bin", "Step", "Gene", "Feature", "Region", "Orientation"]
    )

    out = apply_multiple_testing_corrections(out, alpha=ALPHA)

    out.to_csv(OUTPUT_FILE, index=False)

    print(f"\n✅ Correlazioni salvate in: {OUTPUT_FILE}")
    print(out.head(20))


if __name__ == "__main__":
    main()