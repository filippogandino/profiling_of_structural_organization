#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script che converte i dati di varianza (dalla directory 4_KL_csv) 
in valori di divergenza KL (Kullback-Leibler), organizzati per bande di entropia H.

Procedimento:
1. Legge i file di varianza media dai diversi valori di H (H=0.50 ... H=0.60)
2. Calcola due metriche di divergenza KL con riferimento a H=0.56:
   - gruppo 1 (0.50<H<0.56): somma delle differenze per H = 0.50, 0.51, ..., 0.55 vs H=0.56
   - gruppo 2 (0.57<H<0.60): somma delle differenze per H = 0.57, 0.58, 0.59, 0.60 vs H=0.56
3. Normalizza i valori KL per ogni gene (diviso per il massimo)
4. Organizza i risultati in due cartelle per bande di entropia:
   - H_50_56/: file con divergenza KL per la banda inferiore
   - H_57_60/: file con divergenza KL per la banda superiore

Output: file CSV per gene con 3 colonne (step, KL_value, KL_normalized)
"""
import os
import glob
import pandas as pd
import re

in_dir = r"C:\Users\Chiara Panico\OneDrive - Politecnico di Torino\Desktop\verifica codici\ANALISI 0.45-0.70\importa45-70"
out_dir = r"C:\Users\Chiara Panico\OneDrive - Politecnico di Torino\Desktop\verifica codici\ANALISI 0.45-0.70\KL_45-70"
os.makedirs(out_dir, exist_ok=True)

# crea sottocartelle per i bin
out1_dir = os.path.join(out_dir, "H_50_56")
out2_dir = os.path.join(out_dir, "H_57_60")
os.makedirs(out1_dir, exist_ok=True)
os.makedirs(out2_dir, exist_ok=True)

# ===================== HELPERS =====================

def normalize_kl_values(df, value_column):
    """
    Normalizza i valori KL dividendo per il massimo.
    Restituisce il dataframe con la colonna 'KL_normalized' aggiunta.
    """
    kl_values = pd.to_numeric(df[value_column], errors="coerce").dropna()
    
    if len(kl_values) == 0:
        return df, None
    
    kl_max = kl_values.max()
    if kl_max == 0:
        return df, None
    
    kl_normalized = (kl_values / kl_max).values
    
    df["KL_normalized"] = pd.NA
    df.loc[kl_values.index, "KL_normalized"] = kl_normalized
    
    return df, kl_max

ref = "H=0.56"
group1 = ["H=0.50", "H=0.51", "H=0.52", "H=0.53", "H=0.54", "H=0.55"]
group2 = ["H=0.57", "H=0.58", "H=0.59", "H=0.60"]

gene_pat = re.compile(r"media_(.+?)_sequence_varianze\.csv")

for path in glob.glob(os.path.join(in_dir, "*.csv")):
    fname = os.path.basename(path)
    m = gene_pat.match(fname)
    if not m:
        print(f"Skip: nome file non riconosciuto {fname}")
        continue

    gene = m.group(1)
    df = pd.read_csv(path)

    missing = [c for c in ([ref] + group1 + group2) if c not in df.columns]
    if missing:
        print(f"Skip {fname}: mancano colonne {missing}")
        continue

    # ----- primo file -----
    out1 = pd.DataFrame({
        "step": df["step"],
        "0.50<H<0.56": sum(df[c] - df[ref] for c in group1)
    })

    # ----- secondo file -----
    out2 = pd.DataFrame({
        "step": df["step"],
        "0.57<H<0.60": sum(df[c] - df[ref] for c in group2)
    })
    
    #Normalizza i valori KL
    out1, out1_max = normalize_kl_values(out1, "0.50<H<0.56")
    out2, out2_max = normalize_kl_values(out2, "0.57<H<0.60")

    out1_path = os.path.join(out1_dir, f"{gene}_10.csv")
    out2_path = os.path.join(out2_dir, f"{gene}_10.csv")

    out1.to_csv(out1_path, index=False)
    out2.to_csv(out2_path, index=False)

    msg1 = f" max={out1_max:.6e}" if out1_max else " (skip normalizzazione)"
    msg2 = f" max={out2_max:.6e}" if out2_max else " (skip normalizzazione)"
    
    print("Creati:")
    print(f"  {out1_path}{msg1}")
    print(f"  {out2_path}{msg2}")