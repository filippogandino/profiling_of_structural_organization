#questo codice prende i file di Filippo, modifica il formato di righe e colonne e li converte in csv per poterli usare nel passo successivo (per calcolare KL come somma di differenze)


import re
import glob
import os
import pandas as pd

# --- configurazione ---
input_pattern = r"C:\Users\Chiara Panico\Dropbox\Chiara Panico\GENI\varianza\media_*_sequence_varianze.txt"
out_dir = r"C:\Users\Chiara Panico\OneDrive - Politecnico di Torino\Desktop\verifica codici\4_KL_verifica"
os.makedirs(out_dir, exist_ok=True)

# cattura il gene dal nome file: media_<GENE>_sequence_varianze.txt
gene_pat = re.compile(r"^media_(.+?)_sequence_varianze\.txt$")

# cattura: H, step, valore
line_pat = re.compile(r"H=(\d+)\s+step=(\d+):\s+([0-9.eE+-]+)")

for path in glob.glob(input_pattern):
    fname = os.path.basename(path)
    m = gene_pat.match(fname)
    if not m:
        print(f"Skip (nome non riconosciuto): {fname}")
        continue

    gene = m.group(1)

    rows = []
    with open(path, "r") as f:
        for line in f:
            matches = line_pat.findall(line)
            if not matches:
                continue

            row = {}
            for H, step, value in matches:
                row["step"] = int(step)
                col = f"H={int(H)/100:.2f}"  # H=50 -> H=0.50
                row[col] = float(value)

            rows.append(row)

    if not rows:
        print(f"Nessun dato trovato in: {fname}")
        continue

    df = pd.DataFrame(rows).sort_values("step")

    # assicura ordine colonne: step poi H=0.50..H=0.60 (se presenti)
    h_cols = sorted([c for c in df.columns if c.startswith("H=")],
                    key=lambda x: float(x.split("=")[1]))
    df = df[["step"] + h_cols]

    out_csv = os.path.join(out_dir, f"media_{gene}_sequence_varianze.csv")
    df.to_csv(out_csv, index=False)
    print(f"Creato: {out_csv}")
