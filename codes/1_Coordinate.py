#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GRCh38 • Per ogni gene in genes.txt:
  - Regulatory features nella finestra [1000 upstream] + [gene] + [1000 downstream]
    con risoluzione della CLASSE (Promoter, Enhancer, CTCF binding site, Open chromatin region, EMAR, ...)
  - Esone-set MANE (MANE Select; fallback: MANE Plus Clinical; poi Canonico) = EXON_MANE
  - Unione esoni di tutti i trascritti (merge) = EXON_UNION
  - Conversione in coordinate locali 0-based half-open (inizio gene = 1000)

Output: tuttoforse.csv
  colonne chiave:
    subset ∈ {REGULATORY, EXON_MANE, EXON_UNION}
    feature_id (per REGULATORY = ENSR..., per EXON_* lasciato vuoto)
    feature_type (grezza: 'regulatory' o 'exon')
    feature_label (umana: Promoter/Enhancer/CTCF/... oppure 'Exon')
    feature_start, feature_end (genomiche 1-based inclusive)
    local_start0, local_end0_exclusive (locali 0-based half-open)
    mane_transcript_id (per EXON_MANE)
"""

import csv
import sys
import time
from functools import lru_cache
from typing import Dict, Any, List, Optional, Tuple

import requests

# ======================= CONFIG =======================
SPECIES = "homo_sapiens"          # Ensembl default = GRCh38
GENE_FILE = "C:/Users/Chiara Panico/Dropbox/Chiara Panico/GENI/ANALISI_GENI/dati_geni/0_genes.txt"           # uno per riga: simbolo o ENSG
OUT_CSV = "1_Coordinate_geni.csv"
UPSTREAM = 1000
DOWNSTREAM = 1000
SLEEP_BETWEEN_GENES = 0.2
# ======================================================

BASE = "https://rest.ensembl.org"
HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json",
    "User-Agent": "gene-features-mane-union/1.3 (contact: you@example.com)"
}

# Etichette umane (include sinonimi comuni)
TYPE_MAP = {
    "promoter": "Promoter",
    "promoter_flanking_region": "Promoter flanking region",
    "enhancer": "Enhancer",
    "ctcf_binding_site": "CTCF binding site",
    "ctcf": "CTCF binding site",
    "tf_binding_site": "TF binding site",
    "open_chromatin_region": "Open chromatin region",
    "open chromatin": "Open chromatin region",
    "tss": "Transcription start site (TSS)",
    "tss_flanking_region": "TSS flanking region",
    "regulatory_feature": "Regulatory feature",
    "epigenetically modified accessible region": "EMAR",
    "emar": "EMAR",
    "exon": "Exon"
}

def pretty_type(t: Optional[str]) -> Optional[str]:
    if not t:
        return None
    key = str(t).strip().lower()
    return TYPE_MAP.get(key, key.replace("_", " ").title())

# ------------------ HTTP helper ------------------
SESSION = requests.Session()
SESSION.headers.update(HEADERS)

def _get(url: str, params: Optional[Dict[str, Any]] = None, sleep_s: float = 0.1):
    r = SESSION.get(url, params=params, timeout=120)
    if r.status_code == 429:
        retry_after = float(r.headers.get("Retry-After", "1.0"))
        time.sleep(retry_after)
        r = SESSION.get(url, params=params, timeout=120)
    r.raise_for_status()
    time.sleep(sleep_s)
    return r.json()

# un wrapper identico (usato nella parte "patch")
def _get_json(url: str, params: Optional[Dict[str, Any]] = None, sleep_s: float = 0.1):
    r = SESSION.get(url, params=params, timeout=120)
    if r.status_code == 429:
        retry_after = float(r.headers.get("Retry-After", "1.0"))
        time.sleep(retry_after)
        r = SESSION.get(url, params=params, timeout=120)
    r.raise_for_status()
    time.sleep(sleep_s)
    return r.json()
# -------------------------------------------------

# ------------------ Ensembl queries ------------------
def lookup_gene(symbol_or_id: str, species: str = SPECIES) -> Optional[Dict[str, Any]]:
    """
    /lookup/symbol o /lookup/id con expand=1 per ottenere Transcripts->Exons.
    """
    if symbol_or_id.upper().startswith("ENSG"):
        url = f"{BASE}/lookup/id/{symbol_or_id}"
    else:
        url = f"{BASE}/lookup/symbol/{species}/{symbol_or_id}"
    try:
        return _get(url, params={"expand": 1})
    except requests.HTTPError as e:
        sys.stderr.write(f"[WARN] Gene non trovato per '{symbol_or_id}': {e}\n")
        return None

def get_regulatory_in_region(chr_name: str, start: int, end: int, species: str = SPECIES) -> List[Dict[str, Any]]:
    """
    Regulatory features via /overlap/region (grezzo: 'feature_type' = 'regulatory', 'id' = ENSR...).
    """
    region = f"{chr_name}:{start}-{end}"
    url = f"{BASE}/overlap/region/{species}/{region}"
    try:
        data = _get(url, params={"feature": "regulatory"})
        return data if isinstance(data, list) else []
    except requests.HTTPError:
        return []
# -----------------------------------------------------

# ============ PATCH: risoluzione label delle regulatory ============
@lru_cache(maxsize=100000)
def fetch_regulatory_label(reg_id: str, species: str = SPECIES) -> Optional[str]:
    """
    Data una feature ENSR, prova a ottenere la sua CLASSE (Promoter/Enhancer/CTCF/…)
    dall'endpoint di dettaglio di Ensembl Regulation.
    Ritorna una stringa leggibile o None se non disponibile.
    """
    if not reg_id:
        return None

    # endpoint principale (Ensembl REST Regulation)
    url = f"{BASE}/regulatory/species/{species}/id/{reg_id}"
    try:
        data = _get_json(url)
    except requests.HTTPError:
        # fallback: prova overlap per ID (alcune release lo risolvono)
        try:
            data2 = _get_json(f"{BASE}/overlap/id/{reg_id}", params={"feature": "regulatory"})
            if isinstance(data2, list) and data2:
                d0 = data2[0]
                label = d0.get("description") or d0.get("feature_class") or d0.get("feature_set")
                if label:
                    return str(label)
        except requests.HTTPError:
            return None
        return None

    # Cerca campi utili: 'feature_class', 'description', 'type', 'label', 'class'
    if isinstance(data, dict):
        for key in ("feature_class", "description", "type", "label", "class"):
            val = data.get(key)
            if val:
                return str(val)

        # a volte c'è 'attributes': [{name,value}, ...]
        attrs = data.get("attributes")
        if isinstance(attrs, list):
            for a in attrs:
                if str(a.get("name", "")).lower() in ("class", "feature_class", "type"):
                    v = a.get("value")
                    if v:
                        return str(v)

    return None

def normalize_regulatory_rows_with_labels(reg_feats: List[Dict[str, Any]]) -> List[Tuple[int,int,str,str,str]]:
    """
    Estrae (start, end, feature_id, human_label, raw_type) per le REGULATORY.
    - feature_id: ENSR…
    - human_label: Promoter/Enhancer/CTCF/... se recuperabile, altrimenti 'Regulatory feature'
    - raw_type: tipicamente 'regulatory'
    """
    out = []
    for f in reg_feats:
        s = int(f.get("start", 0)); e = int(f.get("end", 0))
        raw = f.get("feature_type", "regulatory")
        rid = f.get("id") or f.get("ID") or ""
        resolved = fetch_regulatory_label(rid)
        if resolved:
            human = pretty_type(resolved) or "Regulatory feature"
        else:
            human = "Regulatory feature"
        out.append((s, e, rid, human, raw))
    return out
# =====================================================

# ------------------ EXON helpers ------------------
def extract_all_exons_from_gene(gene_info: Dict[str, Any]) -> List[Tuple[int,int]]:
    """
    Unisce gli esoni da TUTTI i trascritti (coordinate distinte esatte).
    """
    exs = []
    seen = set()
    for tr in gene_info.get("Transcript", []) or []:
        for ex in tr.get("Exon", []) or []:
            s = int(ex.get("start")); e = int(ex.get("end"))
            key = (s, e)
            if key in seen:
                continue
            seen.add(key)
            exs.append((s, e))
    exs.sort(key=lambda x: (x[0], x[1]))
    return exs

def pick_mane_transcript(gene_info: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """
    Sceglie il trascritto MANE:
      1) MANE Select
      2) MANE Plus Clinical
      3) Canonico Ensembl
      4) Fallback: quello con più esoni (poi più lungo)
    """
    trs = gene_info.get("Transcript", []) or []
    if not trs:
        return None

    def has_tag(tr, prefix: str) -> bool:
        return any((t or "").lower().startswith(prefix) for t in (tr.get("tags") or []))

    mane_select = [t for t in trs if has_tag(t, "mane select")]
    if mane_select:
        mane_select.sort(key=lambda t: (len(t.get("Exon", []) or []),
                                        int(t.get("end",0))-int(t.get("start",0))), reverse=True)
        return mane_select[0]

    mane_pc = [t for t in trs if has_tag(t, "mane plus clinical")]
    if mane_pc:
        mane_pc.sort(key=lambda t: (len(t.get("Exon", []) or []),
                                    int(t.get("end",0))-int(t.get("start",0))), reverse=True)
        return mane_pc[0]

    canonical = [t for t in trs if int(t.get("is_canonical", 0)) == 1]
    if canonical:
        canonical.sort(key=lambda t: (len(t.get("Exon", []) or []),
                                      int(t.get("end",0))-int(t.get("start",0))), reverse=True)
        return canonical[0]

    trs.sort(key=lambda t: (len(t.get("Exon", []) or []),
                            int(t.get("end",0))-int(t.get("start",0))), reverse=True)
    return trs[0]

def extract_exons_from_transcript(tr: Dict[str, Any]) -> List[Tuple[int,int]]:
    """Ritorna (start,end) 1-based inclusive per il trascritto scelto."""
    ex = [(int(e["start"]), int(e["end"])) for e in (tr.get("Exon") or [])]
    ex.sort(key=lambda x: (x[0], x[1]))
    return ex

def merge_intervals(intervals: List[Tuple[int,int]]) -> List[Tuple[int,int]]:
    """Merge 1-based inclusive: fonde intervalli sovrapposti o adiacenti (+1)."""
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        ms, me = merged[-1]
        if s <= me + 1:
            merged[-1] = (ms, max(me, e))
        else:
            merged.append((s, e))
    return merged
# -----------------------------------------------------

# ------------------ Row builder ------------------
def add_rows(rows: List[Dict[str, Any]],
             gene_symbol: str,
             gene_info: Dict[str, Any],
             items: List[Tuple[int,int]],
             subset: str,
             feature_label: str,
             feature_type_raw: str,
             extra: Optional[Dict[str, Any]] = None):
    """
    Aggiunge righe al CSV per ciascun intervallo (start,end).
    Calcola le coordinate locali 0-based half-open rispetto alla finestra estesa.
    """
    gene_start = int(gene_info["start"])
    gene_end   = int(gene_info["end"])
    chr_name   = gene_info["seq_region_name"]
    window_start = gene_start - UPSTREAM

    for s1, e1 in items:
        local_start0 = s1 - window_start
        local_end0_exclusive = e1 - window_start
        row = {
            "subset": subset,                       # REGULATORY / EXON_MANE / EXON_UNION
            "gene_symbol": gene_symbol,
            "gene_id": gene_info.get("id"),
            "gene_chr": chr_name,
            "gene_start": gene_start,
            "gene_end": gene_end,
            "gene_strand": gene_info.get("strand"),

            "feature_id": "",
            "feature_type": feature_type_raw,       # 'exon' o 'regulatory'
            "feature_label": feature_label,         # 'Exon', 'Promoter', 'Enhancer', ...
            "feature_chr": chr_name,
            "feature_start": s1,                    # genomiche 1-based inclusive
            "feature_end": e1,
            "feature_strand": gene_info.get("strand"),

            "local_start0": local_start0,           # 0-based
            "local_end0_exclusive": local_end0_exclusive,  # 0-based half-open
            "logic_name": subset.lower(),
            "mane_transcript_id": ""
        }
        if extra:
            row.update(extra)
        rows.append(row)

# -----------------------------------------------------

def main():
    # carica lista geni
    try:
        with open(GENE_FILE, "r", encoding="utf-8") as fh:
            gene_list = [line.strip() for line in fh if line.strip()]
    except FileNotFoundError:
        sys.stderr.write(f"❌ File {GENE_FILE} non trovato.\n")
        sys.exit(1)

    all_rows: List[Dict[str, Any]] = []

    for sym in gene_list:
        print(f"🔍 Gene: {sym}  [GRCh38]")
        gene = lookup_gene(sym)
        if not gene:
            continue

        chr_name = gene["seq_region_name"]
        gstart = int(gene["start"]); gend = int(gene["end"])
        ext_start = gstart - UPSTREAM
        ext_end   = gend + DOWNSTREAM

        # 1) REGULATORY (finestra estesa) con label risolta
        regs = get_regulatory_in_region(chr_name, ext_start, ext_end)
        reg_norm = normalize_regulatory_rows_with_labels(regs)
        for s, e, rid, lbl, raw in reg_norm:
            add_rows(
                all_rows, sym, gene, [(s, e)],
                subset="REGULATORY",
                feature_label=lbl,            # es. 'Promoter', 'Enhancer', 'CTCF binding site', ...
                feature_type_raw=raw,         # 'regulatory'
                extra={"feature_id": rid}     # salva l'ENSR id nel CSV
            )

        # 2) EXON_MANE (MANE Select > MANE Plus Clinical > Canonico > fallback)
        tr = pick_mane_transcript(gene)
        exon_mane = extract_exons_from_transcript(tr) if tr else []
        add_rows(
            all_rows, sym, gene, exon_mane,
            subset="EXON_MANE",
            feature_label="Exon",
            feature_type_raw="exon",
            extra={"mane_transcript_id": (tr.get("id") if tr else "")}
        )

        # 3) EXON_UNION (merge di TUTTI i trascritti)
        ex_all = extract_all_exons_from_gene(gene)
        exon_union = merge_intervals(ex_all)
        add_rows(
            all_rows, sym, gene, exon_union,
            subset="EXON_UNION",
            feature_label="Exon",
            feature_type_raw="exon"
        )

        print(f"   ↳ Regulatory: {len(reg_norm)} | EXON_MANE: {len(exon_mane)} | EXON_UNION: {len(exon_union)}")
        time.sleep(SLEEP_BETWEEN_GENES)

    # scrivi CSV unico
    fieldnames = [
        "subset",
        "gene_symbol","gene_id","gene_chr","gene_start","gene_end","gene_strand",
        "feature_id","feature_type","feature_label",
        "feature_chr","feature_start","feature_end","feature_strand",
        "local_start0","local_end0_exclusive","logic_name",
        "mane_transcript_id"
    ]
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(all_rows)

    print(f"✅ Scritto {len(all_rows)} righe in {OUT_CSV}")
    print(f"ℹ️ Coordinate locali: 0-based half-open; inizio gene = indice {UPSTREAM} (UPSTREAM={UPSTREAM}).")


if __name__ == "__main__":
    main()
