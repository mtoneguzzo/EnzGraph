#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, csv, sys, time, requests, re
from urllib.parse import quote

HEADER_MAP = {
    "Entry": "acc",
    "Entry Name": "entry",
    "Protein names": "protein",
    "Gene Names (primary)": "gene",
    "Organism": "organism",
    "Organism (ID)": "taxid",
    "Length": "length",
    "EC number": "ec_list",
    "Pfam": "pfam_list",
    "InterPro": "interpro_list",
    "PDB": "pdb_list",
    "Gene Ontology (molecular function)": "go_mf",
    "Gene Ontology (biological process)": "go_bp",
    "Subcellular location [CC]": "subcellular",
    "Catalytic activity": "catalytic",
}

def parse_args():
    ap = argparse.ArgumentParser(description="Enrich UniProt TSV with reviewed + CAZy/ESTHER xrefs (JSON API)")
    ap.add_argument("inp", help="Input TSV WITH HEADER (UniProt export)")
    ap.add_argument("out", help="Output TSV (same columns + enrichment)")
    ap.add_argument("--batch", type=int, default=200)
    ap.add_argument("--sleep", type=float, default=0.25)
    ap.add_argument("--retries", type=int, default=5)
    return ap.parse_args()

def normalize_header(h):
    h2 = h.strip()
    if h2.lower() == "ntry":  # in case header got truncated
        return "Entry"
    return h2

def read_with_header(inp):
    import csv
    with open(inp, newline='') as f:
        reader = csv.DictReader(f, delimiter="\t")
        norm_fieldnames = [normalize_header(h) for h in reader.fieldnames]
        to_internal = {}
        for raw in norm_fieldnames:
            to_internal[raw] = HEADER_MAP.get(raw, raw)
        rows = []
        for r in reader:
            nr = {}
            for raw_key, val in r.items():
                nk = to_internal.get(normalize_header(raw_key), raw_key)
                nr[nk] = val
            rows.append(nr)
        return rows, list(to_internal.values())

def safe_prop(props, key):
    for kv in props or []:
        if kv.get("key") == key:
            return kv.get("value","")
    return ""

def uniprot_fetch_json_batch(accs, retries=5, sleep=0.25):
    """
    Query UniProtKB search JSON and extract:
      - reviewed (bool via entryType)
      - xref_cazy (string; CAZy ids or family names joined by ';')
      - xref_esther (string; ESTHER ids/family names joined by ';')
    """
    if not accs:
        return {}

    # Build OR query on accession; request a generous page size
    q = " OR ".join(f"accession:{a}" for a in accs)
    url = f"https://rest.uniprot.org/uniprotkb/search?query={quote(q)}&format=json&size=500"

    attempt = 0
    while True:
        attempt += 1
        r = requests.get(url, timeout=60)
        if r.status_code == 200:
            break
        if r.status_code in (429, 502, 503, 504) and attempt <= retries:
            time.sleep(min(5.0, sleep * attempt * 2))
            continue
        r.raise_for_status()

    data = r.json()
    results = data.get("results", [])
    out = {}
    for rec in results:
        acc = rec.get("primaryAccession")
        if not acc:
            continue

        entry_type = rec.get("entryType","")
        reviewed = "reviewed" if entry_type.lower().startswith("uniprotkb reviewed") else "unreviewed"

        # Scan cross-references
        xrefs = rec.get("uniProtKBCrossReferences", []) or []
        cazy_vals, esther_vals = [], []
        for x in xrefs:
            db = x.get("database","")
            xid = x.get("id","")
            props = x.get("properties", [])
            if db == "CAZy":
                # pick FamilyName if present, else the id
                fam = safe_prop(props, "FamilyName")
                cazy_vals.append(fam or xid)
            elif db == "ESTHER":
                fam = safe_prop(props, "FamilyName")
                esther_vals.append(fam or xid)

        # deduplicate while preserving order
        def uniq(seq):
            seen=set(); outl=[]
            for t in seq:
                if t and t not in seen:
                    seen.add(t); outl.append(t)
            return outl

        out[acc] = {
            "reviewed": reviewed,
            "xref_cazy": "; ".join(uniq(cazy_vals)),
            "xref_esther": "; ".join(uniq(esther_vals)),
        }
    return out

def main():
    args = parse_args()
    rows, colnames = read_with_header(args.inp)

    accs = [r.get("acc","") for r in rows if r.get("acc")]
    seen, uniq = set(), []
    for a in accs:
        if a not in seen:
            seen.add(a); uniq.append(a)

    enrich = {}
    for i in range(0, len(uniq), args.batch):
        batch = uniq[i:i+args.batch]
        enrich.update(uniprot_fetch_json_batch(batch, retries=args.retries, sleep=args.sleep))
        time.sleep(args.sleep)

    out_cols = colnames + ["reviewed","xref_cazy","xref_esther"]
    with open(args.out, "w", newline='') as g:
        w = csv.DictWriter(g, delimiter="\t", fieldnames=out_cols)
        w.writeheader()
        for r in rows:
            e = enrich.get(r.get("acc",""), {})
            r2 = dict(r)
            r2["reviewed"] = e.get("reviewed","")
            r2["xref_cazy"] = e.get("xref_cazy","")
            r2["xref_esther"] = e.get("xref_esther","")
            w.writerow(r2)

if __name__ == "__main__":
    main()