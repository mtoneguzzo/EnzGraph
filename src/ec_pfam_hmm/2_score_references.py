#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, csv, sys, re
from collections import Counter

# ----------------------------- Header mapping ---------------------------------

HEADER_MAP = {
    "Entry": "acc",
    "Entry Name": "entry",
    "Protein names": "protein",
    "Gene Names (primary)": "gene",
    "Gene Names": "gene",              # added for robustness
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
    # enrichment columns from step 1:
    "reviewed": "reviewed",
    "Reviewed": "reviewed",            # added for robustness
    "xref_cazy": "xref_cazy",
    "xref_esther": "xref_esther",
}

# ------------------------------- CLI ------------------------------------------

def parse_args():
    ap = argparse.ArgumentParser(
        description="Score UniProt records into Gold/Silver/Bronze (header-aware, lenient when no PFAMs learned)"
    )
    ap.add_argument("inp", help="TSV from step 1 (enriched) or raw (with header)")
    ap.add_argument("--target-ecs", nargs="+", required=True,
                    help="ECs of interest, e.g. 3.1.1.74 3.1.1.101 (exact; use multiple for sets)")
    ap.add_argument("--min-abs", type=int, default=2,
                    help="Min occurrences for Pfam to be considered 'expected' (default: 2)")
    ap.add_argument("--min-frac", type=float, default=0.10,
                    help="Min fraction for Pfam to be considered 'expected' (default: 0.10)")
    ap.add_argument("--write-full", default="",
                    help="Optional TSV to write all rows + Tier appended")
    ap.add_argument("--no-esther-proxy", action="store_true",
                    help="Disable ESTHER proxy (InterPro/Pfam) fallback when no ESTHER xref")
    return ap.parse_args()

# ------------------------------ IO helpers ------------------------------------

def normalize_header(h):
    h2 = (h or "").strip()
    if h2.lower() == "ntry":  # fix occasional truncated header
        return "Entry"
    return h2

def read_with_header(inp):
    with open(inp, newline='') as f:
        reader = csv.DictReader(f, delimiter="\t")
        norm_fieldnames = [normalize_header(h) for h in reader.fieldnames]
        to_internal = {raw: HEADER_MAP.get(raw, raw) for raw in norm_fieldnames}
        rows = []
        for r in reader:
            nr = {}
            for raw_key, val in r.items():
                nk = to_internal.get(normalize_header(raw_key), raw_key)
                nr[nk] = val
            rows.append(nr)
        return rows, list(to_internal.values())

# -------------------------------- Utils ---------------------------------------

def parse_semilist(x):
    x = (x or "").strip()
    # tolerate trailing semicolons/spaces/commas
    return [t for t in re.split(r"[;,\s]+", x) if t]

def has_secreted_eco269(subcell):
    s = subcell or ""
    return ("Secreted" in s) and ("ECO:0000269" in s)

def is_reviewed(row):
    # robust detection regardless of exact header
    return ((row.get("reviewed") or "").strip().lower() == "reviewed")

# ---------------------- Catalytic / evidence parsing --------------------------

def scan_catalytic(catalytic, target_ecs):
    """
    Parse CATALYTIC ACTIVITY to detect, per EC in target_ecs:
      - presence of RHEA id
      - presence of ECO:0000269
    Also set eco_any269 if 269 appears anywhere in the catalytic text.
    """
    out = {ec: {"rhea": False, "eco269": False, "eco_any269": False} for ec in target_ecs}
    if not catalytic:
        return out
    blocks = catalytic.split("CATALYTIC ACTIVITY:")
    for b in blocks:
        ecs_in_block = set(re.findall(r"EC=([\d\.]+)", b))
        rhea_here = bool(re.search(r"RHEA:\d+", b))
        ecos = re.findall(r"ECO:\d{7}", b)
        eco269_here = "ECO:0000269" in ecos
        for ec in (ecs_in_block & target_ecs):
            out[ec]["rhea"]   = out[ec]["rhea"] or rhea_here
            out[ec]["eco269"] = out[ec]["eco269"] or eco269_here
        if eco269_here:
            for ec in target_ecs:
                out[ec]["eco_any269"] = True
    return out

def high_conf_for_ec(row, ec):
    """
    High-confidence seed:
      - (ECO:0000269 and RHEA) for that EC, OR
      - reviewed AND (ECO:0000269 or RHEA) for that EC
    """
    reviewed_ok = is_reviewed(row)
    cat = scan_catalytic(row.get("catalytic",""), {ec})
    eco269 = cat[ec]["eco269"]; rhea = cat[ec]["rhea"]
    return (eco269 and rhea) or (reviewed_ok and (eco269 or rhea))

# -------------------- Learn expected PFAMs per EC -----------------------------

def derive_expected_pfams(rows, target_ecs, min_abs=2, min_frac=0.10):
    seeds = {ec: [] for ec in target_ecs}
    for r in rows:
        ecs = set(parse_semilist(r.get("ec_list","")))
        pfams = set(parse_semilist(r.get("pfam_list","")))
        for ec in (ecs & target_ecs):
            if high_conf_for_ec(r, ec):
                seeds[ec].extend(pfams)
    pfam_expected = {ec: set() for ec in target_ecs}
    for ec, lst in seeds.items():
        if not lst:
            continue
        cnt = Counter(lst); total = sum(cnt.values())
        keep = {pf for pf, n in cnt.items() if n >= min_abs and (n/total) >= min_frac}
        pfam_expected[ec] = keep
    return pfam_expected

# ---- CAZy normalization and ESTHER proxy helpers (fallback family signal) ----

def norm_token(s):
    """Uppercase + strip non-alphanumerics: 'CE-5'->'CE5', 'GH 13'->'GH13'."""
    return re.sub(r"[^A-Z0-9]+", "", (s or "").upper())

def esther_proxy(row):
    """
    Proxy for α/β-hydrolase / cutinase families when no ESTHER xref is present.
    """
    iprs = set(parse_semilist(row.get("interpro_list","")))
    pfams = set(parse_semilist(row.get("pfam_list","")))
    return (
        ("IPR029058" in iprs) or                 # Alpha/beta hydrolase superfamily
        ("PF01083" in pfams) or                  # Cutinase (fungal/bacterial)
        ("PF12740" in pfams)                     # Cutinase-like
    )

def cazy_esther_ok(row, use_esther_proxy=True):
    cx = row.get("xref_cazy","") or ""
    ex = row.get("xref_esther","") or ""
    toks = [norm_token(t) for t in re.split(r"[;,]\s*", cx) if t.strip()]
    # CAZy families indicating carbohydrate-active enzymes (GH, PL, CE, CBM)
    has_poly_cazy = any(t.startswith(("GH","PL","CE","CBM")) for t in toks)
    has_esther = bool(ex.strip())
    if use_esther_proxy:
        has_esther = has_esther or esther_proxy(row)
    return has_poly_cazy or has_esther

# ------------------------------ Tiering ---------------------------------------

def assign_tier(row, target_ecs, pfam_expected, use_esther_proxy=True):
    ecs   = set(parse_semilist(row.get("ec_list","")))
    pfams = set(parse_semilist(row.get("pfam_list","")))
    pdbs  = set(parse_semilist(row.get("pdb_list","")))
    subcell = row.get("subcellular","") or ""
    catalytic = row.get("catalytic","") or ""

    # filter by ECs of interest
    if not (ecs & target_ecs):
        return "skip"

    # Pfam expected? (lenient: if nothing was learned for any target EC, do not block)
    pfam_ok = False
    has_any_expected_defined = False
    for tec in (ecs & target_ecs):
        expected = pfam_expected.get(tec, set())
        if expected:
            has_any_expected_defined = True
            if pfams & expected:
                pfam_ok = True
                break

    # Fallback: CAZy/ESTHER (or ESTHER proxy) as acceptable family signal
    if not pfam_ok and cazy_esther_ok(row, use_esther_proxy=use_esther_proxy):
        pfam_ok = True

    # If no expectations were learned for any matching EC → do not penalize:
    if not has_any_expected_defined:
        pfam_ok = True

    # Catalytic evidence
    cat = scan_catalytic(catalytic, target_ecs)
    has_rhea_target = any(cat[ec]["rhea"] for ec in target_ecs)
    has_eco269_target = any(cat[ec]["eco269"] for ec in target_ecs)
    eco_any269 = any(cat[ec]["eco_any269"] for ec in target_ecs)
    has_pdb = len(pdbs) > 0
    secreted_exp = has_secreted_eco269(subcell)

    # Tiers
    is_gold = pfam_ok and has_eco269_target and (has_rhea_target or has_pdb or is_reviewed(row))
    # Allow silver if catalytic evidence exists even if family expectations were not met (pfam_ok),
    # but keep pfam_ok requirement to avoid too many silvers if you prefer stricter behavior.
    is_silver = pfam_ok and not is_gold and (has_rhea_target or has_pdb or eco_any269 or secreted_exp or is_reviewed(row))

    return "gold" if is_gold else ("silver" if is_silver else "bronze")

# -------------------------------- Main ----------------------------------------

def main():
    args = parse_args()
    target_ecs = set(args.target_ecs)

    rows, _cols = read_with_header(args.inp)

    # Learn expected Pfams from your own high-confidence seeds
    pfam_expected = derive_expected_pfams(rows, target_ecs, args.min_abs, args.min_frac)
    sys.stderr.write("Learned expected Pfams per EC (non-empty only):\n")
    any_learned = False
    for ec in sorted(target_ecs):
        vals = sorted(pfam_expected.get(ec, []))
        if vals:
            any_learned = True
            sys.stderr.write(f"  {ec}: {vals}\n")
    if not any_learned:
        sys.stderr.write("  (none learned; tiering will NOT penalize family expectations)\n")

    # Print compact table
    print("acc\torganism\tlen\tECs\tPfam\tPDB?\tSecreted(269)?\tTier")
    full_rows = []
    for r in rows:
        t = assign_tier(
            r, target_ecs, pfam_expected,
            use_esther_proxy=(not args.no_esther_proxy)
        )
        if t == "skip":
            continue
        pfams = ",".join(parse_semilist(r.get("pfam_list","")))
        pdb_flag = "yes" if parse_semilist(r.get("pdb_list","")) else "no"
        sec269 = "yes" if has_secreted_eco269(r.get("subcellular","")) else "no"
        print(f'{r.get("acc","")}\t{r.get("organism","")}\t{r.get("length","")}\t'
              f'{r.get("ec_list","")}\t{pfams}\t{pdb_flag}\t{sec269}\t{t}')
        if args.write_full:
            r2 = dict(r); r2["Tier"] = t; full_rows.append(r2)

    # Optional full dump with Tier appended
    if args.write_full and full_rows:
        fieldnames = list(full_rows[0].keys())
        with open(args.write_full, "w", newline='') as g:
            w = csv.DictWriter(g, delimiter="\t", fieldnames=fieldnames)
            w.writeheader()
            for r in full_rows:
                w.writerow(r)

if __name__ == "__main__":
    main()