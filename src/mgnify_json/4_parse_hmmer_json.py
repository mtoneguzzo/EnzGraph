#!/usr/bin/env python3
import json
import sys
import csv

# ---- helpers ----
def _join_pairs(pairs):
    """Pairs like [["ERP1","url1"], ["ERP2","url2"]] -> 'ERP1|url1;ERP2|url2'."""
    try:
        return ";".join([f"{a}|{b}" for a, b in pairs])
    except Exception:
        return ""

def _join_list(x):
    """List of scalars -> 'a;b;c'. If not a list, best-effort str()."""
    if isinstance(x, list):
        return ";".join(str(v) for v in x)
    return str(x) if x is not None else ""

def _safe_float(x):
    try:
        return float(x)
    except Exception:
        return None

def _safe_int(x):
    try:
        return int(x)
    except Exception:
        return None

def _get(d, k, default=None):
    return d.get(k, default) if isinstance(d, dict) else default

def flatten_hit(hit):
    """
    Yield one row per domain for a given hit.
    If no domains, yield a single row with n_domains=0.
    """
    rows = []

    # ---- hit-level fields ----
    hit_acc   = _get(hit, "acc", "")
    hit_name  = _get(hit, "name", "")
    hit_score = _get(hit, "score", "")
    hit_bias  = _get(hit, "bias", "")
    hit_pval  = _get(hit, "pvalue", "")
    hit_eval  = _get(hit, "evalue", "")
    ndom      = _get(hit, "ndom", _get(hit, "nreported", 0))
    extlink   = _get(hit, "extlink", "")

    studies    = _join_pairs(_get(hit, "studies", []))
    assemblies = _join_pairs(_get(hit, "assemblies", []))
    biomes     = _join_list(_get(hit, "biome", []))
    samples    = _join_list(_get(hit, "samples", []))

    domains = _get(hit, "domains", [])
    if not isinstance(domains, list) or len(domains) == 0:
        # Emit one empty-domain row so the hit isn't lost
        rows.append({
            "hit_acc": hit_acc,
            "hit_name": hit_name,
            "hit_score": hit_score,
            "hit_bias": hit_bias,
            "hit_pvalue": hit_pval,
            "hit_evalue": hit_eval,
            "n_domains": _safe_int(ndom) or 0,

            "dom_bitscore": "",
            "dom_bias": "",
            "dom_ievalue": "",
            "dom_cevalue": "",
            "dom_aliL": "",
            "dom_aliId": "",
            "dom_aliSim": "",
            "dom_hmm_from": "",
            "dom_hmm_to": "",
            "dom_hmm_acc": "",
            "dom_hmm_name": "",
            "dom_sq_from": "",
            "dom_sq_to": "",
            "dom_model_cov": "",
            "dom_seq_cov": "",

            "studies": studies,
            "assemblies": assemblies,
            "biomes": biomes,
            "samples": samples,
            "extlink": extlink
        })
        return rows

    for d in domains:
        bitscore = _get(d, "bitscore", "")
        dbias    = _get(d, "bias", "")
        ievalue  = _get(d, "ievalue", "")
        cevalue  = _get(d, "cevalue", "")
        aliL     = _get(d, "aliL", "")
        aliId    = _get(d, "aliId", "")
        aliSim   = _get(d, "aliSim", "")

        hmm_from = _get(d, "alihmmfrom", _get(d, "alihmm_from", ""))
        hmm_to   = _get(d, "alihmmto", _get(d, "alihmm_to", ""))
        hmm_acc  = _get(d, "alihmmacc", _get(d, "alihmm_acc", ""))
        hmm_name = _get(d, "alihmmname", _get(d, "alihmm_name", ""))

        sq_from  = _get(d, "alisqfrom", _get(d, "alisq_from", ""))
        sq_to    = _get(d, "alisqto", _get(d, "alisq_to", ""))

        # Coverage (best effort — JSON does not carry full model length)
        try:
            hmm_span = (int(hmm_to) - int(hmm_from) + 1) if hmm_from and hmm_to else None
        except Exception:
            hmm_span = None
        try:
            seq_span = (int(sq_to) - int(sq_from) + 1) if sq_from and sq_to else None
        except Exception:
            seq_span = None

        # We expose spans directly as coverage surrogates (0..1 would require model length)
        dom_model_cov = hmm_span
        dom_seq_cov   = seq_span

        rows.append({
            "hit_acc": hit_acc,
            "hit_name": hit_name,
            "hit_score": hit_score,
            "hit_bias": hit_bias,
            "hit_pvalue": hit_pval,
            "hit_evalue": hit_eval,
            "n_domains": _safe_int(ndom) or 0,

            "dom_bitscore": bitscore,
            "dom_bias": dbias,
            "dom_ievalue": ievalue,
            "dom_cevalue": cevalue,
            "dom_aliL": aliL,
            "dom_aliId": aliId,
            "dom_aliSim": aliSim,
            "dom_hmm_from": hmm_from,
            "dom_hmm_to": hmm_to,
            "dom_hmm_acc": hmm_acc,
            "dom_hmm_name": hmm_name,
            "dom_sq_from": sq_from,
            "dom_sq_to": sq_to,
            "dom_model_cov": dom_model_cov,
            "dom_seq_cov": dom_seq_cov,

            "studies": studies,
            "assemblies": assemblies,
            "biomes": biomes,
            "samples": samples,
            "extlink": extlink
        })
    return rows

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.json> <output.csv>", file=sys.stderr)
        sys.exit(2)

    in_json, out_csv = sys.argv[1], sys.argv[2]

    with open(in_json, "r") as f:
        data = json.load(f)

    # Defensive extraction of hits
    hits = []
    if isinstance(data, dict):
        res = _get(data, "results", {})
        hits = _get(res, "hits", [])
    elif isinstance(data, list):
        # Very rare: results already a list of hits
        hits = data

    if not isinstance(hits, list):
        print("[ERROR] JSON does not contain a hits list under results.hits", file=sys.stderr)
        sys.exit(1)

    rows = []
    for h in hits:
        if not isinstance(h, dict):
            # Skip unexpected types defensively
            continue
        rows.extend(flatten_hit(h))

    # Write CSV
    fieldnames = [
        "hit_acc","hit_name","hit_score","hit_bias","hit_pvalue","hit_evalue","n_domains",
        "dom_bitscore","dom_bias","dom_ievalue","dom_cevalue","dom_aliL","dom_aliId","dom_aliSim",
        "dom_hmm_from","dom_hmm_to","dom_hmm_acc","dom_hmm_name","dom_sq_from","dom_sq_to",
        "dom_model_cov","dom_seq_cov",
        "studies","assemblies","biomes","samples","extlink"
    ]
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

if __name__ == "__main__":
    main()