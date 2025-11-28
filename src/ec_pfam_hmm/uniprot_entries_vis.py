#!/usr/bin/env python3
import os
import re
import sys
import json
import argparse
import textwrap
import subprocess
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")  # render plots to files (no GUI)
import matplotlib.pyplot as plt

def parse_args():
    p = argparse.ArgumentParser(
        description="Summarize UniProt reviewed entries for an EC number; visualize PFAM/organisms/length and optionally download Pfam HMMs."
    )
    p.add_argument("--ec", default=os.environ.get("EC", "3.1.1.74"),
                   help="EC number (default: env EC or 3.1.1.74)")
    p.add_argument("--workdir", default=os.environ.get("WORKDIR", "."),
                   help="Working directory containing EC<EC>/ (default: current dir)")
    p.add_argument("--top-n-pfams", type=int, default=20,
                   help="Top N PFAMs to display (default: 20)")
    p.add_argument("--download-hmms", action="store_true",
                   help="Download Pfam HMMs found in the set")
    p.add_argument("--no-download-hmms", dest="download_hmms", action="store_false")
    p.set_defaults(download_hmms=True)
    p.add_argument("--hmmpress", action="store_true",
                   help="Run hmmpress on the combined HMM file after download")
    return p.parse_args()

def parse_pfams(x):
    if pd.isna(x):
        return []
    return re.findall(r"PF\d{5}", str(x))

def wrap_labels(labels, width=30):
    return ["\n".join(textwrap.wrap(str(x), width=width)) for x in labels]

def main():
    args = parse_args()

    ec = args.ec
    workdir = os.path.abspath(args.workdir)
    indir = os.path.join(workdir, f"EC{ec}")
    tsv = os.path.join(indir, f"EC{ec}.reviewed.tsv")
    fa  = os.path.join(indir, f"EC{ec}.reviewed.fa")
    hmm_outdir = os.path.join(indir, "pfam_hmms")

    if not os.path.exists(tsv):
        sys.exit(f"[error] TSV not found: {tsv}\nRun the UniProt fetch step first.")

    # Load TSV
    t = pd.read_csv(tsv, sep="\t")
    t = t.rename(columns=lambda c: c.strip())

    # Detect Pfam column
    pfam_col_candidates = [c for c in t.columns if "pfam" in c.lower()]
    pfam_col = pfam_col_candidates[0] if pfam_col_candidates else None
    if pfam_col is None:
        print("[warn] No Pfam column found in TSV; PFAM plots will be skipped.")
        t["Pfam_list"] = [[]]*len(t)
    else:
        t["Pfam_list"] = t[pfam_col].apply(parse_pfams)

    print(f"[info] EC={ec}")
    print(f"[info] Sequences: {len(t)}")
    print(f"[info] Columns: {list(t.columns)}")
    os.makedirs(indir, exist_ok=True)

    # ===== Visual 1: PFAM distribution =====
    pf = t["Pfam_list"].explode().dropna()
    pfam_counts = None
    if not pf.empty:
        pfam_counts = pf.value_counts()
        top = pfam_counts.head(args.top_n_pfams)
        plt.figure(figsize=(8,4))
        top.plot(kind="bar")
        plt.xticks(rotation=45, ha="right")
        plt.ylabel("Count")
        plt.title(f"Top Pfam domains in reviewed EC {ec}")
        plt.tight_layout()
        pf_png = os.path.join(indir, f"EC{ec}.pfam_top.png")
        plt.savefig(pf_png, dpi=150)
        plt.close()
        print(f"[save] {pf_png}")

    # ===== Visual 2: Organism distribution (top 20) =====
    org_col = [c for c in t.columns if "organism name" in c.lower()]
    if org_col:
        org_counts = t[org_col[0]].value_counts().head(20)
        if not org_counts.empty:
            plt.figure(figsize=(8,4))
            org_counts.plot(kind="bar")
            plt.xticks(rotation=45, ha="right")
            plt.ylabel("Entries")
            plt.title(f"Top organisms in EC {ec} (reviewed)")
            plt.tight_layout()
            org_png = os.path.join(indir, f"EC{ec}.organisms_top.png")
            plt.savefig(org_png, dpi=150)
            plt.close()
            print(f"[save] {org_png}")

    # ===== Visual 3: Length distribution =====
    length_col = None
    for c in t.columns:
        if c.lower() == "length":
            length_col = c
            break

    if length_col:
        plt.figure(figsize=(6,3.5))
        t[length_col].plot(kind="hist", bins=30)
        plt.xlabel("Protein length (aa)")
        plt.ylabel("Count")
        plt.title(f"Length distribution — EC {ec} reviewed")
        plt.tight_layout()
        len_png = os.path.join(indir, f"EC{ec}.length_hist.png")
        plt.savefig(len_png, dpi=150)
        plt.close()
        print(f"[save] {len_png}")

    # # ===== Optional: download Pfam HMMs =====
    # combined_path = None
    # if args.download_hmms and pfam_counts is not None:
    #     uniq_pfams = sorted(set(pf.dropna().tolist()))
    #     if uniq_pfams:
    #         os.makedirs(hmm_outdir, exist_ok=True)
    #         print(f"[info] Downloading {len(uniq_pfams)} Pfam HMMs → {hmm_outdir}")
    #         for pfid in uniq_pfams:
    #             url = f"https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/current/individual_families/{pfid}.hmm.gz"
    #             out_gz = os.path.join(hmm_outdir, f"{pfid}.hmm.gz")
    #             out_hmm = os.path.join(hmm_outdir, f"{pfid}.hmm")
    #             if os.path.exists(out_hmm):
    #                 continue
    #             # download
    #             subprocess.run(["curl", "-sL", url, "-o", out_gz], check=False)
    #             if os.path.exists(out_gz) and os.path.getsize(out_gz) > 0:
    #                 subprocess.run(["gunzip", "-f", out_gz], check=False)
    #             else:
    #                 print(f"[warn] Could not fetch {pfid} (empty download).")
    #                 try: os.remove(out_gz)
    #                 except: pass

    #         # Concatenate into one HMM file
    #         combined_path = os.path.join(hmm_outdir, f"EC{ec}.Pfam-A.subset.hmm")
    #         with open(combined_path, "w") as w:
    #             for pfid in uniq_pfams:
    #                 p = os.path.join(hmm_outdir, f"{pfid}.hmm")
    #                 if os.path.exists(p):
    #                     with open(p) as f:
    #                         w.write(f.read().rstrip() + "\n")
    #         print(f"[save] Combined HMMs: {combined_path}")

    #         if args.hmmpress and os.path.exists(combined_path):
    #             print("[info] Running hmmpress ...")
    #             subprocess.run(["hmmpress", combined_path], check=False)

    # ===== Save a small JSON summary =====
    summary = {
        "ec": ec,
        "n_sequences": int(len(t)),
        "has_pfam": bool(pfam_counts is not None),
        "top_pfams": pfam_counts.head(20).to_dict() if pfam_counts is not None else {},
        "plots": {
            "pfam_top": os.path.basename(pf_png) if pfam_counts is not None else None,
            "organisms_top": os.path.basename(org_png) if org_col else None,
            "length_hist": os.path.basename(len_png) if length_col else None,
        },
        #"combined_hmm": os.path.relpath(combined_path, indir) if combined_path else None,
        "indir": indir
    }
    with open(os.path.join(indir, f"EC{ec}.summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[save] {os.path.join(indir, f'EC{ec}.summary.json')}")

if __name__ == "__main__":
    main()