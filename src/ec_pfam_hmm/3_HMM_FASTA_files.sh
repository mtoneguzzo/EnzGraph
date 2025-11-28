#!/usr/bin/env bash
set -u

# 3_HMM_FASTA_files.min.sh — fetch only per-PFAM HMMs + one concatenated subset
# Usage:
#   bash 3_HMM_FASTA_files.min.sh <TSV> <OUTPREFIX> <PFAM_DB>
# Example:
#   TSV="$DATA_DIR/EC$EC/EC$EC.enriched.tsv"
#   OUTPREFIX="$DATA_DIR/EC$EC/EC$EC"
#   PFAM_DB="$PFAM_A_DIR/Pfam-A.hmm"
# Add brew’s bin for Apple Silicon macOS
export PATH="/opt/homebrew/bin:$PATH"

which hmmfetch || echo "no hmmfetch"
which hmmpress || echo "no hmmpress"


TSV="${1:?TSV path missing}"
OUTPREFIX="${2:?OUTPREFIX missing}"   # e.g. /path/EC3.1.1.74/EC3.1.1.74
PFAM_DB="${3:?Pfam-A.hmm path missing}"

OUTDIR="$(dirname "$OUTPREFIX")"
BASE="$(basename "$OUTPREFIX")"
PFAM_OUT_DIR="${OUTDIR}/pfam_hmms"
mkdir -p "$PFAM_OUT_DIR"

log(){ printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
die(){ printf '[ERROR] %s\n' "$*" >&2; exit 2; }

[[ -f "$TSV" ]]     || die "TSV not found: $TSV"
[[ -f "$PFAM_DB" ]] || die "Pfam DB not found: $PFAM_DB"

# --- helpers ------------------------------------------------------------------
col_idx() {
  # return 1-based col index by exact header match (after trimming)
  local name="$1" file="$2"
  awk -F'\t' -v col="$name" '
    NR==1{
      for(i=1;i<=NF;i++){
        h=$i
        gsub(/\r/,"",h); sub(/^\xEF\xBB\xBF/,"",h); gsub(/^[ \t]+|[ \t]+$/,"",h)
        if(h==col){ print i; exit }
      }
    }' "$file"
}

pick_col() {
  # try a list of aliases; print first index found
  local file="$1"; shift
  local name idx
  for name in "$@"; do
    idx="$(col_idx "$name" "$file" || true)"
    if [[ -n "${idx:-}" ]]; then
      echo "$idx"
      return 0
    fi
  done
  return 1
}

dump_headers() {
  awk -F'\t' 'NR==1{
    for(i=1;i<=NF;i++){
      h=$i; gsub(/\r/,"",h); sub(/^\xEF\xBB\xBF/,"",h); gsub(/^[ \t]+|[ \t]+$/,"",h)
      printf("Header %02d: [%s]\n",i,h)
    }}' "$TSV" >&2
}

require_bin() { command -v "$1" >/dev/null 2>&1 || die "Missing binary: $1"; }

require_bin hmmfetch



# --- detect PFAM column only ---------------------------------------------------
PFAM_COL="$(pick_col "$TSV" "pfam_list" "Pfam" "pfam" "Pfam(s)" || true)"
[[ -n "${PFAM_COL:-}" ]] || { dump_headers; die "No PFAM column (pfam_list/Pfam/pfam)."; }


# --- collect bare PFAM IDs to a temp file -------------------------------------
PF_TMP="$(mktemp)"
trap 'rm -f "$PF_TMP" "$PF_TMP.ver" "$PF_TMP.cat"' EXIT

awk -F'\t' -v P="$PFAM_COL" 'NR>1 {print $P}' "$TSV" \
  | tr ';, ' '\n' \
  | sed -E 's/[^A-Za-z0-9\.]//g' \
  | awk 'NF' \
  | grep -E '^PF[0-9]{5}(\.[0-9]+)?$' \
  | sed -E 's/^([^\.]+).*$/\1/' \
  | sort -u > "$PF_TMP"

[[ -s "$PF_TMP" ]] || die "No PFAM IDs detected after parsing the PFAM column."
log "PFAM IDs: $(wc -l < "$PF_TMP")"

# --- map PFxxxxx -> latest versioned ACC present in DB ------------------------
# writes to $PF_TMP.ver (no permanent file)
> "$PF_TMP.ver"
while read -r PF; do
  ver=$(grep -n "ACC[[:space:]]\+$PF\.[0-9]\+" "$PFAM_DB" \
        | sed -E 's/^.*ACC[[:space:]]+//' \
        | awk '{print $1}' \
        | sort -V | tail -n1 || true)
  if [[ -n "$ver" ]]; then
    echo "$ver" >> "$PF_TMP.ver"
  else
    log "WARN: No versioned ACC found for $PF in Pfam DB; skipping."
  fi
done < "$PF_TMP"
sort -u -o "$PF_TMP.ver" "$PF_TMP.ver"

# --- fetch per-PFAM HMMs and build concatenated subset ------------------------
SUBSET_HMM="${OUTPREFIX}.Pfam-A.subset.hmm"
: > "$SUBSET_HMM"

n_ok=0
while read -r PFV; do
  out_one="${PFAM_OUT_DIR}/${PFV}.hmm"   # keep version in filename
  if hmmfetch "$PFAM_DB" "$PFV" > "$out_one" 2>/dev/null; then
    cat "$out_one" >> "$SUBSET_HMM"
    ((n_ok++))
  else
    log "WARN: missing $PFV in $PFAM_DB"
  fi
done < "$PF_TMP.ver"

[[ -s "$SUBSET_HMM" ]] || die "Subset HMM is empty; nothing fetched."
log "Fetched $n_ok HMM(s)."
log "Wrote:
  - per-PFAM HMMs in: $PFAM_OUT_DIR/
  - concatenated subset: $SUBSET_HMM"

# No hmmpress, no FASTA, no side files left behind