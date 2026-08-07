#!/usr/bin/env bash
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
#
# Build a decoy aware Salmon index
#
# Requires: salmon, gffread
#
# Usage:
#         REF_DIR=/path/to/refdata-gex-GRCh38-2020-A
#         bash bulk/01_build_index.sh
#
#         FORCE=1 bash bulk/01_build_index.sh    # rebuild over a finished index
# ==============================================================================
set -euo pipefail

# WSL Ubuntu 24.04.3 LTS
# conda activate rrms

REF_DIR="${REF_DIR:-bulk/reference/refdata-gex-GRCh38-2020-A}"
OUT_DIR="${OUT_DIR:-bulk/reference/salmon_index}"
THREADS="${THREADS:-8}"
KMER="${KMER:-31}"          # 31 is good for reads >= 75 bp; these are 150 bp

GENOME="$REF_DIR/fasta/genome.fa"
GTF="$REF_DIR/genes/genes.gtf"

mkdir -p "$OUT_DIR"
VER="$OUT_DIR/versions.txt"

#
if [[ -s "$OUT_DIR/idx/info.json" && "${FORCE:-0}" != "1" ]]; then
  echo "index already present at $OUT_DIR/idx - nothing to do."
  echo "  versions.txt left untouched; it records how that index was built."
  echo "  To rebuild:  rm -rf $OUT_DIR/idx     (or run with FORCE=1)"
  echo
  exit 0
fi

for f in "$GENOME" "$GTF"; do
  [[ -s "$f" ]] || { echo "missing: $f" >&2
    echo "Point REF_DIR at an unpacked refdata-gex-GRCh38-2020-A." >&2; exit 1; }
done

# Any .part file is the debris of an interrupted run. Never reuse it.
rm -f "$OUT_DIR"/*.part

: > "$VER"
{ echo "built: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host : $(uname -a)"
  echo "ref  : $REF_DIR"
  echo "kmer : $KMER"
  echo "--- salmon ---";  salmon --version 2>&1 || true
  echo "--- gffread ---"; gffread --version 2>&1 || true
} >> "$VER"
echo "versions -> $VER"

# ---- preflight --------------
if [[ "${PREFLIGHT:-1}" == "1" ]]; then
  echo "preflight: toy index ..."
  PF=$(mktemp -d)
  printf '>t1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n' >  "$PF/toy.fa"
  printf '>t2\nTTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA\n' >> "$PF/toy.fa"
  if ! salmon index -t "$PF/toy.fa" -i "$PF/idx" -k 11 -p 1 > "$PF/log" 2>&1; then
    echo "PREFLIGHT FAILED - salmon cannot build even a toy index here." >&2
    echo "--- salmon output ---" >&2; tail -20 "$PF/log" >&2
    echo "--- environment ---" >&2
    echo "salmon  : $(command -v salmon)" >&2
    echo "version : $(salmon --version 2>&1)" >&2
    echo "LANG    : ${LANG-<unset>}   LC_ALL: ${LC_ALL-<unset>}   LOCPATH: ${LOCPATH-<unset>}" >&2
    echo "locales : $(locale -a 2>/dev/null | tr '\n' ' ')" >&2
    rm -rf "$PF"; exit 1
  fi
  rm -rf "$PF"
  echo "preflight: OK"
fi

if [[ -s "$OUT_DIR/transcripts.fa" && -n "$(tail -c 1 "$OUT_DIR/transcripts.fa")" ]]; then
  echo "transcripts.fa does not end in a newline - treating as truncated" >&2
  mv -f "$OUT_DIR/transcripts.fa" "$OUT_DIR/transcripts.fa.suspect"
  [[ -e "$OUT_DIR/gentrome.fa.gz" ]] && mv -f "$OUT_DIR/gentrome.fa.gz" "$OUT_DIR/gentrome.fa.gz.suspect"
  echo "  renamed to *.suspect; both will be rebuilt" >&2
fi
if [[ -s "$OUT_DIR/gentrome.fa.gz" ]]; then
  echo "checking gentrome.fa.gz integrity ..."
  if ! gzip -t "$OUT_DIR/gentrome.fa.gz" 2>/dev/null; then
    echo "gentrome.fa.gz fails gzip -t - treating as truncated" >&2
    mv -f "$OUT_DIR/gentrome.fa.gz" "$OUT_DIR/gentrome.fa.gz.suspect"
    echo "  renamed to gentrome.fa.gz.suspect; it will be rebuilt" >&2
  fi
fi

# ---- transcript FASTA from the SAME GTF the counts will be summarized on -----
if [[ ! -s "$OUT_DIR/transcripts.fa" ]]; then
  echo "extracting transcripts ..."
  gffread -w "$OUT_DIR/transcripts.fa.part" -g "$GENOME" "$GTF"
  mv -f "$OUT_DIR/transcripts.fa.part" "$OUT_DIR/transcripts.fa"
fi

# ---- decoys: every genome sequence name -------------------------------------
if [[ ! -s "$OUT_DIR/decoys.txt" ]]; then
  echo "collecting decoy names ..."
  grep '^>' "$GENOME" | cut -d ' ' -f 1 | sed 's/^>//' > "$OUT_DIR/decoys.txt.part"
  mv -f "$OUT_DIR/decoys.txt.part" "$OUT_DIR/decoys.txt"
fi

# ---- gentrome: transcripts FIRST, then the genome ---------------------------
if [[ ! -s "$OUT_DIR/gentrome.fa.gz" ]]; then
  echo "assembling gentrome ..."
  cat "$OUT_DIR/transcripts.fa" "$GENOME" | gzip -c > "$OUT_DIR/gentrome.fa.gz.part"
  mv -f "$OUT_DIR/gentrome.fa.gz.part" "$OUT_DIR/gentrome.fa.gz"
fi

echo "indexing (this needs ~16 GB RAM and a while) ..."
salmon index \
  -t "$OUT_DIR/gentrome.fa.gz" \
  -d "$OUT_DIR/decoys.txt" \
  -i "$OUT_DIR/idx" \
  -k "$KMER" \
  -p "$THREADS" \
  --gencode 2>&1 | tee -a "$VER"

echo
echo "index -> $OUT_DIR/idx"
echo
