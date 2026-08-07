#!/usr/bin/env bash
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
#
# CD19 Bulk RNA-seq Salmon quantification
#
# Re-runnable: note, if a quant.sf already exists, it's effectively cached
#
# Usage:  JOBS=3 THREADS=6 bash bulk/02_quantify.sh
#         BENCH=4 bash bulk/02_quantify.sh
# ==============================================================================
set -uo pipefail

### WSL Ubuntu 24.04.3 LTS
# conda activate rrms
# cd /mnt/c/Users/devin/Desktop/rrms
# BENCH=6 JOBS=6 THREADS=3 bash bulk/02_quantify.sh
# nohup env JOBS=6 THREADS=3 bash bulk/02_quantify.sh > ~/02_quantify.out 2>&1 & echo "PID $!"

FASTQ_DIR="${FASTQ_DIR:-bulk/fastq}"
INDEX="${INDEX:-bulk/reference/salmon_index/idx}"
OUT_DIR="${OUT_DIR:-zenodo/bulk/counts}"
JOBS="${JOBS:-3}"          # concurrent samples
THREADS="${THREADS:-6}"    # threads per sample
N_BOOT="${N_BOOT:-0}"
CHECKSUM="${CHECKSUM:-1}"
BENCH="${BENCH:-0}"        # >0: quantify only this many samples, then stop

[[ -d "$INDEX" ]] || { echo "no index at $INDEX" >&2; exit 1; }
[[ -d "$FASTQ_DIR" ]] || { echo "no $FASTQ_DIR" >&2; exit 1; }
command -v salmon >/dev/null || { echo "salmon not on PATH" >&2; exit 1; }
mkdir -p "$OUT_DIR"

MANIFEST="$OUT_DIR/manifest.tsv"
LOG="$OUT_DIR/02_quantify.log"
PARTS="$OUT_DIR/.manifest.d"
rm -rf "$PARTS"; mkdir -p "$PARTS"
: > "$LOG"

NCORE=$(nproc 2>/dev/null || echo 0)
SALMON_V=$(salmon --version 2>&1 | tr -d '\n')
{
  echo "salmon    : $SALMON_V"
  echo "cores     : $NCORE     JOBS=$JOBS x THREADS=$THREADS = $((JOBS * THREADS))"
  echo "bootstraps: $N_BOOT    checksums: $CHECKSUM"
  echo "index     : $INDEX"
  echo "fastq     : $FASTQ_DIR"
  echo "out       : $OUT_DIR"
} | tee -a "$LOG"
if [[ "$NCORE" -gt 0 && $((JOBS * THREADS)) -gt "$NCORE" ]]; then
  echo "WARNING: JOBS x THREADS ($((JOBS * THREADS))) exceeds $NCORE cores ... they will contend." | tee -a "$LOG"
fi


shopt -s nullglob
FQ_FILES=( "$FASTQ_DIR"/*.fq.gz "$FASTQ_DIR"/*.fastq.gz )
shopt -u nullglob
if (( ${#FQ_FILES[@]} == 0 )); then
  echo "no FASTQ files in $FASTQ_DIR (looked for *.fq.gz and *.fastq.gz)" | tee -a "$LOG" >&2
  exit 1
fi
mapfile -t SAMPLES < <(printf '%s\n' "${FQ_FILES[@]##*/}" \
  | sed -E 's/_S[0-9]+_L[0-9]{3}_R[12]_001\.f(ast)?q\.gz$//; s/_(R?[12])\.f(ast)?q\.gz$//' \
  | sort -u)
echo "fastq files: ${#FQ_FILES[@]}    samples: ${#SAMPLES[@]}" | tee -a "$LOG"


# Resolve read pairs 
TASKS="$OUT_DIR/.tasks.tsv"; : > "$TASKS"
n_skip=0
for SID in "${SAMPLES[@]}"; do
  R1=""; R2=""
  for pat in "_R1.fastq.gz:_R2.fastq.gz" "_R1.fq.gz:_R2.fq.gz" \
             "_1.fastq.gz:_2.fastq.gz"   "_1.fq.gz:_2.fq.gz"; do
    a="${pat%%:*}"; b="${pat##*:}"
    if [[ -s "$FASTQ_DIR/${SID}${a}" && -s "$FASTQ_DIR/${SID}${b}" ]]; then
      R1="$FASTQ_DIR/${SID}${a}"; R2="$FASTQ_DIR/${SID}${b}"; break
    fi
  done
  # bcl2fastq style: <Sample_ID>_S#_L00#_R1_001.fastq.gz. Exactly one match is
  # required ; several means the sample spans lanes and must be merged first
  if [[ -z "$R1" ]]; then
    shopt -s nullglob
    cand=( "$FASTQ_DIR/${SID}"_S[0-9]*_L[0-9][0-9][0-9]_R1_001.f*q.gz )
    shopt -u nullglob
    if (( ${#cand[@]} == 1 )); then
      R1="${cand[0]}"; R2="${R1/_R1_001./_R2_001.}"
      [[ -s "$R2" ]] || { echo "SKIP $SID : R1 present but no matching R2" | tee -a "$LOG"; R1=""; }
    elif (( ${#cand[@]} > 1 )); then
      echo "SKIP $SID : ${#cand[@]} R1 files match (spans lanes?); merge first" | tee -a "$LOG"
    fi
  fi
  if [[ -z "$R1" ]]; then
    echo "SKIP $SID ; no read pair found in $FASTQ_DIR" | tee -a "$LOG"
    printf '%s\tNA\tNA\tNA\tNA\t%s\tNA\tmissing_fastq\n' "$SID" "$SALMON_V" > "$PARTS/$SID.tsv"
    n_skip=$((n_skip + 1)); continue
  fi
  printf '%s\t%s\t%s\n' "$SID" "$R1" "$R2" >> "$TASKS"
done
n_task=$(wc -l < "$TASKS")
echo "resolved: $n_task read pairs, $n_skip without one" | tee -a "$LOG"

if [[ "$BENCH" -gt 0 ]]; then
  head -n "$BENCH" "$TASKS" > "$TASKS.bench" && mv "$TASKS.bench" "$TASKS"
  n_task=$(wc -l < "$TASKS")
  echo "BENCH=$BENCH : quantifying only $n_task sample(s) to time the settings." | tee -a "$LOG"
fi

# ---- single sample -------------------------------------------------------------
quant_one() {
  local SID="$1" R1="$2" R2="$3"
  local d="$OUT_DIR/$SID" t0 t1 secs status
  if [[ -s "$d/quant.sf" ]]; then
    echo "have $SID ; skipping (delete $d to redo)"
    status=cached; secs=0
  else
    t0=$(date +%s)
    local boot=(); [[ "$N_BOOT" -gt 0 ]] && boot=(--numBootstraps "$N_BOOT")
    if salmon quant -i "$INDEX" -l A -1 "$R1" -2 "$R2" \
         -p "$THREADS" --gcBias --seqBias "${boot[@]}" \
         -o "$d" > "$PARTS/$SID.salmon.log" 2>&1; then
      status=ok
    else
      status=failed
      echo "FAILED $SID ; see $PARTS/$SID.salmon.log"
    fi
    t1=$(date +%s); secs=$((t1 - t0))
    echo "done $SID  ${secs}s  $status"
  fi
  local m1=NA m2=NA
  if [[ "$CHECKSUM" == "1" && "$status" != "failed" ]]; then
    if command -v md5sum >/dev/null; then m1=$(md5sum "$R1" | cut -d' ' -f1); m2=$(md5sum "$R2" | cut -d' ' -f1)
    else m1=$(shasum -a 256 "$R1" | cut -d' ' -f1); m2=$(shasum -a 256 "$R2" | cut -d' ' -f1); fi
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$SID" "$R1" "$R2" "$m1" "$m2" "$SALMON_V" "$secs" "$status" > "$PARTS/$SID.tsv"
}
export -f quant_one
export OUT_DIR INDEX THREADS N_BOOT CHECKSUM PARTS SALMON_V

# ---- run --------------------------------------------------------------------
echo | tee -a "$LOG"
echo "quantifying $n_task sample(s), $JOBS at a time ..." | tee -a "$LOG"
T0=$(date +%s)
# -I{} keeps one line per invocation
xargs -a "$TASKS" -d '\n' -P "$JOBS" -I{} bash -c \
  'IFS=$'"'"'\t'"'"' read -r s a b <<< "{}"; quant_one "$s" "$a" "$b"' 2>&1 | tee -a "$LOG"
T1=$(date +%s)

printf 'Sample_ID\tR1\tR2\tR1_md5\tR2_md5\tsalmon_version\tseconds\tstatus\n' > "$MANIFEST"
for SID in "${SAMPLES[@]}"; do [[ -s "$PARTS/$SID.tsv" ]] && cat "$PARTS/$SID.tsv" >> "$MANIFEST"; done
rm -f "$TASKS"

n_ok=$(awk -F'\t' 'NR>1 && ($8=="ok"||$8=="cached")' "$MANIFEST" | wc -l)
n_bad=$(awk -F'\t' 'NR>1 && $8!="ok" && $8!="cached"' "$MANIFEST" | wc -l)
{
  echo
  echo "quantified: $n_ok    skipped/failed: $n_bad"
  echo "wall clock: $(( (T1-T0)/60 )) min $(( (T1-T0)%60 )) s"
  awk -F'\t' 'NR>1 && $7 ~ /^[0-9]+$/ && $7>0 {s+=$7; n++} END{if(n) printf "per-sample CPU time: mean %.1f min over %d samples\n", s/n/60, n}' "$MANIFEST"
  echo "manifest -> $MANIFEST"
} | tee -a "$LOG"
