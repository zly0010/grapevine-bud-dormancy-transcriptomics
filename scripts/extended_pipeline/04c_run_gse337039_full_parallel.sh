#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "Usage: $0 <work_dir> <manifest.tsv> <mapping.tsv> <results_dir>" >&2
  exit 2
fi

WORK_DIR="$(realpath "$1")"
MANIFEST="$(realpath "$2")"
MAPPING="$(realpath "$3")"
RESULTS_DIR="$4"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DOWNLOAD_JOBS="${DOWNLOAD_JOBS:-1}"
ALIGN_JOBS="${ALIGN_JOBS:-2}"
SAMPLE_THREADS="${SAMPLE_THREADS:-3}"

if [[ "$ALIGN_JOBS" != "2" ]]; then
  echo "This runner currently requires ALIGN_JOBS=2" >&2
  exit 3
fi

SRA_DIR="$WORK_DIR/star_sra"
COUNT_DIR="$WORK_DIR/star_counts"
LOG_DIR="$WORK_DIR/logs_star"
META_DIR="$WORK_DIR/metadata_star"
mkdir -p "$SRA_DIR" "$COUNT_DIR" "$LOG_DIR" "$META_DIR" "$RESULTS_DIR"

download_run() {
  local run="$1"
  local count_file="$COUNT_DIR/$run.forward.txt"
  local run_dir="$SRA_DIR/$run"
  local sra_file="$run_dir/$run.sra"
  local part_file="$sra_file.part"
  local url="https://sra-pub-run-odp.s3.amazonaws.com/sra/$run/$run"
  if [[ -s "$count_file" ]]; then
    return 0
  fi
  mkdir -p "$run_dir"
  if [[ -s "$sra_file" ]]; then
    if vdb-validate "$sra_file" > "$LOG_DIR/${run}_validate_existing.log" 2>&1; then
      return 0
    fi
    rm -f "$sra_file"
  fi
  echo "DOWNLOAD_START $run $(date --iso-8601=seconds)"
  if ! aria2c --continue=true --max-connection-per-server=16 --split=16 \
    --min-split-size=1M --file-allocation=none --max-tries=5 --retry-wait=5 \
    --dir "$run_dir" --out "$run.sra.part" "$url" \
    > "$LOG_DIR/${run}_download.log" 2>&1; then
    echo "Download failed: $run" >&2
    return 1
  fi
  mv "$part_file" "$sra_file"
  if ! vdb-validate "$sra_file" > "$LOG_DIR/${run}_validate.log" 2>&1; then
    rm -f "$sra_file"
    echo "SRA validation failed: $run" >&2
    return 1
  fi
  rm -f "$run_dir/$run.sra.tmp" "$run_dir/$run.sra.lock" "$run_dir/$run.sra.prf"
  echo "DOWNLOAD_DONE $run $(date --iso-8601=seconds)"
}
export -f download_run
export SRA_DIR COUNT_DIR LOG_DIR

echo "Downloading missing SRA files with $DOWNLOAD_JOBS workers"
tail -n +2 "$MANIFEST" | cut -f1 | tr -d '"\r' | grep -E '^SRR[0-9]+$' | \
  xargs -r -n1 -P "$DOWNLOAD_JOBS" bash -c 'download_run "$1"' _

PART1="$META_DIR/GSE337039_manifest_part1.tsv"
PART2="$META_DIR/GSE337039_manifest_part2.tsv"
awk -v p1="$PART1" -v p2="$PART2" '
  NR == 1 { print > p1; print > p2; next }
  NR % 2 == 0 { print > p1; next }
  { print > p2 }
' "$MANIFEST"

echo "Running two STAR/featureCounts workers"
set +e
RUN_FINAL_QC=0 RUN_LABEL=part1 SAMPLE_THREADS="$SAMPLE_THREADS" \
  bash "$SCRIPT_DIR/04b_rescue_gse337039_star.sh" "$WORK_DIR" "$PART1" \
  > "$LOG_DIR/worker_part1.log" 2>&1 &
PID1=$!
RUN_FINAL_QC=0 RUN_LABEL=part2 SAMPLE_THREADS="$SAMPLE_THREADS" \
  bash "$SCRIPT_DIR/04b_rescue_gse337039_star.sh" "$WORK_DIR" "$PART2" \
  > "$LOG_DIR/worker_part2.log" 2>&1 &
PID2=$!
wait "$PID1"; STATUS1=$?
wait "$PID2"; STATUS2=$?
set -e
if [[ "$STATUS1" -ne 0 || "$STATUS2" -ne 0 ]]; then
  echo "STAR workers failed: part1=$STATUS1 part2=$STATUS2" >&2
  exit 4
fi

FORWARD_COUNT=$(find "$COUNT_DIR" -maxdepth 1 -type f -name '*.forward.txt' | wc -l)
if [[ "$FORWARD_COUNT" -ne 60 ]]; then
  echo "Expected 60 forward count files, found $FORWARD_COUNT" >&2
  exit 5
fi

multiqc "$WORK_DIR/qc/fastp_star" "$WORK_DIR/star_align" "$COUNT_DIR" \
  -o "$WORK_DIR/qc/multiqc_star" --force 2>&1 | tee "$LOG_DIR/multiqc.log"

Rscript "$SCRIPT_DIR/05b_gse337039_featurecounts_deseq2.R" \
  "$WORK_DIR" "$MANIFEST" "$MAPPING" "$RESULTS_DIR" auto \
  2>&1 | tee "$LOG_DIR/deseq2_featurecounts.log"

echo "GSE337039_FULL_PARALLEL_OK"
