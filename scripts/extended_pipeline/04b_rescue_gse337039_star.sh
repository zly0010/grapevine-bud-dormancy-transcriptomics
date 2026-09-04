#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <work_dir> <GSE337039_sample_manifest.tsv>" >&2
  exit 2
fi

WORK_DIR="$(realpath "$1")"
MANIFEST="$(realpath "$2")"
THREADS="${THREADS:-6}"
SAMPLE_THREADS="${SAMPLE_THREADS:-4}"
RUN_FINAL_QC="${RUN_FINAL_QC:-1}"
RUN_LABEL="${RUN_LABEL:-main}"

REF_DIR="$WORK_DIR/reference/ncbi_12X_refseq"
INDEX_DIR="$REF_DIR/star_index"
SRA_DIR="$WORK_DIR/star_sra"
FASTQ_DIR="$WORK_DIR/star_fastq"
CLEAN_DIR="$WORK_DIR/star_clean_fastq"
QC_DIR="$WORK_DIR/qc/fastp_star"
ALIGN_DIR="$WORK_DIR/star_align"
COUNT_DIR="$WORK_DIR/star_counts"
LOG_DIR="$WORK_DIR/logs_star"
META_DIR="$WORK_DIR/metadata_star"
mkdir -p "$REF_DIR" "$INDEX_DIR" "$SRA_DIR" "$FASTQ_DIR" "$CLEAN_DIR" \
  "$QC_DIR" "$ALIGN_DIR" "$COUNT_DIR" "$LOG_DIR" "$META_DIR"

for tool in prefetch fasterq-dump vdb-validate fastp STAR featureCounts samtools pigz curl aria2c; do
  command -v "$tool" >/dev/null || { echo "Missing required tool: $tool" >&2; exit 3; }
done

BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/745/GCF_000003745.3_12X"
download() {
  local url="$1"
  local output="$2"
  if [[ ! -s "$output" ]]; then
    curl --fail --location --retry 5 --retry-all-errors --continue-at - "$url" --output "$output"
  fi
}

download "$BASE_URL/GCF_000003745.3_12X_genomic.fna.gz" "$REF_DIR/genome.fa.gz"
download "$BASE_URL/GCF_000003745.3_12X_genomic.gtf.gz" "$REF_DIR/genes.gtf.gz"
[[ -s "$REF_DIR/genome.fa" ]] || pigz -dc "$REF_DIR/genome.fa.gz" > "$REF_DIR/genome.fa"
[[ -s "$REF_DIR/genes.gtf" ]] || pigz -dc "$REF_DIR/genes.gtf.gz" > "$REF_DIR/genes.gtf"

if [[ ! -s "$INDEX_DIR/Genome" || ! -s "$INDEX_DIR/SA" || ! -s "$INDEX_DIR/SAindex" || ! -s "$INDEX_DIR/genomeParameters.txt" ]]; then
  # A complete STAR index needs all four files; remove partial output from an interrupted build.
  find "$INDEX_DIR" -mindepth 1 -maxdepth 1 -delete
  STAR \
    --runMode genomeGenerate \
    --runThreadN "$THREADS" \
    --genomeDir "$INDEX_DIR" \
    --genomeFastaFiles "$REF_DIR/genome.fa" \
    --sjdbGTFfile "$REF_DIR/genes.gtf" \
    --sjdbOverhang 85 \
    --genomeSAindexNbases 13 \
    --genomeSAsparseD 2 \
    --limitGenomeGenerateRAM 9000000000 \
    2>&1 | tee "$LOG_DIR/star_index.log"
fi

cp "$MANIFEST" "$META_DIR/GSE337039_sample_manifest_${RUN_LABEL}.tsv"
sha256sum "$MANIFEST" "$REF_DIR/genome.fa.gz" "$REF_DIR/genes.gtf.gz" \
  > "$META_DIR/input_sha256_${RUN_LABEL}.tsv"

tail -n +2 "$MANIFEST" | cut -f1 | tr -d '"\r' | while read -r RUN; do
  [[ -n "$RUN" ]] || continue
  if [[ -s "$COUNT_DIR/$RUN.unstranded.txt" && -s "$COUNT_DIR/$RUN.forward.txt" ]]; then
    echo "SKIP completed $RUN"
    continue
  fi

  echo "START $RUN $(date --iso-8601=seconds)"
  SRA_FILE="$SRA_DIR/$RUN/$RUN.sra"
  if [[ -s "$SRA_FILE" ]] && ! vdb-validate "$SRA_FILE" \
    > "$LOG_DIR/${RUN}_validate_existing.log" 2>&1; then
    rm -f "$SRA_FILE"
  fi
  if [[ ! -s "$SRA_FILE" ]]; then
    mkdir -p "$SRA_DIR/$RUN"
    AWS_URL="https://sra-pub-run-odp.s3.amazonaws.com/sra/$RUN/$RUN"
    PART_FILE="$SRA_FILE.part"
    if aria2c --continue=true --max-connection-per-server=16 --split=16 \
      --min-split-size=1M --file-allocation=none --max-tries=5 --retry-wait=5 \
      --dir "$SRA_DIR/$RUN" --out "$RUN.sra.part" "$AWS_URL" \
      2>&1 | tee "$LOG_DIR/${RUN}_download.log"; then
      mv "$PART_FILE" "$SRA_FILE"
    fi
    if [[ ! -s "$SRA_FILE" ]] || ! vdb-validate "$SRA_FILE" \
      2>&1 | tee "$LOG_DIR/${RUN}_validate.log"; then
      rm -f "$SRA_FILE" "$PART_FILE"
      prefetch "$RUN" --output-directory "$SRA_DIR" --max-size 10G \
        2>&1 | tee "$LOG_DIR/${RUN}_prefetch.log"
    fi
    rm -f "$SRA_DIR/$RUN/$RUN.sra.tmp" "$SRA_DIR/$RUN/$RUN.sra.lock" \
      "$SRA_DIR/$RUN/$RUN.sra.prf"
  fi
  [[ -s "$SRA_FILE" ]] || { echo "Missing prefetched SRA: $SRA_FILE" >&2; exit 4; }

  if [[ ! -s "$FASTQ_DIR/$RUN.fastq" ]]; then
    fasterq-dump "$SRA_FILE" --threads "$SAMPLE_THREADS" --outdir "$FASTQ_DIR" \
      2>&1 | tee "$LOG_DIR/${RUN}_fasterq.log"
  fi
  [[ -s "$FASTQ_DIR/$RUN.fastq" ]] || { echo "Missing single-end FASTQ for $RUN" >&2; exit 5; }
  [[ ! -e "$FASTQ_DIR/${RUN}_2.fastq" ]] || { echo "Unexpected paired-end output" >&2; exit 6; }

  if [[ ! -s "$CLEAN_DIR/$RUN.clean.fastq.gz" ]]; then
    fastp \
      --in1 "$FASTQ_DIR/$RUN.fastq" \
      --out1 "$CLEAN_DIR/$RUN.clean.fastq.gz" \
      --thread "$SAMPLE_THREADS" \
      --qualified_quality_phred 20 \
      --unqualified_percent_limit 40 \
      --length_required 30 \
      --html "$QC_DIR/$RUN.html" \
      --json "$QC_DIR/$RUN.json" \
      2>&1 | tee "$LOG_DIR/${RUN}_fastp.log"
  fi

  PREFIX="$ALIGN_DIR/$RUN."
  TMP_STAR_DIR="${TMPDIR:-/tmp}/GSE337039_STAR_${RUN}_$$"
  trap 'rm -rf "$TMP_STAR_DIR"' EXIT
  STAR \
    --runThreadN "$SAMPLE_THREADS" \
    --genomeDir "$INDEX_DIR" \
    --readFilesIn "$CLEAN_DIR/$RUN.clean.fastq.gz" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$PREFIX" \
    --outTmpDir "$TMP_STAR_DIR" \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts \
    --outSAMattributes NH HI AS nM \
    2>&1 | tee "$LOG_DIR/${RUN}_star.log"
  rm -rf "$TMP_STAR_DIR"
  trap - EXIT

  BAM="$PREFIX"'Aligned.sortedByCoord.out.bam'
  [[ -s "$BAM" ]] || { echo "Missing STAR BAM for $RUN" >&2; exit 7; }
  featureCounts -T "$SAMPLE_THREADS" -a "$REF_DIR/genes.gtf" -t exon -g gene_id \
    -s 0 -o "$COUNT_DIR/$RUN.unstranded.txt" "$BAM" \
    2>&1 | tee "$LOG_DIR/${RUN}_featureCounts_unstranded.log"
  featureCounts -T "$SAMPLE_THREADS" -a "$REF_DIR/genes.gtf" -t exon -g gene_id \
    -s 1 -o "$COUNT_DIR/$RUN.forward.txt" "$BAM" \
    2>&1 | tee "$LOG_DIR/${RUN}_featureCounts_forward.log"

  [[ -s "$COUNT_DIR/$RUN.forward.txt" ]] || { echo "featureCounts failed for $RUN" >&2; exit 8; }
  rm -f "$BAM" "$FASTQ_DIR/$RUN.fastq" "$CLEAN_DIR/$RUN.clean.fastq.gz" "$SRA_FILE"
  rmdir "$SRA_DIR/$RUN" 2>/dev/null || true
  echo "DONE $RUN $(date --iso-8601=seconds)"
done

if [[ "$RUN_FINAL_QC" == "1" ]]; then
  multiqc "$QC_DIR" "$ALIGN_DIR" "$COUNT_DIR" -o "$WORK_DIR/qc/multiqc_star" --force \
    2>&1 | tee "$LOG_DIR/multiqc.log"

  {
    echo -e "tool\tversion"
    echo -e "sra_tools\t$(prefetch --version 2>&1 | head -n1)"
    echo -e "fastp\t$(fastp --version 2>&1 | head -n1)"
    echo -e "STAR\t$(STAR --version 2>&1 | head -n1)"
    echo -e "featureCounts\t$(featureCounts -v 2>&1 | head -n1)"
    echo -e "samtools\t$(samtools --version 2>&1 | head -n1)"
  } > "$META_DIR/software_versions.tsv"
fi

echo "GSE337039_STAR_COUNTING_OK"
