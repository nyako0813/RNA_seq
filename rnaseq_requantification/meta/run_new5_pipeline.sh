#!/usr/bin/env bash
# Process the 5 additional SRA runs with the SAME method as the existing
# 6-sample pipeline (see meta/align_sample.sh, meta/run_featurecounts.sh):
#   fasterq-dump -> bwa mem -> samtools sort/index/flagstat -> featureCounts
#
# Deviation from the original per-sample scripts (functionally equivalent,
# documented in the final report): fastq are aligned directly in
# uncompressed form (bwa mem reads plain or gzipped fastq identically via
# htslib), instead of gzip-compressing first, to avoid ~5h of single-
# threaded gzip time on ~9GB/file uncompressed fastq. Raw fastq/sra are
# deleted after each sample's BAM+index+flagstat are confirmed present, to
# control disk usage, matching how the original run only left
# SRR20281611's BAM behind.
set -uo pipefail

ROOT=/home/nyako/rnaseq_requantification
FASTQ_DIR="$ROOT/fastq"
BAM_DIR="$ROOT/bam"
LOG_DIR="$ROOT/logs_new5"
GENOME="$ROOT/genome/NC_003552.1.fna"
MM="$ROOT/tools/micromamba"
export MAMBA_ROOT_PREFIX="$ROOT/tools/mamba_root"
TMP_FQ="$FASTQ_DIR/tmp_fasterq"

mkdir -p "$BAM_DIR" "$LOG_DIR" "$TMP_FQ"

SAMPLES=(SRR20281612 SRR20281613 SRR20650029 SRR20650030 SRR20650031)

for SRR in "${SAMPLES[@]}"; do
  SRA="$FASTQ_DIR/$SRR/$SRR.sra"
  R1="$FASTQ_DIR/${SRR}_1.fastq"
  R2="$FASTQ_DIR/${SRR}_2.fastq"
  OUT_SORTED="$BAM_DIR/${SRR}.sorted.bam"

  echo "$(date '+%F %T') [$SRR] === start ===" | tee -a "$LOG_DIR/progress.log"

  # 1) fasterq-dump (skip if fastq already present from earlier manual test run)
  if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    if [[ ! -f "$SRA" ]]; then
      echo "$(date '+%F %T') [$SRR] SRA file missing, skipping" | tee -a "$LOG_DIR/progress.log"
      continue
    fi
    echo "$(date '+%F %T') [$SRR] fasterq-dump started" | tee -a "$LOG_DIR/progress.log"
    "$MM" run -n rnaseq fasterq-dump -O "$FASTQ_DIR" -t "$TMP_FQ" -e 4 -p --split-3 "$SRA" \
      > "$LOG_DIR/${SRR}.fasterq.log" 2>&1
    rc=$?
    if [[ $rc -ne 0 || ! -f "$R1" || ! -f "$R2" ]]; then
      echo "$(date '+%F %T') [$SRR] FAILED fasterq-dump (rc=$rc)" | tee -a "$LOG_DIR/progress.log"
      continue
    fi
    echo "$(date '+%F %T') [$SRR] fasterq-dump done" | tee -a "$LOG_DIR/progress.log"
  else
    echo "$(date '+%F %T') [$SRR] fastq already present, skipping fasterq-dump" | tee -a "$LOG_DIR/progress.log"
  fi

  # 2) bwa mem -> sort -> index -> flagstat (same params as meta/align_sample.sh)
  if [[ -f "$OUT_SORTED" && -f "${OUT_SORTED}.bai" ]]; then
    echo "$(date '+%F %T') [$SRR] already aligned, skipping" | tee -a "$LOG_DIR/progress.log"
  else
    echo "$(date '+%F %T') [$SRR] bwa mem started" | tee -a "$LOG_DIR/progress.log"
    "$MM" run -n rnaseq bash -c "bwa mem -t 16 '$GENOME' '$R1' '$R2' 2> '$LOG_DIR/${SRR}.bwa.log' | samtools sort -@ 16 -o '$OUT_SORTED' -"
    rc=$?
    if [[ $rc -ne 0 ]]; then
      echo "$(date '+%F %T') [$SRR] FAILED bwa/sort (rc=$rc)" | tee -a "$LOG_DIR/progress.log"
      continue
    fi
    "$MM" run -n rnaseq samtools index "$OUT_SORTED"
    "$MM" run -n rnaseq samtools flagstat "$OUT_SORTED" > "$BAM_DIR/${SRR}.flagstat.txt"
    echo "$(date '+%F %T') [$SRR] alignment done" | tee -a "$LOG_DIR/progress.log"
    cat "$BAM_DIR/${SRR}.flagstat.txt" | tee -a "$LOG_DIR/progress.log"
  fi

  # 3) cleanup raw fastq/sra to control disk usage, only once BAM+index exist
  if [[ -f "$OUT_SORTED" && -f "${OUT_SORTED}.bai" ]]; then
    rm -f "$R1" "$R2" "$SRA"
    rmdir "$FASTQ_DIR/$SRR" 2>/dev/null || true
    echo "$(date '+%F %T') [$SRR] cleaned up raw fastq/sra" | tee -a "$LOG_DIR/progress.log"
  fi

  echo "$(date '+%F %T') [$SRR] === end ===" | tee -a "$LOG_DIR/progress.log"
done

echo "$(date '+%F %T') ALL 5 SAMPLES PROCESSED" | tee -a "$LOG_DIR/progress.log"
