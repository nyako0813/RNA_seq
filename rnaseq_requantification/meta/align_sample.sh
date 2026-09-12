#!/usr/bin/env bash
# Align one paired-end RNA-seq sample against the M. acetivorans C2A genome
# with BWA-mem, sort, and index. Usage: align_sample.sh <SRR_ID>
set -euo pipefail

SRR="$1"
THREADS=16
FASTQ_DIR=/home/nyako/rnaseq_requantification/fastq
BAM_DIR=/home/nyako/rnaseq_requantification/bam
GENOME=/home/nyako/rnaseq_requantification/genome/NC_003552.1.fna
MM=/home/nyako/rnaseq_requantification/tools/micromamba
export MAMBA_ROOT_PREFIX=/home/nyako/rnaseq_requantification/tools/mamba_root

mkdir -p "$BAM_DIR"

R1="$FASTQ_DIR/${SRR}_1.fastq.gz"
R2="$FASTQ_DIR/${SRR}_2.fastq.gz"
OUT_SORTED="$BAM_DIR/${SRR}.sorted.bam"

if [[ ! -f "$R1" || ! -f "$R2" ]]; then
  echo "[$SRR] FASTQ files missing, skipping" >&2
  exit 1
fi

if [[ -f "$OUT_SORTED" && -f "${OUT_SORTED}.bai" ]]; then
  echo "[$SRR] already aligned, skipping"
  exit 0
fi

echo "[$SRR] bwa mem started $(date)"
"$MM" run -n rnaseq bash -c "bwa mem -t $THREADS '$GENOME' '$R1' '$R2' 2> '$BAM_DIR/${SRR}.bwa.log' | samtools sort -@ $THREADS -o '$OUT_SORTED' -"
"$MM" run -n rnaseq samtools index "$OUT_SORTED"
echo "[$SRR] done $(date)"
"$MM" run -n rnaseq samtools flagstat "$OUT_SORTED" > "$BAM_DIR/${SRR}.flagstat.txt"
cat "$BAM_DIR/${SRR}.flagstat.txt"
