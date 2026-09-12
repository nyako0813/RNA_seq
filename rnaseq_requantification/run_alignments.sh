#!/bin/bash
set -uo pipefail
cd /home/nyako/rnaseq_requantification

export MAMBA_ROOT_PREFIX=/home/nyako/rnaseq_requantification/tools/mamba_root
MM=/home/nyako/rnaseq_requantification/tools/micromamba

GENOME=/home/nyako/rnaseq_requantification/genome/NC_003552.1.fna
FASTQ_DIR=/home/nyako/rnaseq_requantification/fastq
BAM_DIR=/home/nyako/rnaseq_requantification/bam
mkdir -p "$BAM_DIR"

THREADS_BWA=4
THREADS_SORT=1
MAX_PARALLEL=3

# The 6 samples actually downloaded (resolved from SRX16315152-155 /
# SRX16314424/429 via ENA -- see meta/SRX*.tsv for the mapping).
samples=(SRR20281608 SRR20281609 SRR20281610 SRR20281611 SRR20280855 SRR20280860)

run_one() {
  local srr="$1"
  local r1="$FASTQ_DIR/${srr}_1.fastq.gz"
  local r2="$FASTQ_DIR/${srr}_2.fastq.gz"
  local out="$BAM_DIR/${srr}.sorted.bam"

  if [[ -f "$out" && -f "${out}.bai" ]]; then
    echo "$(date '+%F %T'): ${srr} already done, skipping" >> progress.log
    return 0
  fi

  "$MM" run -n rnaseq bash -c "bwa mem -t $THREADS_BWA '$GENOME' '$r1' '$r2' 2> '$BAM_DIR/${srr}.align.log' | samtools sort -@ $THREADS_SORT -o '$out' -"
  if [[ $? -ne 0 ]]; then
    echo "$(date '+%F %T'): ${srr} FAILED (bwa/sort)" >> progress.log
    return 1
  fi
  "$MM" run -n rnaseq samtools index "$out"
  "$MM" run -n rnaseq samtools flagstat "$out" > "$BAM_DIR/${srr}.flagstat.txt"
  echo "$(date '+%F %T'): ${srr} done" >> progress.log
}

running=0
for srr in "${samples[@]}"; do
  run_one "$srr" &
  running=$((running+1))
  if [ "$running" -ge "$MAX_PARALLEL" ]; then
    wait -n
    running=$((running-1))
  fi
done
wait
echo "ALL DONE: $(date '+%F %T')" >> progress.log
