#!/usr/bin/env bash
# Run featureCounts on all sorted BAMs at once (consistent library-size
# normalization context across samples), using locus_tag as the gene
# identifier (present on all 4770 "gene" features in the GFF3, unlike
# old_locus_tag which is missing on 549 of them, incl. rnpB).
set -euo pipefail

BAM_DIR=/home/nyako/rnaseq_requantification/bam
COUNTS_DIR=/home/nyako/rnaseq_requantification/counts
GFF=/home/nyako/projects/ProteinHunter/data/databases/target/methanosarcina_acetivorans/ncbi_dataset/data/GCF_000007345.1/genomic.gff
MM=/home/nyako/rnaseq_requantification/tools/micromamba
export MAMBA_ROOT_PREFIX=/home/nyako/rnaseq_requantification/tools/mamba_root

mkdir -p "$COUNTS_DIR"

BAMS=$(ls "$BAM_DIR"/*.sorted.bam)
echo "BAMs: $BAMS"

"$MM" run -n rnaseq featureCounts -p --countReadPairs -T 16 -t gene -g locus_tag \
  -a "$GFF" -o "$COUNTS_DIR/counts_locus_tag.txt" $BAMS \
  2>&1 | tee "$COUNTS_DIR/featureCounts.log"

echo "=== featureCounts summary ==="
cat "$COUNTS_DIR/counts_locus_tag.txt.summary"
