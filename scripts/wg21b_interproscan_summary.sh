#!/bin/bash
# wg21b: SLURM wrapper for wg21b_interproscan_summary.py ().
# Per set (hap1_all, hap1_nonTE, hap2_all): % genes with >= 1 InterPro entry / Pfam / GO / any member-DB hit, and per-gene
# best Pfam/PANTHER product (lowest e-value) for a later NCBI submission. Definitions: header of the .py.
# Inputs : ~/scratch/california_oaks/paperA_runs/wg21/hap{1,2}/chunk_NN.{tsv,ids} (20 per hap), wg19 prot_hap{1,2}.fa,
#          cg05 hap1 GFF3, wg19 Qwis_hap2.helixer.gff3, results/genes/gene_te_overlap.tsv
# Outputs: results/annotation/interproscan_summary.tsv, results/annotation/interproscan_gene_products_hap{1,2}.tsv
#          SMOKE=1: hap1 only, from the wg21 smoke chunk (wg10 prot.fa, 50 proteins) -> ~/scratch/california_oaks/paperA_runs/wg21/smoke/
# Usage  : sbatch --export=ALL,SMOKE=1 --dependency=afterok:<wg21 smoke> scripts/wg21b_interproscan_summary.sh
#          sbatch --dependency=afterok:<wg21 array>:<wg21b smoke> scripts/wg21b_interproscan_summary.sh
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0:30:00
#SBATCH --job-name=wg21b_ipssum
#SBATCH --output=/home/sashanik/scratch/california_oaks/paperA_runs/logs/wg21b_%j.log
set -euo pipefail
SMOKE=${SMOKE:-0}
WG=/home/sashanik/oak-genomics-wislizeni/california_oaks/wislizeni_genome
PY=/home/sashanik/projects/def-santang3/sashanik/software/venv_plotting/bin/python
GFF1=/home/sashanik/scratch/california_oaks/cg05_helixer/Qwislizeni/Qwislizeni.helixer.gff3
W19=/home/sashanik/scratch/california_oaks/paperA_runs/wg19
W=/home/sashanik/scratch/california_oaks/paperA_runs/wg21
TEO=$WG/results/genes/gene_te_overlap.tsv
echo "[$(date '+%F %T')] wg21b SMOKE=$SMOKE"
if [ "$SMOKE" = 1 ]; then
    O=$W/smoke; mkdir -p "$O"
    "$PY" "$WG/scripts/wg21b_interproscan_summary.py" "$O" 1 "$TEO" \
        hap1 "$GFF1" /home/sashanik/scratch/california_oaks/wislizeni_reasm/final/wg10_proteome/prot.fa "$W/smoke" 1
else
    O=$WG/results/annotation; mkdir -p "$O"
    "$PY" "$WG/scripts/wg21b_interproscan_summary.py" "$O" 0 "$TEO" \
        hap1 "$GFF1" "$W19/prot_hap1.fa" "$W/hap1" 20 \
        hap2 "$W19/Qwis_hap2.helixer.gff3" "$W19/prot_hap2.fa" "$W/hap2" 20
fi
echo "[$(date '+%F %T')] DONE wg21b"
