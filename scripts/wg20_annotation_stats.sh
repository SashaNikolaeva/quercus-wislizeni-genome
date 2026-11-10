#!/bin/bash
# wg20: gene-model statistics for hap1_all, hap1_nonTE
# and hap2_all. Statistic definitions: header of
# wg20_annotation_stats.py (introns = gaps between consecutive exons of one mRNA). Fails unless gene counts are
# 36,682 / 32,788 (= gene_te_summary.tsv) / 36,582 and n_introns = n_exons - n_mRNA.
# Inputs : cg05 hap1 GFF3; ~/scratch/california_oaks/paperA_runs/wg19/Qwis_hap2.helixer.gff3 (from wg19);
#          results/genes/gene_te_overlap.tsv, results/genes/gene_te_summary.tsv
# Output : results/annotation/annotation_stats.tsv
#          SMOKE=1: chr12 only, from the wg19 smoke GFF3s -> ~/scratch/california_oaks/paperA_runs/wg20/smoke/annotation_stats.tsv
# Usage  : sbatch --export=ALL,SMOKE=1 --dependency=afterok:<wg19 smoke> scripts/wg20_annotation_stats.sh
#          sbatch --dependency=afterok:<wg19 full>:<wg20 smoke> scripts/wg20_annotation_stats.sh
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0:30:00
#SBATCH --job-name=wg20_annstats
#SBATCH --output=/home/sashanik/scratch/california_oaks/paperA_runs/logs/wg20_%j.log
set -euo pipefail
SMOKE=${SMOKE:-0}
WG=/home/sashanik/oak-genomics-wislizeni/california_oaks/wislizeni_genome
PY=/home/sashanik/projects/def-santang3/sashanik/software/venv_plotting/bin/python
W19=/home/sashanik/scratch/california_oaks/paperA_runs/wg19
W=/home/sashanik/scratch/california_oaks/paperA_runs/wg20
if [ "$SMOKE" = 1 ]; then
    O=$W/smoke; GFF1=$W19/smoke/hap1.helixer.gff3; GFF2=$W19/smoke/Qwis_hap2.helixer.gff3
else
    O=$WG/results/annotation; GFF1=/home/sashanik/scratch/california_oaks/cg05_helixer/Qwislizeni/Qwislizeni.helixer.gff3
    GFF2=$W19/Qwis_hap2.helixer.gff3
fi
mkdir -p "$O"
echo "[$(date '+%F %T')] wg20 SMOKE=$SMOKE"
"$PY" "$WG/scripts/wg20_annotation_stats.py" "$GFF1" "$GFF2" "$WG/results/genes/gene_te_overlap.tsv" \
    "$WG/results/genes/gene_te_summary.tsv" "$O/annotation_stats.tsv" "$SMOKE"
cat "$O/annotation_stats.tsv"
echo "[$(date '+%F %T')] DONE wg20"
