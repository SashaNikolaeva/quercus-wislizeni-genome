#!/bin/bash
# wg24b: SLURM wrapper for wg24b_busco_dup_check.py. For every
# protein-mode odb12 Duplicated BUSCO: n copies, copies within 3 genes on one chromosome, TE-derived (hap1), and
# copy length / median length of Single hits (and / odb12 refseq median). Tests the fixed prediction of explanation (a)
# 'artifact: fragments or TE-derived models' vs (b) 'real paralogs'. Definitions: header of the .py.
# Inputs : ~/scratch/california_oaks/paperA_runs/wg24/protein_hap{1,2}/full_table.tsv (wg24 tasks 3-4),
#          wg19 prot_hap{1,2}.fa, cg05 hap1 GFF3, wg19 Qwis_hap2.helixer.gff3, results/genes/gene_te_overlap.tsv,
#          software/compleasm_lineages/embryophyta_odb12/refseq_db.faa.gz
# Outputs: results/annotation/busco_dup_check.tsv, results/annotation/busco_dup_check_summary.tsv
# Usage  : sbatch --dependency=afterok:<wg24 tasks 3-4> scripts/wg24b_busco_dup_check.sh
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0:30:00
#SBATCH --job-name=wg24b_dupchk
#SBATCH --output=/home/sashanik/scratch/california_oaks/paperA_runs/logs/wg24b_%j.log
set -euo pipefail
WG=/home/sashanik/oak-genomics-wislizeni/california_oaks/wislizeni_genome; R=$WG/results/annotation; mkdir -p "$R"
PY=/home/sashanik/projects/def-santang3/sashanik/software/venv_plotting/bin/python
W24=/home/sashanik/scratch/california_oaks/paperA_runs/wg24
W19=/home/sashanik/scratch/california_oaks/paperA_runs/wg19
echo "[$(date '+%F %T')] wg24b"
"$PY" "$WG/scripts/wg24b_busco_dup_check.py" "$R/busco_dup_check.tsv" "$R/busco_dup_check_summary.tsv" \
    "$WG/results/genes/gene_te_overlap.tsv" \
    /home/sashanik/projects/def-santang3/sashanik/software/compleasm_lineages/embryophyta_odb12/refseq_db.faa.gz \
    hap1 "$W24/protein_hap1/full_table.tsv" "$W19/prot_hap1.fa" /home/sashanik/scratch/california_oaks/cg05_helixer/Qwislizeni/Qwislizeni.helixer.gff3 \
    hap2 "$W24/protein_hap2/full_table.tsv" "$W19/prot_hap2.fa" "$W19/Qwis_hap2.helixer.gff3"
echo "[$(date '+%F %T')] DONE wg24b"
