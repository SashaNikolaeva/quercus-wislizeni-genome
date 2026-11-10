#!/bin/bash
# wg23: rRNA genes (18S, 5.8S, 28S, 5S) with barrnap 0.9 on the WHOLE assembly
# per haplotype (array task 1 = hap1, 2 = hap2). rDNA arrays may sit on unplaced contigs, so wg22b reports chromosomes
# and unplaced sequence separately. barrnap/0.9 exists only under StdEnv/2020, so it runs in a subshell.
# Command: barrnap --kingdom euk --threads 8 <fasta>  (defaults: --reject 0.25 = predictions shorter than 25% of the
#   expected length are dropped; --lencutoff 0.8 = shorter than 80% of expected is labelled partial; --evalue 1e-6).
# barrnap 0.9 writes the partial label into the product attribute ("product=18S ribosomal RNA (partial)"), not into Name
# (Name=18S_rRNA etc.); wg22b therefore calls a hit partial when the attributes contain "(partial)".
# Assembled copies are collapsed arrays, not copy number (k-mer estimate in results/rdna/).
# Inputs: /home/sashanik/projects/def-santang3/sashanik/california_oaks/genomes/Qwislizeni_reasm_hap{1,2}.fa
# Outputs: ~/scratch/california_oaks/paperA_runs/wg23/hap{1,2}.barrnap.gff3
#   SMOKE=1: hap2 chr12 only -> wg23/smoke/hap2_chr12.barrnap.gff3; never results/
# Usage: sbatch --export=ALL,SMOKE=1 --array=2 --time=1:00:00 scripts/wg23_barrnap.sh   (smoke)
#        sbatch --array=1-2 --dependency=afterok:<smoke> scripts/wg23_barrnap.sh        (full)
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --job-name=wg23_barrnap
#SBATCH --output=/home/sashanik/scratch/california_oaks/paperA_runs/logs/wg23_%A_%a.log
set -euo pipefail
SMOKE=${SMOKE:-0}
GD=/home/sashanik/projects/def-santang3/sashanik/california_oaks/genomes
W=/home/sashanik/scratch/california_oaks/paperA_runs/wg23
hap=${SLURM_ARRAY_TASK_ID:?array task id required}
G=$GD/Qwislizeni_reasm_hap$hap.fa
export TMPDIR=${SLURM_TMPDIR:-/tmp}   # barrnap uses File::Temp
if [ "$SMOKE" = 1 ]; then
    O=$W/smoke; mkdir -p "$O"
    c=$(printf 'Qwis_hap%d_chr12' "$hap")
    ( module load StdEnv/2023 samtools/1.22.1 && samtools faidx "$G" "$c" > "$SLURM_TMPDIR/$c.fa" )
    IN=$SLURM_TMPDIR/$c.fa; OUT=$O/hap${hap}_chr12.barrnap.gff3
else
    O=$W; mkdir -p "$O"; IN=$G; OUT=$O/hap$hap.barrnap.gff3
fi
rm -f "$OUT"
( module load StdEnv/2020 barrnap/0.9 && barrnap --version 2>&1 && \
  barrnap --kingdom euk --threads 8 "$IN" > "$OUT.tmp" )
mv "$OUT.tmp" "$OUT"
echo "barrnap hits by Name (partial = '(partial)' in attributes):"
awk -F'\t' '!/^#/ && NF>=9 {split($9,a,";"); n=a[1]; p=($9 ~ /\(partial\)/) ? "partial" : "full"; k[n" "p]++}
            END{for (x in k) print x, k[x]}' "$OUT" | sort
echo "DONE wg23 hap$hap"
