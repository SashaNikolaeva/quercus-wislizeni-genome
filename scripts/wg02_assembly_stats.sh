#!/bin/bash
# wg02: like-for-like assembly statistics for the paper's Table 2 / v1-vs-v2 table.
# (1) compleasm 0.2.10 with embryophyta_odb10 (the lineage the v1 paper's BUSCO 5.0.0 used; n = 1,614) on v2 hap1/hap2
#     AND on the v1 CCGP hap1/hap2, so all four numbers come from one tool. odb12 numbers are in rp10.
# (2) gfastats 1.3.11 with expected genome size 739,879,015 bp (GenomeScope2 point estimate, rp08 plot; was 741,000,159 = top of its range until 2026-10-22) for contig/scaffold NG50, the
#     quality code (Rhie et al. 2021) and gaps; whole assemblies, as in v1.
# Output: ~/scratch/california_oaks/wislizeni_reasm/wg02/, summary wg02_summary.tsv
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --array=1-4
#SBATCH --job-name=wg02_stats
#SBATCH --output=/home/sashanik/scratch/california_oaks/logs/wg02_stats_%A_%a.log
set -euo pipefail
G=/home/sashanik/projects/def-santang3/sashanik/california_oaks/genomes
NAMES=(v2_hap1 v2_hap2 v1_hap1 v1_hap2)
FILES=($G/Qwislizeni_reasm_hap1.fa $G/Qwislizeni_reasm_hap2.fa $G/Qwislizeni.fa $G/Qwislizeni_hap2.fa)
i=$((SLURM_ARRAY_TASK_ID-1)); N=${NAMES[$i]}; FA=${FILES[$i]}
ENV=/home/sashanik/projects/def-santang3/sashanik/software/miniforge3/envs
LIB=/home/sashanik/projects/def-santang3/sashanik/software/compleasm_lineages
export PATH="$ENV/compleasm/bin:$PATH"   # compleasm looks up miniprot on PATH (6699576 failed without it)
O=/home/sashanik/scratch/california_oaks/wislizeni_reasm/wg02; mkdir -p "$O"; cd "$O"
echo "[$(date '+%F %T')] $N $FA"
"$ENV/asm/bin/gfastats" "$FA" 739879015 --tabular > "$N.gfastats.tsv"
[ -s "$N.compleasm_odb10/summary.txt" ] || "$ENV/compleasm/bin/compleasm" run -a "$FA" -o "$N.compleasm_odb10" -l embryophyta --odb odb10 -L "$LIB" -t 16 > "$N.compleasm.log" 2>&1
cat "$N.compleasm_odb10/summary.txt"
echo "[$(date '+%F %T')] DONE $N"
