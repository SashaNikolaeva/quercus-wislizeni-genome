#!/bin/bash
# wg24: compleasm 0.2.10, embryophyta_odb12 (N = 2,026), 
# (genome mode on the FINAL scaffolds; rp10's odb12 numbers came from hifiasm contigs) and items 31/35 support (protein mode).
# Tasks: 1 genome hap1, 2 genome hap2 (Qwislizeni_reasm_hap{1,2}.fa)   -> submit now, --array=1-2
#        3 protein hap1, 4 protein hap2 (wg19 prot_hap{1,2}.fa)        -> separate submission, --array=3-4 afterok wg19
# Genome-mode check (reported, not enforced): S and D within +-0.2 percentage points of rp10 contigs
#   (hap1 S 97.88 / D 2.07; hap2 S 97.58 / D 2.42). A larger difference is reported, not explained away.
# No smoke.
# Outputs: results/annotation/compleasm_odb12_{genome,protein}_hap{1,2}.txt (summary.txt copies)
#          ~/scratch/california_oaks/paperA_runs/wg24/{genome,protein}_hap<h>/ (genome: embryophyta_odb12/full_table_busco_format.tsv
#          stays here for wg27; protein: full_table.tsv for wg24b)
# Usage  : sbatch --array=1-2 scripts/wg24_compleasm_odb12.sh
#          sbatch --array=3-4 --dependency=afterok:<wg19 full> scripts/wg24_compleasm_odb12.sh
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --job-name=wg24_cmpl12
#SBATCH --output=/home/sashanik/scratch/california_oaks/paperA_runs/logs/wg24_%A_%a.log
set -euo pipefail
t=${SLURM_ARRAY_TASK_ID:?run as an array task 1-4}
WG=/home/sashanik/oak-genomics-wislizeni/california_oaks/wislizeni_genome; R=$WG/results/annotation; mkdir -p "$R"
G=/home/sashanik/projects/def-santang3/sashanik/california_oaks/genomes
ENV=/home/sashanik/projects/def-santang3/sashanik/software/miniforge3/envs
LIB=/home/sashanik/projects/def-santang3/sashanik/software/compleasm_lineages
export PATH="$ENV/compleasm/bin:$PATH"   # compleasm looks up miniprot/hmmsearch on PATH (6699576 failed without it)
W=/home/sashanik/scratch/california_oaks/paperA_runs/wg24; mkdir -p "$W"; cd "$W"
h=$(( (t - 1) % 2 + 1 ))
if [ "$t" -le 2 ]; then MODE=genome; else MODE=protein; fi
D=${MODE}_hap$h
echo "[$(date '+%F %T')] wg24 task $t $MODE hap$h"
rm -rf "$D"
if [ "$MODE" = genome ]; then
    "$ENV/compleasm/bin/compleasm" run -a "$G/Qwislizeni_reasm_hap$h.fa" -o "$D" -l embryophyta --odb odb12 -L "$LIB" -t 16 > "$D.log" 2>&1
    [ -s "$D/embryophyta_odb12/full_table_busco_format.tsv" ] || { echo "ERROR: no full_table_busco_format.tsv"; exit 1; }
else
    P=/home/sashanik/scratch/california_oaks/paperA_runs/wg19/prot_hap$h.fa
    [ -s "$P" ] || { echo "ERROR: missing $P"; exit 1; }
    "$ENV/compleasm/bin/compleasm" protein -p "$P" -o "$D" -l embryophyta --odb odb12 -L "$LIB" -t 16 > "$D.log" 2>&1
    [ -s "$D/full_table.tsv" ] || { echo "ERROR: no full_table.tsv"; exit 1; }
fi
grep -q "embryophyta_odb12" "$D/summary.txt" || { echo "ERROR: summary is not odb12"; exit 1; }
cp "$D/summary.txt" "$R/compleasm_odb12_${MODE}_hap$h.txt"; cat "$D/summary.txt"
if [ "$MODE" = genome ]; then
    # rp10 contig-level odb12 numbers, |difference| > 0.2 pp is reported
    if [ "$h" = 1 ]; then RS=97.88; RD=2.07; else RS=97.58; RD=2.42; fi
    awk -v rs=$RS -v rd=$RD -v h=$h '/^S:/{sub(/^S:/,"");sub(/%.*/,"");s=$0} /^D:/{sub(/^D:/,"");sub(/%.*/,"");d=$0}
        END{ds=s-rs; dd=d-rd; flag=(ds>0.2||ds<-0.2||dd>0.2||dd<-0.2)?"REPORT: differs by >0.2 pp":"within 0.2 pp";
            printf "hap%s scaffolds S %.2f D %.2f vs rp10 contigs S %.2f D %.2f (dS %+.2f dD %+.2f): %s\n",h,s,d,rs,rd,ds,dd,flag}' "$D/summary.txt"
fi
echo "[$(date '+%F %T')] DONE wg24 task $t"
