#!/bin/bash
# wg21: functional annotation with InterProScan 5.73-104.0
# (StdEnv/2023 module; loads interproscan_data/5.73-104.0), member DBs Pfam, PANTHER, Gene3D, SUPERFAMILY, SMART, CDD,
# GO terms, no lookup service (-dp). eggNOG dropped.
# Chunking (reproducible): proteins of prot_hap<h>.fa sorted by ID (LC_ALL=C sort), protein with 0-based sorted index i
#   goes to chunk (i mod 20) + 1. Array task t: hap = 1 + (t > 20), chunk = ((t - 1) mod 20) + 1.
# Inputs : ~/scratch/california_oaks/paperA_runs/wg19/prot_hap{1,2}.fa (wg19)
#          SMOKE=1: the first 50 records of ~/scratch/california_oaks/wislizeni_reasm/final/wg10_proteome/prot.fa (hap1)
# Outputs: ~/scratch/california_oaks/paperA_runs/wg21/hap<h>/chunk_<kk>.{tsv,ids}  (raw TSVs stay here; wg21b summarises)
#          SMOKE=1: ~/scratch/california_oaks/paperA_runs/wg21/smoke/chunk_smoke.{tsv,ids}; smoke must contain Pfam rows
#          with GO terms and PANTHER rows with GO terms, else exit 1.
# Checks : chunk TSV non-empty; distinct proteins in the TSV <= proteins in the chunk.
# Temp   : $SLURM_TMPDIR/ips (InterProScan writes many files; scratch file quota).
# Usage  : sbatch --export=ALL,SMOKE=1 --array=1 --time=1:00:00 scripts/wg21_interproscan.sh
#          sbatch --dependency=afterok:<wg19 full>:<wg21 smoke> scripts/wg21_interproscan.sh
#          (if PANTHER pushes a task over 8 h: rerun only that task, --array=<t> --time=16:00:00, nothing else changed)
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=8:00:00
#SBATCH --array=1-40
#SBATCH --job-name=wg21_ips
#SBATCH --output=/home/sashanik/scratch/california_oaks/paperA_runs/logs/wg21_%A_%a.log
set -euo pipefail
SMOKE=${SMOKE:-0}
W=/home/sashanik/scratch/california_oaks/paperA_runs/wg21
T=${SLURM_TMPDIR:?no SLURM_TMPDIR}
t=${SLURM_ARRAY_TASK_ID:-1}
if [ "$SMOKE" = 1 ]; then
    O=$W/smoke; mkdir -p "$O"; TAG=chunk_smoke
    awk '/^>/{n++} n<=50' /home/sashanik/scratch/california_oaks/wislizeni_reasm/final/wg10_proteome/prot.fa > "$T/chunk.fa"
else
    h=$(( t > 20 ? 2 : 1 )); k=$(( (t - 1) % 20 + 1 )); TAG=$(printf 'chunk_%02d' "$k")
    IN=/home/sashanik/scratch/california_oaks/paperA_runs/wg19/prot_hap$h.fa
    O=$W/hap$h; mkdir -p "$O"
    # linearise -> sort by ID -> take every 20th starting at k-1
    awk '/^>/{if(n)print n"\t"s; n=substr($1,2); s=""; next}{s=s $0} END{if(n)print n"\t"s}' "$IN" \
        | LC_ALL=C sort -t$'\t' -k1,1 \
        | awk -F'\t' -v k="$k" '(NR-1)%20==k-1{print ">"$1"\n"$2}' > "$T/chunk.fa"
fi
grep '^>' "$T/chunk.fa" | sed 's/^>//' > "$O/$TAG.ids"
N=$(wc -l < "$O/$TAG.ids")
echo "[$(date '+%F %T')] wg21 SMOKE=$SMOKE task=$t $TAG proteins=$N"
[ "$N" -gt 0 ] || { echo "ERROR: empty chunk"; exit 1; }

module load StdEnv/2023 interproscan/5.73-104.0
mkdir -p "$T/ips"; cd "$T"
interproscan.sh -i chunk.fa -f tsv -dp -goterms -appl Pfam,PANTHER,Gene3D,SUPERFAMILY,SMART,CDD \
    -cpu 16 -T "$T/ips" -o "$T/chunk.tsv"
[ -s "$T/chunk.tsv" ] || { echo "ERROR: empty InterProScan TSV"; exit 1; }
NP=$(cut -f1 "$T/chunk.tsv" | sort -u | wc -l)
echo "proteins with >=1 hit: $NP of $N"
[ "$NP" -le "$N" ] || { echo "ERROR: more distinct proteins in output than input"; exit 1; }
cut -f4 "$T/chunk.tsv" | sort | uniq -c
cp "$T/chunk.tsv" "$O/$TAG.tsv.tmp" && mv "$O/$TAG.tsv.tmp" "$O/$TAG.tsv"

if [ "$SMOKE" = 1 ]; then
    # GO column = 14 (with -goterms); '-' or empty means none
    PF=$(awk -F'\t' '$4=="Pfam" && $14!="" && $14!="-"' "$O/$TAG.tsv" | wc -l)
    PA=$(awk -F'\t' '$4=="PANTHER" && $14!="" && $14!="-"' "$O/$TAG.tsv" | wc -l)
    echo "smoke: Pfam rows with GO = $PF, PANTHER rows with GO = $PA"
    [ "$PF" -gt 0 ] && [ "$PA" -gt 0 ] || { echo "SMOKE FAILED: need Pfam and PANTHER rows with GO terms"; exit 1; }
fi
echo "[$(date '+%F %T')] DONE wg21 $TAG"
