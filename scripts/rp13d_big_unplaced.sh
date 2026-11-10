#!/bin/bash
# rp13d: placement + origin of the large (>= 1 Mb) unplaced hap1 sequences of the new wislizeni assembly.
# minimap2 2.28 -x asm20 of those sequences vs the old CCGP hap1 + hap2, Qagrifolia, Qparvula, Qrubra_full, Qlobata;
# then Omni-C contact counts and coverage summary (rp13d_big_unplaced.py; method in its docstring).
# Output: wislizeni_reasm/final/rp13d_hap1_big_unplaced.tsv
# Usage: sbatch scripts/rp13d_big_unplaced.sh
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --job-name=rp13d_unplaced
#SBATCH --output=/home/sashanik/scratch/california_oaks/logs/rp13d_unplaced_%j.log
set -euo pipefail
W=/home/sashanik/scratch/california_oaks/wislizeni_reasm; F=$W/final; D=$F/rp13d; mkdir -p "$D"
G=/home/sashanik/projects/def-santang3/sashanik/california_oaks/genomes
REPO=/home/sashanik/oak-genomics/california_oaks/wislizeni_genome
PYH=/home/sashanik/projects/def-santang3/sashanik/software/miniforge3/envs/hic/bin/python
module load StdEnv/2023 minimap2/2.28 samtools/1.22.1
cd "$D"
samtools faidx "$F/Qwis_hap1.fa" $(awk '$1~/_unplaced_/ && $2>=1000000{print $1}' "$F/Qwis_hap1.fa.fai") > big_unplaced.fa
grep -c ">" big_unplaced.fa
for g in Qwislizeni Qwislizeni_hap2 Qagrifolia Qparvula Qrubra_full Qlobata; do
    [ -s "$g.paf" ] || minimap2 -x asm20 -t 16 "$G/$g.fa" big_unplaced.fa > "$g.paf"
done
"$PYH" "$REPO/scripts/rp13d_big_unplaced.py" "$F" "$W/yahs" hap1 "$F/rp13d_hap1_big_unplaced.tsv" "$D"/*.paf
echo "DONE rp13d"
