#!/bin/bash
# rp13c: explain the hap1 vs hap2 size difference of the new wislizeni assembly.
# Aligns each haplotype's unplaced sequences (minimap2 2.28 -x asm5, -c for identities) to its own 12 chromosomes and to
# the other haplotype's 12 chromosomes, then runs rp13c_hap_sizes.py (classes and thresholds in that file's docstring).
# Output: wislizeni_reasm/final/rp13c_{chromosomes,unplaced,unplaced_summary}.tsv
# Usage: sbatch scripts/rp13c_hap_sizes.sh    (needs rp13 names: final/Qwis_hap{1,2}.{fa,names.tsv})
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --job-name=rp13c_hapsize
#SBATCH --output=/home/sashanik/scratch/california_oaks/logs/rp13c_hapsize_%j.log
set -euo pipefail
W=/home/sashanik/scratch/california_oaks/wislizeni_reasm; F=$W/final
REPO=/home/sashanik/oak-genomics/california_oaks/wislizeni_genome
PYH=/home/sashanik/projects/def-santang3/sashanik/software/miniforge3/envs/hic/bin/python
module load StdEnv/2023 minimap2/2.28 samtools/1.22.1
cd "$F"
for h in hap1 hap2; do
    samtools faidx "Qwis_$h.fa" $(awk '$1~/_chr/{print $1}' "Qwis_$h.fa.fai") > "rp13c_${h}_chroms.fa"
    samtools faidx "Qwis_$h.fa" $(awk '$1~/_unplaced_/{print $1}' "Qwis_$h.fa.fai") > "rp13c_${h}_unpl.fa"
done
for h in hap1 hap2; do
    o=$([ $h = hap1 ] && echo hap2 || echo hap1)
    minimap2 -x asm5 -c -t 16 "rp13c_${h}_chroms.fa" "rp13c_${h}_unpl.fa" > "rp13c_${h}_unpl_vs_own.paf"
    minimap2 -x asm5 -c -t 16 "rp13c_${o}_chroms.fa" "rp13c_${h}_unpl.fa" > "rp13c_${h}_unpl_vs_other.paf"
done
"$PYH" "$REPO/scripts/rp13c_hap_sizes.py" "$F" "$W/hifiasm" "$W/yahs" "$F/rp13c"
rm -f rp13c_hap?_chroms.fa rp13c_hap?_unpl.fa
echo "DONE rp13c"
