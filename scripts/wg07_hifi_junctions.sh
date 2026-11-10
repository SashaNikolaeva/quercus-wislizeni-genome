#!/bin/bash
# wg07: HiFi read support for the wislizeni inversions (hap1-vs-hap2 INV >= 100 kb, the chr06
# complex, and the chr05 inversion both haplotypes carry). In the genome report only step 1 is used: the HiFi mapping to
# the diploid assembly (dip.bam), whose read depth is used by wg25 (organelle classification) and wg27 (BlobToolKit).
# 1. HiFi SRR19243399 (the hifiasm input, ~47x) -> minimap2 2.28 -ax map-hifi to the DIPLOID reference hap1 + hap2 (all
#    sequences, final/Qwis_hap{1,2}.fa) -> samtools sort (node-local -T) -> dip.bam; samtools coverage -q 20 per sequence.
# 2. wg07_hifi_junctions.py: spanning-read traces across every junction zone, reads supporting the other arrangement,
#    depth inside each block; verdict per candidate (rules in the .py header).
# Output: ~/scratch/california_oaks/wislizeni_reasm/hifi_junctions/ (dip.bam, coverage, traces);
#         tables -> results/wislizeni_reasm/hifi_junctions/, figure -> figures/wislizeni_reasm/hifi_junctions/
# Usage: sbatch scripts/wg07_hifi_junctions.sh
#        SMOKE=1 sbatch --export=ALL --time=1:00:00 -c 16 --mem=48G scripts/wg07_hifi_junctions.sh   (first 20,000 reads,
#        own dir smoke/: tests the code path, verdicts meaningless)
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=20:00:00
#SBATCH --tmp=300G
#SBATCH --job-name=wg07_hifi_junc
#SBATCH --output=/home/sashanik/scratch/california_oaks/logs/wg07_hifi_junc_%j.log
set -euo pipefail
REPO=/home/sashanik/oak-genomics-wislizeni/california_oaks/wislizeni_genome
WR=/home/sashanik/scratch/california_oaks/wislizeni_reasm; F=$WR/final
HIFI=$WR/raw/SRR19243399_subreads.fastq.gz
SYRI=/home/sashanik/scratch/california_oaks/cg03_syri/Qwislizeni_hap2/syri.out
PY=/home/sashanik/projects/def-santang3/sashanik/software/venv_plotting/bin/python
O=$WR/hifi_junctions; RES=$REPO/results/wislizeni_reasm/hifi_junctions; FIG=$REPO/figures/wislizeni_reasm/hifi_junctions
if [ "${SMOKE:-0}" = 1 ]; then O=$O/smoke; RES=$O/res; FIG=$O/fig; fi
mkdir -p "$O"; cd "$O"
module load StdEnv/2023 minimap2/2.28 samtools/1.22.1
T=${SLURM_CPUS_PER_TASK:-32}
ts() { date '+%F %T'; }
[ -s dip.fa ] || { cat "$F/Qwis_hap1.fa" "$F/Qwis_hap2.fa" > dip.fa; samtools faidx dip.fa; }
if [ ! -s dip.bam.bai ]; then
    if [ "${SMOKE:-0}" = 1 ]; then READS=(<(zcat "$HIFI" | awk 'NR<=80000')); else READS=("$HIFI"); fi
    echo "[$(ts)] minimap2 map-hifi on $(hostname), $T threads"
    minimap2 -ax map-hifi -t "$T" dip.fa "${READS[@]}" 2> minimap2.log \
        | samtools sort -@ 8 -m 3G -T "$SLURM_TMPDIR/wg07sort" -o dip.bam -
    samtools index -@ 8 dip.bam
fi
echo "[$(ts)] samtools coverage"
[ -s coverage_q20.tsv ] || samtools coverage -q 20 dip.bam > coverage_q20.tsv
echo "[$(ts)] junction test"
$PY "$REPO/scripts/wg07_hifi_junctions.py" dip.bam "$SYRI" coverage_q20.tsv "$RES" "$FIG" "$O"
column -t "$RES/inversion_verdicts.tsv"
echo "[$(ts)] DONE wg07"
