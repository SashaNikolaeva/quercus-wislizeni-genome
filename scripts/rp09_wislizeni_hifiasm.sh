#!/bin/bash
# rp09: hifiasm 0.25.0 in Hi-C mode, HiFi SRR19243399 + paired Omni-C
# (SRR19243397, SRR19243398 _1/_2; lists from rp08). Default purging; --telo-m CCCTAAA = plant telomere (TTTAGGG)n. *.bin files are kept so reruns with other
# parameters skip overlap/error correction.  GFA -> FASTA for hap1, hap2, primary; gfastats summary per FASTA.
# Output: ~/scratch/california_oaks/wislizeni_reasm/hifiasm/wis.hic.{hap1,hap2,p_ctg}.fa (+ .gfastats.tsv), rp09.ok
# Usage: sbatch --dependency=afterok:<rp08 job> scripts/rp09_wislizeni_hifiasm.sh
#SBATCH --account=def-santang3_cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=rp09_wis_hifiasm
#SBATCH --output=/home/sashanik/scratch/california_oaks/logs/rp09_wis_hifiasm_%j.log
set -euo pipefail
W=/home/sashanik/scratch/california_oaks/wislizeni_reasm; QC=$W/qc; H=$W/hifiasm; mkdir -p "$H"
[ -f "$QC/rp08.ok" ] || { echo "ERROR: rp08 read QC did not pass"; exit 1; }
source /home/sashanik/projects/def-santang3/sashanik/software/miniforge3/etc/profile.d/conda.sh; conda activate asm
ts() { date '+%F %T'; }
echo "[$(ts)] hifiasm $(hifiasm --version) on $(hostname)"
cd "$H"
hifiasm -o wis -t 64 --telo-m CCCTAAA --h1 "$(cat "$QC/omnic_r1.txt")" --h2 "$(cat "$QC/omnic_r2.txt")" \
    "$W/raw/SRR19243399_subreads.fastq.gz" 2> hifiasm.log || { tail -30 hifiasm.log; exit 1; }
tail -25 hifiasm.log
for g in hap1.p_ctg hap2.p_ctg p_ctg; do
    awk '/^S/{print ">"$2"\n"$3}' "wis.hic.$g.gfa" > "wis.hic.${g%.p_ctg}.fa"
    gfastats "wis.hic.${g%.p_ctg}.fa" --tabular > "wis.hic.${g%.p_ctg}.gfastats.tsv"
    echo "== ${g%.p_ctg}"; grep -E 'of contigs|Contig N50|Total contig length|Largest contig' "wis.hic.${g%.p_ctg}.gfastats.tsv"
done
touch rp09.ok
echo "[$(ts)] DONE rp09"
