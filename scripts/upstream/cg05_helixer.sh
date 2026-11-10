#!/bin/bash
# cg05: Helixer gene annotation, same model on every genome.
# Helixer v0.3.7 (container software/containers/helixer_v0.3.7.sif), model land_plant_v0.3_a_0080.h5, with the
# recommended land-plant windowing: --subsequence-length 64152 --overlap-offset 32076 --overlap-core-length 48114.
# Input: the 12 named chromosomes ONLY (unplaced sequence dropped for every genome alike, so retained haplotigs and
# debris in unplaced contigs cannot inflate gene-family counts; the cost is losing the few real genes on unplaced contigs).
# Output: ~/scratch/california_oaks/cg05_helixer/<G>/<G>.helixer.gff3 (+ chromosome FASTA), cg05.ok
# Usage: sbatch --array=1%1 scripts/cg05_helixer.sh   (pilot), then --array=2-11
#SBATCH --account=def-santang3_gpu
#SBATCH --gpus-per-node=a100:1
#SBATCH --cpus-per-task=12
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --tmp=200G
#SBATCH --job-name=cg05_helixer
#SBATCH --output=/home/sashanik/scratch/california_oaks/logs/cg05_helixer_%A_%a.log
set -euo pipefail
GEN=(Qwislizeni Qagrifolia Qparvula Qrubra_full Qlobata Qengelmannii Qberberidifolia Qtomentella Qcerris Qilex Cmollissima)
G=${GEN[$((SLURM_ARRAY_TASK_ID-1))]}
CA=/home/sashanik/projects/def-santang3/sashanik/california_oaks
SRC=$CA/genomes_named/$G.fa; [ "$G" = Qwislizeni ] && SRC=$CA/genomes/Qwislizeni_reasm_hap1.fa
SW=/home/sashanik/projects/def-santang3/sashanik/software
SIF=$SW/containers/helixer_v0.3.7.sif; MODEL=$SW/helixer_models/land_plant/land_plant_v0.3_a_0080.h5
W=/home/sashanik/scratch/california_oaks/cg05_helixer/$G; mkdir -p "$W"; cd "$W"
[ -f cg05.ok ] && { echo "$G done"; exit 0; }
module load StdEnv/2023 samtools/1.22.1 apptainer/1.4.5
samtools faidx "$SRC" $(awk '$1~/_chr[0-9][0-9]$/{print $1}' "$SRC.fai") > "$G.chr.fa"; samtools faidx "$G.chr.fa"
echo "[$(date '+%F %T')] $G: $(wc -l < "$G.chr.fa.fai") chromosomes, $(awk '{s+=$2} END{printf "%.1f Mb", s/1e6}' "$G.chr.fa.fai") on $(hostname)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
T=${SLURM_TMPDIR:-/tmp}/helixer_$G; mkdir -p "$T"
apptainer exec --nv -B /lustre06,/lustre07,/home/sashanik,"$T" "$SIF" Helixer.py \
    --fasta-path "$W/$G.chr.fa" --lineage land_plant --model-filepath "$MODEL" \
    --subsequence-length 64152 --overlap-offset 32076 --overlap-core-length 48114 \
    --species "$G" --temporary-dir "$T" --gff-output-path "$W/$G.helixer.gff3" > helixer.log 2>&1 || { tail -30 helixer.log; exit 1; }
echo "genes: $(awk '$3=="gene"' "$G.helixer.gff3" | wc -l)  mRNAs: $(awk '$3=="mRNA"' "$G.helixer.gff3" | wc -l)"
touch cg05.ok
echo "[$(date '+%F %T')] DONE cg05 $G"
