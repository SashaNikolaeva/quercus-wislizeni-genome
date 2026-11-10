#!/bin/bash
# wg01: Helixer gene annotation of oak haplotypes (array, one haplotype per task), 12 chromosomes only, with the same
# container, model and land-plant windowing as upstream/cg05_helixer.sh. In the genome report only the Q. wislizeni
# haplotype 2 output (Qwislizeni_h2) is used; haplotype 1 is annotated by cg05.
# Input: /project/.../oak_sv_barrier/genomes/named/<H>.chr.fa
# Output: ~/scratch/california_oaks/wg01_helixer/<H>/<H>.helixer.gff3, wg01.ok
# Usage: sbatch --array=1-9 scripts/wg01_helixer_haps.sh
#SBATCH --account=def-santang3_gpu
#SBATCH --gpus-per-node=a100:1
#SBATCH --cpus-per-task=12
#SBATCH --mem=120G
#SBATCH --time=8:00:00
#SBATCH --tmp=200G
#SBATCH --job-name=wg01_helixer
#SBATCH --output=/home/sashanik/scratch/california_oaks/logs/wg01_helixer_%A_%a.log
set -euo pipefail
HAPS=(Qwislizeni_h1 Qwislizeni_h2 Qparvula_h1 Qparvula_h2 Qagrifolia_h1 Qagrifolia_h2 Qrubra Qilex_h1_red Qilex_h2_red)
H=${HAPS[$((SLURM_ARRAY_TASK_ID-1))]}
SRC=/project/def-santang3/sashanik/oak_sv_barrier/genomes/named/$H.chr.fa
SW=/home/sashanik/projects/def-santang3/sashanik/software
SIF=$SW/containers/helixer_v0.3.7.sif; MODEL=$SW/helixer_models/land_plant/land_plant_v0.3_a_0080.h5
W=/home/sashanik/scratch/california_oaks/wg01_helixer/$H; mkdir -p "$W"; cd "$W"
[ -f wg01.ok ] && { echo "$H done"; exit 0; }
module load StdEnv/2023 samtools/1.22.1 apptainer/1.4.5
n=$(wc -l < "$SRC.fai"); [ "$n" -eq 12 ] || { echo "ERROR: $SRC has $n sequences, expected 12"; exit 1; }
echo "[$(date '+%F %T')] $H: $n chromosomes, $(awk '{s+=$2} END{printf "%.1f Mb", s/1e6}' "$SRC.fai") on $(hostname)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
T=$SLURM_TMPDIR/helixer_$H; mkdir -p "$T"
apptainer exec --nv -B /lustre06,/lustre07,/home/sashanik,/project,"$T" "$SIF" Helixer.py \
    --fasta-path "$SRC" --lineage land_plant --model-filepath "$MODEL" \
    --subsequence-length 64152 --overlap-offset 32076 --overlap-core-length 48114 \
    --species "$H" --temporary-dir "$T" --gff-output-path "$W/$H.helixer.gff3" > helixer.log 2>&1 || { tail -30 helixer.log; exit 1; }
echo "genes: $(awk '$3=="gene"' "$H.helixer.gff3" | wc -l)  mRNAs: $(awk '$3=="mRNA"' "$H.helixer.gff3" | wc -l)"
touch wg01.ok
echo "[$(date '+%F %T')] DONE wg01 $H"
