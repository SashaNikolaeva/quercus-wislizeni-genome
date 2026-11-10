#!/bin/bash
# wg26b: SLURM wrapper for wg26b_organelle_check.py (definitions, equations and the
# fixed plastid/mito success criteria are in that file's header). Summarises the three oatk runs of wg26
# (-c 150/250/400) per organelle: n sequences, circular flag, total length, GFA segments/length, depth CV along the path,
# identity + aligned fraction vs NC_062091.1 (plastid) / CP129458.1 (mito) with minimap2 -cx asm5, LSC/IR/SSC.
# Inputs: ~/scratch/california_oaks/paperA_runs/wg26/c{150,250,400}/wg26_c*.{pltd,mito}.{ctg.fasta,gfa}
# Output: results/organelle/oatk_summary.tsv
# Usage: sbatch --dependency=afterok:<wg26 array> scripts/wg26b_organelle_check.sh
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=0:30:00
#SBATCH --job-name=wg26b_orgcheck
#SBATCH --output=/home/sashanik/scratch/california_oaks/paperA_runs/logs/wg26b_%j.log
set -euo pipefail
module load StdEnv/2023 minimap2/2.28
WG=/home/sashanik/oak-genomics-wislizeni/california_oaks/wislizeni_genome
R=$WG/results/organelle; mkdir -p "$R"
/home/sashanik/projects/def-santang3/sashanik/software/venv_plotting/bin/python "$WG/scripts/wg26b_organelle_check.py" \
    --workdir /home/sashanik/scratch/california_oaks/paperA_runs/wg26 --cvals c150:150,c250:250,c400:400 \
    --out "$R/oatk_summary.tsv" --tmp "${SLURM_TMPDIR:?run under SLURM}"
echo "DONE wg26b"
