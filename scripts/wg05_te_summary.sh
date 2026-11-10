#!/bin/bash
# wg05: TE section numbers and figures for the paper (see wg05_te_summary.py for methods); reads cg04 EDTA + cg05 Helixer.
# Output: results/te/, figures/te/
# Usage: sbatch [--export=ALL,MU=..,MU_ALT=..] scripts/wg05_te_summary.sh
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --job-name=wg05_te
#SBATCH --output=/home/sashanik/scratch/california_oaks/logs/wg05_te_%j.log
set -euo pipefail
cd ~/oak-genomics-wislizeni/california_oaks/wislizeni_genome
/home/sashanik/projects/def-santang3/sashanik/software/venv_plotting/bin/python scripts/wg05_te_summary.py results/te figures/te "${MU:-1.3e-8}" "${MU_ALT:-1.3e-8}"
