#!/bin/bash
# wg04b: NCBI FCS-adaptor (adaptor/vector screen, eukaryote mode) on
# assembly v2 hap1 and hap2 (whole assemblies, as they would be submitted). Container fcs-adaptor.sif (wg04a).
# Output: ~/scratch/california_oaks/wislizeni_reasm/wg04_fcs/adaptor_<hap>/ (fcs_adaptor_report.txt; empty report = clean)
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --array=1-2
#SBATCH --job-name=wg04b_fcs_adaptor
#SBATCH --output=/home/sashanik/scratch/california_oaks/logs/wg04b_fcs_adaptor_%A_%a.log
set -euo pipefail
HAP=hap$SLURM_ARRAY_TASK_ID
F=/home/sashanik/projects/def-santang3/sashanik/software/fcs
FA=/home/sashanik/projects/def-santang3/sashanik/california_oaks/genomes/Qwislizeni_reasm_$HAP.fa
O=/home/sashanik/scratch/california_oaks/wislizeni_reasm/wg04_fcs/adaptor_$HAP; mkdir -p "$O"
module load StdEnv/2023 apptainer/1.4.5
export APPTAINER_BIND=/lustre06,/lustre07,/home/sashanik
cp "$FA" "$SLURM_TMPDIR/Qwis_$HAP.fa"   # the runner binds the FASTA's directory into the container
cd "$SLURM_TMPDIR"
"$F/run_fcsadaptor.sh" --fasta-input "$SLURM_TMPDIR/Qwis_$HAP.fa" --output-dir "$O" --euk \
    --container-engine singularity --image "$F/fcs-adaptor.sif"
echo "report lines (excluding header): $(grep -vc '^#' "$O/fcs_adaptor_report.txt" || true)"
head -20 "$O/fcs_adaptor_report.txt"
