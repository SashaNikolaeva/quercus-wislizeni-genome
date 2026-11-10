#!/bin/bash
# wg26: de novo organelle assembly from the HiFi reads with oatk 1.0
# (conda env software/miniforge3/envs/oatk, bioconda oatk-1.0-h577a1d6_1) and the OatkDB HMM profiles
# (software/OatkDB, git commit 75e8db0ac4a7d508a9a518d900876003ceb70737, 2026-01-08; db v20230921,
# embryophyta_mito.fam + embryophyta_pltd.fam).
# Command: oatk -k 1001 -c <c> -t 32 -m embryophyta_mito.fam -p embryophyta_pltd.fam -o wg26_c<c>
#   array task 1/2/3 -> -c 150/250/400 (minimum k-mer coverage; all well above the ~47x nuclear coverage, so nuclear
#   k-mers drop out and organelle k-mers remain).
# Success criteria:
#   plastid (check only): one circular sequence of 150-170 kb, >= 99.5% identity over >= 95% of NC_062091.1
#     (either SSC orientation accepted);
#   mito (deposit candidate): in >= 1 -c value, one circular sequence of 300-800 kb, depth CV < 0.3, >= 80% of
#     CP129458.1 aligned; if several pass use -c 250 if it passes, else the lowest passing -c.
#   If no -c passes: "mitochondrial graph of X segments / Y kb" in the supplement, not deposited.
# Inputs: ~/scratch/california_oaks/wislizeni_reasm/raw/SRR19243399_subreads.fastq.gz (HiFi, 23.7 GB).
# Outputs: ~/scratch/california_oaks/paperA_runs/wg26/c<c>/wg26_c<c>.* (oatk GFA/FASTA/BED/annotation); nothing in
#   results/ (wg26b writes results/organelle/oatk_summary.tsv).
# Smoke (SMOKE=1, no array): first 400,000 FASTQ lines (100,000 reads), -c 20; tests the code path only; output in
#   <work>/smoke/; then runs wg26b on the smoke output (into smoke/, non-fatal) to test that code path too.
# Usage: SMOKE=1 sbatch --cpus-per-task=8 --mem=16G --time=1:00:00 scripts/wg26_oatk.sh            # smoke
#        sbatch --array=1-3 --dependency=afterok:<smoke> scripts/wg26_oatk.sh                       # full
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --job-name=wg26_oatk
#SBATCH --output=/home/sashanik/scratch/california_oaks/paperA_runs/logs/wg26_%A_%a.log
set -euo pipefail
SW=/home/sashanik/projects/def-santang3/sashanik/software
ENV=$SW/miniforge3/envs/oatk
DB=$SW/OatkDB/v20230921
FQ=/home/sashanik/scratch/california_oaks/wislizeni_reasm/raw/SRR19243399_subreads.fastq.gz
W=/home/sashanik/scratch/california_oaks/paperA_runs/wg26
T=${SLURM_TMPDIR:?run under SLURM}
TH=${SLURM_CPUS_PER_TASK:-32}
export PATH="$ENV/bin:$PATH"   # oatk calls nhmmscan from PATH (hmmer in the same env)
SMOKE=${SMOKE:-0}
if [ "$SMOKE" = 1 ]; then
    C=20; OUT=$W/smoke; mkdir -p "$OUT"
    # 400,000 lines = 100,000 FASTQ records; awk exits early (no `| head` under pipefail)
    zcat "$FQ" | awk 'NR<=400000{print; next}{exit}' > "$T/smoke.fq" || true
    n=$(wc -l < "$T/smoke.fq"); [ "$n" = 400000 ] || { echo "smoke FASTQ has $n lines, expected 400000" >&2; exit 1; }
    IN=$T/smoke.fq
else
    CS=(150 250 400); C=${CS[$((SLURM_ARRAY_TASK_ID-1))]}; OUT=$W/c$C; mkdir -p "$OUT"; IN=$FQ
fi
cd "$OUT"
echo "oatk $(oatk --version 2>&1) -c $C -t $TH input $IN"
mkdir -p "$T/oatk_tmp"
oatk -k 1001 -c "$C" -t "$TH" -m "$DB/embryophyta_mito.fam" -p "$DB/embryophyta_pltd.fam" \
     -T "$T/oatk_tmp" -o "wg26_c$C" "$IN" > "wg26_c$C.oatk.log" 2>&1
ls -l "$OUT"
if [ "$SMOKE" = 1 ]; then
    module load StdEnv/2023 minimap2/2.28
    "$SW/venv_plotting/bin/python" /home/sashanik/oak-genomics-wislizeni/california_oaks/wislizeni_genome/scripts/wg26b_organelle_check.py \
        --workdir "$W" --cvals smoke:20 --out "$OUT/oatk_summary_smoke.tsv" --tmp "$T" \
        || echo "WARNING: wg26b on the smoke output failed (code path to fix before wg26b runs)"
fi
echo "DONE wg26 c=$C (SMOKE=$SMOKE)"
