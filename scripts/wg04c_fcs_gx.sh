#!/bin/bash
# wg04c: NCBI FCS-GX foreign contamination screen of assembly v2 hap1 and
# hap2 (whole assemblies), tax-id 67830 (Quercus wislizeni, NCBI taxonomy). Database 'all' (wg04a, release 2023-01-24)
# is copied to /dev/shm first, as NCBI recommends (GX needs the ~499 GB db in RAM to run in about an hour per genome).
# Large-memory node (Narval cpularge, 2 TB). Output: ~/scratch/.../wg04_fcs/gx_<hap>/ (*.fcs_gx_report.txt, *.taxonomy.rpt)
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=48
#SBATCH --mem=750G
#SBATCH --time=02:00:00   # 6713997: db copy 5 min, hap1 screen ~5 min
#SBATCH --job-name=wg04c_fcs_gx
#SBATCH --output=/home/sashanik/scratch/california_oaks/logs/wg04c_fcs_gx_%j.log
set -euo pipefail
F=/home/sashanik/projects/def-santang3/sashanik/software/fcs
G=/home/sashanik/projects/def-santang3/sashanik/california_oaks/genomes
O=/home/sashanik/scratch/california_oaks/wislizeni_reasm/wg04_fcs; mkdir -p "$O"
[ -f "$F/gxdb/download.ok" ] || { echo "ERROR: gxdb download not verified (wg04a)"; exit 1; }
module load StdEnv/2023 python/3.11 apptainer/1.4.5
SHM=/dev/shm/gxdb_$SLURM_JOB_ID; mkdir -p "$SHM"
trap 'rm -rf "$SHM"' EXIT
echo "[$(date '+%F %T')] copying gxdb to $SHM"; cp "$F"/gxdb/all.* "$SHM/"; echo "[$(date '+%F %T')] copied"
export FCS_DEFAULT_IMAGE="$F/fcs-gx.sif" APPTAINER_BIND=/lustre06,/lustre07,/home/sashanik,/dev/shm
for HAP in hap1 hap2; do
    ls "$O/gx_$HAP"/*.fcs_gx_report.txt >/dev/null 2>&1 && { echo "$HAP already screened, skip"; continue; }
    echo "[$(date '+%F %T')] screen $HAP"
    python3 "$F/fcs.py" screen genome --fasta "$G/Qwislizeni_reasm_$HAP.fa" --out-dir "$O/gx_$HAP" \
        --gx-db "$SHM/all" --tax-id 67830
    R=$(ls "$O/gx_$HAP"/*.fcs_gx_report.txt); echo "$HAP actions: $(grep -vc '^#' "$R" || true)"; { grep -v "^#" "$R" || true; } | awk "NR<=30"   # not head: SIGPIPE under pipefail; || true: a clean report has no non-# lines (6713997 died on that)
done
echo "[$(date '+%F %T')] DONE"
