#!/bin/bash
# cg04: per-genome EDTA, the first stage of panEDTA.
# EDTA 2.2.2 (container software/containers/edta_2.2.2.sif, biocontainers 2.2.2--hdfd78af_1, pandas 2.2.3): EDTA.pl --species others --step all --sensitive 0 --anno 1.
# Input: the 12 named chromosomes only (as in cg05), with headers shortened to chr01..chr12 because EDTA needs short
# sequence names (map back via <G>.chrnames.tsv). The panEDTA stage (one consistent library, re-annotation of every genome,
# optionally with Helixer CDS) is a later script (cg04b) once all 11 EDTA runs finish.
# Output: ~/scratch/california_oaks/cg04_edta/<G>/ (EDTA outputs: *.mod.EDTA.TElib.fa, *.EDTA.TEanno.gff3, *.EDTA.TEanno.sum, ...), cg04.ok
# --cleanenv is REQUIRED: without it the host's exported `which` function (CC: --tty-only) breaks BusyBox `which` in the
# container, and EDTA/LTR_retriever then randomly fail to find cd-hit-est / rmblastn (array 6211789: LTR step 0 bp, one hang).
# Runs in node-local $SLURM_TMPDIR (2026-10-21): array 6395525 died in the LINE step because RepeatModeler round 5 writes
# ~283K batch files per genome (wislizeni: 282,872), and 4 concurrent tasks filled the /scratch 1M file quota on
# 2026-10-20 (RepeatModeler crashed after its 16 h all-by-other step; EDTA's `-pa` retry then overwrote repeatmodeler.log).
# Now: <W> is copied to local disk without its LINE/ dir, EDTA runs there (default --overwrite 0, so finished LTR/TIR/
# Helitron/SINE modules are reused), and results are copied back WITHOUT the RepeatModeler round-* dirs: on exit, on
# failure, and 30 min before the time limit (USR1).
# Usage: sbatch --array=0 -c 8 --time=6:00:00 scripts/cg04_edta.sh   (test: first 8 Mb of wislizeni chr01; LTR step must give > 0 bp)
#        sbatch --array=1-11%4 scripts/cg04_edta.sh
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=4-00:00:00
#SBATCH --tmp=150G
#SBATCH --signal=B:USR1@1800
#SBATCH --job-name=cg04_edta
#SBATCH --output=/home/sashanik/scratch/california_oaks/logs/cg04_edta_%A_%a.log
set -euo pipefail
GEN=(Qwislizeni Qagrifolia Qparvula Qrubra_full Qlobata Qengelmannii Qberberidifolia Qtomentella Qcerris Qilex Cmollissima)
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then G=EDTAtest; else G=${GEN[$((SLURM_ARRAY_TASK_ID-1))]}; fi
CA=/home/sashanik/projects/def-santang3/sashanik/california_oaks
SRC=$CA/genomes_named/$G.fa; [ "$G" = Qwislizeni ] || [ "$G" = EDTAtest ] && SRC=$CA/genomes/Qwislizeni_reasm_hap1.fa
SIF=/home/sashanik/projects/def-santang3/sashanik/software/containers/edta_2.2.2.sif   # NOT 2.3.0: its pandas 3.0.1 breaks TIR-Learner (x[0] KeyError, test 6286794)
W=/home/sashanik/scratch/california_oaks/cg04_edta/$G; mkdir -p "$W"
[ -f "$W/cg04.ok" ] && { echo "$G done"; exit 0; }
L=$SLURM_TMPDIR/$G; mkdir -p "$L"
copy_back() { rsync -a --exclude 'RM_*/round-*' --exclude 'RM_*/*.n??' "$L/" "$W/" || echo "WARNING: copy-back failed"; }
trap 'echo "[$(date "+%F %T")] USR1: time limit near, copying back"; copy_back; exit 1' USR1
trap 'copy_back' EXIT
rsync -a --exclude 'LINE/' --exclude 'edta.log' --exclude 'cg04.ok' "$W/" "$L/"
cd "$L"
module load StdEnv/2023 samtools/1.22.1 apptainer/1.4.5
if [ "$G" = EDTAtest ] && [ ! -s "$G.fa" ]; then
    samtools faidx "$SRC" Qwis_hap1_chr01:1-8000000 | sed -E 's/^>.*/>chr01/' > "$G.fa"; printf "chr01\tQwis_hap1_chr01:1-8000000\n" > "$G.chrnames.tsv"
fi
if [ ! -s "$G.fa" ]; then
    samtools faidx "$SRC" $(awk '$1~/_chr[0-9][0-9]$/{print $1}' "$SRC.fai") | sed -E 's/^>(.*_(chr[0-9][0-9]))$/>\2/' > "$G.fa"
    awk '$1~/_chr[0-9][0-9]$/{n=$1; sub(/.*_/,"",n); print n"\t"$1}' "$SRC.fai" > "$G.chrnames.tsv"
fi
echo "[$(date '+%F %T')] $G: $(grep -c '>' "$G.fa") chromosomes on $(hostname), local dir $L ($(df -h "$SLURM_TMPDIR" | awk 'NR==2{print $4}') free)"
echo "reused modules: $(ls -d "$G.fa.mod.EDTA.raw"/*.raw.fa 2>/dev/null | xargs -r -n1 basename | tr '\n' ' ')"
# EDTA in the background + wait, so the USR1 trap can fire while it runs
apptainer exec --cleanenv -B /lustre06,/lustre07,/home/sashanik,"$SLURM_TMPDIR" "$SIF" EDTA.pl --genome "$G.fa" --species others --step all \
    --sensitive 0 --anno 1 --threads "${SLURM_CPUS_PER_TASK:-32}" > edta.log 2>&1 &
wait $! || { tail -40 edta.log; exit 1; }
grep -E "LTR result file has 0 bp|LINE result file has 0 bp|is not found in the .* path|engine is not installed|^Error: |RepeatModeler did not run correctly" edta.log && { echo "ERROR: EDTA dependency/LTR/LINE failure (see edta.log)"; exit 1; }
ls -la "$G.fa.mod.EDTA.TElib.fa" "$G.fa.mod.EDTA.TEanno.sum" 2>/dev/null || echo "WARNING: expected EDTA outputs missing"
{ grep -A25 -i "Repeat Classes" "$G.fa.mod.EDTA.TEanno.sum" 2>/dev/null || true; } | awk 'NR<=30'
[ -s "$G.fa.mod.EDTA.TEanno.gff3" ] && [ -s "$G.fa.mod.EDTA.TElib.fa" ] || { echo "ERROR: EDTA outputs missing; not marking done"; exit 1; }
touch cg04.ok
echo "[$(date '+%F %T')] DONE cg04 $G"
