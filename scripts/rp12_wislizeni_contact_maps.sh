#!/bin/bash
# rp12: Omni-C contact maps of the YaHS scaffolds, one haplotype per array task.
#  1. final-scaffold map: juicer pre (yahs 1.2.2) -> sort -> juicer_tools 1.22.01 pre -> <hap>_final.hic   (view in Juicebox)
#  2. assembly-mode map:  juicer pre -a -> juicer_tools pre -> <hap>_JBAT.hic + .assembly  (manual curation in Juicebox JBAT)
#  3. PNG: 12 largest scaffolds, 1 Mb bins, log10(contacts+1) (scripts/rp12b_contact_png.py), from the step-1 text file
# Commands follow the YaHS README. MAPQ filter is juicer pre's default (10); the BAM was already MAPQ >= 30 (pairtools).
# Output: ~/scratch/california_oaks/wislizeni_reasm/yahs/<hap>/maps/
# Usage: sbatch --array=1-2 scripts/rp12_wislizeni_contact_maps.sh
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --tmp=200G
#SBATCH --job-name=rp12_wis_maps
#SBATCH --output=/home/sashanik/scratch/california_oaks/logs/rp12_wis_maps_%A_%a.log
set -euo pipefail
HAP=hap${SLURM_ARRAY_TASK_ID:?array task 1 or 2}
Y=/home/sashanik/scratch/california_oaks/wislizeni_reasm/yahs/$HAP; M=$Y/maps; mkdir -p "$M"
JT=/home/sashanik/projects/def-santang3/sashanik/software/juicer_tools/juicer_tools_1.22.01.jar
REPO=/home/sashanik/oak-genomics/california_oaks/wislizeni_genome
[ -f "$Y/rp11.ok" ] || { echo "ERROR: rp11 not finished for $HAP"; exit 1; }
module load java/17.0.6
source /home/sashanik/projects/def-santang3/sashanik/software/miniforge3/etc/profile.d/conda.sh; conda activate hic
export JAVA_TOOL_OPTIONS=-Xmx80g
ts() { date '+%F %T'; }
T=${SLURM_TMPDIR:-/tmp}/rp12_$HAP.$$; mkdir -p "$T"; trap 'rm -rf "$T"' EXIT
P=${HAP}_yahs; cd "$M"

echo "[$(ts)] final-scaffold map"
juicer pre "$Y/$P.bin" "$Y/${P}_scaffolds_final.agp" "$Y/$HAP.fa.fai" 2> juicer_pre_final.log \
  | sort -k2,2d -k6,6d -T "$T" --parallel=16 -S48G | awk 'NF' > "$T/final.txt"
cut -f1,2 "$Y/${P}_scaffolds_final.fa.fai" > "$HAP.final.chrom.sizes"
java -jar "$JT" pre --threads 16 "$T/final.txt" "${HAP}_final.hic" "$HAP.final.chrom.sizes" > juicer_tools_final.log 2>&1

echo "[$(ts)] PNG"
python "$REPO/scripts/rp12b_contact_png.py" "$T/final.txt" "$HAP.final.chrom.sizes" "${HAP}_contacts_top12_1Mb.png" "$HAP"

echo "[$(ts)] assembly-mode (JBAT) map"
juicer pre -a -o "${HAP}_JBAT" "$Y/$P.bin" "$Y/${P}_scaffolds_final.agp" "$Y/$HAP.fa.fai" > juicer_pre_JBAT.log 2>&1
java -jar "$JT" pre --threads 16 "${HAP}_JBAT.txt" "${HAP}_JBAT.hic" \
     <(grep PRE_C_SIZE juicer_pre_JBAT.log | awk '{print $2" "$3}') > juicer_tools_JBAT.log 2>&1
rm -f "${HAP}_JBAT.txt"   # large; regenerate with juicer pre -a if needed
ls -la "$M"
touch "$M/rp12.ok"
echo "[$(ts)] DONE rp12 $HAP"
