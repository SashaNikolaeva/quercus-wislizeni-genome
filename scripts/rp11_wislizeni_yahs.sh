#!/bin/bash
# rp11: Omni-C scaffolding with YaHS, one haplotype per array task (1 = hap1, 2 = hap2).
# Omni-C = paired SRR19243397 + SRR19243398 _1/_2 (lists from rp08 qc/omnic_r{1,2}.txt), streamed with zcat.
# Alignment follows Dovetail's Omni-C recommendations (parameters + rationale in WISLIZENI_ASSEMBLY_PARAMETERS.md step 4):
#   bwa-mem2 -5SP -T0 | pairtools parse (MAPQ 30, 5unique, gap 30) | sort | dedup | split -> name-sorted BAM -> yahs
# yahs: defaults (no -e: Omni-C has no restriction site; contig error correction left on).
# Submitted in parallel with rp10; scaffolds are only used if the rp10 gates pass.
# Output: ~/scratch/california_oaks/wislizeni_reasm/yahs/<hap>/ (<hap>_yahs_scaffolds_final.{fa,agp}, .bin, dedup.stats,
#         mapped.pairs.gz, <hap>.omnic.nsort.bam, gfastats, rp11.ok)
# Usage: sbatch --array=1-2 scripts/rp11_wislizeni_yahs.sh
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=16:00:00
#SBATCH --tmp=300G
#SBATCH --job-name=rp11_wis_yahs
#SBATCH --output=/home/sashanik/scratch/california_oaks/logs/rp11_wis_yahs_%A_%a.log
set -euo pipefail
HAP=hap${SLURM_ARRAY_TASK_ID:?array task 1 or 2}
W=/home/sashanik/scratch/california_oaks/wislizeni_reasm; QC=$W/qc; O=$W/yahs/$HAP; mkdir -p "$O"
[ -f "$W/hifiasm/rp09.ok" ] || { echo "ERROR: rp09 hifiasm not finished"; exit 1; }
[ -f "$O/rp11.ok" ] && { echo "$HAP already done"; exit 0; }
source /home/sashanik/projects/def-santang3/sashanik/software/miniforge3/etc/profile.d/conda.sh; conda activate hic
ts() { date '+%F %T'; }
T=${SLURM_TMPDIR:-/tmp}/rp11_$HAP.$$; mkdir -p "$T"; trap 'rm -rf "$T"' EXIT
IFS=, read -r -a R1 < "$QC/omnic_r1.txt"; IFS=, read -r -a R2 < "$QC/omnic_r2.txt"
echo "=== $HAP on $(hostname): R1 ${R1[*]} | R2 ${R2[*]}"

cd "$O"
FA=$HAP.fa
[ -s "$FA" ] || cp "$W/hifiasm/wis.hic.$HAP.fa" "$FA"
echo "[$(ts)] index"
[ -s "$FA.bwt.2bit.64" ] || bwa-mem2 index "$FA" > index.log 2>&1
samtools faidx "$FA"; cut -f1,2 "$FA.fai" > "$HAP.genome"

echo "[$(ts)] bwa-mem2 -5SP -T0 | pairtools parse | sort | dedup | split"
bwa-mem2 mem -5SP -T0 -t 32 "$FA" <(zcat "${R1[@]}") <(zcat "${R2[@]}") 2> bwa.log \
 | pairtools parse --min-mapq 30 --walks-policy 5unique --max-inter-align-gap 30 \
       --nproc-in 8 --nproc-out 8 --chroms-path "$HAP.genome" \
 | pairtools sort --nproc 16 --tmpdir "$T" \
 | pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats dedup.stats \
 | pairtools split --nproc-in 8 --nproc-out 8 --output-pairs mapped.pairs.gz --output-sam - \
 | samtools sort -n -@ 16 -m 2G -T "$T/s" -o "$HAP.omnic.nsort.bam" -
echo "== Omni-C library QC (Dovetail guidance: cis >1 kb >= 20% of unique pairs)"
awk -F'\t' '{v[$1]=$2} END{u=v["total_nodups"]; printf "total %d  mapped %d  dups %d (%.1f%% of mapped)  unique %d\n",v["total"],v["total_mapped"],v["total_dups"],100*v["total_dups"]/v["total_mapped"],u;
  printf "cis %.1f%%  cis>1kb %.1f%%  cis>10kb %.1f%%  trans %.1f%% (of unique)\n",100*v["cis"]/u,100*v["cis_1kb+"]/u,100*v["cis_10kb+"]/u,100*v["trans"]/u}' dedup.stats | tee omnic_qc.txt

echo "[$(ts)] yahs $(yahs --version)"
yahs "$FA" "$HAP.omnic.nsort.bam" -o "${HAP}_yahs" > yahs.log 2>&1 || { tail -30 yahs.log; exit 1; }
F=${HAP}_yahs_scaffolds_final.fa
conda activate --stack /home/sashanik/projects/def-santang3/sashanik/software/miniforge3/envs/asm
gfastats "$F" --tabular > "${HAP}_yahs.gfastats.tsv"
grep -E '^# scaffolds|Total scaffold length|Scaffold N50|Largest scaffold|^# contigs|Contig N50|^# gaps' "${HAP}_yahs.gfastats.tsv"
samtools faidx "$F"; sort -k2,2nr "$F.fai" | awk 'NR<=12{s+=$2} {t+=$2} END{printf "largest 12 scaffolds: %.1f of %.1f Mb (%.1f%%)\n",s/1e6,t/1e6,100*s/t}' | tee -a omnic_qc.txt
touch rp11.ok
echo "[$(ts)] DONE rp11 $HAP"
