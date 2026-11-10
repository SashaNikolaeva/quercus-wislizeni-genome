#!/bin/bash
# wg22: tRNA genes of the wislizeni haplotypes with tRNAscan-SE 2.0.12, eukaryotic
# Infernal mode, on the 12 chromosomes only (same scope as Helixer; unplaced sequence is not annotated).
# One array task per haplotype x chromosome: task t -> hap = 1 + (t > 12), chr = ((t - 1) mod 12) + 1.
# The chromosome is extracted with samtools faidx into $SLURM_TMPDIR.
# Command: tRNAscan-SE -E -H --thread 4 --detail -Q -o chrNN.out -f chrNN.ss -s chrNN.iso -m chrNN.stats -b chrNN.bed
# (the original spec wrote "-o# -f# -s# -m# -b#"; explicit file names are used instead so that no -p prefix is needed,
#    and -Q (--forceow) so a rerun never stops at an overwrite prompt; neither changes the predictions)
# then: EukHighConfidenceFilter -i chrNN.out -s chrNN.ss -o hc/ -p chrNN   (tag mode, no -r: filtered hits stay in the
#   file with a Note tag; "high confidence set" in Note = high-confidence tRNA). Default cutoffs (from the 2.0.12 script,
#   read 2026-11-03): cmscore1 (-c1) 50, ssscore1 (-m1) 10, isoscore1 (-e1) 70, isoscore2 (-e2) 70, isomaxscore2 (-x) 95,
#   anticodon count cutoff 40 (hard-coded).
# Flags confirmed from the installed EukHighConfidenceFilter (Getopt: result|i, ss|s, output|o, prefix|p, remove|r,
#   cmscore1|c1, ssscore1|m1, isoscore1|e1, isoscore2|e2, isomaxscore2|x, help|h); the smoke run also prints -h to its log.
# NOTE: the anticodon/isotype tertiary filter uses counts over its input, so a per-chromosome filter is not the same as a
#   genome-wide one. The per-chromosome filter here is a check; wg22b re-runs the filter once per haplotype on the
#   concatenated 12-chromosome output, and that haplotype-level set is the reported high-confidence count.
# Inputs: /home/sashanik/projects/def-santang3/sashanik/california_oaks/genomes/Qwislizeni_reasm_hap{1,2}.fa (+ .fai)
# Outputs: ~/scratch/california_oaks/paperA_runs/wg22/hap{1,2}/Qwis_hap{1,2}_chrNN.{out,ss,iso,stats,bed} and hc/*
#   SMOKE=1: first 5 Mb of hap2 chr12 (renamed to Qwis_hap2_chr12, coordinates unchanged) -> wg22/smoke/; never results/
# Usage: sbatch --export=ALL,SMOKE=1 --array=1 --time=1:00:00 scripts/wg22_trnascan.sh   (smoke)
#        sbatch --array=1-24 --dependency=afterok:<smoke> scripts/wg22_trnascan.sh         (full)
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=8:00:00
#SBATCH --job-name=wg22_trnascan
#SBATCH --output=/home/sashanik/scratch/california_oaks/paperA_runs/logs/wg22_%A_%a.log
set -euo pipefail
SMOKE=${SMOKE:-0}
GD=/home/sashanik/projects/def-santang3/sashanik/california_oaks/genomes
W=/home/sashanik/scratch/california_oaks/paperA_runs/wg22
module load StdEnv/2023 gcc/12.3 trnascan-se/2.0.12 samtools/1.22.1
export TMPDIR=${SLURM_TMPDIR:-/tmp}   # tRNAscan-SE temp_dir (/tmp in the conf) is overridden by TMPDIR
t=${SLURM_ARRAY_TASK_ID:?array task id required}
if [ "$SMOKE" = 1 ]; then
    hap=2; chr=12; O=$W/smoke
else
    hap=$(( 1 + (t > 12) )); chr=$(( (t - 1) % 12 + 1 )); O=$W/hap$hap
fi
c=$(printf 'Qwis_hap%d_chr%02d' "$hap" "$chr")
G=$GD/Qwislizeni_reasm_hap$hap.fa
mkdir -p "$O/hc"
grep -qP "^${c}\t" "$G.fai" || { echo "ERROR: $c not in $G.fai"; exit 1; }
FA=$SLURM_TMPDIR/$c.fa
if [ "$SMOKE" = 1 ]; then
    samtools faidx "$G" "$c:1-5000000" | sed "1s/.*/>$c/" > "$FA"
    echo "== smoke: tRNAscan-SE -h ==";            tRNAscan-SE -h > "$O/tRNAscan-SE_help.txt" 2>&1 || true; cat "$O/tRNAscan-SE_help.txt"
    echo "== smoke: EukHighConfidenceFilter -h =="; EukHighConfidenceFilter -h > "$O/EukHighConfidenceFilter_help.txt" 2>&1 || true
    cat "$O/EukHighConfidenceFilter_help.txt"
else
    samtools faidx "$G" "$c" > "$FA"
fi
echo "task $t hap $hap $c length $(awk '!/^>/{n+=length($0)} END{print n}' "$FA")"
cd "$O"
rm -f "$c".{out,ss,iso,stats,bed}
# -H (--breakdown) writes the HMM and 2'Str score columns that EukHighConfidenceFilter reads. Fix 2026-11-03: without
# -H the filter saw hmm_score/ss_score as undefined and returned 0 high-confidence tRNAs (smoke 8444915).
tRNAscan-SE -E -H --thread 4 --detail -Q -o "$c.out" -f "$c.ss" -s "$c.iso" -m "$c.stats" -b "$c.bed" "$FA"
n=$(awk '!/^(Sequence|Name|-----)/ && NF' "$c.out" | wc -l)
echo "tRNAscan-SE hits: $n"
if [ "$n" -gt 0 ]; then
    # guard: the filter needs the HMM and 2'Str columns, named in the first header line ("Sequence ... HMM 2'Str ...")
    # (fix 2026-11-03: the first guard looked at the "Name" line and failed smoke 8445761 although -H had worked)
    awk -F'\t' '/^Sequence/{f=1; if ($0 !~ /HMM/ || $0 !~ /Str/) {print "ERROR: no HMM/2Str columns in .out (needs -H)"; exit 1}}
                END{if (!f) {print "ERROR: no Sequence header line in .out"; exit 1}}' "$c.out"
    EukHighConfidenceFilter -i "$c.out" -s "$c.ss" -o "$O/hc" -p "$c"
    echo "per-chromosome high confidence: $(grep -c 'high confidence set' "$O/hc/$c.out" || true)"
else
    echo "zero tRNAscan-SE hits on $c; filter skipped (a zero is reported as a zero)"
fi
echo "DONE wg22 $c"
