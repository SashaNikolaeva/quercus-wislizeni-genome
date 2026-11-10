#!/bin/bash
# rp13f: why are 2 of 48 wislizeni chromosome ends uncapped (hap1 chr07 left, hap2 chr01 right)?
# 1. telomere-repeat scan along the two chromosomes + unplaced contigs with telomeric ends (rp13f_telomere_scan.py)
# 2. which sequences cover the missing ends, from the existing rp13/rp13c PAFs (no new alignments):
#    hap1 unplaced vs hap2 chromosomes (asm5 -c, blocks >= 5 kb) on hap2 chr07 0-0.6 Mb, and hap2 sequences on hap1 chr01 >= 52 Mb.
# Login node OK (one chromosome scan + awk on PAFs). Output: results/wislizeni_reasm/rp13f_missing_telomeres.txt
set -euo pipefail
F=/home/sashanik/scratch/california_oaks/wislizeni_reasm/final; REPO=$(cd "$(dirname "$0")/.." && pwd)
PY=/home/sashanik/projects/def-santang3/sashanik/software/miniforge3/envs/hic/bin/python
{ "$PY" -I "$REPO/scripts/rp13f_telomere_scan.py" "$F/Qwis_hap1.fa" Qwis_hap1_chr07
  "$PY" -I "$REPO/scripts/rp13f_telomere_scan.py" "$F/Qwis_hap2.fa" Qwis_hap2_chr01
  echo "== hap2 chr07 vs hap1 chr07, first 1.5 Mb of hap1"
  awk -F'\t' '$1=="Qwis_hap2_chr07" && $6=="Qwis_hap1_chr07" && $8<1500000 && $11>=5000{printf "hap2 %.2f-%.2f -> hap1 %.2f-%.2f\n",$3/1e6,$4/1e6,$8/1e6,$9/1e6}' "$F/Qwis_hap2_vs_new_hap1.paf"
  echo "== hap1 unplaced on hap2 chr07 0-0.6 Mb"
  awk -F'\t' '$6=="Qwis_hap2_chr07" && $8<600000 && $11>=5000{printf "%s (%d kb) %d-%d -> hap2 chr07 %.3f-%.3f Mb ident %.3f\n",$1,$2/1e3,$3,$4,$8/1e6,$9/1e6,$10/$11}' "$F/rp13c_hap1_unpl_vs_other.paf" | sort -t'>' -k2
  echo "== hap2 sequences on hap1 chr01 >= 52 Mb"
  awk -F'\t' '$6=="Qwis_hap1_chr01" && $9>52000000 && $11>=5000{printf "%s %.2f-%.2f -> hap1 %.2f-%.2f\n",$1,$3/1e6,$4/1e6,$8/1e6,$9/1e6}' "$F/Qwis_hap2_vs_new_hap1.paf" | sort; } > "$REPO/results/wislizeni_reasm/rp13f_missing_telomeres.txt"
cat "$REPO/results/wislizeni_reasm/rp13f_missing_telomeres.txt"
