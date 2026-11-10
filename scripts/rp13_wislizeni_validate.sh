#!/bin/bash
# rp13: validation + chromosome naming of the YaHS scaffolds (hap1, hap2).
#  1. minimap2 -x asm20 scaffolds vs lobata12 -> name/orient the 12 largest as Qwis_<hap>_chrNN (lobata numbering),
#     rest Qwis_<hap>_unplaced_NNNN (scripts/rp13_synteny.py name)
#  2. named assembly vs Qparvula, Qagrifolia, old CCGP Qwislizeni hap1 + hap2, lobata12; hap1 vs hap2:
#     per-chromosome best partner / % collinear / % minus strand, broken-old-join table, dotplots (rp13_synteny.py report)
#  3. telomere repeat copies (CCCTAAA/TTTAGGG) in the 20 kb at each chromosome end (rp13_synteny.py telo)
#  4. before/after table (gfastats + compleasm/Merqury from rp10)
# Gates (both haps): 12 distinct lobata chromosomes assigned; each chromosome >= 30 Mb; >= 80% of each chromosome's
# lobata-aligned bp on its assigned chromosome. Pass -> copy to project genomes/Qwislizeni_reasm_<hap>.fa, rp13.pass.
# Output: ~/scratch/california_oaks/wislizeni_reasm/final/
# Usage: sbatch scripts/rp13_wislizeni_validate.sh
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=16:00:00
#SBATCH --job-name=rp13_wis_val
#SBATCH --output=/home/sashanik/scratch/california_oaks/logs/rp13_wis_val_%j.log
set -euo pipefail
W=/home/sashanik/scratch/california_oaks/wislizeni_reasm; Y=$W/yahs; F=$W/final; mkdir -p "$F"
P=/home/sashanik/projects/def-santang3/sashanik/california_oaks; G=$P/genomes; LOB=$P/phylogeny/ref/lobata12.fa
REPO=/home/sashanik/oak-genomics/california_oaks/wislizeni_genome; PY=$REPO/scripts/rp13_synteny.py
for h in hap1 hap2; do [ -f "$Y/$h/rp11.ok" ] || { echo "ERROR: rp11 not finished for $h"; exit 1; }; done
module load StdEnv/2023 minimap2/2.28
source /home/sashanik/projects/def-santang3/sashanik/software/miniforge3/etc/profile.d/conda.sh; conda activate hic
conda activate --stack /home/sashanik/projects/def-santang3/sashanik/software/miniforge3/envs/asm
ts() { date '+%F %T'; }
PYH=/home/sashanik/projects/def-santang3/sashanik/software/miniforge3/envs/hic/bin/python   # has numpy + matplotlib
mm() { [ -s "$3" ] && return 0; minimap2 -x asm20 -t 32 "$1" "$2" > "$3.tmp" 2>> "$F/minimap2.log"; mv "$3.tmp" "$3"; }   # target query out
cd "$F"

for h in hap1 hap2; do
    echo "[$(ts)] $h: align to lobata12 and name"
    mm "$LOB" "$Y/$h/${h}_yahs_scaffolds_final.fa" "${h}_yahs_vs_lobata12.paf"
    "$PYH" "$PY" name "$Y/$h/${h}_yahs_scaffolds_final.fa" "${h}_yahs_vs_lobata12.paf" "$LOB.fai" "Qwis_$h" "Qwis_$h.fa" "Qwis_$h.names.tsv"
    samtools faidx "Qwis_$h.fa"
done

echo "[$(ts)] synteny"
for h in hap1 hap2; do
    for ref in lobata12:$LOB Qparvula:$G/Qparvula.fa Qagrifolia:$G/Qagrifolia.fa old_wis_hap1:$G/Qwislizeni.fa old_wis_hap2:$G/Qwislizeni_hap2.fa; do
        lab=${ref%%:*}; fa=${ref#*:}
        mm "$fa" "Qwis_$h.fa" "Qwis_${h}_vs_$lab.paf"
        "$PYH" "$PY" report "Qwis_$h.fa.fai" "$lab" "Qwis_${h}_vs_$lab.paf" "Qwis_${h}_vs_$lab"
    done
done
mm Qwis_hap1.fa Qwis_hap2.fa Qwis_hap2_vs_new_hap1.paf
"$PYH" "$PY" report Qwis_hap2.fa.fai new_hap1 Qwis_hap2_vs_new_hap1.paf Qwis_hap2_vs_new_hap1

echo "[$(ts)] telomeres"
for h in hap1 hap2; do "$PYH" "$PY" telo "Qwis_$h.fa" "Qwis_$h.telomeres.tsv"; done

echo "[$(ts)] summaries"
for h in hap1 hap2; do
    echo "== $h names"; awk 'NR<=13' "Qwis_$h.names.tsv" | column -t   # was column|head: SIGPIPE killed job 6141602
    echo "== $h synteny (best partner, % of aligned bp on it, % minus, aligned/length)"
    cat Qwis_${h}_vs_*.chrom.tsv | awk -F'\t' 'NR==1 || $1!="chrom"' | column -t
    echo "== $h old CCGP scaffolds >= 1 Mb split over >= 2 new chromosomes (>= 10% each)"
    awk -F'\t' 'FNR>1 && $4>=2' Qwis_${h}_vs_old_wis_*.target.tsv | column -t
    echo "== $h telomeres (20 kb ends, capped = >= 100 copies)"; column -t "Qwis_$h.telomeres.tsv"
    awk -F'\t' 'NR>1{c+=$5+$6} END{print "capped ends: " c " of " 2*(NR-1)}' "Qwis_$h.telomeres.tsv"
done | tee summary.txt

echo "[$(ts)] before/after"
{ printf "assembly\tsequences\ttotal_Mb\tlargest12_Mb\tcontig_N50_Mb\tscaffold_N50_Mb\tgaps\n"
  for a in old_hap1:$G/Qwislizeni.fa old_hap2:$G/Qwislizeni_hap2.fa new_hap1:Qwis_hap1.fa new_hap2:Qwis_hap2.fa; do
      lab=${a%%:*}; fa=${a#*:}
      t12=$(sort -k2,2nr "$fa.fai" | awk 'NR<=12{s+=$2} END{printf "%.1f", s/1e6}')
      gfastats "$fa" --tabular | awk -F'\t' -v l=$lab -v t=$t12 '{v[$1]=$2} END{printf "%s\t%d\t%.1f\t%s\t%.1f\t%.1f\t%d\n",l,v["# scaffolds"],v["Total scaffold length"]/1e6,t,v["Contig N50"]/1e6,v["Scaffold N50"]/1e6,v["# gaps in scaffolds"]}'
  done; } > before_after.tsv
column -t before_after.tsv | tee -a summary.txt

echo "[$(ts)] gates"
{ printf "hap\tcheck\tvalue\tresult\n"
  for h in hap1 hap2; do
      n=$(awk -F'\t' 'NR>1 && $1~/_chr/{print $4}' "Qwis_$h.names.tsv" | sort -u | wc -l)
      printf "%s\tdistinct_lobata_chromosomes\t%d\t%s\n" $h $n $([ "$n" -eq 12 ] && echo PASS || echo FAIL)
      awk -F'\t' -v h=$h 'NR>1 && $1~/_chr/{if($3<3e7)s++; if($6<0.8)f++} END{printf "%s\tchromosomes_lt_30Mb\t%d\t%s\n%s\tchromosomes_lt_80pct_on_assigned\t%d\t%s\n",h,s+0,(s?"FAIL":"PASS"),h,f+0,(f?"FAIL":"PASS")}' "Qwis_$h.names.tsv"
  done; } > rp13_gates.tsv
column -t rp13_gates.tsv | tee -a summary.txt
if grep -q FAIL rp13_gates.tsv; then
    echo "GATES: at least one FAIL; assemblies NOT copied to genomes/"
else
    for h in hap1 hap2; do cp "Qwis_$h.fa" "$G/Qwislizeni_reasm_$h.fa"; samtools faidx "$G/Qwislizeni_reasm_$h.fa"; done
    touch rp13.pass; echo "GATES: all PASS; copied to $G/Qwislizeni_reasm_hap{1,2}.fa"
fi
echo "[$(ts)] DONE rp13"
