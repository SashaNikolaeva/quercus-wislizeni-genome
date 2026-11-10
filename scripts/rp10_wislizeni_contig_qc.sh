#!/bin/bash
# rp10: contig QC of the hifiasm Hi-C assemblies (hap1, hap2, primary p_ctg).
#  - gfastats --tabular (from rp09) + small-contig tally (<100 kb, >=1 Mb)
#  - compleasm 0.2.10 (miniprot) with embryophyta_odb12 and eudicotyledons_odb12 (no Fagales lineage in odb12);
#    lineages pre-downloaded to project software/compleasm_lineages (compute nodes have no internet)
#  - Merqury on the hifi_k21.meryl db from rp08: hap1+hap2 as a pair (QV, k-mer completeness, spectra-cn), p_ctg alone
# Gates (WISLIZENI_ASSEMBLY_PARAMETERS.md step 3): each hap 700-800 Mb, contig N50 > 10 Mb; compleasm complete >= 97%
# and duplicated < 5% (embryophyta); Merqury QV > 50. Results -> qc/rp10_gates.tsv; qc/rp10.pass only if all pass.
# A failed gate does not fail the job (exit 0) so the full report is always written.
# Usage: sbatch scripts/rp10_wislizeni_contig_qc.sh
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --job-name=rp10_wis_qc
#SBATCH --output=/home/sashanik/scratch/california_oaks/logs/rp10_wis_qc_%j.log
set -euo pipefail
W=/home/sashanik/scratch/california_oaks/wislizeni_reasm; H=$W/hifiasm; QC=$W/qc/rp10; mkdir -p "$QC"
LIB=/home/sashanik/projects/def-santang3/sashanik/software/compleasm_lineages
CONDA=/home/sashanik/projects/def-santang3/sashanik/software/miniforge3/etc/profile.d/conda.sh
[ -f "$H/rp09.ok" ] || { echo "ERROR: rp09 hifiasm not finished"; exit 1; }
source "$CONDA"
ts() { date '+%F %T'; }
ASM="hap1 hap2 p_ctg"

echo "[$(ts)] contig size classes"
for a in $ASM; do
    awk -v a=$a '/^>/{if(l)print l; l=0; next}{l+=length($0)} END{print l}' "$H/wis.hic.$a.fa" |
      awk -v a=$a '{n++; s+=$1} $1<100000{sn++; ss+=$1} $1>=1000000{bn++; bs+=$1}
        END{printf "%s\tcontigs %d\t%.1f Mb\t<100kb: %d (%.1f Mb)\t>=1Mb: %d (%.1f Mb)\n",a,n,s/1e6,sn,ss/1e6,bn,bs/1e6}'
done | tee "$QC/size_classes.txt"

echo "[$(ts)] compleasm"
conda activate compleasm
for a in $ASM; do for lin in embryophyta eudicotyledons; do
    o=$QC/compleasm_${a}_$lin
    [ -s "$o/summary.txt" ] || compleasm run -a "$H/wis.hic.$a.fa" -o "$o" -l $lin -L "$LIB" -t 16
    echo "== $a $lin"; cat "$o/summary.txt"
done; done
conda deactivate

echo "[$(ts)] Merqury"
conda activate asm
export MERQURY=${MERQURY:-$CONDA_PREFIX/share/merqury}
cd "$QC"
[ -s wis.qv ] || merqury.sh "$W/hifi_k21.meryl" "$H/wis.hic.hap1.fa" "$H/wis.hic.hap2.fa" wis > merqury_wis.log 2>&1
[ -s wis_p.qv ] || merqury.sh "$W/hifi_k21.meryl" "$H/wis.hic.p_ctg.fa" wis_p > merqury_wis_p.log 2>&1
echo "== QV (asm, kmers only in asm, kmers in both, QV, error rate)"; cat wis.qv wis_p.qv
echo "== k-mer completeness (asm, all, found, total, %)"; cat wis.completeness.stats wis_p.completeness.stats

echo "[$(ts)] gates"
# compleasm summary lines look like "S:95.10%, 1234"; C = S + D
pct() { awk -v k="$2" -F'[:%,]' '$1==k{print $2}' "$1"; }
{ printf "assembly\tcheck\tvalue\tthreshold\tresult\n"
  for a in hap1 hap2; do
      read -r len n50 < <(awk -F'\t' '$1=="Total contig length"{l=$2} $1=="Contig N50"{n=$2} END{print l, n}' "$H/wis.hic.$a.gfastats.tsv")
      awk -v a=$a -v v=$len 'BEGIN{printf "%s\tlength_bp\t%d\t700-800 Mb\t%s\n",a,v,(v>=7e8&&v<=8e8)?"PASS":"FAIL"}'
      awk -v a=$a -v v=$n50 'BEGIN{printf "%s\tcontig_N50_bp\t%d\t>10 Mb\t%s\n",a,v,(v>1e7)?"PASS":"FAIL"}'
      s=$(pct "$QC/compleasm_${a}_embryophyta/summary.txt" S); d=$(pct "$QC/compleasm_${a}_embryophyta/summary.txt" D)
      awk -v a=$a -v s=$s -v d=$d 'BEGIN{printf "%s\tcompleasm_embryophyta_complete_pct\t%.2f\t>=97\t%s\n",a,s+d,(s+d>=97)?"PASS":"FAIL";
                                        printf "%s\tcompleasm_embryophyta_duplicated_pct\t%.2f\t<5\t%s\n",a,d,(d<5)?"PASS":"FAIL"}'
      q=$(awk -v f="wis.hic.$a" '$1==f{print $4}' wis.qv)
      awk -v a=$a -v q=$q 'BEGIN{printf "%s\tmerqury_QV\t%.2f\t>50\t%s\n",a,q,(q>50)?"PASS":"FAIL"}'
  done; } > rp10_gates.tsv
column -t rp10_gates.tsv
if grep -q FAIL rp10_gates.tsv; then echo "GATES: at least one FAIL (see rp10_gates.tsv)"; else touch rp10.pass; echo "GATES: all PASS"; fi
echo "[$(ts)] DONE rp10"
