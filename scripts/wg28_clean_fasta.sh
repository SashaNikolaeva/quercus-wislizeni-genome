#!/bin/bash
# wg28:
# writes NEW cleaned FASTAs without the wg25 organelle candidate_remove contigs; the original assemblies are never modified.
# Then recomputes checklist items 11 (assembly size), 12 (contig/scaffold N50/NG50) and 13 (% in chromosomes) on the cleaned
# FASTAs with gfastats 1.3.11 and the same expected genome size as wg02 (739,879,015 bp, GenomeScope2 point estimate).
# Inputs : genomes/Qwislizeni_reasm_hap{1,2}.fa (+ .fai), results/organelle/organelle_removal_candidates_hap{1,2}.txt
# Output : genomes/Qwislizeni_reasm_hap{1,2}.clean.fa (+ .fai), results/assembly_clean/{hapN.gfastats.tsv, clean_summary.tsv}
# Checks (exit 1): every candidate ID exists in the original; none is a chromosome; cleaned seq count = original - n_candidates;
#                  cleaned bp = original bp - candidate bp; all 12 chromosomes byte-identical (seqkit sum) to the original.
# Usage  : sbatch scripts/wg28_clean_fasta.sh
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --job-name=wg28_clean
#SBATCH --output=/home/sashanik/scratch/california_oaks/paperA_runs/logs/wg28_%j.log
set -euo pipefail
WG=/home/sashanik/oak-genomics-wislizeni/california_oaks/wislizeni_genome
G=/home/sashanik/projects/def-santang3/sashanik/california_oaks/genomes
B=/home/sashanik/projects/def-santang3/sashanik/software/miniforge3/envs/asm/bin
O=$WG/results/assembly_clean; mkdir -p "$O"
GSIZE=739879015
printf "hap\tseqs_orig\tbp_orig\tn_removed\tbp_removed\tseqs_clean\tbp_clean\tbp_chr\tpct_in_chr_clean\n" > "$O/clean_summary.tsv.tmp"
for h in 1 2; do
    IN=$G/Qwislizeni_reasm_hap$h.fa; OUT=$G/Qwislizeni_reasm_hap$h.clean.fa
    L=$WG/results/organelle/organelle_removal_candidates_hap$h.txt
    echo "[$(date '+%F %T')] hap$h: $(grep -c . "$L") candidates"
    awk -F'\t' 'NR==FNR{fai[$1]=$2; next} NF{ if(!($1 in fai)){print "ERROR missing " $1; exit 1}
                 if($1 ~ /_chr[0-9]+$/){print "ERROR chromosome in list " $1; exit 1} }' "$IN.fai" "$L"
    "$B/seqkit" grep -v -n -r -f <(awk 'NF{print "^" $1 "$"}' "$L") "$IN" > "$OUT.tmp"
    mv "$OUT.tmp" "$OUT"; "$B/samtools" faidx "$OUT"
    read so bo <<< "$(awk '{n++; s+=$2} END{print n, s}' "$IN.fai")"
    read nr br <<< "$(awk -F'\t' 'NR==FNR{fai[$1]=$2; next} NF{n++; s+=fai[$1]} END{print n, s}' "$IN.fai" "$L")"
    read sc bc <<< "$(awk '{n++; s+=$2} END{print n, s}' "$OUT.fai")"
    bchr=$(awk '$1 ~ /_chr[0-9]+$/{s+=$2} END{print s}' "$OUT.fai")
    [ "$sc" -eq $((so - nr)) ] || { echo "ERROR hap$h seqs $sc != $so - $nr"; exit 1; }
    [ "$bc" -eq $((bo - br)) ] || { echo "ERROR hap$h bp $bc != $bo - $br"; exit 1; }
    d=$(diff <("$B/seqkit" grep -r -p '_chr[0-9]+$' "$IN" | "$B/seqkit" sum) \
             <("$B/seqkit" grep -r -p '_chr[0-9]+$' "$OUT" | "$B/seqkit" sum) | awk 'END{print NR}')
    [ "$d" -eq 0 ] || { echo "ERROR hap$h chromosome checksum differs"; exit 1; }
    "$B/gfastats" "$OUT" $GSIZE --tabular > "$O/hap$h.gfastats.tsv"
    awk -v h=hap$h -v so=$so -v bo=$bo -v nr=$nr -v br=$br -v sc=$sc -v bc=$bc -v bchr=$bchr \
        'BEGIN{printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\n", h, so, bo, nr, br, sc, bc, bchr, 100*bchr/bc}' >> "$O/clean_summary.tsv.tmp"
done
mv "$O/clean_summary.tsv.tmp" "$O/clean_summary.tsv"
cat "$O/clean_summary.tsv"; grep -iE 'Total scaffold length|Scaffold N50|Scaffold NG50|Contig N50|Contig NG50|# scaffolds|# contigs' "$O"/hap*.gfastats.tsv
echo "[$(date '+%F %T')] DONE"
