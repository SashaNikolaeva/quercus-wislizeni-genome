#!/bin/bash
# wg22b: summary of tRNAs (wg22) and rRNAs (wg23), with the NUPT/NUMT flag.
# Step 1 (haplotype-level tRNA filter): concatenate the 12 per-chromosome tRNAscan-SE .out / .ss files of each haplotype
#   and run EukHighConfidenceFilter once on them (tag mode, default cutoffs, see wg22 header). The filter's anticodon
#   count cutoff (40) and isotype score cutoffs are computed over its input, so the haplotype-level set is the reported
#   one; the per-chromosome sets from wg22 are kept as a check column.
# Step 2 (NUPT/NUMT blocks): minimap2 -cx asm20 -N 1000 -p 0.1 -t 8, queries = NC_062091.1 (Q. wislizeni plastid) and
#   CP129458.1 (Q. variabilis mito), target = that haplotype's 12 chromosomes. Blocks kept if PAF col11 (alignment block
#   length) >= 500 bp. Primary and secondary alignments are both kept (-N 1000 is meant to find every copy).
# Step 3 (wg22b_ncrna_summary.py): a tRNA is plastid-like (mito-like) if >= 50% of its length overlaps the merged
#   plastid (mito) blocks; organelle-like = plastid-like OR mito-like. See the .py header for every column definition.
# Inputs: wg22/hap{1,2}/*.{out,ss}, wg23/hap{1,2}.barrnap.gff3, assemblies G (.fai), wg09_organelle/{NC_062091,CP129458}.fa,
#   results/rdna/rdna_kmer_copy_number_summary.tsv
# Outputs: results/annotation/trna_summary.tsv, results/annotation/rrna_summary.tsv, results/annotation/ncrna_checks.tsv;
#   work: ~/scratch/california_oaks/paperA_runs/wg22b/ (filtered .out, PAFs, per-tRNA flag table)
#   SMOKE=1: hap2 only, inputs from wg22/smoke + wg23/smoke, target = hap2 chr12 -> wg22b/smoke/; never results/
# Usage: sbatch --export=ALL,SMOKE=1 --time=0:30:00 --dependency=afterok:<wg22 smoke>:<wg23 smoke> scripts/wg22b_ncrna_summary.sh
#        sbatch --dependency=afterok:<wg22>:<wg23> scripts/wg22b_ncrna_summary.sh
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --job-name=wg22b_ncrna
#SBATCH --output=/home/sashanik/scratch/california_oaks/paperA_runs/logs/wg22b_%j.log
set -euo pipefail
SMOKE=${SMOKE:-0}
WG=/home/sashanik/oak-genomics-wislizeni/california_oaks/wislizeni_genome
PY=/home/sashanik/projects/def-santang3/sashanik/software/venv_plotting/bin/python
GD=/home/sashanik/projects/def-santang3/sashanik/california_oaks/genomes
ORG=/home/sashanik/scratch/california_oaks/wislizeni_reasm/final/wg09_organelle
P=/home/sashanik/scratch/california_oaks/paperA_runs
module load StdEnv/2023 gcc/12.3 trnascan-se/2.0.12 samtools/1.22.1 minimap2/2.28
export TMPDIR=${SLURM_TMPDIR:-/tmp}
if [ "$SMOKE" = 1 ]; then
    HAPS="2"; O=$P/wg22b/smoke; R=$O
else
    HAPS="1 2"; O=$P/wg22b; R=$WG/results/annotation
fi
mkdir -p "$O" "$R"
for hap in $HAPS; do
    if [ "$SMOKE" = 1 ]; then
        T=$P/wg22/smoke; CHRS="12"; BARR=$P/wg23/smoke/hap2_chr12.barrnap.gff3
    else
        T=$P/wg22/hap$hap; CHRS="01 02 03 04 05 06 07 08 09 10 11 12"; BARR=$P/wg23/hap$hap.barrnap.gff3
    fi
    [ -s "$BARR" ] || { echo "ERROR: missing $BARR"; exit 1; }
    : > "$O/hap$hap.trnascan.out"; : > "$O/hap$hap.trnascan.ss"; : > "$O/hap$hap.chrlist"
    for c in $CHRS; do
        id=Qwis_hap${hap}_chr$c
        for x in out ss; do [ -e "$T/$id.$x" ] || { echo "ERROR: missing $T/$id.$x"; exit 1; }; done
        cat "$T/$id.out" >> "$O/hap$hap.trnascan.out"; cat "$T/$id.ss" >> "$O/hap$hap.trnascan.ss"
        echo "$id" >> "$O/hap$hap.chrlist"
    done
    # step 1: haplotype-level high-confidence filter
    rm -rf "$O/hc_hap$hap"
    EukHighConfidenceFilter -i "$O/hap$hap.trnascan.out" -s "$O/hap$hap.trnascan.ss" -o "$O/hc_hap$hap" -p "hap$hap"
    # per-chromosome filter outputs from wg22 (check column)
    : > "$O/hap$hap.hc_perchr.out"
    for c in $CHRS; do f=$T/hc/Qwis_hap${hap}_chr$c.out; if [ -e "$f" ]; then cat "$f" >> "$O/hap$hap.hc_perchr.out"; fi; done
    # step 2: organelle blocks on the 12 chromosomes
    xargs samtools faidx "$GD/Qwislizeni_reasm_hap$hap.fa" < "$O/hap$hap.chrlist" > "$SLURM_TMPDIR/chr_hap$hap.fa"
    for q in NC_062091 CP129458; do
        minimap2 -cx asm20 -N 1000 -p 0.1 -t 8 "$SLURM_TMPDIR/chr_hap$hap.fa" "$ORG/$q.fa" > "$O/hap${hap}_$q.paf"
        echo "hap$hap $q PAF lines: $(wc -l < "$O/hap${hap}_$q.paf")"
    done
    cp "$BARR" "$O/hap$hap.barrnap.gff3"
done
$PY "$WG/scripts/wg22b_ncrna_summary.py" --work "$O" --haps $HAPS --out "$R" \
    --rdna "$WG/results/rdna/rdna_kmer_copy_number_summary.tsv" --fai-dir "$GD"
echo "DONE wg22b"
