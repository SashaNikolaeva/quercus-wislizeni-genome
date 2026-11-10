#!/bin/bash
# wg09: organelle screen of the unplaced contigs of both haplotypes (hap1 648 contigs 34.8 Mb,
# hap2 250 contigs 13.2 Mb). minimap2 -x asm20 against the Q. wislizeni plastome NC_062091.1 (161 kb) and the closest
# complete oak mitogenome on NCBI, Q. variabilis CP129458.1 (420 kb; no red-oak mitogenome exists), both downloaded
# 2026-10-23 into <O>. Per contig: bp covered by plastid / mito alignments (>= 1 kb blocks, merged), fraction of contig.
# A contig is called organellar when >= 50% of it is covered (plastid or mito, whichever is larger).
# Output: results/organelle/unplaced_organelle.tsv (per contig), results/organelle/unplaced_organelle_summary.tsv
# Usage: sbatch scripts/wg09_organelle_screen.sh
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --job-name=wg09_organelle
#SBATCH --output=/home/sashanik/scratch/california_oaks/logs/wg09_organelle_%j.log
set -euo pipefail
module load StdEnv/2023 minimap2/2.28 samtools/1.22.1
R=~/oak-genomics-wislizeni/california_oaks/wislizeni_genome/results/organelle; mkdir -p "$R"
G=/home/sashanik/projects/def-santang3/sashanik/california_oaks/genomes
O=/home/sashanik/scratch/california_oaks/wislizeni_reasm/final/wg09_organelle
cd "$O"
printf "hap\tcontig\tlength\tplastid_bp\tmito_bp\tplastid_frac\tmito_frac\tcall\n" > "$R/unplaced_organelle.tsv"
for h in 1 2; do
    F=$G/Qwislizeni_reasm_hap$h.fa
    awk '$1!~/_chr[0-9][0-9]$/{print $1}' "$F.fai" > unplaced_hap$h.ids
    samtools faidx "$F" -r unplaced_hap$h.ids > unplaced_hap$h.fa
    for o in NC_062091 CP129458; do
        minimap2 -x asm20 -t 8 --secondary=no "$o.fa" unplaced_hap$h.fa > hap${h}_$o.paf 2>/dev/null
    done
    # union of query intervals (>= 1 kb aligned blocks) per contig and organelle
    for o in NC_062091 CP129458; do
        awk '$11>=1000{print $1"\t"$3"\t"$4}' hap${h}_$o.paf | sort -k1,1 -k2,2n \
          | awk -v o=$o '{if($1!=c || $2>e){if(c!="")s[c]+=e-b; c=$1; b=$2; e=$3} else if($3>e) e=$3} END{if(c!="")s[c]+=e-b; for(k in s) print k"\t"o"\t"s[k]}'
    done > hap${h}_cov.tsv
    awk -v h=hap$h 'FNR==NR{if($2=="NC_062091")p[$1]=$3; else m[$1]=$3; next} $1!~/_chr[0-9][0-9]$/{
        pf=(p[$1]+0)/$2; mf=(m[$1]+0)/$2; call=(pf>=0.5||mf>=0.5)?(pf>=mf?"plastid":"mito"):"nuclear";
        printf "%s\t%s\t%d\t%d\t%d\t%.3f\t%.3f\t%s\n", h, $1, $2, p[$1]+0, m[$1]+0, pf, mf, call}' hap${h}_cov.tsv "$F.fai" >> "$R/unplaced_organelle.tsv"
done
awk -F'\t' 'NR>1{k=$1"\t"$8; n[k]++; L[k]+=$3; P[$1]+=$4; M[$1]+=$5} END{
    print "hap\tcall\tcontigs\tMb"; for(k in n) printf "%s\t%d\t%.3f\n", k, n[k], L[k]/1e6
    for(h in P) printf "%s\tall_plastid_aligned_bp\t\t%.3f\n%s\tall_mito_aligned_bp\t\t%.3f\n", h, P[h]/1e6, h, M[h]/1e6}' \
    "$R/unplaced_organelle.tsv" | sort > "$R/unplaced_organelle_summary.tsv"
cat "$R/unplaced_organelle_summary.tsv"
echo "DONE wg09"
