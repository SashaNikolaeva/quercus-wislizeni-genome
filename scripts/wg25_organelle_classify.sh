#!/bin/bash
# wg25: classify the 898 unplaced contigs (hap1 648 / 34.81 Mb, hap2 250 /
# 13.25 Mb) as organelle-derived or nuclear. Writes a CANDIDATE LIST ONLY; the assemblies are NOT modified
#.
# Step 1: minimap2 -cx asm20 --secondary=no of each haplotype's unplaced FASTA (wg09) against NC_062091.1 (Q. wislizeni
#   plastome) + CP129458.1 (Q. variabilis mitogenome) concatenated, so plastid and mito compete for each region.
#   -c is needed for base-level identity (wg09 PAF had none). Blocks >= 1 kb; identity = PAF col10/col11.
#   plastid_cov = merged fraction covered by NC_062091.1 blocks with identity >= 0.99 (same species);
#   mito_cov    = same for CP129458.1 at identity >= 0.95 (cross-species, section Cerris); org_cov = union of the two.
# Step 2: samtools coverage of hifi_junctions/dip.bam with NO -q (default flag filter: primary, non-duplicate,
# MAPQ >= 0). At MAPQ >= 20 organelle-like contigs read 0.3x, because reads on near-identical copies get
#   MAPQ 0. depth_ratio = contig meandepth / median meandepth of that haplotype's 12 chromosomes.
# Step 3 classes: candidate_remove: org_cov >= 0.90 AND depth_ratio >= 2;
# organelle_lowdepth_flag: org_cov >= 0.90 AND depth_ratio < 2; mixed_keep (possible
#   NUPT/NUMT-bearing): 0.50 <= org_cov < 0.90; nuclear: org_cov < 0.50.
#   Caveat: depth is split across redundant copies, which biases depth_ratio LOW; hence lowdepth_flag contigs are kept and listed.
# Success (else exit 1, outputs stay in the work dir, nothing copied to results/): every one of the 898 contigs
#   classified exactly once; lengths sum to 34.81 / 13.25 Mb (+-0.01); chromosome median MAPQ>=0 depth per haplotype in
#   15-30x (hap1 chr01 MAPQ>=20 meandepth is 20.85x in coverage_q20.tsv).
# Inputs: ~/scratch/california_oaks/wislizeni_reasm/final/wg09_organelle/{unplaced_hap1.fa,unplaced_hap2.fa,
#   NC_062091.fa,CP129458.fa}; ~/scratch/california_oaks/wislizeni_reasm/hifi_junctions/{dip.bam,coverage_q20.tsv};
#   results/organelle/unplaced_organelle.tsv (wg09 calls, carried along); assembly .fai (read-only, for lengths).
# Outputs: results/organelle/{dip_coverage_q0.tsv, organelle_contig_classes.tsv, organelle_classes_summary.tsv,
#   organelle_removal_candidates_hap1.txt, organelle_removal_candidates_hap2.txt}; PAFs in the work dir.
# Smoke (SMOKE=1): hap2 only, 20 contigs (first 10 wg09 'nuclear' + first 10 wg09 plastid/mito calls), depth for those
#   contigs only (samtools coverage -r per contig) plus chr:1-1000000 of the 12 hap2 chromosomes (code path for the
#   median only); writes to <work>/smoke/, never results/.
# Usage: SMOKE=1 sbatch --time=0:30:00 scripts/wg25_organelle_classify.sh   # smoke
#        sbatch --dependency=afterok:<smoke> scripts/wg25_organelle_classify.sh  # full
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=3:00:00
#SBATCH --job-name=wg25_organelle
#SBATCH --output=/home/sashanik/scratch/california_oaks/paperA_runs/logs/wg25_%j.log
set -euo pipefail
module load StdEnv/2023 minimap2/2.28 samtools/1.22.1
PY=/home/sashanik/projects/def-santang3/sashanik/software/venv_plotting/bin/python
WG=/home/sashanik/oak-genomics-wislizeni/california_oaks/wislizeni_genome
G=/home/sashanik/projects/def-santang3/sashanik/california_oaks/genomes
O=/home/sashanik/scratch/california_oaks/wislizeni_reasm/final/wg09_organelle
J=/home/sashanik/scratch/california_oaks/wislizeni_reasm/hifi_junctions
W=/home/sashanik/scratch/california_oaks/paperA_runs/wg25
SMOKE=${SMOKE:-0}
if [ "$SMOKE" = 1 ]; then OUT=$W/smoke; HAPS="2"; R=$OUT; else OUT=$W; HAPS="1 2"; R=$WG/results/organelle; fi
mkdir -p "$OUT" "$R"
T=${SLURM_TMPDIR:?run under SLURM}
cd "$OUT"
cat "$O/NC_062091.fa" "$O/CP129458.fa" > "$T/organelles.fa"
HAPARGS=()
for h in $HAPS; do
    if [ "$SMOKE" = 1 ]; then
        awk -F'\t' -v h=hap$h '$1==h && $8=="nuclear"' "$WG/results/organelle/unplaced_organelle.tsv" | awk 'NR<=10{print $2}' > ids_hap$h.txt
        awk -F'\t' -v h=hap$h '$1==h && $8!="nuclear"' "$WG/results/organelle/unplaced_organelle.tsv" | awk 'NR<=10{print $2}' >> ids_hap$h.txt
        samtools faidx "$G/Qwislizeni_reasm_hap$h.fa" -r ids_hap$h.txt > "$T/query_hap$h.fa"
        Q=$T/query_hap$h.fa
    else
        Q=$O/unplaced_hap$h.fa
        awk '$1!~/_chr[0-9][0-9]$/{print $1}' "$G/Qwislizeni_reasm_hap$h.fa.fai" > ids_hap$h.txt
        n_fa=$(grep -c '>' "$Q"); n_ids=$(wc -l < ids_hap$h.txt)
        [ "$n_fa" = "$n_ids" ] || { echo "hap$h: unplaced FASTA has $n_fa seqs, assembly has $n_ids unplaced" >&2; exit 1; }
    fi
    awk 'FNR==NR{k[$1]=1; next} ($1 in k){print $1"\t"$2}' ids_hap$h.txt "$G/Qwislizeni_reasm_hap$h.fa.fai" > lengths_hap$h.tsv
    minimap2 -cx asm20 --secondary=no -t 8 "$T/organelles.fa" "$Q" > hap${h}_organelles.paf 2> hap${h}_minimap2.log
    HAPARGS+=(--hap "hap$h:$OUT/lengths_hap$h.tsv:$OUT/hap${h}_organelles.paf")
done
# Step 2: MAPQ >= 0 depth (no -q)
if [ "$SMOKE" = 1 ]; then
    { printf "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
      for c in $(cat ids_hap2.txt); do samtools coverage -r "$c" "$J/dip.bam" | awk 'NR>1'; done
      for i in 01 02 03 04 05 06 07 08 09 10 11 12; do samtools coverage -r "Qwis_hap2_chr$i:1-1000000" "$J/dip.bam" | awk 'NR>1'; done
    } > "$OUT/dip_coverage_q0.tsv"
else
    samtools coverage "$J/dip.bam" > "$OUT/dip_coverage_q0.tsv"
fi
# Step 3: classes + success checks (exit 1 on failure, before anything is copied to results/)
"$PY" "$WG/scripts/wg25_organelle_classify.py" "${HAPARGS[@]}" --cov_q0 "$OUT/dip_coverage_q0.tsv" \
    --cov_q20 "$J/coverage_q20.tsv" --wg09 "$WG/results/organelle/unplaced_organelle.tsv" --outdir "$OUT" \
    $( [ "$SMOKE" = 1 ] && echo --smoke )
if [ "$SMOKE" != 1 ]; then
    cp dip_coverage_q0.tsv organelle_contig_classes.tsv organelle_classes_summary.tsv \
       organelle_removal_candidates_hap1.txt organelle_removal_candidates_hap2.txt "$R/"
fi
echo "DONE wg25 (SMOKE=$SMOKE) -> $R"
