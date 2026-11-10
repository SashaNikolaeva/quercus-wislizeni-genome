#!/bin/bash
# wg19: annotation prep for .
# Step 1  identity gate: for each of the 24 chromosomes, md5 of the UPPERCASED Helixer-input chromosome vs the md5 of the
#         final-assembly chromosome of the same number, forward and reverse-complement (wg19_annotation_prep.py identity).
# GATE: every row forward-equal, otherwise exit 1 and write nothing further.
# Step 2  hap2 rename: sed 's/Qwislizeni_h2_chr/Qwis_hap2_chr/g' on the wg01 hap2 GFF3 -> <work>/Qwis_hap2.helixer.gff3;
#         the set of GFF seqids must equal exactly the 12 assembly chromosome names (chr12 only in smoke).
# Step 3  proteins: gffread -y (StdEnv/2020 gffread/0.12.3, subshell) against the assembly FASTA (hap2) and
#         cg05 Qwislizeni.chr.fa (hap1); then strip one terminal '.'/'*', drop proteins with an internal stop (logged),
#         write <work>/prot_hap{1,2}.fa. Raw protein count must equal the GFF mRNA count; full mode also checks
#         36,682 (hap1, = wg10 log "36682 raw") and 36,582 (hap2) mRNAs.
# Inputs : hap1 /home/sashanik/scratch/california_oaks/cg05_helixer/Qwislizeni/{Qwislizeni.helixer.gff3,Qwislizeni.chr.fa}
#          hap2 /home/sashanik/scratch/california_oaks/wg01_helixer/Qwislizeni_h2/Qwislizeni_h2.helixer.gff3
#               /project/def-santang3/sashanik/oak_sv_barrier/genomes/named/Qwislizeni_h2.chr.fa
#          assemblies /home/sashanik/projects/def-santang3/sashanik/california_oaks/genomes/Qwislizeni_reasm_hap{1,2}.fa
# Outputs: results/annotation/wg19_seq_identity_check.tsv (hap chr len_helixer len_asm md5_fwd_equal md5_rc_equal)
#          ~/scratch/california_oaks/paperA_runs/wg19/{Qwis_hap2.helixer.gff3,prot_hap1.fa,prot_hap2.fa,wg19_protein_counts.tsv}
#          SMOKE=1: chr12 only, everything (incl. the identity TSV) in ~/scratch/california_oaks/paperA_runs/wg19/smoke/
# Usage  : sbatch --export=ALL,SMOKE=1 scripts/wg19_annotation_prep.sh   (smoke)
#          sbatch --dependency=afterok:<smoke> scripts/wg19_annotation_prep.sh
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0:30:00
#SBATCH --job-name=wg19_annprep
#SBATCH --output=/home/sashanik/scratch/california_oaks/paperA_runs/logs/wg19_%j.log
set -euo pipefail
SMOKE=${SMOKE:-0}
WG=/home/sashanik/oak-genomics-wislizeni/california_oaks/wislizeni_genome
PY=/home/sashanik/projects/def-santang3/sashanik/software/venv_plotting/bin/python
G=/home/sashanik/projects/def-santang3/sashanik/california_oaks/genomes
H1=/home/sashanik/scratch/california_oaks/cg05_helixer/Qwislizeni
H1GFF=$H1/Qwislizeni.helixer.gff3; H1FA=$H1/Qwislizeni.chr.fa
H2GFF=/home/sashanik/scratch/california_oaks/wg01_helixer/Qwislizeni_h2/Qwislizeni_h2.helixer.gff3
H2FA=/project/def-santang3/sashanik/oak_sv_barrier/genomes/named/Qwislizeni_h2.chr.fa
W=/home/sashanik/scratch/california_oaks/paperA_runs/wg19
if [ "$SMOKE" = 1 ]; then O=$W/smoke; R=$O; CHRS="12"; else O=$W; R=$WG/results/annotation; CHRS="all"; fi
mkdir -p "$O" "$R"; cd "$O"
echo "[$(date '+%F %T')] wg19 SMOKE=$SMOKE out=$O results=$R"

# ---- step 1: identity gate (python exits 1 if any row is not forward-equal) ----
"$PY" "$WG/scripts/wg19_annotation_prep.py" identity "$R/wg19_seq_identity_check.tsv" "$CHRS" \
    hap1 "$H1FA" "$G/Qwislizeni_reasm_hap1.fa" hap2 "$H2FA" "$G/Qwislizeni_reasm_hap2.fa"
echo "[$(date '+%F %T')] identity gate PASSED"; cat "$R/wg19_seq_identity_check.tsv"

# ---- step 2: hap2 rename (+ chr12 subset in smoke; hap1 GFF used as is in full mode) ----
if [ "$SMOKE" = 1 ]; then
    awk -F'\t' '/^#/ || $1=="Qwis_hap1_chr12"' "$H1GFF" > hap1.helixer.gff3; GFF1=$O/hap1.helixer.gff3
    sed 's/Qwislizeni_h2_chr/Qwis_hap2_chr/g' "$H2GFF" | awk -F'\t' '/^#/ || $1=="Qwis_hap2_chr12"' > Qwis_hap2.helixer.gff3
    EXPECT2="Qwis_hap2_chr12"
else
    GFF1=$H1GFF
    sed 's/Qwislizeni_h2_chr/Qwis_hap2_chr/g' "$H2GFF" > Qwis_hap2.helixer.gff3
    EXPECT2=$(awk '$1 ~ /^Qwis_hap2_chr[0-9][0-9]$/{print $1}' "$G/Qwislizeni_reasm_hap2.fa.fai" | sort | tr '\n' ' ')
fi
GOT2=$(awk -F'\t' '!/^#/ && NF>=9{print $1}' Qwis_hap2.helixer.gff3 | sort -u | tr '\n' ' ')
if [ "$(echo $GOT2)" != "$(echo $EXPECT2)" ]; then
    echo "ERROR: renamed hap2 seqids [$GOT2] != assembly chromosomes [$EXPECT2]"; exit 1; fi
echo "renamed hap2 seqids OK: $GOT2"

# ---- step 3: proteins ----
( module load StdEnv/2020 gffread/0.12.3 >/dev/null 2>&1 && \
  gffread -y prot_raw_hap1.fa -g "$H1FA" "$GFF1" && \
  gffread -y prot_raw_hap2.fa -g "$G/Qwislizeni_reasm_hap2.fa" Qwis_hap2.helixer.gff3 )
"$PY" "$WG/scripts/wg19_annotation_prep.py" proteins "$O/wg19_protein_counts.tsv" "$SMOKE" \
    hap1 "$GFF1" prot_raw_hap1.fa prot_hap1.fa hap2 Qwis_hap2.helixer.gff3 prot_raw_hap2.fa prot_hap2.fa
cat "$O/wg19_protein_counts.tsv"
rm -f prot_raw_hap1.fa prot_raw_hap2.fa
echo "[$(date '+%F %T')] DONE wg19"
