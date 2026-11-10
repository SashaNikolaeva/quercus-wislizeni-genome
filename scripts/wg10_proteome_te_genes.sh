#!/bin/bash
# wg10: (1) completeness of the Helixer proteome of wislizeni hap1 (cg05, 12 chromosomes):
# compleasm 0.2.10 protein mode, embryophyta_odb10 (same lineage as wg02's genome-mode numbers); proteins from gffread -y
# (gffread only exists under StdEnv/2020, loaded in a subshell). (2) gene models whose CDS is > 50% covered by EDTA TE
# annotation (cg04 Qwislizeni, uppercase input), see wg10_proteome_te_genes.py. Helixer ran on the unmasked genome, so
# TE-derived models inflate the gene count; if > ~10% are TE-derived the filtered count becomes the headline.
# Output: results/genes/proteome_compleasm_odb10.txt, results/genes/gene_te_overlap.tsv, results/genes/gene_te_summary.tsv
# Usage: sbatch scripts/wg10_proteome_te_genes.sh
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --job-name=wg10_proteome
#SBATCH --output=/home/sashanik/scratch/california_oaks/logs/wg10_proteome_%j.log
set -euo pipefail
WG=~/oak-genomics-wislizeni/california_oaks/wislizeni_genome; R=$WG/results/genes; mkdir -p "$R"
H=/home/sashanik/scratch/california_oaks/cg05_helixer/Qwislizeni
TE=/home/sashanik/scratch/california_oaks/cg04_edta/Qwislizeni/Qwislizeni.fa.mod.EDTA.TEanno.gff3
ENV=/home/sashanik/projects/def-santang3/sashanik/software/miniforge3/envs
LIB=/home/sashanik/projects/def-santang3/sashanik/software/compleasm_lineages
O=/home/sashanik/scratch/california_oaks/wislizeni_reasm/final/wg10_proteome; mkdir -p "$O"; cd "$O"
( module load StdEnv/2020 gffread/0.12.3 >/dev/null 2>&1 && gffread -y prot_raw.fa -g "$H/Qwislizeni.chr.fa" "$H/Qwislizeni.helixer.gff3" )
# drop terminal stop '.' and skip proteins with internal stops (memory: '.' chars hang/confuse homology tools)
awk '/^>/{if(n&&s!~/\./)print n"\n"s; n=$0; s=""; next}{s=s $0} END{if(n&&s!~/\./)print n"\n"s}' \
    <(sed -E '/^>/!s/\.$//' prot_raw.fa) > prot.fa
echo "proteins: $(grep -c '>' prot_raw.fa) raw, $(grep -c '>' prot.fa) without internal stops"
export PATH="$ENV/compleasm/bin:$PATH"
rm -rf compleasm_prot
"$ENV/compleasm/bin/compleasm" protein -p prot.fa -l embryophyta --odb odb10 -L "$LIB" -t 16 -o compleasm_prot > compleasm_prot.log 2>&1
cp compleasm_prot/summary.txt "$R/proteome_compleasm_odb10.txt"; cat "$R/proteome_compleasm_odb10.txt"
/home/sashanik/projects/def-santang3/sashanik/software/venv_plotting/bin/python "$WG/scripts/wg10_proteome_te_genes.py" \
    "$H/Qwislizeni.helixer.gff3" "$TE" "$R/gene_te_overlap.tsv" "$R/gene_te_summary.tsv"
echo "DONE wg10"
