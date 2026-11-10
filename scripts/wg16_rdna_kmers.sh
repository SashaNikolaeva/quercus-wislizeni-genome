#!/bin/bash
# wg16: rDNA copy number of the reference tree from HiFi k-mers (hifi_k21.meryl),
# independent of assembly collapse. Assembled copies are only a floor.
# Query sequences = wislizeni's own rDNA copies from hap1 (c08 blastn hits of A. thaliana 18S X16077 / 25S X52320; 5S
# from Rfam RF00001 cmsearch): 3 full-length copies each of 18S (X16077 78-1901), 25S (X52320 848-4230), 5S (full model),
# plus 5 single-copy control regions (2 kb, gene-rich, collinear hap1-hap2 windows) to check the scale.
# meryl-lookup -wig-count -> per-k-mer count (k = 21, canonical). Copy number (diploid, both haplotypes together) =
# median k-mer count / kmercov, kmercov = 22.79 = GenomeScope2 per-haploid-copy k-mer coverage (wislizeni_model.txt);
# a single-copy homozygous k-mer is expected at 2 x 22.79 = 45.6. Haploid estimate = diploid / 2.
# Output: results/rdna/{rdna_kmer_copy_number.tsv,README.md}; work in ~/scratch/california_oaks/paperA_runs/wg16
# Usage: sbatch scripts/wg16_rdna_kmers.sh
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=2:00:00
#SBATCH --job-name=wg16_rdna
#SBATCH --output=/home/sashanik/scratch/california_oaks/paperA_runs/logs/wg16_%j.log
set -euo pipefail
WG=/home/sashanik/oak-genomics-wislizeni/california_oaks/wislizeni_genome
PY=/home/sashanik/projects/def-santang3/sashanik/software/venv_plotting/bin/python
W=/home/sashanik/scratch/california_oaks/paperA_runs/wg16; mkdir -p "$W"; cd "$W"
module load StdEnv/2023 samtools/1.22.1
ML=/home/sashanik/projects/def-santang3/sashanik/software/miniforge3/envs/asm/bin/meryl-lookup
$PY "$WG/scripts/wg16_rdna_kmers.py" select "$W"
for f in q_*.fa; do
    $ML -wig-count -sequence "$f" -output "${f%.fa}.wig" -mers /home/sashanik/scratch/california_oaks/wislizeni_reasm/hifi_k21.meryl
done
$PY "$WG/scripts/wg16_rdna_kmers.py" count "$W" "$WG/results/rdna"
echo "DONE wg16"
