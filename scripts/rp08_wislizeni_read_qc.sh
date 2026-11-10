#!/bin/bash
# rp08: QC the raw reads before hifiasm.
#  - seqkit stats -a on every file -> qc/read_stats.tsv
#  - HiFi gate: mean read quality >= 20 and read N50 >= 10 kb (the ENA file is named "_subreads"; headers were stripped,
#    qualities spot-checked as CCS-like ~Q93). Fail = stop: subreads need ccs first.
#  - Omni-C gate: R1/R2 read counts must match for each run. SRR19243398's bare file is characterised (read count,
#    /1 vs /2 suffixes, overlap of its read IDs with _1) and NOT used: only properly paired _1/_2 go to hifiasm.
#  - meryl k=21 on HiFi (also Merqury's db for step 3) -> histogram -> GenomeScope2 (-p 2): size, heterozygosity.
# Writes qc/omnic_r1.txt, qc/omnic_r2.txt (comma lists for hifiasm --h1/--h2) and qc/rp08.ok on success.
# Usage: sbatch scripts/rp08_wislizeni_read_qc.sh
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --job-name=rp08_wis_qc
#SBATCH --output=/home/sashanik/scratch/california_oaks/logs/rp08_wis_qc_%j.log
set -euo pipefail
W=/home/sashanik/scratch/california_oaks/wislizeni_reasm; RAW=$W/raw; QC=$W/qc; mkdir -p "$QC"
source /home/sashanik/projects/def-santang3/sashanik/software/miniforge3/etc/profile.d/conda.sh; conda activate asm
ts() { date '+%F %T'; }
HIFI=$RAW/SRR19243399_subreads.fastq.gz

echo "[$(ts)] seqkit stats"
seqkit stats -a -T -j 16 "$RAW"/*.fastq.gz > "$QC/read_stats.tsv"
column -t "$QC/read_stats.tsv"

# HiFi gate (seqkit columns: file format type num_seqs sum_len min_len avg_len max_len Q1 Q2 Q3 sum_gap N50 ... AvgQual)
read -r N50 AQ < <(awk -F'\t' -v f="$HIFI" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} $1==f{print $h["N50"], $h["AvgQual"]}' "$QC/read_stats.tsv")
echo "HiFi: read N50 $N50, mean quality $AQ"
awk -v n="$N50" -v q="$AQ" 'BEGIN{exit !(n >= 10000 && q >= 20)}' || { echo "FAIL HiFi gate (not CCS-like); stop"; exit 1; }

# Omni-C pairing gate
R1L=""; R2L=""
for run in SRR19243397 SRR19243398; do
    n1=$(awk -F'\t' -v f="$RAW/${run}_1.fastq.gz" '$1==f{print $4}' "$QC/read_stats.tsv")
    n2=$(awk -F'\t' -v f="$RAW/${run}_2.fastq.gz" '$1==f{print $4}' "$QC/read_stats.tsv")
    echo "$run: R1 $n1 reads, R2 $n2 reads"
    [ "$n1" = "$n2" ] || { echo "FAIL $run R1/R2 counts differ"; exit 1; }
    R1L+="${R1L:+,}$RAW/${run}_1.fastq.gz"; R2L+="${R2L:+,}$RAW/${run}_2.fastq.gz"
done
B=$RAW/SRR19243398.fastq.gz
echo "bare file $B read-name suffixes:"; zcat "$B" | awk 'NR%4==1{print $NF ~ /\/1$/ ? "/1" : ($NF ~ /\/2$/ ? "/2" : "other")}' | sort | uniq -c
# head closes the pipe early -> zcat/awk get SIGPIPE (141); pipefail off for this line only (job 5868141 died here)
set +o pipefail; zcat "$B" | awk 'NR%4==1{sub(/^@/,"",$1); print $1}' | head -200000 | sort > "$QC/bare_ids.txt"; set -o pipefail
zcat "$RAW/SRR19243398_1.fastq.gz" | awk 'NR%4==1{sub(/^@/,"",$1); print $1}' | sort > "$QC/r1_ids.txt"
echo "bare-file IDs (first 200k) also present in _1: $(comm -12 "$QC/bare_ids.txt" "$QC/r1_ids.txt" | wc -l)"
rm -f "$QC/bare_ids.txt" "$QC/r1_ids.txt"
echo "$R1L" > "$QC/omnic_r1.txt"; echo "$R2L" > "$QC/omnic_r2.txt"

echo "[$(ts)] meryl k=21"
[ -d "$W/hifi_k21.meryl" ] || meryl count k=21 threads=16 memory=56 "$HIFI" output "$W/hifi_k21.meryl"
meryl histogram "$W/hifi_k21.meryl" > "$QC/hifi_k21.histo"
echo "[$(ts)] GenomeScope2"
genomescope2 -i "$QC/hifi_k21.histo" -o "$QC/genomescope2" -k 21 -p 2 -n wislizeni
cat "$QC/genomescope2/wislizeni_summary.txt"
touch "$QC/rp08.ok"
echo "[$(ts)] DONE rp08"
