#!/bin/bash
# rp13e: taxonomy of the large unplaced hap1 contigs.
# 20 evenly spaced 2 kb chunks of Qwis_hap1_unplaced_0001 (4.76 Mb) and 10 of _0002 (1.66 Mb) (rp13e_blast/chunks_*.fa,
# made with samtools faidx), megablast vs NCBI nt REMOTE (blast+ 2.17.0, -remote). This must run on the LOGIN node:
# compute nodes have no internet. One query at a time, to stay within NCBI's usage limits.
# Output: rp13e_blast/<contig>.tsv (outfmt 6 + staxids sscinames scomnames sskingdoms), top 5 hits per chunk.
# Usage: bash scripts/rp13e_blast_unplaced.sh   (login node, in the background)
set -euo pipefail
D=/home/sashanik/scratch/california_oaks/wislizeni_reasm/final/rp13e_blast
module load StdEnv/2023 blast+/2.17.0
cd "$D"
for c in 0001 0002; do
    : > "unplaced_$c.tsv"
    awk '/^>/{if(s)print n"\t"s; n=substr($0,2); s=""; next}{s=s $0} END{print n"\t"s}' "chunks_$c.fa" | while IFS=$'\t' read -r n s; do
        printf ">%s\n%s\n" "$n" "$s" > q.fa
        blastn -task megablast -remote -db nt -query q.fa -max_target_seqs 5 -evalue 1e-20 \
            -outfmt "6 qseqid sseqid pident length evalue bitscore staxids sscinames scomnames sskingdoms" >> "unplaced_$c.tsv" 2>> blast.err \
            || echo "$n FAILED" >> blast.err
        sleep 10
    done
done
rm -f q.fa
echo "DONE rp13e"
