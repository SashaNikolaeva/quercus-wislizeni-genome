#!/bin/bash
# wg27: BlobToolKit snail + blob plots of the final
# assemblies, with and without the wg25 organelle removal candidates. No new mapping, no diamond.
# Env: software/miniforge3/envs/blobtk = blobtoolkit 4.5.1 (pip; `blobtools` CLI) + blobtk 0.8.0 (bioconda), python 3.11.
# Fields added to the BlobDir (blobtools create computes gc, length, ncount):
#   <lineage>_busco : wg24 genome-mode compleasm odb12 full_table_busco_format.tsv on the final scaffolds, header line 1
#                     rewritten so the BUSCO parser accepts it (see wg27_blobtoolkit.py; Missing rows dropped)
#   depth_q0        : samtools coverage meandepth of dip.bam, MAPQ >= 0 (results/organelle/dip_coverage_q0.tsv, wg25)
#   gx_taxon/gx_taxid/gx_div : FCS-GX taxonomy.rpt rows (split at gaps, seqid~start..end) aggregated to the seqid,
#                     top tax-id by summed score-1, else "no-hit"
#   wg25_class      : chromosome / nuclear / mixed_keep / organelle_lowdepth_flag / candidate_remove
# Plots (blobtk plot): snail (BUSCO + length/GC), blob (x = gc, y = depth_q0, z = length, colour = wg25_class);
#   each also as _minus_candidates = filter wg25_class:inv=candidate_remove (the contigs listed in
#   results/organelle/organelle_removal_candidates_hap<h>.txt).
# Fallback: if blobtk cannot render the blob plot, figures/blobtools/gc_depth_hap<h>.png (matplotlib, venv_plotting):
#   GC fraction vs log10 MAPQ>=0 depth, coloured by wg25 class, labelled "blob-style".
# Success: a snail plot per haplotype whose BUSCO numbers equal wg24 genome-mode odb12 to 0.1 pp; the log
#   prints the S/D/F % recomputed from the table, the snail JSON, and the wg24 summary for this comparison.
# Inputs: G = ~/projects/def-santang3/sashanik/california_oaks/genomes/Qwislizeni_reasm_hap{1,2}.fa;
#   wg24 table (~/scratch/california_oaks/paperA_runs/wg24/**/embryophyta_odb12/full_table_busco_format.tsv with
#   'genome' and 'hap<h>' in the path; override with BUSCO_H1 / BUSCO_H2); results/organelle/{dip_coverage_q0.tsv,
#   organelle_contig_classes.tsv}; ~/scratch/california_oaks/wislizeni_reasm/wg04_fcs/gx_hap<h>/*.taxonomy.rpt.
# Outputs: figures/blobtools/{snail,blob}_hap<h>{,_minus_candidates}.{svg,png} (or gc_depth_hap<h>.png fallback);
#   BlobDirs + snail JSON in ~/scratch/california_oaks/paperA_runs/wg27/.
# Smoke (SMOKE=1): hap2, 12 chromosomes + 30 unplaced contigs (the 20 wg25-smoke contigs + the first 10 other unplaced;
#   the 10 have no depth/class, filled with 0 / NA), wg25-smoke depth/classes, and - because wg24 has not run - the
#   wg02 odb10 genome-mode table of the same final hap2 scaffolds (code path only); writes to <work>/smoke/ only.
# Usage: SMOKE=1 sbatch --time=1:00:00 --dependency=afterok:<wg25 smoke> scripts/wg27_blobtoolkit.sh       # smoke
#        sbatch --dependency=afterok:<wg25 full>:<wg24 tasks 1-2>:<wg27 smoke> scripts/wg27_blobtoolkit.sh    # full
#SBATCH --account=def-santang3_cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH --job-name=wg27_blobtk
#SBATCH --output=/home/sashanik/scratch/california_oaks/paperA_runs/logs/wg27_%j.log
set -euo pipefail
module load StdEnv/2023 samtools/1.22.1
SW=/home/sashanik/projects/def-santang3/sashanik/software
ENV=$SW/miniforge3/envs/blobtk
PY=$SW/venv_plotting/bin/python
WG=/home/sashanik/oak-genomics-wislizeni/california_oaks/wislizeni_genome
G=/home/sashanik/projects/def-santang3/sashanik/california_oaks/genomes
GX=/home/sashanik/scratch/california_oaks/wislizeni_reasm/wg04_fcs
PR=/home/sashanik/scratch/california_oaks/paperA_runs
W=$PR/wg27
T=${SLURM_TMPDIR:?run under SLURM}
SMOKE=${SMOKE:-0}
if [ "$SMOKE" = 1 ]; then
    HAPS="2"; OUT=$W/smoke; FIG=$OUT; ORG=$PR/wg25/smoke
else
    HAPS="1 2"; OUT=$W; FIG=$WG/figures/blobtools; ORG=$WG/results/organelle
fi
mkdir -p "$OUT" "$FIG"
cd "$OUT"
export PATH="$ENV/bin:$PATH"
echo "blobtoolkit $("$ENV/bin/python" -c 'import importlib.metadata as m; print(m.version("blobtoolkit"))'), $(blobtk --version)"
for h in $HAPS; do
    H=hap$h; D=$OUT/blobdir_$H; rm -rf "$D" "$OUT/prep_$H"; mkdir -p "$OUT/prep_$H"
    if [ "$SMOKE" = 1 ]; then
        { awk '$1~/_chr[0-9][0-9]$/{print $1}' "$G/Qwislizeni_reasm_hap$h.fa.fai"
          cat "$ORG/ids_hap$h.txt"
          awk 'FNR==NR{k[$1]=1; next} $1!~/_chr[0-9][0-9]$/ && !($1 in k){print $1}' "$ORG/ids_hap$h.txt" "$G/Qwislizeni_reasm_hap$h.fa.fai" | awk 'NR<=10'
        } > "$OUT/prep_$H/ids.txt"
        samtools faidx "$G/Qwislizeni_reasm_hap$h.fa" -r "$OUT/prep_$H/ids.txt" > "$T/$H.fa"
        FA=$T/$H.fa; IDS=(--ids "$OUT/prep_$H/ids.txt" --smoke)
        B=/home/sashanik/scratch/california_oaks/wislizeni_reasm/wg02/v2_hap$h.compleasm_odb10/embryophyta_odb10/full_table_busco_format.tsv
    else
        FA=$G/Qwislizeni_reasm_hap$h.fa; IDS=()
        var=BUSCO_H$h; B=${!var:-}
        if [ -z "$B" ]; then
            mapfile -t BS < <(find "$PR/wg24" -path "*embryophyta_odb12*" -name full_table_busco_format.tsv | grep -i genome | grep -E "hap$h([^0-9]|$)" || true)
            [ "${#BS[@]}" = 1 ] || { echo "found ${#BS[@]} wg24 genome tables for $H: ${BS[*]:-none}; set BUSCO_H$h" >&2; exit 1; }
            B=${BS[0]}
        fi
    fi
    echo "$H: FASTA $FA; BUSCO $B"
    "$PY" "$WG/scripts/wg27_blobtoolkit.py" prep --fasta "$FA" --cov "$ORG/dip_coverage_q0.tsv" \
        --classes "$ORG/organelle_contig_classes.tsv" --gx "$(ls "$GX"/gx_$H/*.taxonomy.rpt)" --busco "$B" \
        --outdir "$OUT/prep_$H" "${IDS[@]}"
    blobtools create --fasta "$FA" "$D"
    blobtools add --busco "$OUT/prep_$H/busco.tsv" --text "$OUT/prep_$H/fields.tsv" --text-header "$D"
    # 2026-11-05: blobtk 0.8.0 filters on wg25_class return only the candidate contigs whatever the syntax
    # (":inv=candidate_remove" in 8577172 and an explicit ":keys=" keep-list in 8594923 both gave 81/35 sequences).
    # No blobtk filter any more: the _minus_candidates plots come from their own BlobDir built from the wg28 cleaned
    # FASTA (Qwislizeni_reasm_hap<h>.clean.fa = the final assembly), checked against the count expected.
    VARS=(""); [ "$SMOKE" = 1 ] || VARS+=(_minus_candidates)
    for v in "${VARS[@]}"; do
        DD=$D
        if [ -n "$v" ]; then
            CF=$G/Qwislizeni_reasm_hap$h.clean.fa; DD=$OUT/blobdir_${H}_clean; P=$OUT/prep_${H}_clean
            rm -rf "$DD" "$P"; mkdir -p "$P"; cut -f1 "$CF.fai" > "$P/ids.txt"
            "$PY" "$WG/scripts/wg27_blobtoolkit.py" prep --fasta "$CF" --cov "$ORG/dip_coverage_q0.tsv" \
                --classes "$ORG/organelle_contig_classes.tsv" --gx "$(ls "$GX"/gx_$H/*.taxonomy.rpt)" --busco "$B" \
                --outdir "$P" --ids "$P/ids.txt"
            blobtools create --fasta "$CF" "$DD"
            blobtools add --busco "$P/busco.tsv" --text "$P/fields.tsv" --text-header "$DD"
        fi
        blobtk plot -d "$DD" -v snail -o "$FIG/snail_$H$v.svg" -o "$FIG/snail_$H$v.png" -o "$OUT/snail_$H$v.json"
        if [ -n "$v" ]; then
            nall=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["scaffold_count"])' "$OUT/snail_$H.json")
            nmin=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["scaffold_count"])' "$OUT/snail_$H$v.json")
            ncand=$(awk -F'\t' -v h="$H" '$1==h && $11=="candidate_remove"' "$ORG/organelle_contig_classes.tsv" | awk 'END{print NR}')
            echo "$H: snail sequences all=$nall minus_candidates=$nmin candidates=$ncand"
            [ "$nmin" -eq $((nall - ncand)) ] || { echo "ERROR $H: minus_candidates count $nmin != $nall - $ncand"; exit 1; }
        fi
        if ! blobtk plot -d "$DD" -v blob -x gc -y depth_q0 -z length -c wg25_class \
                -o "$FIG/blob_$H$v.svg" -o "$FIG/blob_$H$v.png"; then
            echo "WARNING: blobtk blob plot failed for $H$v; using the matplotlib fallback"
            BLOBFAIL=1
        fi
    done
    echo "snail JSON busco entries ($H):"; grep -o '"busco[^]]*]' "$OUT/snail_$H.json" || grep -o -i 'busco.\{0,200\}' "$OUT/snail_$H.json" | awk 'NR<=5' || true
    S=$WG/results/annotation/compleasm_odb12_genome_$H.txt; [ -f "$S" ] && { echo "wg24 summary $S:"; cat "$S"; }
    if [ "${BLOBFAIL:-0}" = 1 ]; then
        "$PY" "$WG/scripts/wg27_blobtoolkit.py" fallback --gc "$OUT/prep_$H/gc.tsv" --hap "$H" --out "$FIG/gc_depth_$H.png"
        BLOBFAIL=0
    fi
done
ls -l "$FIG"
echo "DONE wg27 (SMOKE=$SMOKE)"
