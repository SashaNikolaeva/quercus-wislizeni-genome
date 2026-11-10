#!/bin/bash
# rp07: download the raw reads behind the Q. wislizeni CCGP reference (GCA_023055345.1, isolate QUWIS_CCGP_1,
# BioSample SAMN26368104) for reassembly: PacBio HiFi SRR19243399 + Dovetail Omni-C SRR19243397/SRR19243398 (ENA).
# Parallel byte ranges: ENA gave ~1 MB/s per connection, so every file is split into CHUNK-byte
# ranges fetched by P concurrent curl -r; each chunk resumes by appending; chunks are concatenated, md5-checked
# against the ENA manifest, moved to raw/ and appended to verified.txt.
# Login node (downloads are the only login-node exception). Resumable: rerun the same command.
# Usage: nohup bash scripts/rp07_download_wislizeni_raw.sh >> $OUT/download.log 2>&1 &
set -uo pipefail
MANIFEST=${MANIFEST:-$HOME/oak-genomics/california_oaks/wislizeni_genome/data/wislizeni_ccgp_raw_reads.tsv}
OUT=${OUT:-$HOME/scratch/california_oaks/wislizeni_reasm}
P=${P:-12}; CHUNK=${CHUNK:-1073741824}
mkdir -p "$OUT/raw" "$OUT/part"; touch "$OUT/verified.txt"
N=$(($(wc -l < "$MANIFEST") - 1))
echo "$(date '+%F %T') START $(grep -c . "$OUT/verified.txt")/$N verified; P=$P chunk=$CHUNK"

chunk() {  # file idx start end url
    local f=$1 i=$2 s=$3 e=$4 url=$5 c have want
    c=$OUT/part/$f.c$(printf %03d "$i"); want=$((e - s + 1))
    for try in 1 2 3 4 5 6 7 8; do
        have=$(stat -c %s "$c" 2>/dev/null || echo 0)
        [ "$have" -eq "$want" ] && return 0
        [ "$have" -gt "$want" ] && { rm -f "$c"; have=0; }
        curl -sS -L --retry 3 --speed-limit 20000 --speed-time 120 -r "$((s + have))-$e" "$url" >> "$c"
    done
    echo "$(date '+%F %T') CHUNK FAIL $f c$i"
}
export -f chunk; export OUT

tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r run md5 bytes url; do
    f=$(basename "$url"); grep -qx "$f" "$OUT/verified.txt" && continue
    for ((i = 0, s = 0; s < bytes; i++, s += CHUNK)); do
        e=$((s + CHUNK - 1)); [ $e -ge "$bytes" ] && e=$((bytes - 1))
        echo "$f $i $s $e $url"
    done
done | xargs -P "$P" -n 5 bash -c 'chunk "$@"' _

tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r run md5 bytes url; do
    f=$(basename "$url"); grep -qx "$f" "$OUT/verified.txt" && continue
    cat "$OUT/part/$f".c* > "$OUT/part/$f"
    if [[ $(stat -c %s "$OUT/part/$f") -eq $bytes && $(md5sum "$OUT/part/$f" | cut -d' ' -f1) == "$md5" ]]; then
        mv "$OUT/part/$f" "$OUT/raw/$f"; rm -f "$OUT/part/$f".c*; echo "$f" >> "$OUT/verified.txt"
        echo "$(date '+%F %T') OK $f"
    else
        rm -f "$OUT/part/$f"; echo "$(date '+%F %T') md5/size mismatch $f (chunks kept; rerun)"
    fi
done
n=$(grep -c . "$OUT/verified.txt")
echo "$(date '+%F %T') END $n/$N verified"
[[ $n -eq $N ]] && echo WISLIZENI_RAW_DONE
