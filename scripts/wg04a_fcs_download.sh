#!/bin/bash
# wg04a: download NCBI FCS (foreign contamination screen) for assembly v2.
# LOGIN-NODE download (compute nodes have no internet), FCS-GX database 'all' (~499 GB, release
# 2023-01-24 per its manifest) + the fcs-gx / fcs-adaptor containers + fcs.py / run_fcsadaptor.sh runners.
# S3 ran at ~0.5 MB/s per stream, so the files come from NCBI FTP; the two big files (all.gxs 177 GB, all.gxi 321 GB) are fetched as 2 GiB
# byte ranges, $PAR at a time; each part's size is checked (retries), parts are appended in order and deleted, and every
# file is md5-checked against the manifest. Resumable: finished parts / files are skipped.
# Run: nohup bash scripts/wg04a_fcs_download.sh > <fcs dir>/download.log 2>&1 &
set -uo pipefail
D=/home/sashanik/projects/def-santang3/sashanik/software/fcs; mkdir -p "$D/gxdb/parts"; cd "$D"
PAR=${PAR:-4}; CHUNK=$((2*1024*1024*1024))
REL=https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/latest
for f in fcs-gx.sif fcs-adaptor.sif; do [ -s "$f" ] || curl -sS -o "$f" "$REL/$f" || echo "FAIL $f"; done
[ -s fcs.py ] || curl -sSL -o fcs.py https://github.com/ncbi/fcs/raw/main/dist/fcs.py
[ -s run_fcsadaptor.sh ] || { curl -sSL -o run_fcsadaptor.sh https://github.com/ncbi/fcs/raw/main/dist/run_fcsadaptor.sh; chmod +x run_fcsadaptor.sh; }
M=https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest   # NCBI FTP: 73 MB/s per stream vs 0.5 MB/s from S3 (tested 2026-10-21)
[ -s gxdb/all.manifest ] || curl -sS -o gxdb/all.manifest "$M/all.manifest"
python3 - "$D/gxdb/all.manifest" > gxdb/files.tsv <<'PY'
import json, sys
for f in json.load(open(sys.argv[1]))["fileDetails"]: print(f["fileName"], f["hashValue"], f["fileSize"], sep="\t")
PY
echo "[$(date '+%F %T')] start (PAR=$PAR, chunk 2 GiB)"
get_part() {   # name idx start end
    local n=$1 i=$2 s=$3 e=$4 p="$D/gxdb/parts/$1.$2" want=$(( $4 - $3 + 1 ))
    [ -f "$p.ok" ] && return 0
    for t in 1 2 3 4 5 6 7 8; do
        curl -sS --fail -r "$s-$e" -o "$p" "$M/$n" && [ "$(stat -c %s "$p")" -eq "$want" ] && { touch "$p.ok"; return 0; }
        echo "retry $t $n part $i"; sleep 30
    done
    echo "PART FAIL $n $i"; return 1
}
export -f get_part; export D M
while IFS=$'\t' read -r n h s; do
    if [ "$s" -lt 1000000000 ]; then
        [ "$(stat -c %s "gxdb/$n" 2>/dev/null || echo 0)" -eq "$s" ] || curl -sS -o "gxdb/$n" "$M/$n"
        continue
    fi
    [ -f "gxdb/$n.done" ] && continue
    np=$(( (s + CHUNK - 1) / CHUNK ))
    for ((i=0; i<np; i++)); do st=$((i*CHUNK)); en=$(( st + CHUNK - 1 )); [ $en -ge $s ] && en=$((s-1)); echo "$n $i $st $en"; done
done < gxdb/files.tsv > gxdb/parts/jobs.txt
echo "$(wc -l < gxdb/parts/jobs.txt) parts to fetch"
xargs -P "$PAR" -L 1 bash -c 'get_part "$@"' _ < gxdb/parts/jobs.txt
echo "[$(date '+%F %T')] parts fetched; assembling"
while IFS=$'\t' read -r n h s; do
    [ "$s" -lt 1000000000 ] && continue; [ -f "gxdb/$n.done" ] && continue
    np=$(( (s + CHUNK - 1) / CHUNK )); miss=0
    for ((i=0; i<np; i++)); do [ -f "gxdb/parts/$n.$i.ok" ] || miss=1; done
    [ $miss -eq 1 ] && { echo "ERROR: $n has missing parts; rerun this script"; continue; }
    : > "gxdb/$n"
    for ((i=0; i<np; i++)); do cat "gxdb/parts/$n.$i" >> "gxdb/$n" && rm -f "gxdb/parts/$n.$i" "gxdb/parts/$n.$i.ok"; done
    [ "$(stat -c %s "gxdb/$n")" -eq "$s" ] && touch "gxdb/$n.done" || echo "ERROR: $n size mismatch"
done < gxdb/files.tsv
echo "[$(date '+%F %T')] md5 check"
ok=1; while IFS=$'\t' read -r n h s; do got=$(md5sum "gxdb/$n" | cut -d' ' -f1); [ "$got" = "$h" ] && echo "OK $n" || { echo "MD5 FAIL $n"; ok=0; }; done < gxdb/files.tsv
[ $ok -eq 1 ] && touch gxdb/download.ok && echo "[$(date '+%F %T')] ALL OK"
