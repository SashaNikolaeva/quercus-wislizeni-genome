#!/usr/bin/env python
"""rp13d: where do the large unplaced hap1 sequences of the new wislizeni assembly belong?

For each unplaced sequence >= 1 Mb (rp13c: Qwis_hap1_unplaced_0001 4.76 Mb, _0002 1.66 Mb; both single-copy depth,
genome-like GC, < 12% aligned to either haplotype's chromosomes):
 1. Omni-C contacts (yahs/<hap>/mapped.pairs.gz, contig coordinates, MAPQ >= 30 unique deduplicated pairs) between the
    contig(s) of the unplaced sequence and each named chromosome; reported as raw counts and per Mb of chromosome.
    A strong excess on one chromosome = a placement candidate.
 2. Coverage by other genomes (PAFs from rp13d.sh, minimap2 asm20, blocks >= 1 kb): fraction of the sequence covered and
    the best target sequence. Present in other oaks = real oak sequence; absent everywhere = contamination candidate.
Usage: rp13d_big_unplaced.py <final_dir> <yahs_dir> <hap> <out.tsv> <paf> [<paf> ...]
"""
import sys
import gzip
from collections import defaultdict

fin, ydir, hap, out = sys.argv[1:5]
pafs = sys.argv[5:]
names = [l.rstrip("\n").split("\t") for l in open(f"{fin}/Qwis_{hap}.names.tsv")][1:]
yahs2new = {r[1]: r[0] for r in names}
length = {r[0]: int(r[2]) for r in names}
big = [r[0] for r in names if "_unplaced_" in r[0] and int(r[2]) >= 1_000_000]
ctg2new = {}
for l in open(f"{ydir}/{hap}/{hap}_yahs_scaffolds_final.agp"):
    f = l.rstrip("\n").split("\t")
    if len(f) > 5 and f[4] == "W" and f[0] in yahs2new:
        ctg2new[f[5]] = yahs2new[f[0]]
chroms = [r[0] for r in names if "_chr" in r[0]]
cnt = defaultdict(lambda: defaultdict(int))
with gzip.open(f"{ydir}/{hap}/mapped.pairs.gz", "rt") as fh:
    for l in fh:
        if l[0] == "#":
            continue
        f = l.split("\t", 6)
        a, b = ctg2new.get(f[1]), ctg2new.get(f[3])
        if a is None or b is None or a == b:
            continue
        if a in big:
            cnt[a][b] += 1
        if b in big:
            cnt[b][a] += 1
cov = defaultdict(dict)
for p in pafs:
    lab = p.rsplit("/", 1)[-1].replace(".paf", "")
    iv, tgt = defaultdict(list), defaultdict(lambda: defaultdict(int))
    for l in open(p):
        f = l.split("\t")
        if f[0] in big and int(f[10]) >= 1000:
            iv[f[0]].append((int(f[2]), int(f[3])))
            tgt[f[0]][f[5]] += int(f[3]) - int(f[2])
    for q in big:
        ivs = sorted(iv[q])
        tot, cs, ce = 0, -1, -1
        for s, e in ivs:
            if s > ce:
                tot += max(0, ce - cs)
                cs, ce = s, e
            else:
                ce = max(ce, e)
        tot += max(0, ce - cs)
        best = max(tgt[q], key=tgt[q].get) if tgt[q] else "-"
        cov[q][lab] = (tot / length[q], best)
with open(out, "w") as fo:
    for q in big:
        fo.write(f"## {q} ({length[q]/1e6:.2f} Mb)\n")
        tot = sum(cnt[q].values())
        fo.write("# Omni-C contacts to chromosomes (count, per Mb of chromosome, % of all its inter-sequence contacts)\n")
        for c in sorted(chroms, key=lambda c: -cnt[q][c] / length[c]):
            fo.write(f"{q}\tcontacts\t{c}\t{cnt[q][c]}\t{cnt[q][c] / (length[c] / 1e6):.1f}\t{100 * cnt[q][c] / tot if tot else 0:.1f}\n")
        fo.write(f"{q}\tcontacts\tunplaced_other\t{sum(v for k, v in cnt[q].items() if k not in chroms)}\t-\t-\n")
        fo.write("# coverage by other genomes (fraction covered by asm20 blocks >= 1 kb, best target)\n")
        for lab, (c, b) in cov[q].items():
            fo.write(f"{q}\tcoverage\t{lab}\t{c:.3f}\t{b}\n")
print(open(out).read())
