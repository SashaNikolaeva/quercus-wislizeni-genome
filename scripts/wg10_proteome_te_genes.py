#!/usr/bin/env python
"""wg10: fraction of each Helixer gene model's CDS covered by EDTA TE annotation (wislizeni hap1, 12 chromosomes).

Per mRNA: CDS bp (union of its CDS features) and CDS bp overlapping the union of EDTA TEanno features, counted twice:
all EDTA features, and without 'repeat_fragment' (unclassified fragments). Per gene: the mRNA with the most CDS bp.
A gene is 'TE-derived' when > 50% of its CDS bp is TE-covered.
Names: Helixer uses Qwis_hap1_chrNN, EDTA chrNN (same sequence, cg04 input = reasm hap1 chromosomes).
Usage: wg10_proteome_te_genes.py <helixer.gff3> <EDTA.TEanno.gff3> <out_genes.tsv> <out_summary.tsv>
"""
import sys, re
from collections import defaultdict
import numpy as np

GFF, TE, OUT, SUM = sys.argv[1:5]

def merge(iv):
    iv.sort(); out = []
    for s, e in iv:
        if out and s <= out[-1][1]: out[-1][1] = max(out[-1][1], e)
        else: out.append([s, e])
    return np.array(out, dtype=np.int64).reshape(-1, 2)

te_all, te_cls = defaultdict(list), defaultdict(list)
for l in open(TE):
    if l.startswith("#"): continue
    t = l.split("\t")
    if len(t) < 9: continue
    s, e = int(t[3]) - 1, int(t[4])
    te_all[t[0]].append((s, e))
    if t[2] != "repeat_fragment": te_cls[t[0]].append((s, e))
te_all = {c: merge(v) for c, v in te_all.items()}; te_cls = {c: merge(v) for c, v in te_cls.items()}

def covered(c, iv, te):
    """bp of merged intervals iv covered by merged TE intervals on chrom c"""
    if c not in te: return 0
    T = te[c]; n = 0
    for s, e in iv:
        i = np.searchsorted(T[:, 1], s, side="right")
        while i < len(T) and T[i, 0] < e:
            n += min(e, T[i, 1]) - max(s, T[i, 0]); i += 1
    return n

cds, parent, chrom = defaultdict(list), {}, {}
for l in open(GFF):
    if l.startswith("#"): continue
    t = l.rstrip("\n").split("\t")
    if len(t) < 9: continue
    a = dict(x.split("=", 1) for x in t[8].split(";") if "=" in x)
    if t[2] == "mRNA": parent[a["ID"]] = a["Parent"]
    elif t[2] == "CDS":
        m = a["Parent"]; cds[m].append((int(t[3]) - 1, int(t[4]))); chrom[m] = re.sub(r"^.*_(chr\d\d)$", r"\1", t[0])

best = {}
for m, iv in cds.items():
    iv = merge(iv); L = int((iv[:, 1] - iv[:, 0]).sum())
    g = parent.get(m, m)
    if g not in best or L > best[g][1]:
        best[g] = (m, L, covered(chrom[m], iv, te_all), covered(chrom[m], iv, te_cls), chrom[m])

with open(OUT, "w") as f:
    f.write("gene\tmrna\tchrom\tcds_bp\tte_bp_all\tte_bp_classified\tfrac_all\tfrac_classified\n")
    for g, (m, L, a, b, c) in best.items():
        f.write(f"{g}\t{m}\t{c}\t{L}\t{a}\t{b}\t{a/L:.3f}\t{b/L:.3f}\n")
fa = np.array([v[2] / v[1] for v in best.values()]); fc = np.array([v[3] / v[1] for v in best.values()])
n = len(fa)
with open(SUM, "w") as f:
    f.write("measure\tvalue\n")
    f.write(f"genes\t{n}\n")
    for k, x in (("all_EDTA", fa), ("classified_TE_only", fc)):
        for thr in (0.5, 0.8):
            f.write(f"genes_cds_gt{int(thr*100)}pct_TE_{k}\t{int((x > thr).sum())} ({100*(x > thr).mean():.1f}%)\n")
        f.write(f"genes_after_removing_gt50pct_TE_{k}\t{int((x <= 0.5).sum())}\n")
print(open(SUM).read())
