#!/usr/bin/env python
"""wg20: gene-model statistics of the Helixer annotations (12 chromosomes per haplotype).

Usage: wg20_annotation_stats.py <hap1.gff3> <hap2_renamed.gff3> <gene_te_overlap.tsv> <gene_te_summary.tsv> <out.tsv> <smoke 0|1>

Sets: hap1_all (cg05 GFF3), hap1_nonTE (hap1 genes with frac_all <= 0.5 in gene_te_overlap.tsv; TE-derived =
frac_all > 0.5: > 50% of the gene's CDS bp covered by any EDTA feature), hap2_all (wg19 renamed GFF3).
Definitions (all lengths in bp, GFF coordinates 1-based inclusive, length = end - start + 1):
  genes, mRNA            = number of 'gene' / 'mRNA' features in the set (mRNA of a gene belong to the gene's set)
  exons_per_mRNA         = number of exon features / number of mRNA
  CDS_per_mRNA           = number of CDS features / number of mRNA
  gene_length            = gene feature span (end - start + 1); mean and median over genes
  transcript_length      = sum of the mRNA's exon lengths (spliced length); mean/median over mRNA
  CDS_length             = sum of the mRNA's CDS feature lengths; mean/median over mRNA
  exon_length            = each exon feature; mean/median over all exons
  intron_length          = gaps between consecutive exons of one mRNA after sorting by start:
                           start(exon i+1) - end(exon i) - 1; mean/median over all introns
                           (check: n_introns == n_exons - n_mRNA)
  pct_mono_exonic        = 100 * mRNA with exactly 1 exon / mRNA
  pct_with_5UTR / 3UTR   = 100 * mRNA with >= 1 five_prime_UTR / three_prime_UTR feature / mRNA
Checks (full mode, exit 1 on failure): genes hap1_all 36682, hap1_nonTE 32788 and equal to
'genes_after_removing_gt50pct_TE_all_EDTA' in gene_te_summary.tsv, hap2_all 36582; every hap1 gene present in
gene_te_overlap.tsv; intron count identity in every set (also in smoke); no intron length <= 0.
Output: TSV, rows = statistic, columns = sets (written as <out>.FAILED if a check fails).
"""
import sys
from collections import defaultdict
import numpy as np

G1, G2, TEO, TES, OUT, SMOKE = sys.argv[1:7]
smoke = SMOKE == "1"


def attr(a, key):
    for kv in a.strip().split(";"):
        if kv.startswith(key + "="):
            return kv[len(key) + 1:]
    return None


def read_gff(path):
    genes = {}                      # gene id -> length
    mrna_gene = {}                  # mRNA id -> gene id
    exons = defaultdict(list); cds = defaultdict(list); u5 = set(); u3 = set()
    with open(path) as f:
        for l in f:
            if l.startswith("#"):
                continue
            t = l.rstrip("\n").split("\t")
            if len(t) < 9:
                continue
            ty, s, e = t[2], int(t[3]), int(t[4])
            if ty == "gene":
                genes[attr(t[8], "ID")] = e - s + 1
            elif ty == "mRNA":
                mrna_gene[attr(t[8], "ID")] = attr(t[8], "Parent")
            elif ty == "exon":
                exons[attr(t[8], "Parent")].append((s, e))
            elif ty == "CDS":
                cds[attr(t[8], "Parent")].append((s, e))
            elif ty == "five_prime_UTR":
                u5.add(attr(t[8], "Parent"))
            elif ty == "three_prime_UTR":
                u3.add(attr(t[8], "Parent"))
    return genes, mrna_gene, exons, cds, u5, u3


def stats(name, gset, gff):
    genes, mrna_gene, exons, cds, u5, u3 = gff
    mr = [m for m, g in mrna_gene.items() if g in gset]
    glen = np.array([genes[g] for g in gset])
    tlen, clen, elen, ilen = [], [], [], []
    n_exon = n_cds = mono = w5 = w3 = 0
    for m in mr:
        ex = sorted(exons.get(m, []))
        n_exon += len(ex); n_cds += len(cds.get(m, []))
        elen += [e - s + 1 for s, e in ex]
        tlen.append(sum(e - s + 1 for s, e in ex))
        clen.append(sum(e - s + 1 for s, e in cds.get(m, [])))
        ilen += [ex[i + 1][0] - ex[i][1] - 1 for i in range(len(ex) - 1)]
        mono += len(ex) == 1; w5 += m in u5; w3 += m in u3
    nm = len(mr)
    ok = True
    if len(ilen) != n_exon - nm:
        print("ERROR %s: introns %d != exons %d - mRNA %d" % (name, len(ilen), n_exon, nm)); ok = False
    if ilen and min(ilen) <= 0:
        print("ERROR %s: %d introns of length <= 0" % (name, sum(1 for x in ilen if x <= 0))); ok = False
    f = lambda a: (float(np.mean(a)), float(np.median(a))) if len(a) else (float("nan"), float("nan"))
    r = {"genes": len(gset), "mRNA": nm, "exons": n_exon, "introns": len(ilen), "CDS_features": n_cds,
         "exons_per_mRNA": n_exon / nm, "CDS_per_mRNA": n_cds / nm}
    for lab, a in (("gene_length", glen), ("transcript_length", tlen), ("CDS_length", clen),
                   ("exon_length", elen), ("intron_length", ilen)):
        r[lab + "_mean"], r[lab + "_median"] = f(a)
    r["pct_mono_exonic"] = 100.0 * mono / nm
    r["pct_with_5UTR"] = 100.0 * w5 / nm
    r["pct_with_3UTR"] = 100.0 * w3 / nm
    print("%s: genes %d mRNA %d exons %d introns %d" % (name, len(gset), nm, n_exon, len(ilen)))
    return r, ok


te = {}
with open(TEO) as f:
    hdr = f.readline().rstrip("\n").split("\t")
    # exact fraction from bp columns: the rounded frac_all column (3 dp) turns 0.5004 into "0.500" (bug 2026-11-05,
    # job 8445437 counted 32789 instead of wg10's exact 32788)
    gi, ci, ai = hdr.index("gene"), hdr.index("cds_bp"), hdr.index("te_bp_all")
    for l in f:
        t = l.rstrip("\n").split("\t")
        te[t[gi]] = int(t[ai]) / int(t[ci]) > 0.5
summ = {}
with open(TES) as f:
    for l in f:
        t = l.rstrip("\n").split("\t")
        if len(t) >= 2:
            summ[t[0]] = t[1]

g1 = read_gff(G1); g2 = read_gff(G2)
ok = True
missing = [g for g in g1[0] if g not in te]
if missing:
    print("ERROR: %d hap1 genes not in gene_te_overlap.tsv (e.g. %s)" % (len(missing), missing[:3])); ok = False
sets = {"hap1_all": (set(g1[0]), g1),
        "hap1_nonTE": ({g for g in g1[0] if g in te and not te[g]}, g1),
        "hap2_all": (set(g2[0]), g2)}
res = {}
for n, (gs, gff) in sets.items():
    res[n], o = stats(n, gs, gff); ok &= o

if not smoke:
    exp = {"hap1_all": 36682, "hap1_nonTE": 32788, "hap2_all": 36582}
    for n, v in exp.items():
        if res[n]["genes"] != v:
            print("ERROR: %s genes %d != expected %d" % (n, res[n]["genes"], v)); ok = False
    ref = int(summ.get("genes_after_removing_gt50pct_TE_all_EDTA", "-1"))
    if res["hap1_nonTE"]["genes"] != ref:
        print("ERROR: hap1_nonTE %d != gene_te_summary.tsv %d" % (res["hap1_nonTE"]["genes"], ref)); ok = False

names = list(sets)
if not ok:
    OUT = OUT + ".FAILED"           # never leave a failed table under its final name
with open(OUT, "w") as o:
    o.write("statistic\t" + "\t".join(names) + "\n")
    for k in res[names[0]]:
        vals = [res[n][k] for n in names]
        o.write(k + "\t" + "\t".join(str(v) if isinstance(v, int) else "%.2f" % v for v in vals) + "\n")
print("wrote", OUT)
if not ok:
    sys.exit(1)
