#!/usr/bin/env python
"""wg24b: saved check of protein-mode compleasm odb12 'Duplicated' BUSCOs .

Competing explanations, fixed before the test:
  (a) artifact: fragments or TE-derived gene models are scored above the cutoff;
  (b) real paralogs.
  Prediction of (a): duplicated copies are enriched, relative to single-copy BUSCOs, for TE-derived models
  (hap1 only, gene_te_overlap.tsv frac_all > 0.5) or for length ratios < 0.5. Reported whichever way it falls.

Usage: wg24b_busco_dup_check.py <out.tsv> <out_summary.tsv> <gene_te_overlap.tsv> <refseq_db.faa.gz>
           hap1 <full_table.tsv> <prot.fa> <gff3>  hap2 <full_table.tsv> <prot.fa> <gff3>
Definitions:
  copies       = rows of full_table.tsv with Status 'Duplicated' for that BUSCO (one row per protein)
  n_copies     = number of distinct proteins (mRNA) among them
  within_3_genes = 'yes' if at least one pair of copies lies on the same chromosome with |gene index difference| <= 3,
                 gene index = rank of the gene by start coordinate among the genes of that chromosome in the GFF3
  TE_derived   = per copy, hap1 only (frac_all > 0.5), 'NA' for hap2
  len_ratio_single = copy protein length / median length of the proteins that are 'Single' for the same BUSCO in the
                 protein-mode tables of hap1 and hap2 pooled (the spec's 'median length of the BUSCO's complete
                 single-copy hits'; NA if the BUSCO is Single in neither haplotype)
  len_ratio_refseq = copy protein length / median length of the BUSCO's sequences in the odb12 refseq_db.faa.gz
                 (reference independent of our proteome, so the same ratio can be computed for single-copy BUSCOs)
  Protein length = residues in prot_hap<h>.fa (stop removed).
Summary (out_summary.tsv), per haplotype, copies of Duplicated BUSCOs vs proteins of Single BUSCOs:
  n, n TE-derived (hap1), n len_ratio_refseq < 0.5, and a two-sided Fisher exact test of each 2x2 table;
  number of Duplicated BUSCOs with within_3_genes = yes.
"""
import sys, gzip, re
from collections import defaultdict
import numpy as np
from scipy.stats import fisher_exact

OUT, SUM, TEO, REF = sys.argv[1:5]
rest = sys.argv[5:]
haps = [rest[i:i + 4] for i in range(0, len(rest), 4)]
core = lambda b: re.sub(r"at\d+$", "", b)

te = {}
with open(TEO) as f:
    h = f.readline().rstrip("\n").split("\t"); gi, mi, fi = h.index("gene"), h.index("mrna"), h.index("frac_all")
    for l in f:
        t = l.rstrip("\n").split("\t"); te[t[gi]] = int(t[h.index("te_bp_all")]) / int(t[h.index("cds_bp")]) > 0.5  # exact, not rounded frac_all (fix 2026-11-05)

reflen = defaultdict(list)
with gzip.open(REF, "rt") as f:
    name, n = None, 0
    for l in f:
        if l.startswith(">"):
            if name: reflen[name].append(n)
            name, n = core(l[1:].split("_")[0]), 0
        else:
            n += len(l.strip())
    if name: reflen[name].append(n)
refmed = {b: float(np.median(v)) for b, v in reflen.items()}


def read_hap(ft, prot, gff):
    plen, name = {}, None
    for l in open(prot):
        if l.startswith(">"): name = l[1:].split()[0]; plen[name] = 0
        else: plen[name] += len(l.strip())
    genes = defaultdict(list); m2g = {}; gchr = {}
    for l in open(gff):
        if l.startswith("#"): continue
        t = l.rstrip("\n").split("\t")
        if len(t) < 9: continue
        a = dict(kv.split("=", 1) for kv in t[8].strip().split(";") if "=" in kv)
        if t[2] == "gene": genes[t[0]].append((int(t[3]), a["ID"])); gchr[a["ID"]] = t[0]
        elif t[2] == "mRNA": m2g[a["ID"]] = a["Parent"]
    gidx = {}
    for c, v in genes.items():
        for i, (_, g) in enumerate(sorted(v)): gidx[g] = i
    status = defaultdict(lambda: defaultdict(set))   # status -> busco -> proteins
    for l in open(ft):
        if l.startswith("#"): continue
        t = l.rstrip("\n").split("\t")
        if len(t) < 3 or t[2] in ("", "-", "None"): continue
        status[t[1]][core(t[0])].add(t[2])
    return plen, m2g, gchr, gidx, status


data = {h: read_hap(*a) for h, *a in haps}
single_len = defaultdict(list)
for h, (plen, m2g, gchr, gidx, st) in data.items():
    for b, ps in st["Single"].items():
        single_len[b] += [plen[p] for p in ps if p in plen]
single_med = {b: float(np.median(v)) for b, v in single_len.items() if v}

fmt = lambda x: "NA" if x is None else "%.3f" % x
with open(OUT, "w") as o, open(SUM, "w") as s:
    o.write("hap\tbusco\tn_copies\twithin_3_genes\tprotein\tgene\tchrom\tgene_index\tprotein_len\tTE_derived\t"
            "median_len_single\tlen_ratio_single\tmedian_len_refseq\tlen_ratio_refseq\n")
    s.write("hap\tgroup\tn_busco\tn_proteins\tn_TE_derived\tn_len_ratio_refseq_lt_0.5\tn_busco_within_3_genes\t"
            "fisher_p_TE_dup_vs_single\tfisher_p_short_dup_vs_single\n")
    for h, (plen, m2g, gchr, gidx, st) in data.items():
        cnt = {}
        for grp in ("Duplicated", "Single"):
            n_te = n_short = n_prot = n_adj = 0
            for b, ps in sorted(st[grp].items()):
                ps = sorted(ps)
                gs = [m2g.get(p) for p in ps]
                adj = False
                for i in range(len(gs)):
                    for j in range(i + 1, len(gs)):
                        gi_, gj = gs[i], gs[j]
                        if gi_ and gj and gchr.get(gi_) == gchr.get(gj) and abs(gidx[gi_] - gidx[gj]) <= 3:
                            adj = True
                n_adj += adj
                for p, g in zip(ps, gs):
                    L = plen.get(p)
                    tder = ("yes" if te.get(g) else "no") if h == "hap1" else "NA"
                    rs = L / single_med[b] if (L and b in single_med) else None
                    rr = L / refmed[b] if (L and b in refmed) else None
                    n_prot += 1; n_te += tder == "yes"; n_short += (rr is not None and rr < 0.5)
                    if grp == "Duplicated":
                        o.write("\t".join(map(str, [h, b, len(ps), "yes" if adj else "no", p, g, gchr.get(g, "NA"),
                                gidx.get(g, "NA"), L if L else "NA", tder, fmt(single_med.get(b)), fmt(rs),
                                fmt(refmed.get(b)), fmt(rr)])) + "\n")
            cnt[grp] = (len(st[grp]), n_prot, n_te, n_short, n_adj)
        for grp in ("Duplicated", "Single"):
            nb, npr, nte, nsh, nadj = cnt[grp]
            if grp == "Duplicated":
                d, sg = cnt["Duplicated"], cnt["Single"]
                pte = fisher_exact([[d[2], d[1] - d[2]], [sg[2], sg[1] - sg[2]]])[1] if h == "hap1" else None
                psh = fisher_exact([[d[3], d[1] - d[3]], [sg[3], sg[1] - sg[3]]])[1]
                ptxt = ["NA" if pte is None else "%.3g" % pte, "%.3g" % psh]
            else:
                ptxt = ["", ""]
            s.write("\t".join(map(str, [h, grp, nb, npr, nte if h == "hap1" else "NA", nsh,
                                        nadj if grp == "Duplicated" else "NA"] + ptxt)) + "\n")
print(open(SUM).read())
