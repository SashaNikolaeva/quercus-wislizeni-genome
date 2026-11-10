#!/usr/bin/env python
"""rp13f helper: telomere-repeat (CCCTAAA + TTTAGGG) copies per 10 kb window along named chromosomes (windows >= 50 copies
reported, to spot missing / internal telomere arrays), plus unplaced sequences with >= 100 copies in their first or last
20 kb (candidate detached chromosome ends). Usage: rp13f_telomere_scan.py <fasta> <chrom[,chrom...]>"""
import sys
fa, targets = sys.argv[1], sys.argv[2].split(",")
W = 10000
def seqs(path):
    n, b = None, []
    for l in open(path):
        if l[0] == ">":
            if n: yield n, "".join(b)
            n, b = l[1:].split()[0], []
        else: b.append(l.strip())
    if n: yield n, "".join(b)
for n, s in seqs(fa):
    s = s.upper()
    if n in targets:
        hits = [(i, s[i:i+W].count("CCCTAAA") + s[i:i+W].count("TTTAGGG")) for i in range(0, len(s), W)]
        top = [(i, c) for i, c in hits if c >= 50]
        print(f"{n} len {len(s)/1e6:.2f} Mb; 10kb windows with >=50 telomere-repeat copies: " + (", ".join(f"{i/1e6:.2f}Mb:{c}" for i, c in top) or "none"))
    elif "_unplaced_" in n:
        l = s[:20000].count("CCCTAAA") + s[:20000].count("TTTAGGG"); r = s[-20000:].count("CCCTAAA") + s[-20000:].count("TTTAGGG")
        if max(l, r) >= 100: print(f"  unplaced with telomeric end: {n} len {len(s)/1e3:.0f} kb left {l} right {r}")
