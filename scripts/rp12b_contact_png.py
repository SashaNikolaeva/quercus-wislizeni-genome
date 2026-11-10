#!/usr/bin/env python
"""rp12b: whole-genome Omni-C contact heatmap of the 12 largest scaffolds from a `juicer pre` (non -a) text file.

Input lines (juicer short format): str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 (whitespace separated).
Bins of 1 Mb, symmetric matrix, colour = log10(contacts + 1). Scaffold boundaries drawn as lines.
Usage: rp12b_contact_png.py <juicer_pre.txt> <chrom.sizes> <out.png> <title>
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BIN = 1_000_000
txt, sizes_f, out, title = sys.argv[1:5]
sizes = [l.split() for l in open(sizes_f)]
sizes = sorted(((n, int(s)) for n, s in sizes), key=lambda x: -x[1])[:12]
offset, edges, tot = {}, [], 0
for n, s in sizes:
    offset[n] = tot
    tot += s // BIN + 1
    edges.append(tot)
mat = np.zeros((tot, tot), dtype=np.float64)
with open(txt) as fh:
    for line in fh:
        f = line.split()
        c1, c2 = f[1], f[5]
        if c1 not in offset or c2 not in offset:
            continue
        i = offset[c1] + int(f[2]) // BIN
        j = offset[c2] + int(f[6]) // BIN
        mat[i, j] += 1
        if i != j:
            mat[j, i] += 1
fig, ax = plt.subplots(figsize=(10, 10), dpi=150)
ax.imshow(np.log10(mat + 1), cmap="Reds", interpolation="none", origin="upper")
for e in edges[:-1]:
    ax.axhline(e - 0.5, color="0.4", lw=0.4)
    ax.axvline(e - 0.5, color="0.4", lw=0.4)
mids = [offset[n] + (s // BIN + 1) / 2 for n, s in sizes]
ax.set_xticks(mids)
ax.set_xticklabels([n for n, _ in sizes], rotation=90, fontsize=6)
ax.set_yticks(mids)
ax.set_yticklabels([n for n, _ in sizes], fontsize=6)
ax.set_title(f"{title}: Omni-C contacts, 12 largest YaHS scaffolds, 1 Mb bins, log10(n+1)")
fig.tight_layout()
fig.savefig(out)
print(f"wrote {out}: {tot} bins, {int(mat.sum())} contacts counted")
