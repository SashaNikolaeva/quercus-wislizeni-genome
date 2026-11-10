#!/usr/bin/env python
"""wg03: Figure 2 of the paper.
(1) per haplotype, whole-genome Omni-C contact map of the 12 named chromosomes in chr01..chr12 order, 1 Mb bins,
    log10(n+1). Pairs = yahs/<hap>/mapped.pairs.gz (contig coordinates; unique MAPQ >= 30 deduplicated pairs);
    contig -> YaHS scaffold (AGP) -> named chromosome (rp13 names.tsv; '-' = reverse-complemented).
(2) composite: A GenomeScope2 profile (rp08), B Merqury spectra-cn of hap1+hap2 (rp10), C/D hap1/hap2 contact maps.
Usage: wg03_fig2.py <wislizeni_reasm dir> <repo figures/wislizeni_reasm dir> <out dir>
"""
import sys, gzip
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

W, FIGIN, OUT = sys.argv[1:4]
BIN = 1_000_000

def contact_matrix(hap):
    rows = [l.rstrip("\n").split("\t") for l in list(open(f"{W}/final/Qwis_{hap}.names.tsv"))[1:]]
    chroms = sorted((r for r in rows if "_chr" in r[0]), key=lambda r: r[0])
    scaf = {r[1]: (r[0], int(r[2]), r[4]) for r in chroms}          # yahs scaffold -> (name, length, orientation)
    comp = {}                                                         # contig -> (chrom, scaf start, ctg start, ctg end, ori)
    for l in open(f"{W}/yahs/{hap}/{hap}_yahs_scaffolds_final.agp"):
        f = l.rstrip("\n").split("\t")
        if f[0] in scaf and f[4] == "W":
            comp[f[5]] = (f[0], int(f[1]), int(f[6]), int(f[7]), f[8])
    names = [r[0] for r in chroms]; L = {r[0]: int(r[2]) for r in chroms}
    off, tot = {}, 0
    for n in names:
        off[n] = tot; tot += L[n] // BIN + 1
    def to_bin(ctg, pos):
        sc, s0, cs, ce, o = comp[ctg]
        n, Ls, so = scaf[sc]
        sp = s0 + (pos - cs) if o == "+" else s0 + (ce - pos)
        cp = sp if so == "+" else Ls - sp + 1
        return off[n] + cp // BIN
    M = np.zeros((tot, tot))
    with gzip.open(f"{W}/yahs/{hap}/mapped.pairs.gz", "rt") as fh:
        for l in fh:
            if l[0] == "#": continue
            f = l.split("\t", 6)
            if f[1] not in comp or f[3] not in comp: continue
            i, j = to_bin(f[1], int(f[2])), to_bin(f[3], int(f[4]))
            M[i, j] += 1
            if i != j: M[j, i] += 1
    edges = [off[n] for n in names] + [tot]
    return M, names, edges

def draw_contacts(ax, M, names, edges, title):
    ax.imshow(np.log10(M + 1), cmap="Reds", interpolation="none", origin="upper")
    for e in edges[1:-1]:
        ax.axhline(e - 0.5, color="0.45", lw=0.3); ax.axvline(e - 0.5, color="0.45", lw=0.3)
    mids = [(edges[k] + edges[k + 1]) / 2 for k in range(len(names))]
    lab = [n.split("_")[-1].replace("chr", "") for n in names]
    ax.set_xticks(mids); ax.set_xticklabels(lab, fontsize=6)
    ax.set_yticks(mids); ax.set_yticklabels(lab, fontsize=6)
    ax.tick_params(length=0)
    ax.set_title(title, fontsize=9)

mats = {}
for hap in ("hap1", "hap2"):
    M, names, edges = contact_matrix(hap)
    mats[hap] = (M, names, edges)
    np.save(f"{OUT}/contacts_1Mb_{hap}.npy", M)
    fig, ax = plt.subplots(figsize=(6, 6), dpi=200)
    draw_contacts(ax, M, names, edges, f"Q. wislizeni v2 {hap}: Omni-C contacts, 1 Mb bins, log10(n+1)")
    fig.tight_layout(); fig.savefig(f"{OUT}/contacts_chr_1Mb_{hap}.png"); plt.close(fig)
    print(hap, int(M.sum()), "contacts on chromosomes;", f"{100*np.trace(M)/M.sum():.1f}% in diagonal bins")

fig = plt.figure(figsize=(10, 10.5), dpi=200)
panels = [("A", f"{FIGIN}/genomescope2_k21_linear.png"), ("B", f"{FIGIN}/merqury_hap1+hap2_spectra-cn.png")]
for k, (lab, png) in enumerate(panels):
    ax = fig.add_subplot(2, 2, k + 1); ax.imshow(mpimg.imread(png)); ax.axis("off")
    ax.text(-0.02, 1.02, lab, transform=ax.transAxes, fontsize=16, fontweight="bold", va="bottom")
for k, (lab, hap) in enumerate((("C", "hap1"), ("D", "hap2"))):
    ax = fig.add_subplot(2, 2, k + 3)
    draw_contacts(ax, *mats[hap], f"{hap} (chr01-chr12)")
    ax.text(-0.02, 1.06, lab, transform=ax.transAxes, fontsize=16, fontweight="bold", va="bottom")
fig.tight_layout(); fig.savefig(f"{OUT}/fig2_assembly_qc.png"); fig.savefig(f"{OUT}/fig2_assembly_qc.pdf")
print("wrote", f"{OUT}/fig2_assembly_qc.png")
