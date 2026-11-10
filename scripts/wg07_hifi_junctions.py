#!/usr/bin/env python
"""wg07: are the wislizeni inversions supported by HiFi reads?

HiFi reads (SRR19243399, ~47x) mapped to the DIPLOID reference hap1 + hap2 (all sequences), so each read lands on the
haplotype it matches best (heterozygosity ~1%: a 14 kb read carries ~100 differences).
Candidates: every hap1-vs-hap2 INV >= 100 kb in SyRI (cg03 Qwislizeni_hap2, asm5) plus the chr06 INVTR (the chr06
complex), and the chr05 inversion that both haplotypes carry vs the other oaks (coordinates set in this script).
Per candidate and haplotype, at each junction:
  zone     = the gap between the inverted block and the nearest SyRI SYN block (hap1 vs hap2) +- 5 kb, capped at 500 kb;
             chr05 (collinear between haplotypes) = breakpoint +- 300 kb
  spanning = number of primary reads covering a position with >= 2 kb aligned on each side, every 1 kb across the zone;
             min over the zone (MAPQ >= 0, and MAPQ >= 20). A join the reads do not support (orientation error, or a
             repeat longer than the reads) gives a dip to ~0
  alt      = reads on this haplotype that support the OTHER arrangement: one segment in the left zone and one in the right
             zone (primary + supplementary, from the SA tag), on opposite strands
  depth    = mean MAPQ >= 20 depth inside the block / median chromosome depth of that haplotype (1 = haploid; ~2 = the
             other haplotype collapsed onto this one)
Verdict:
  heterozygous (hap1 vs hap2): every junction on both haplotypes min spanning >= 5, alt <= 1 on each haplotype, depth
    ratio 0.6-1.4 on both -> "read-supported heterozygous"; any junction with min spanning 0 and alt <= 1 -> "unresolved
    (junction not spanned by reads)"; otherwise "not supported"
  chr05 (homozygous): same rule on both haplotypes -> "read-supported (both haplotypes)"
Usage: wg07_hifi_junctions.py <dip.bam> <syri.out (hap1 vs hap2, asm5)> <coverage_q20.tsv> <res_dir> <fig_dir> <trace_dir>
"""
import sys, os, re, subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BAM, SYRI, COV, RES, FIG, TR = sys.argv[1:7]
for d in (RES, FIG, TR):
    os.makedirs(d, exist_ok=True)
ANCHOR, STEP, PAD, CAP, MINLEN = 2000, 1000, 5000, 500000, 100000
H = {"hap1": "Qwis_hap1_", "hap2": "Qwis_hap2_"}
CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def ref_len(cigar):
    return sum(int(n) for n, op in CIG.findall(cigar) if op in "MDN=X")


# SyRI hap1 vs hap2: SYN blocks (for junction zones) and the candidate inversions
syn = {}
cand = []
for l in open(SYRI):
    t = l.rstrip("\n").split("\t")
    if len(t) < 11:
        continue
    c, ty = t[0], t[10]
    if ty not in ("SYN", "INV", "INVTR"):
        continue
    s, e = sorted((int(t[1]), int(t[2])))
    qs, qe = sorted((int(t[6]), int(t[7])))
    if ty == "SYN":
        syn.setdefault(c, []).append((s, e, qs, qe))
    elif e - s + 1 >= MINLEN and t[5] == c:
        cand.append(dict(id=f"{c}_{ty}_{s/1e6:.2f}", chr=c, type=ty, h1=(s, e), h2=(qs, qe)))


def zones(c, s, e, k):
    """left and right junction zones on haplotype k (0 = hap1 = SyRI ref, 1 = hap2 = SyRI qry)"""
    blocks = [(b[2 * k], b[2 * k + 1]) for b in syn.get(c, [])]
    left = max([b[1] for b in blocks if b[1] <= s], default=s - CAP)
    right = min([b[0] for b in blocks if b[0] >= e], default=e + CAP)
    capped = (s - left > CAP) or (right - e > CAP)
    left, right = max(left, s - CAP), min(right, e + CAP)
    return (left - PAD, s + PAD), (e - PAD, right + PAD), capped


for C in cand:
    C["zones"] = {h: zones(C["chr"], *C[h], k) for k, h in enumerate(("h1", "h2"))}
hom = dict(id="chr05_hom_50.64", chr="chr05", type="INV_both_haps", h1=(50640000, 52425000), h2=(50335000, 52115000))
hom["zones"] = {h: ((hom[h][0] - 300000, hom[h][0] + 300000), (hom[h][1] - 300000, hom[h][1] + 300000), False)
                for h in ("h1", "h2")}
cand.append(hom)


def view(region, flags):
    p = subprocess.run(["samtools", "view", "-F", flags, BAM, region], capture_output=True, text=True, check=True)
    return [l.split("\t") for l in p.stdout.splitlines()]


def spanning(ctg, z):
    """min / median number of primary reads spanning each 1 kb position of zone z (>= ANCHOR bp aligned each side)"""
    a, b = z
    grid = np.arange(a, b + 1, STEP)
    out, mq = {}, []
    recs = view(f"{ctg}:{max(1, a - 60000)}-{b + 60000}", "0x904")
    for q in (0, 20):
        d = np.zeros(len(grid) + 1)
        for r in recs:
            if int(r[4]) < q:
                continue
            st = int(r[3])
            en = st + ref_len(r[5]) - 1
            lo, hi = np.searchsorted(grid, st + ANCHOR), np.searchsorted(grid, en - ANCHOR, side="right")
            if hi > lo:
                d[lo] += 1
                d[hi] -= 1
        tr = np.cumsum(d)[:-1]
        out[q] = tr
    mq = [int(r[4]) for r in recs if int(r[3]) <= b and int(r[3]) + ref_len(r[5]) >= a]
    return grid, out, (np.mean(mq) if mq else np.nan)


def alt_reads(ctg, zl, zr):
    """reads with one segment in each junction zone on opposite strands (support the other arrangement)"""
    segs = {}
    for z in (zl, zr):
        for r in view(f"{ctg}:{max(1, z[0])}-{z[1]}", "0x104"):
            name, flag = r[0], int(r[1])
            S = segs.setdefault(name, set())
            S.add((r[2], int(r[3]), int(r[3]) + ref_len(r[5]) - 1, "-" if flag & 16 else "+"))
            for f in r[11:]:
                if f.startswith("SA:Z:"):
                    for sa in f[5:].strip(";").split(";"):
                        rn, pos, strand, cig = sa.split(",")[:4]
                        S.add((rn, int(pos), int(pos) + ref_len(cig) - 1, strand))
    n = 0
    for S in segs.values():
        L = [x for x in S if x[0] == ctg and x[1] <= zl[1] and x[2] >= zl[0]]
        R = [x for x in S if x[0] == ctg and x[1] <= zr[1] and x[2] >= zr[0]]
        if any(x[3] != y[3] for x in L for y in R if x != y):
            n += 1
    return n


def block_depth(ctg, s, e):
    p = subprocess.run(["samtools", "coverage", "-q", "20", "-r", f"{ctg}:{s}-{e}", BAM], capture_output=True, text=True,
                       check=True)
    return float(p.stdout.splitlines()[1].split("\t")[6])


base = {}
cov = [l.rstrip("\n").split("\t") for l in open(COV) if not l.startswith("#")]
for h, pre in H.items():
    v = [float(t[6]) for t in cov if re.fullmatch(pre + r"chr\d\d", t[0])]
    base[h] = float(np.median(v))
print(f"median chromosome depth (MAPQ>=20): hap1 {base['hap1']:.1f}, hap2 {base['hap2']:.1f}")

rows, ver, traces = [], [], []
for C in cand:
    res = {}
    for hk, h in (("h1", "hap1"), ("h2", "hap2")):
        ctg = H[h] + C["chr"]
        zl, zr, capped = C["zones"][hk]
        a = alt_reads(ctg, zl, zr)
        dr = block_depth(ctg, *C[hk]) / base[h]
        for side, z in (("left", zl), ("right", zr)):
            g, tr, mq = spanning(ctg, z)
            i = int(np.argmin(tr[0]))
            rows.append([C["id"], h, side, ctg, z[0], z[1], z[1] - z[0], int(capped), int(tr[0].min()), int(g[i]),
                         int(np.median(tr[0])), int(tr[20].min()), f"{mq:.1f}", a, f"{dr:.2f}"])
            for x, y0, y20 in zip(g, tr[0], tr[20]):
                traces.append([C["id"], h, side, int(x), int(y0), int(y20)])
            res.setdefault(h, []).append((int(tr[0].min()), a, dr))
    mins = [m for h in res for m, _, _ in res[h]]
    alts = [res[h][0][1] for h in res]
    drs = [res[h][0][2] for h in res]
    if min(mins) >= 5 and max(alts) <= 1 and all(0.6 <= d <= 1.4 for d in drs):
        v = "read-supported (both haplotypes)" if C["type"] == "INV_both_haps" else "read-supported heterozygous"
    elif min(mins) == 0 and max(alts) <= 1:
        v = "unresolved (junction not spanned by reads)"
    else:
        v = "not supported"
    ver.append([C["id"], C["type"], C["chr"], *C["h1"], *C["h2"], min(mins), *alts, *(f"{d:.2f}" for d in drs), v])
    print(C["id"], v)

with open(f"{RES}/junction_support.tsv", "w") as f:
    f.write("candidate\thap\tjunction\tcontig\tzone_start\tzone_end\tzone_bp\tzone_capped\tmin_spanning_q0\tmin_at\t"
            "median_spanning_q0\tmin_spanning_q20\tmean_mapq\talt_arrangement_reads\tblock_depth_ratio\n")
    for r in rows:
        f.write("\t".join(map(str, r)) + "\n")
with open(f"{RES}/inversion_verdicts.tsv", "w") as f:
    f.write("candidate\ttype\tchr\thap1_start\thap1_end\thap2_start\thap2_end\tmin_spanning_all_junctions\talt_hap1\t"
            "alt_hap2\tdepth_ratio_hap1\tdepth_ratio_hap2\tverdict\n")
    for r in ver:
        f.write("\t".join(map(str, r)) + "\n")
with open(f"{TR}/spanning_traces.tsv", "w") as f:
    f.write("candidate\thap\tjunction\tpos\tspanning_q0\tspanning_q20\n")
    for r in traces:
        f.write("\t".join(map(str, r)) + "\n")

# one panel per candidate x haplotype: spanning traces at both junctions
n = len(cand)
fig, axs = plt.subplots(n, 4, figsize=(14, 1.9 * n), squeeze=False)
for i, C in enumerate(cand):
    for j, (h, side) in enumerate((("hap1", "left"), ("hap1", "right"), ("hap2", "left"), ("hap2", "right"))):
        ax = axs[i, j]
        sel = [t for t in traces if t[0] == C["id"] and t[1] == h and t[2] == side]
        x = np.array([t[3] for t in sel]) / 1e6
        ax.plot(x, [t[4] for t in sel], color="black", lw=1, label="MAPQ>=0")
        ax.plot(x, [t[5] for t in sel], color="#888888", lw=1, label="MAPQ>=20")
        ax.axhline(5, color="#c0392b", lw=0.6, ls="--")
        ax.set_ylim(bottom=0)
        ax.tick_params(labelsize=6)
        ax.set_title(f"{C['id']} {h} {side}", fontsize=7, loc="left")
        if i == 0 and j == 0:
            ax.legend(fontsize=6, frameon=False)
fig.supxlabel("position (Mb, own haplotype)", fontsize=8)
fig.supylabel("HiFi reads spanning (>= 2 kb each side)", fontsize=8)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{FIG}/junction_traces.{ext}", dpi=150)
print(f"written {RES}/inversion_verdicts.tsv, junction_support.tsv; {FIG}/junction_traces.png")
