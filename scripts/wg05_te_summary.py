#!/usr/bin/env python
"""wg05: transposable elements of Q. wislizeni hap1 for the genome paper, from the per-genome EDTA 2.2.2 run (cg04).

Inputs (cg04 = ~/scratch/california_oaks/cg04_edta/<G>/, the 12 chromosomes renamed chr01..chr12; cg05 = Helixer):
  <G>.fa.mod.EDTA.TEanno.sum          superfamily content (EDTA's own RepeatMasker-style summary, % of the 12 chromosomes)
  <G>.fa.mod.EDTA.TEanno.gff3         all TE features (structural + homology); TE coverage = union of all features except
                                      target_site_duplication
  <G>.fa.mod.EDTA.raw/LTR/<G>.fa.mod.pass.list   intact LTR retrotransposons (LTR_retriever) with 5'/3' LTR identity
  cg05 Qwislizeni.helixer.gff3        genes (hap1, same chromosomes; names Qwis_hap1_chrNN -> chrNN)
Genomes: wislizeni (new hap1) and the three other red oaks (agrifolia, parvula, rubra). agrifolia, parvula and rubra
were rerun in cg04 with upper-case input sequence. These are
per-genome libraries, so TE content differs between genomes partly by library; intact-LTR identities are structural and
library-independent.

LTR age: K = -3/4 ln(1 - 4/3 d), d = 1 - identity (Jukes-Cantor); T = K / (2 mu), mu = argv[3] (default 1.3e-8
substitutions/site/year, LTR_retriever default, rice). Ages at mu_alt = argv[4] (default 1.3e-8) are written as an extra
column in ltr_ages.tsv (sensitivity); ages scale as 1/mu.
Breakpoints (hap1 coordinates): chr01 1,082,415 and 2,905,222 (hap1-private inversion), chr05 50,640,000 and
52,425,000 (inversion on both haplotypes). For each: TE coverage in a 20 kb and a 100 kb window
centred on it, and its percentile among all non-overlapping windows of the same size on the 12 chromosomes; gene count in
the 100 kb window; TE features overlapping +-2 kb; intact LTRs within +-50 kb.
Outputs: <res>/te_superfamilies.tsv, ltr_ages.tsv (per genome x superfamily), ltr_intact_wislizeni.tsv,
         breakpoints_te.tsv, breakpoints_te_features.tsv, te_gene_density_1Mb.tsv, wg05_summary.txt
         <fig>/te_landscape_hap1.{png,pdf} (TE and gene density per 1 Mb, breakpoints marked),
         <fig>/ltr_ages.{png,pdf} (intact LTR insertion ages, Copia / Gypsy / unknown, 4 genomes)
Usage: wg05_te_summary.py <res_dir> <fig_dir> [mu] [mu_alt]
"""
import sys, os, re, math, bisect
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES, FIG = sys.argv[1:3]
MU = float(sys.argv[3]) if len(sys.argv) > 3 else 1.3e-8
MU_ALT = float(sys.argv[4]) if len(sys.argv) > 4 else 1.3e-8
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)
E = "/home/sashanik/scratch/california_oaks/cg04_edta"
HEL = "/home/sashanik/scratch/california_oaks/cg05_helixer/Qwislizeni/Qwislizeni.helixer.gff3"
GEN = [("Qwislizeni", "Q. wislizeni (this study, hap1)"), ("Qagrifolia", "Q. agrifolia"), ("Qparvula", "Q. parvula"),
       ("Qrubra_full", "Q. rubra")]
BPS = [("chr01 inversion, left", "chr01", 1082415), ("chr01 inversion, right", "chr01", 2905222),
       ("chr05 inversion, left", "chr05", 50640000), ("chr05 inversion, right", "chr05", 52425000)]
out = open(f"{RES}/wg05_summary.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); out.write(s + "\n")

# ---- 1. superfamily content (EDTA summary) ----
def read_sum(g):
    rows, cls = [], None
    for l in open(f"{E}/{g}/{g}.fa.mod.EDTA.TEanno.sum"):
        if l.startswith("Repeat Stats"): break
        m = re.match(r"^(\S.*?)\s+--\s+--", l)
        if m: cls = m.group(1); continue
        m = re.match(r"^\s+(\S+(?: \S+)?)\s+(\d+)\s+(\d+)\s+([\d.]+)%", l)
        if m and m.group(1) != "total interspersed":
            rows.append((cls, m.group(1), int(m.group(2)), int(m.group(3)), float(m.group(4))))
        m = re.match(r"^(repeat_fragment)\s+(\d+)\s+(\d+)\s+([\d.]+)%", l)
        if m: rows.append(("other", "repeat_fragment", int(m.group(2)), int(m.group(3)), float(m.group(4))))
        m = re.match(r"^Total\s+(\d+)\s+(\d+)\s+([\d.]+)%", l)
        if m: tot = (int(m.group(1)), int(m.group(2)), float(m.group(3)))
        m = re.match(r"^Total Length: (\d+)", l)
        if m: glen = int(m.group(1))
    return rows, tot, glen
sf = {g: read_sum(g) for g, _ in GEN}
keys = [(c, s) for c, s, *_ in sf["Qwislizeni"][0]]
with open(f"{RES}/te_superfamilies.tsv", "w") as f:
    f.write("class\tsuperfamily\t" + "\t".join(f"{g}_Mb\t{g}_pct" for g, _ in GEN) + "\n")
    for c, s in keys:
        vals = []
        for g, _ in GEN:
            d = {(r[0], r[1]): r for r in sf[g][0]}
            r = d.get((c, s)); vals += [f"{r[3]/1e6:.1f}", f"{r[4]:.2f}"] if r else ["NA", "NA"]
        f.write(f"{c}\t{s}\t" + "\t".join(vals) + "\n")
    f.write("all\ttotal\t" + "\t".join(f"{sf[g][1][1]/1e6:.1f}\t{sf[g][1][2]:.2f}" for g, _ in GEN) + "\n")
    f.write("all\tchromosome_length\t" + "\t".join(f"{sf[g][2]/1e6:.1f}\tNA" for g, _ in GEN) + "\n")
say("TE content of the 12 chromosomes (EDTA per-genome library):")
for g, n in GEN:
    say(f"  {g}: {sf[g][1][2]:.2f}% of {sf[g][2]/1e6:.1f} Mb ({sf[g][1][1]/1e6:.1f} Mb)")

# ---- 2. intact LTR ages ----
def jc_age(idn, mu=None):
    d = 1 - idn
    k = -0.75 * math.log(1 - 4.0 / 3.0 * d)
    return k / (2 * (mu or MU))
ages = {}
with open(f"{RES}/ltr_ages.tsv", "w") as f:
    f.write("genome\tsuperfamily\tn_intact\tmedian_identity\tmedian_age_Myr\tmean_age_Myr\tpct_younger_1Myr\tpct_identity_1\tmedian_age_Myr_mu_alt\n")
    for g, _ in GEN:
        rows = [l.rstrip("\n").split("\t") for l in open(f"{E}/{g}/{g}.fa.mod.EDTA.raw/LTR/{g}.fa.mod.pass.list") if not l.startswith("#")]
        rec = [(r[9], float(r[7]), r[0]) for r in rows]
        ages[g] = rec
        for s in ("all", "Copia", "Gypsy", "unknown"):
            sub = [x for x in rec if s == "all" or x[0] == s]
            idn = np.array([x[1] for x in sub]); a = np.array([jc_age(x) for x in idn]) / 1e6
            a2 = np.array([jc_age(x, MU_ALT) for x in idn]) / 1e6
            f.write(f"{g}\t{s}\t{len(sub)}\t{np.median(idn):.4f}\t{np.median(a):.2f}\t{a.mean():.2f}\t"
                    f"{100*np.mean(a < 1):.1f}\t{100*np.mean(idn >= 1):.1f}\t{np.median(a2):.2f}\n")
say(f"intact LTR retrotransposons (LTR_retriever pass.list), age with mu = {MU:g} (last column: mu_alt = {MU_ALT:g}):")
for l in list(open(f"{RES}/ltr_ages.tsv"))[1:]:
    say("  " + l.rstrip())
with open(f"{RES}/ltr_intact_wislizeni.tsv", "w") as f:
    f.write("location\tsuperfamily\tidentity\tage_Myr\n")
    for s, idn, loc in ages["Qwislizeni"]:
        f.write(f"{loc}\t{s}\t{idn}\t{jc_age(idn)/1e6:.3f}\n")

fig, axs = plt.subplots(1, 3, figsize=(12, 3.6), sharey=False)
cols = {"Qwislizeni": "#c0392b", "Qagrifolia": "#1b9e77", "Qparvula": "#d95f02", "Qrubra_full": "#7570b3"}
bins = np.arange(0, 6.01, 0.2)
for ax, s in zip(axs, ("Copia", "Gypsy", "unknown")):
    for g, n in GEN:
        a = np.array([jc_age(x[1]) for x in ages[g] if x[0] == s]) / 1e6
        h, _ = np.histogram(np.clip(a, 0, 6), bins=bins)
        ax.plot(bins[:-1] + 0.1, h, color=cols[g], lw=2.2 if g == "Qwislizeni" else 1.4, label=f"{n} (n = {len(a)})")
    ax.set_title(f"{s} intact LTR retrotransposons" if s != "unknown" else "unclassified intact LTR retrotransposons", fontsize=10)
    ax.set_xlabel(f"insertion age (Myr; mu = {MU:g}, last bin = 6+)", fontsize=9); ax.set_ylabel("elements per 0.2 Myr", fontsize=9)
    for k in ("top", "right"): ax.spines[k].set_visible(False)
axs[0].legend(fontsize=7.5, frameon=False)
fig.tight_layout()
for ext in ("png", "pdf"): fig.savefig(f"{FIG}/ltr_ages.{ext}", dpi=200, bbox_inches="tight")
plt.close(fig)

# ---- 3. TE coverage intervals (wislizeni hap1) ----
iv, feats = {}, {}
for l in open(f"{E}/Qwislizeni/Qwislizeni.fa.mod.EDTA.TEanno.gff3"):
    if l.startswith("#"): continue
    t = l.rstrip("\n").split("\t")
    if len(t) < 9 or t[2] == "target_site_duplication": continue
    s, e = int(t[3]) - 1, int(t[4])
    iv.setdefault(t[0], []).append((s, e))
    if t[2] not in ("long_terminal_repeat", "repeat_region"):
        cl = re.search(r"Classification=([^;]+)", t[8]); meth = re.search(r"Method=([^;]+)", t[8])
        feats.setdefault(t[0], []).append((s, e, t[2], cl.group(1) if cl else "", meth.group(1) if meth else ""))
merged = {}
for c, v in iv.items():
    v.sort(); m = []
    for s, e in v:
        if m and s <= m[-1][1]: m[-1][1] = max(m[-1][1], e)
        else: m.append([s, e])
    merged[c] = (np.array([x[0] for x in m]), np.array([x[1] for x in m]))
clen = {}
name = None
for l in open(f"{E}/Qwislizeni/Qwislizeni.fa"):
    if l.startswith(">"): name = l[1:].split()[0]; clen[name] = 0
    else: clen[name] += len(l.rstrip())
def cov(c, a, b):
    S, Ee = merged[c]
    i = max(0, bisect.bisect_right(S, a) - 1); tot = 0
    while i < len(S) and S[i] < b:
        tot += max(0, min(Ee[i], b) - max(S[i], a)); i += 1
    return tot / (b - a)
genes = {}
for l in open(HEL):
    if l.startswith("#"): continue
    t = l.split("\t")
    if len(t) > 3 and t[2] == "gene":
        genes.setdefault(re.sub(r".*_(chr\d\d)$", r"\1", t[0]), []).append((int(t[3]) + int(t[4])) // 2)
for c in genes: genes[c].sort()
def ngenes(c, a, b): return bisect.bisect_left(genes.get(c, []), b) - bisect.bisect_left(genes.get(c, []), a)
chroms = sorted(clen)

# background windows
bg = {}
for w in (20000, 100000):
    bg[w] = np.array([cov(c, s, s + w) for c in chroms for s in range(0, clen[c] - w + 1, w)])
say(f"genome-wide TE coverage of non-overlapping windows: 20 kb median {np.median(bg[20000]):.3f}, 100 kb median {np.median(bg[100000]):.3f}")
gw100 = np.array([ngenes(c, s, s + 100000) for c in chroms for s in range(0, clen[c] - 100000 + 1, 100000)])
with open(f"{RES}/breakpoints_te.tsv", "w") as f, open(f"{RES}/breakpoints_te_features.tsv", "w") as g2:
    f.write("breakpoint\tchrom\tpos\tTE_cov_20kb\tpercentile_20kb\tTE_cov_100kb\tpercentile_100kb\tgenes_100kb\tgenome_median_genes_100kb\tintact_LTR_within_50kb\n")
    g2.write("breakpoint\tchrom\tstart\tend\ttype\tclassification\tmethod\tdistance_to_breakpoint\n")
    for lab, c, p in BPS:
        r = []
        for w in (20000, 100000):
            x = cov(c, p - w // 2, p + w // 2); r += [x, 100 * np.mean(bg[w] < x)]
        ltr = [x for x in ages["Qwislizeni"] if x[2].split(":")[0] == c and
               abs((lambda a, b: (min(a, b) + max(a, b)) // 2)(*map(int, x[2].split(":")[1].split(".."))) - p) <= 50000]
        f.write(f"{lab}\t{c}\t{p}\t{r[0]:.3f}\t{r[1]:.0f}\t{r[2]:.3f}\t{r[3]:.0f}\t{ngenes(c, p-50000, p+50000)}\t"
                f"{np.median(gw100):.0f}\t{len(ltr)}\n")
        for s, e, t, cl, me in sorted(feats.get(c, [])):
            if e >= p - 2000 and s <= p + 2000:
                d = 0 if s <= p <= e else (s - p if s > p else e - p)
                g2.write(f"{lab}\t{c}\t{s+1}\t{e}\t{t}\t{cl}\t{me}\t{d}\n")
say("breakpoints (TE coverage and percentile among genome-wide windows):")
for l in list(open(f"{RES}/breakpoints_te.tsv"))[1:]:
    say("  " + l.rstrip())

# ---- 4. landscape: TE and gene density per 1 Mb ----
W = 1_000_000
with open(f"{RES}/te_gene_density_1Mb.tsv", "w") as f:
    f.write("chrom\tstart\tend\tTE_cov\tgenes\n")
    land = {}
    for c in chroms:
        rows = [(s, min(s + W, clen[c])) for s in range(0, clen[c], W)]
        land[c] = [(s, e, cov(c, s, e), ngenes(c, s, e)) for s, e in rows if e - s >= W // 2]
        for s, e, t, n in land[c]: f.write(f"{c}\t{s}\t{e}\t{t:.4f}\t{n}\n")
allte = np.array([x[2] for c in chroms for x in land[c]]); allg = np.array([x[3] for c in chroms for x in land[c]])
rho = np.corrcoef(np.argsort(np.argsort(allte)), np.argsort(np.argsort(allg)))[0, 1]
say(f"1 Mb windows: TE coverage median {np.median(allte):.3f} (range {allte.min():.3f}-{allte.max():.3f}); "
    f"Spearman rho TE vs genes = {rho:.2f} (n = {len(allte)})")
fig, axs = plt.subplots(6, 2, figsize=(12, 10), sharey=True)
gmax = allg.max()
for ax, c in zip(axs.T.flatten(), chroms):
    x = np.array([(s + e) / 2e6 for s, e, _, _ in land[c]])
    ax.fill_between(x, [t * 100 for _, _, t, _ in land[c]], color="#c0392b", alpha=0.35, lw=0, step="mid")
    ax.plot(x, [t * 100 for _, _, t, _ in land[c]], color="#c0392b", lw=1, drawstyle="steps-mid")
    ax.plot(x, [n / gmax * 100 for *_, n in land[c]], color="#1f4e79", lw=1, drawstyle="steps-mid")
    # breakpoints not drawn since 2026-11-10 (paper does not present the inversions); still tabulated in breakpoints_te.tsv
    ax.set_xlim(0, 93); ax.set_ylim(0, 100)
    ax.text(0.99, 0.92, c, transform=ax.transAxes, ha="right", va="top", fontsize=9, fontweight="bold")
    ax.tick_params(labelsize=7)
    for k in ("top", "right"): ax.spines[k].set_visible(False)
for ax in axs[-1]: ax.set_xlabel("position (Mb)", fontsize=8.5)
fig.text(0.005, 0.5, f"TE coverage (%, red) / genes per Mb (blue, scaled to max {gmax})", rotation=90, va="center", fontsize=9)
fig.suptitle("Q. wislizeni hap1: TE coverage and gene density per 1 Mb (EDTA, Helixer)",
             fontsize=10, x=0.01, ha="left")
fig.tight_layout(rect=(0.015, 0, 1, 0.97))
for ext in ("png", "pdf"): fig.savefig(f"{FIG}/te_landscape_hap1.{ext}", dpi=200, bbox_inches="tight")
out.close()
