#!/usr/bin/env python3
"""wg27 helper. Called by scripts/wg27_blobtoolkit.sh.

prep: per haplotype, build the BlobToolKit inputs
  fields.tsv   identifiers, depth_q0 (samtools coverage meandepth of dip.bam, MAPQ >= 0, from wg25),
               gx_taxon / gx_taxid / gx_div (FCS-GX taxonomy.rpt: rows are split at gaps as `seqid~start..end`;
               aggregated to the seqid by summing score-1 per tax-id-1; top tax-id = max summed score; no rows or
               no tax-id -> "no-hit"), wg25_class (chromosomes -> "chromosome"; unplaced -> wg25 class)
  gc.tsv       identifiers, length, GC fraction = (G+C) / (A+C+G+T) (N excluded), for the fallback plot
  busco.tsv    compleasm full_table_busco_format.tsv rewritten so blobtools' BUSCO parser accepts it: line 1 becomes
               "# BUSCO version is: 5.0.0 (compleasm table, ...)" with no further colon (the parser reads the major version from line 1 and uses
               the v3 layout for < 4); Missing rows (no sequence) dropped; everything else unchanged.
               Logs S/D/F % from the table: S = ids with exactly one Complete row, D = ids with Duplicated rows,
               F = ids with Fragmented rows, as % of N (lineage count from line 2), for the 0.1 pp check vs wg24.
fallback: blob-style plot, GC fraction (x) vs log10 MAPQ>=0 depth (y), point area ~ length, coloured by wg25 class
          (fixed categorical order, validated palette), figures/blobtools/gc_depth_<hap>.png
"""
import argparse
import collections
import math
import re
import sys

CLASSES = ["chromosome", "nuclear", "mixed_keep", "organelle_lowdepth_flag", "candidate_remove"]
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
MARKERS = ["o", "s", "^", "D", "X"]


def fasta_gc(path, keep):
    out, name, cnt = {}, None, None
    def flush():
        if name is not None and (keep is None or name in keep):
            acgt = sum(cnt[b] for b in "ACGT")
            out[name] = (sum(cnt.values()), (cnt["G"] + cnt["C"]) / acgt if acgt else float("nan"))
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                flush()
                name, cnt = line[1:].split()[0], collections.Counter()
            elif keep is None or name in keep:
                cnt.update(line.strip().upper())
    flush()
    return out


def prep(a):
    keep = set(l.strip() for l in open(a.ids)) if a.ids else None
    gc = fasta_gc(a.fasta, keep)
    depth = {}
    for line in open(a.cov):
        if not line.startswith("#"):
            x = line.split("\t")
            depth[x[0]] = float(x[6])
    cls = {}
    for line in open(a.classes):
        if line.startswith("#") or line.startswith("hap\t"):
            continue
        x = line.rstrip("\n").split("\t")
        cls[x[1]] = x[10]
    score, tname, tdiv = collections.defaultdict(lambda: collections.defaultdict(float)), {}, {}
    for line in open(a.gx):
        if line.startswith("#"):
            continue
        x = line.rstrip("\n").split("\t")
        sid = x[0].split("~")[0]
        if len(x) > 10 and x[6] not in ("", "0"):
            score[sid][x[6]] += float(x[10] or 0)
            tname[x[6]] = x[5] or x[6]
            tdiv[x[6]] = x[7] or "NA"
    with open(f"{a.outdir}/fields.tsv", "w") as f, open(f"{a.outdir}/gc.tsv", "w") as g:
        f.write("identifiers\tdepth_q0\tgx_taxon\tgx_taxid\tgx_div\twg25_class\n")
        g.write("identifiers\tlength\tgc\tdepth_q0\twg25_class\n")
        miss = 0
        for sid, (L, gcf) in gc.items():
            c = "chromosome" if re.search(r"_chr\d\d$", sid) else cls.get(sid, "NA")
            if sid not in depth or c == "NA":
                miss += 1
            d = depth.get(sid, 0.0)
            if score[sid]:
                tid = max(score[sid], key=score[sid].get)
                tx, dv = tname[tid].replace(" ", "_"), tdiv[tid]
            else:
                tid, tx, dv = "NA", "no-hit", "no-hit"
            # taxid written as a category ("txid<N>" / "no-hit"): as a bare number blobtools typed the column numeric
            # and crashed on "NA" (smoke 8572764, 2026-11-05)
            f.write(f"{sid}\t{d:.3f}\t{tx}\t{'txid' + tid if tid != 'NA' else 'no-hit'}\t{dv}\t{c}\n")
            g.write(f"{sid}\t{L}\t{gcf:.5f}\t{d:.3f}\t{c}\n")
        print(f"prep: {len(gc)} sequences; {miss} without depth or class")
        if miss and not a.smoke:
            sys.exit("ERROR: sequences without wg25 depth/class")
    if a.busco:
        lines = open(a.busco).read().split("\n")
        n = int(re.search(r"number of BUSCOs:\s*(\d+)", lines[1]).group(1))
        st = collections.defaultdict(list)
        with open(f"{a.outdir}/busco.tsv", "w") as f:
            f.write("# BUSCO version is: 5.0.0 (compleasm table, " + lines[0].lstrip("# ").replace(":", "") + ")\n")  # parser splits on ':' (job 8445487)
            f.write(lines[1] + "\n" + lines[2] + "\n")
            for line in lines[3:]:
                x = line.split("\t")
                if len(x) < 3 or x[1] == "Missing" or not x[2]:
                    continue
                if keep is not None and x[2] not in keep:
                    continue
                st[x[0]].append(x[1])
                f.write(line + "\n")
        S = sum(1 for v in st.values() if v == ["Complete"])
        D = sum(1 for v in st.values() if "Duplicated" in v)
        F = sum(1 for v in st.values() if "Fragmented" in v)
        print(f"busco table: N={n} S={100*S/n:.2f}% D={100*D/n:.2f}% F={100*F/n:.2f}%")


# plain-language legend labels for the paper figure (2026-11-05; internal class names stay in the tables)
LABELS = {"chromosome": "chromosomes",
          "nuclear": "unplaced, nuclear",
          "mixed_keep": "unplaced, 50-90% organelle (kept)",
          "organelle_lowdepth_flag": "unplaced, organelle, < 2x depth (kept)",
          "candidate_remove": "unplaced, organelle, ≥ 2x depth (removed)"}


def fallback(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = [l.rstrip("\n").split("\t") for l in open(a.gc)][1:]
    fig, ax = plt.subplots(figsize=(6.5, 5), dpi=200)
    for c, col, mk in zip(CLASSES, COLORS, MARKERS):
        sel = [r for r in rows if r[4] == c]
        if not sel:
            continue
        x = [float(r[2]) for r in sel]
        y = [math.log10(max(float(r[3]), 0.01)) for r in sel]
        s = [max(6, min(300, int(r[1]) / 2e5)) for r in sel]
        mb = sum(int(r[1]) for r in sel) / 1e6
        ax.scatter(x, y, s=s, c=col, marker=mk, alpha=0.75, edgecolors="white", linewidths=0.5,
                   label=f"{LABELS[c]} (n={len(sel)}, {mb:.2f} Mb)")
    ax.set_xlabel("GC fraction")
    ax.set_ylabel("log10 HiFi read depth")
    ax.set_title(f"Haplotype {a.hap.replace('hap', '')}: GC content vs. read depth")
    ax.grid(color="#e6e5e0", linewidth=0.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(frameon=False, fontsize=7, loc="best")
    fig.tight_layout()
    fig.savefig(a.out)
    print(f"fallback plot -> {a.out}")


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("prep")
    for k in ("fasta", "cov", "classes", "gx", "outdir"):
        p.add_argument("--" + k, required=True)
    p.add_argument("--busco")
    p.add_argument("--ids")
    p.add_argument("--smoke", action="store_true")
    q = sp.add_parser("fallback")
    q.add_argument("--gc", required=True)
    q.add_argument("--hap", required=True)
    q.add_argument("--out", required=True)
    a = ap.parse_args()
    prep(a) if a.cmd == "prep" else fallback(a)


if __name__ == "__main__":
    main()
