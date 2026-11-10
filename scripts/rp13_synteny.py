#!/usr/bin/env python
"""rp13 helper: chromosome naming, synteny summaries, dotplots and telomere counts for the YaHS wislizeni scaffolds.

Modes
  name    <query.fa> <query_vs_lobata12.paf> <lobata12.fa.fai> <prefix> <out.fa> <out_names.tsv>
          The 12 largest scaffolds get the lobata chromosome with the most aligned query bp (PAF blocks >= 5 kb,
          MAPQ >= 1); reverse-complemented if most of those bp are on the minus strand. Chromosome number = order of
          the sequence in lobata12.fa.fai (CM012293.1 = chr01 ... CM012304.1 = chr12). Everything else becomes
          <prefix>_unplaced_NNNN, largest first.
  report  <named.fa.fai> <ref_label> <paf> <out_prefix>
          Per chromosome: best target sequence, % of aligned bp on it, % minus-strand bp, aligned bp / length;
          for each target sequence >= 1 Mb: how its aligned bp splits over the new chromosomes (broken old joins).
          Writes <out_prefix>.chrom.tsv, <out_prefix>.target.tsv, <out_prefix>.png (dotplot, blocks >= 10 kb).
  telo    <named.fa> <out.tsv>
          CCCTAAA + TTTAGGG copies in the first and last 20 kb of every chromosome; an end is "capped" at >= 100 copies.
"""
import sys
from collections import defaultdict

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
MIN_BLOCK_NAME, MIN_BLOCK_PLOT = 5_000, 10_000


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name:
        yield name, "".join(seq)


def paf(path, min_block, min_mapq=1):
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if int(f[10]) < min_block or int(f[11]) < min_mapq:
                continue
            yield dict(q=f[0], ql=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4],
                       t=f[5], tl=int(f[6]), ts=int(f[7]), te=int(f[8]), aln=int(f[10]))


def write_fa(fh, name, seq):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), 80):
        fh.write(seq[i:i + 80] + "\n")


def mode_name(qfa, pafp, lob_fai, prefix, out_fa, out_tsv):
    lob_order = [l.split()[0] for l in open(lob_fai)]
    seqs = dict(read_fasta(qfa))
    top12 = sorted(seqs, key=lambda n: -len(seqs[n]))[:12]
    bp = defaultdict(lambda: defaultdict(int))  # q -> (t, strand) -> bp
    for r in paf(pafp, MIN_BLOCK_NAME):
        if r["q"] in top12:
            bp[r["q"]][(r["t"], r["strand"])] += r["qe"] - r["qs"]
    rows, used = [], set()
    for q in top12:
        per_t = defaultdict(int)
        for (t, s), v in bp[q].items():
            per_t[t] += v
        tot = sum(per_t.values())
        best = max(per_t, key=per_t.get)
        minus = bp[q][(best, "-")]
        frac = per_t[best] / tot if tot else 0
        rows.append((q, best, lob_order.index(best) + 1, "-" if minus > per_t[best] / 2 else "+", frac, tot))
        used.add(best)
    rows.sort(key=lambda r: r[2])
    with open(out_fa, "w") as fo, open(out_tsv, "w") as ft:
        ft.write("new_name\tyahs_name\tlength\tlobata_chr\torientation\tfrac_aligned_bp_on_best\taligned_bp\n")
        for q, best, num, ori, frac, tot in rows:
            n = f"{prefix}_chr{num:02d}"
            s = seqs[q] if ori == "+" else seqs[q].translate(COMP)[::-1]
            write_fa(fo, n, s)
            ft.write(f"{n}\t{q}\t{len(s)}\t{best}\t{ori}\t{frac:.4f}\t{tot}\n")
        rest = sorted((n for n in seqs if n not in top12), key=lambda n: -len(seqs[n]))
        for i, q in enumerate(rest, 1):
            n = f"{prefix}_unplaced_{i:04d}"
            write_fa(fo, n, seqs[q])
            ft.write(f"{n}\t{q}\t{len(seqs[q])}\t-\t+\tNA\tNA\n")
    print(f"named 12 chromosomes; distinct lobata chromosomes used: {len(used)}; unplaced: {len(seqs) - 12}")
    for r in rows:
        print(f"  {r[0]} -> chr{r[2]:02d} ({r[1]}, {r[3]}), {100 * r[4]:.1f}% of aligned bp on it")


def mode_report(fai, label, pafp, outp):
    lens = {l.split()[0]: int(l.split()[1]) for l in open(fai)}
    chroms = [n for n in lens if "_chr" in n]
    q_t = defaultdict(lambda: defaultdict(int))
    q_minus = defaultdict(lambda: defaultdict(int))
    t_q = defaultdict(lambda: defaultdict(int))
    tlen, blocks = {}, []
    for r in paf(pafp, MIN_BLOCK_NAME):
        v = r["qe"] - r["qs"]
        q_t[r["q"]][r["t"]] += v
        if r["strand"] == "-":
            q_minus[r["q"]][r["t"]] += v
        t_q[r["t"]][r["q"]] += r["te"] - r["ts"]
        tlen[r["t"]] = r["tl"]
        if r["aln"] >= MIN_BLOCK_PLOT and r["q"] in chroms:
            blocks.append(r)
    with open(f"{outp}.chrom.tsv", "w") as fo:
        fo.write("chrom\tlength\tref\tbest_target\tpct_aligned_bp_on_best\tpct_minus_on_best\taligned_bp_over_length\n")
        for c in chroms:
            tot = sum(q_t[c].values())
            if not tot:
                fo.write(f"{c}\t{lens[c]}\t{label}\tNA\t0\t0\t0\n")
                continue
            b = max(q_t[c], key=q_t[c].get)
            fo.write(f"{c}\t{lens[c]}\t{label}\t{b}\t{100 * q_t[c][b] / tot:.1f}\t"
                     f"{100 * q_minus[c][b] / q_t[c][b]:.1f}\t{tot / lens[c]:.3f}\n")
    with open(f"{outp}.target.tsv", "w") as fo:
        fo.write("target\ttarget_length\tref\tn_new_chroms_ge10pct\tsplit\n")
        for t in sorted(t_q, key=lambda t: -tlen[t]):
            if tlen[t] < 1_000_000:
                continue
            tot = sum(t_q[t].values())
            parts = sorted(((v / tot, q) for q, v in t_q[t].items() if v / tot >= 0.10), reverse=True)
            fo.write(f"{t}\t{tlen[t]}\t{label}\t{len(parts)}\t" + ",".join(f"{q}:{100 * p:.0f}%" for p, q in parts) + "\n")
    # dotplot: x = new chromosomes, y = target sequences with >= 1 Mb aligned to the chromosomes
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    xoff, x = {}, 0
    for c in chroms:
        xoff[c] = x
        x += lens[c]
    tgt_bp = defaultdict(int)
    for r in blocks:
        tgt_bp[r["t"]] += r["te"] - r["ts"]
    # order targets by the chromosome (and position) they mostly align to
    tpos = {}
    for t in tgt_bp:
        best = max((q for q in t_q[t] if q in xoff), key=lambda q: t_q[t][q], default=None)
        tpos[t] = xoff.get(best, 1e12)
    targets = sorted((t for t, v in tgt_bp.items() if v >= 1_000_000), key=lambda t: tpos[t])
    yoff, y = {}, 0
    for t in targets:
        yoff[t] = y
        y += tlen[t]
    fig, ax = plt.subplots(figsize=(10, 10), dpi=150)
    for r in blocks:
        if r["t"] not in yoff:
            continue
        x0, x1 = xoff[r["q"]] + r["qs"], xoff[r["q"]] + r["qe"]
        y0, y1 = yoff[r["t"]] + r["ts"], yoff[r["t"]] + r["te"]
        if r["strand"] == "-":
            y0, y1 = y1, y0
        ax.plot([x0 / 1e6, x1 / 1e6], [y0 / 1e6, y1 / 1e6], color="#1f5fa8" if r["strand"] == "+" else "#c0392b", lw=0.8)
    for v in xoff.values():
        ax.axvline(v / 1e6, color="0.85", lw=0.4)
    for v in yoff.values():
        ax.axhline(v / 1e6, color="0.85", lw=0.4)
    ax.set_xticks([(xoff[c] + lens[c] / 2) / 1e6 for c in chroms])
    ax.set_xticklabels([c.split("_")[-1] for c in chroms], rotation=90, fontsize=7)
    ax.set_yticks([(yoff[t] + tlen[t] / 2) / 1e6 for t in targets])
    ax.set_yticklabels(targets if len(targets) <= 40 else [""] * len(targets), fontsize=5)
    ax.set_xlim(0, x / 1e6)
    ax.set_ylim(0, y / 1e6)
    ax.set_xlabel("new assembly (Mb)")
    ax.set_ylabel(f"{label} (Mb)")
    ax.set_title(f"{chroms[0].rsplit('_', 1)[0]} vs {label} (minimap2 asm20, blocks >= 10 kb; blue +, red -)")
    fig.tight_layout()
    fig.savefig(f"{outp}.png")
    print(f"{label}: wrote {outp}.chrom.tsv / .target.tsv / .png")


def mode_telo(fa, out):
    W = 20_000
    with open(out, "w") as fo:
        fo.write("chrom\tlength\tleft_copies\tright_copies\tleft_capped\tright_capped\n")
        for n, s in read_fasta(fa):
            if "_chr" not in n:
                continue
            s = s.upper()
            l = s[:W].count("CCCTAAA") + s[:W].count("TTTAGGG")
            r = s[-W:].count("CCCTAAA") + s[-W:].count("TTTAGGG")
            fo.write(f"{n}\t{len(s)}\t{l}\t{r}\t{int(l >= 100)}\t{int(r >= 100)}\n")


if __name__ == "__main__":
    m, a = sys.argv[1], sys.argv[2:]
    {"name": mode_name, "report": mode_report, "telo": mode_telo}[m](*a)
