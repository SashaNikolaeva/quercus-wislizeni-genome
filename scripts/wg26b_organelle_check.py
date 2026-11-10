#!/usr/bin/env python3
"""wg26b: summarise the oatk organelle assemblies of wg26 and apply the fixed
success criteria. Needs minimap2 on PATH (StdEnv/2023 minimap2/2.28). Run by scripts/wg26b_organelle_check.sh.

Per oatk run (-c value) and organelle (pltd, mito), from <dir>/wg26_c<c>.{pltd,mito}.ctg.fasta and .{pltd,mito}.gfa:
  n_seq, circular flags, total_len     FASTA headers `>name<TAB>length=.. wlength=.. nv=.. circular=true|false path=..`
  gfa_segments, gfa_total_len          S lines of the organelle GFA (length = sequence length, else LN:i)
  depth_cv                             along the path of each sequence: path segments s_i (strand stripped), segment
                                       coverage d_i = SC:f tag of the GFA S line (oatk --seq-c-tag default), m_i = number
                                       of times s_i occurs in that path (repeats traversed m times carry m copies),
                                       per-copy depth x_i = d_i / m_i; CV = sqrt(sum w (x - xbar)^2 / sum w) / xbar with
                                       weights w = segment length per occurrence, xbar = weighted mean. NA if <2
                                       occurrences or tags/names missing.
  ref_identity, ref_aligned_frac       minimap2 -cx asm5 of the assembled sequence(s) (query) on NC_062091.1 (pltd) or
                                       CP129458.1 (mito); blocks >= 1 kb; identity = sum col10 / sum col11 over primary
                                       blocks (tp:A:P); aligned fraction = |union of target intervals of primary and
                                       secondary blocks| / reference length (secondaries let both IR copies count).
  LSC, IR, SSC (pltd)                  self-alignment minimap2 -cx asm5 -P seq seq; IR = longest reverse-strand block
                                       >= 10 kb (copies IRa=[qs,qe), IRb=[ts,te)); the two single-copy regions between
                                       the copies on the circle: LSC = longer, SSC = shorter. NA if not found / wraps.
Criteria:
  pltd pass  : n_seq == 1, circular, 150,000 <= len <= 170,000, ref_identity >= 0.995, ref_aligned_frac >= 0.95
  mito pass  : n_seq == 1, circular, 300,000 <= len <= 800,000, depth_cv < 0.3, ref_aligned_frac >= 0.80
  mito choice: -c 250 if it passes, else the lowest passing -c; none -> "mitochondrial graph of X segments / Y kb"
"""
import argparse
import os
import re
import statistics
import subprocess
import sys

REF = {"pltd": ("NC_062091.fa", "NC_062091.1"), "mito": ("CP129458.fa", "CP129458.1")}
REFDIR = "/home/sashanik/scratch/california_oaks/wislizeni_reasm/final/wg09_organelle"


def read_fasta(path):
    seqs, name, hdr = [], None, None
    if not os.path.exists(path):
        return seqs
    buf = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, hdr, "".join(buf)))
                hdr = line[1:].rstrip("\n")
                name = hdr.split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs.append((name, hdr, "".join(buf)))
    return seqs


def read_gfa(path):
    seg = {}
    if not os.path.exists(path):
        return seg
    with open(path) as f:
        for line in f:
            if not line.startswith("S\t"):
                continue
            x = line.rstrip("\n").split("\t")
            tags = dict((t[:2], t[5:]) for t in x[3:] if len(t) > 5)
            L = len(x[2]) if x[2] != "*" else int(tags.get("LN", 0))
            sc = float(tags["SC"]) if "SC" in tags else None
            seg[x[1]] = (L, sc)
    return seg


def merged(iv):
    tot, b0, e0 = 0, None, None
    for b, e in sorted(iv):
        if e0 is None or b > e0:
            if e0 is not None:
                tot += e0 - b0
            b0, e0 = b, e
        elif e > e0:
            e0 = e
    return tot + (e0 - b0 if e0 is not None else 0)


def run_mm(args):
    r = subprocess.run(["minimap2"] + args, capture_output=True, text=True, check=True)
    return [l.split("\t") for l in r.stdout.splitlines() if l]


def depth_cv(hdr, seg):
    m = re.search(r"path=(\S+)", hdr)
    if not m:
        return "NA"
    names = [re.sub(r"[+-]$", "", p) for p in m.group(1).split(",") if p]
    if len(names) < 2 or any(n not in seg or seg[n][1] is None for n in names):
        return "NA"
    mult = {n: names.count(n) for n in names}
    xs = [(seg[n][1] / mult[n], seg[n][0]) for n in names]
    W = sum(w for _, w in xs)
    mu = sum(x * w for x, w in xs) / W
    var = sum(w * (x - mu) ** 2 for x, w in xs) / W
    return "%.4f" % (var ** 0.5 / mu) if mu > 0 else "NA"


def ir_regions(name, seq, tmp):
    fa = os.path.join(tmp, f"self_{name}.fa")
    with open(fa, "w") as f:
        f.write(f">{name}\n{seq}\n")
    best = None
    for x in run_mm(["-cx", "asm5", "-P", fa, fa]):
        if x[4] == "-" and int(x[10]) >= 10000:
            if best is None or int(x[10]) > int(best[10]):
                best = x
    if best is None:
        return "NA", "NA", "NA"
    L = len(seq)
    a = sorted([(int(best[2]), int(best[3])), (int(best[7]), int(best[8]))])
    (a1, b1), (a2, b2) = a
    if b1 > a2:
        return "NA", "NA", "NA"
    r1, r2 = a2 - b1, L - b2 + a1
    return str(max(r1, r2)), str(round(((b1 - a1) + (b2 - a2)) / 2)), str(min(r1, r2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--cvals", required=True, help="comma list of <subdir>:<c>, e.g. c150:150,c250:250,c400:400")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tmp", required=True)
    a = ap.parse_args()
    rows, mito_pass, mito_graph = [], [], {}
    for item in a.cvals.split(","):
        sub, c = item.split(":")
        d = os.path.join(a.workdir, sub)
        for org in ("pltd", "mito"):
            seqs = read_fasta(os.path.join(d, f"wg26_c{c}.{org}.ctg.fasta"))
            seg = read_gfa(os.path.join(d, f"wg26_c{c}.{org}.gfa"))
            circ = []
            for _, hdr, _ in seqs:
                mm = re.search(r"circular=(\S+)", hdr)
                circ.append(mm.group(1) if mm else "NA")
            tot = sum(len(s) for _, _, s in seqs)
            cvs = [depth_cv(h, seg) for _, h, _ in seqs]
            ident, afrac = "NA", "NA"
            lsc = ir = ssc = "NA"
            if seqs:
                q = os.path.join(a.tmp, f"q_{sub}_{org}.fa")
                with open(q, "w") as f:
                    for n, _, s in seqs:
                        f.write(f">{n}\n{s}\n")
                ref_fa = os.path.join(REFDIR, REF[org][0])
                reflen = sum(len(s) for _, _, s in read_fasta(ref_fa))
                paf = [x for x in run_mm(["-cx", "asm5", ref_fa, q]) if int(x[10]) >= 1000]
                prim = [x for x in paf if "tp:A:P" in x[12:]]
                if prim:
                    ident = "%.5f" % (sum(int(x[9]) for x in prim) / sum(int(x[10]) for x in prim))
                afrac = "%.4f" % (merged([(int(x[7]), int(x[8])) for x in paf]) / reflen)
                if org == "pltd" and len(seqs) == 1:
                    lsc, ir, ssc = ir_regions(seqs[0][0], seqs[0][2], a.tmp)
            one = len(seqs) == 1 and circ[0].lower() in ("true", "yes", "1")
            if org == "pltd":
                ok = (one and 150000 <= tot <= 170000 and ident != "NA" and float(ident) >= 0.995
                      and float(afrac) >= 0.95)
            else:
                ok = (one and 300000 <= tot <= 800000 and cvs[0] != "NA" and float(cvs[0]) < 0.3
                      and afrac != "NA" and float(afrac) >= 0.80)
                mito_graph[c] = (len(seg), sum(v[0] for v in seg.values()))
                if ok:
                    mito_pass.append(int(c))
            rows.append([c, org, len(seqs), ",".join(circ) or "NA", tot, len(seg), sum(v[0] for v in seg.values()),
                         ",".join(cvs) or "NA", REF[org][1], ident, afrac, lsc, ir, ssc, "PASS" if ok else "FAIL"])
    with open(a.out, "w") as f:
        f.write("# wg26b: pltd pass = 1 circular seq 150-170 kb, identity >= 0.995 over >= 95% of NC_062091.1; "
                "mito pass = 1 circular seq 300-800 kb, depth_cv < 0.3, >= 80% of CP129458.1 aligned. "
                "depth_cv = length-weighted CV of per-copy segment depth (SC:f / path multiplicity) along the path.\n")
        f.write("c\torganelle\tn_seq\tcircular\ttotal_len\tgfa_segments\tgfa_total_len\tdepth_cv\treference\t"
                "ref_identity\tref_aligned_frac\tLSC\tIR\tSSC\tcriterion\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
        if mito_pass:
            pick = 250 if 250 in mito_pass else min(mito_pass)
            f.write(f"# mito decision: passing -c = {sorted(mito_pass)}; deposit candidate = -c {pick}\n")
        else:
            desc = "; ".join(f"-c {c}: mitochondrial graph of {n} segments / {L / 1000:.1f} kb"
                             for c, (n, L) in mito_graph.items())
            f.write(f"# mito decision: no -c passes; not deposited; report {desc}\n")
    print(open(a.out).read())


if __name__ == "__main__":
    sys.exit(main())
