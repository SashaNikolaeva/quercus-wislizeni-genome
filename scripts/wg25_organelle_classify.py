#!/usr/bin/env python3
"""wg25: classify unplaced contigs as organelle-derived or nuclear. CANDIDATE LIST ONLY;
the assemblies are not modified.

Called by scripts/wg25_organelle_classify.sh. Inputs per haplotype: a list of unplaced contigs with lengths, a PAF of
`minimap2 -cx asm20 --secondary=no` of those contigs against NC_062091.1 (Q. wislizeni plastome) + CP129458.1
(Q. variabilis mitogenome) concatenated (plastid and mito compete for each region), and a `samtools coverage` table
of dip.bam WITHOUT -q (default flag filter: primary, non-duplicate; MAPQ >= 0).

Definitions:
  block           = PAF line with alignment block length (col11) >= 1000 bp
  identity        = PAF col10 / col11 (matching bases / block length incl. gaps; needs -c for base-level col10)
  plastid_cov     = |union of query intervals of blocks on NC_062091.1 with identity >= 0.99| / contig length
                    (same-species reference)
  mito_cov        = |union of query intervals of blocks on CP129458.1 with identity >= 0.95| / contig length
                    (cross-species reference, section Cerris)
  org_cov         = |union of the plastid and mito interval sets above| / contig length
  depth_q0        = samtools coverage meandepth (MAPQ >= 0)
  depth_ratio     = depth_q0 / median(depth_q0 of that haplotype's 12 chromosomes)
  class:
    candidate_remove         org_cov >= 0.90 and depth_ratio >= 2
    organelle_lowdepth_flag org_cov >= 0.90 and depth_ratio < 2 
    mixed_keep               0.50 <= org_cov < 0.90                 (possible NUPT/NUMT-bearing; kept)
    nuclear                  org_cov < 0.50
Caveat: on a diploid reference with many near-identical organelle-derived copies, reads are split across the
redundant copies (and get MAPQ 0), so depth per copy is biased LOW and depth_ratio underestimates organelle origin.
That is why organelle_lowdepth_flag contigs are kept and listed instead of being silently kept.
depth_q20 (from the earlier coverage_q20.tsv) is carried along only to document why MAPQ >= 20 depth was not used.
"""
import argparse
import statistics
import sys

PLASTID, MITO = "NC_062091.1", "CP129458.1"
MIN_BLOCK, ID_PT, ID_MT = 1000, 0.99, 0.95
ORG_HI, ORG_LO, DEPTH_RATIO = 0.90, 0.50, 2.0


def merged_len(iv):
    tot, cur_b, cur_e = 0, None, None
    for b, e in sorted(iv):
        if cur_e is None or b > cur_e:
            if cur_e is not None:
                tot += cur_e - cur_b
            cur_b, cur_e = b, e
        elif e > cur_e:
            cur_e = e
    if cur_e is not None:
        tot += cur_e - cur_b
    return tot


def read_cov(path):
    d = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            x = line.rstrip("\n").split("\t")
            d[x[0]] = float(x[6])
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hap", action="append", required=True,
                    help="hapN:contig_lengths.tsv:alignments.paf (repeat per haplotype)")
    ap.add_argument("--cov_q0", required=True)
    ap.add_argument("--cov_q20", required=True)
    ap.add_argument("--wg09", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()

    q0, q20 = read_cov(a.cov_q0), read_cov(a.cov_q20)
    wg09 = {}
    with open(a.wg09) as f:
        next(f)
        for line in f:
            x = line.rstrip("\n").split("\t")
            wg09[x[1]] = x[7]

    rows, fails = [], []
    for spec in a.hap:
        hap, lens_path, paf = spec.split(":")
        lens = {}
        with open(lens_path) as f:
            for line in f:
                c, L = line.split()[:2]
                if c in lens:
                    fails.append(f"{hap}: contig {c} listed twice")
                lens[c] = int(L)
        iv = {c: {"pt": [], "mt": []} for c in lens}
        with open(paf) as f:
            for line in f:
                x = line.rstrip("\n").split("\t")
                q, qs, qe, t, nm, bl = x[0], int(x[2]), int(x[3]), x[5], int(x[9]), int(x[10])
                if bl < MIN_BLOCK or q not in iv:
                    continue
                ident = nm / bl
                if t == PLASTID and ident >= ID_PT:
                    iv[q]["pt"].append((qs, qe))
                elif t == MITO and ident >= ID_MT:
                    iv[q]["mt"].append((qs, qe))
        chr_depth = [v for k, v in q0.items() if k.startswith(f"Qwis_{hap}_chr")]
        if len(chr_depth) != 12:
            fails.append(f"{hap}: {len(chr_depth)} chromosome depth rows, expected 12")
        med = statistics.median(chr_depth) if chr_depth else float("nan")
        print(f"{hap}: chromosome median MAPQ>=0 meandepth = {med:.2f}x over {len(chr_depth)} chromosomes")
        if not 15 <= med <= 30:  # smoke: chromosome depth is from the first 1 Mb only, so warn instead
            (print if a.smoke else fails.append)(f"{hap}: chromosome median depth {med:.2f}x outside 15-30x")
        for c, L in lens.items():
            pt, mt = merged_len(iv[c]["pt"]), merged_len(iv[c]["mt"])
            org = merged_len(iv[c]["pt"] + iv[c]["mt"])
            pc, mc, oc = pt / L, mt / L, org / L
            if c not in q0:
                fails.append(f"{hap}: no MAPQ>=0 depth for {c}")
                d0 = float("nan")
            else:
                d0 = q0[c]
            dr = d0 / med if med > 0 else float("nan")
            if oc >= ORG_HI:
                cls = "candidate_remove" if dr >= DEPTH_RATIO else "organelle_lowdepth_flag"
            elif oc >= ORG_LO:
                cls = "mixed_keep"
            else:
                cls = "nuclear"
            rows.append((hap, c, L, wg09.get(c, "NA"), pc, mc, oc, d0, q20.get(c, float("nan")), dr, cls))

    seen = [r[1] for r in rows]
    if len(seen) != len(set(seen)):
        fails.append("some contig classified more than once")
    o = a.outdir
    with open(f"{o}/organelle_contig_classes.tsv", "w") as f:
        f.write("# wg25. CANDIDATE LIST ONLY, assemblies unchanged. block>=1kb; identity=PAF col10/col11; "
                "plastid_cov: NC_062091.1 id>=0.99; mito_cov: CP129458.1 id>=0.95; org_cov=union; "
                "depth_ratio=depth_q0/median chr depth_q0 (MAPQ>=0). candidate_remove: org_cov>=0.90 & ratio>=2; "
                "organelle_lowdepth_flag: org_cov>=0.90 & ratio<2; mixed_keep: 0.50<=org_cov<0.90; nuclear: <0.50. "
                "Depth is split across redundant copies, which biases depth_ratio LOW: lowdepth_flag contigs are "
                "kept and listed.\n")
        f.write("hap\tcontig\tlength\twg09_call\tplastid_cov\tmito_cov\torg_cov\tdepth_q0\tdepth_q20\tdepth_ratio\tclass\n")
        for r in rows:
            f.write("%s\t%s\t%d\t%s\t%.4f\t%.4f\t%.4f\t%.3f\t%.3f\t%.3f\t%s\n" % r)
    classes = ["candidate_remove", "organelle_lowdepth_flag", "mixed_keep", "nuclear"]
    haps = [s.split(":")[0] for s in a.hap]
    with open(f"{o}/organelle_classes_summary.tsv", "w") as f:
        f.write("hap\tclass\tcontigs\tMb\n")
        for h in haps:
            for cl in classes + ["all"]:
                sel = [r for r in rows if r[0] == h and (cl == "all" or r[10] == cl)]
                f.write("%s\t%s\t%d\t%.3f\n" % (h, cl, len(sel), sum(r[2] for r in sel) / 1e6))
    for h in haps:
        with open(f"{o}/organelle_removal_candidates_{h}.txt", "w") as f:
            for r in rows:
                if r[0] == h and r[10] == "candidate_remove":
                    f.write(r[1] + "\n")

    # success criteria: 898 contigs (648 + 250) classified once; lengths 34.81 / 13.25 Mb (+-0.01)
    expect = {"hap1": (648, 34.81), "hap2": (250, 13.25)}
    for h in haps:
        sel = [r for r in rows if r[0] == h]
        mb = sum(r[2] for r in sel) / 1e6
        print(f"{h}: {len(sel)} contigs classified, {mb:.3f} Mb")
        if not a.smoke and (len(sel) != expect[h][0] or abs(mb - expect[h][1]) > 0.01):
            fails.append(f"{h}: {len(sel)} contigs / {mb:.3f} Mb, expected {expect[h][0]} / {expect[h][1]}")
    if a.smoke and len(rows) != 20:
        fails.append(f"smoke: {len(rows)} contigs classified, expected 20")
    print(open(f"{o}/organelle_classes_summary.tsv").read())
    if fails:
        print("SUCCESS CRITERIA FAILED:\n  " + "\n  ".join(fails), file=sys.stderr)
        sys.exit(1)
    print("all success criteria passed")


if __name__ == "__main__":
    main()
