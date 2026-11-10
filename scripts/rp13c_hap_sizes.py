#!/usr/bin/env python
"""rp13c: where does the hap1 vs hap2 size difference of the new wislizeni assembly come from?

1. Chromosomes: per lobata chromosome, hap1 vs hap2 length (rp13 names tables).
2. Unplaced sequences, each classified by
   - hifiasm read depth (rd:i: of its contig in wis.hic.<hap>.p_ctg.gfa; the unplaced YaHS scaffold -> contig via the AGP;
     multi-contig scaffolds take the length-weighted mean),
   - GC, and
   - best minimap2 -x asm5 coverage against its OWN haplotype's 12 chromosomes and the OTHER haplotype's 12 chromosomes
     (fraction of the sequence covered by alignments >= 1 kb with identity >= 0.95; identity = matches / block length).
   Classes, applied in this order:
     organelle_or_highcopy  rd >= 3 x median chromosome-contig rd
     dup_of_own_chrom       own coverage >= 0.8                       (false duplication / retained haplotig)
     other_hap_only         other-hap coverage >= 0.8, own < 0.8      (this haplotype's allele, left unplaced)
     partial                own or other coverage 0.2-0.8
     unique_or_repeat       neither >= 0.2                            (satellite, rDNA, or sequence absent from chromosomes)
Usage: rp13c_hap_sizes.py <final_dir> <hifiasm_dir> <yahs_dir> <out_prefix>
   expects <final_dir>/Qwis_<hap>.{names.tsv,fa}, <final_dir>/rp13c_<hap>_unpl_vs_{own,other}.paf
"""
import sys
import re
import statistics
from collections import defaultdict

fin, hdir, ydir, outp = sys.argv[1:5]
HAPS = ("hap1", "hap2")


def names(h):
    rows = [l.rstrip("\n").split("\t") for l in open(f"{fin}/Qwis_{h}.names.tsv")][1:]
    return {r[0]: dict(yahs=r[1], length=int(r[2]), lob=r[3]) for r in rows}


def gfa_rd(h):
    rd, ln = {}, {}
    for l in open(f"{hdir}/wis.hic.{h}.p_ctg.gfa"):
        if l.startswith("S\t"):
            f = l.rstrip("\n").split("\t")
            m = re.search(r"rd:i:(\d+)", l)
            rd[f[1]] = int(m.group(1)) if m else None
            ln[f[1]] = len(f[2])
    return rd, ln


def agp_contigs(h):
    """yahs scaffold -> [(contig, length)]"""
    out = defaultdict(list)
    for l in open(f"{ydir}/{h}/{h}_yahs_scaffolds_final.agp"):
        f = l.rstrip("\n").split("\t")
        if len(f) > 5 and f[4] == "W":
            out[f[0]].append((f[5], int(f[7]) - int(f[6]) + 1))
    return out


def gc_by_seq(fa, keep):
    gc, name, buf = {}, None, []
    def flush():
        if name in keep:
            s = "".join(buf).upper()
            acgt = sum(s.count(c) for c in "ACGT")
            gc[name] = (s.count("G") + s.count("C")) / acgt if acgt else float("nan")
    with open(fa) as fh:
        for l in fh:
            if l.startswith(">"):
                flush()
                name, buf = l[1:].split()[0], []
            else:
                buf.append(l.strip())
    flush()
    return gc


def coverage(paf):
    iv = defaultdict(list)
    for l in open(paf):
        f = l.split("\t")
        blk, mat = int(f[10]), int(f[9])
        if blk >= 1000 and mat / blk >= 0.95:
            iv[f[0]].append((int(f[2]), int(f[3])))
    cov = {}
    for q, ivs in iv.items():
        ivs.sort()
        tot, cs, ce = 0, -1, -1
        for s, e in ivs:
            if s > ce:
                tot += max(0, ce - cs)
                cs, ce = s, e
            else:
                ce = max(ce, e)
        tot += max(0, ce - cs)
        cov[q] = tot
    return cov


N = {h: names(h) for h in HAPS}
# 1. chromosomes
lob2 = {h: {v["lob"]: (k, v["length"]) for k, v in N[h].items() if "_chr" in k} for h in HAPS}
with open(f"{outp}_chromosomes.tsv", "w") as fo:
    fo.write("chrom\tlobata\thap1_Mb\thap2_Mb\thap1_minus_hap2_Mb\n")
    t1 = t2 = 0
    for lob in sorted(lob2["hap1"], key=lambda x: lob2["hap1"][x][0]):
        a, b = lob2["hap1"][lob][1], lob2["hap2"][lob][1]
        t1, t2 = t1 + a, t2 + b
        fo.write(f"{lob2['hap1'][lob][0].split('_')[-1]}\t{lob}\t{a/1e6:.2f}\t{b/1e6:.2f}\t{(a-b)/1e6:+.2f}\n")
    fo.write(f"total\t-\t{t1/1e6:.2f}\t{t2/1e6:.2f}\t{(t1-t2)/1e6:+.2f}\n")

# 2. unplaced
rows = []
for h in HAPS:
    other = "hap2" if h == "hap1" else "hap1"
    rd, cl = gfa_rd(h)
    agp = agp_contigs(h)
    chrom_ctgs = [c for k, v in N[h].items() if "_chr" in k for c, _ in agp[v["yahs"]]]
    med = statistics.median(rd[c] for c in chrom_ctgs if rd.get(c) is not None)
    unpl = {k: v for k, v in N[h].items() if "_unplaced_" in k}
    gc = gc_by_seq(f"{fin}/Qwis_{h}.fa", set(unpl))
    own = coverage(f"{fin}/rp13c_{h}_unpl_vs_own.paf")
    oth = coverage(f"{fin}/rp13c_{h}_unpl_vs_other.paf")
    for k, v in unpl.items():
        parts = agp[v["yahs"]]
        w = sum(L for c, L in parts if rd.get(c) is not None)
        d = sum(rd[c] * L for c, L in parts if rd.get(c) is not None) / w if w else float("nan")
        L = v["length"]
        co, ct = own.get(k, 0) / L, oth.get(k, 0) / L
        if d >= 3 * med:
            cls = "organelle_or_highcopy"
        elif co >= 0.8:
            cls = "dup_of_own_chrom"
        elif ct >= 0.8:
            cls = "other_hap_only"
        elif max(co, ct) >= 0.2:
            cls = "partial"
        else:
            cls = "unique_or_repeat"
        rows.append((h, k, L, len(parts), d, d / med, gc.get(k, float("nan")), co, ct, cls))
    print(f"{h}: median chromosome-contig rd {med}")
with open(f"{outp}_unplaced.tsv", "w") as fo:
    fo.write("hap\tseq\tlength\tn_contigs\trd\trd_over_median\tGC\tcov_own_chroms\tcov_other_hap_chroms\tclass\n")
    for r in rows:
        fo.write("\t".join([r[0], r[1], str(r[2]), str(r[3]), f"{r[4]:.1f}", f"{r[5]:.2f}", f"{r[6]:.3f}",
                            f"{r[7]:.3f}", f"{r[8]:.3f}", r[9]]) + "\n")
summ = defaultdict(lambda: [0, 0])
for r in rows:
    summ[(r[0], r[9])][0] += 1
    summ[(r[0], r[9])][1] += r[2]
with open(f"{outp}_unplaced_summary.tsv", "w") as fo:
    fo.write("hap\tclass\tn\tMb\n")
    for h in HAPS:
        for c in ("organelle_or_highcopy", "dup_of_own_chrom", "other_hap_only", "partial", "unique_or_repeat"):
            n, b = summ[(h, c)]
            fo.write(f"{h}\t{c}\t{n}\t{b/1e6:.2f}\n")
        fo.write(f"{h}\ttotal\t{sum(summ[(h, c)][0] for c in [k[1] for k in summ if k[0] == h])}\t"
                 f"{sum(summ[(h, c)][1] for c in [k[1] for k in summ if k[0] == h])/1e6:.2f}\n")
print(open(f"{outp}_chromosomes.tsv").read())
print(open(f"{outp}_unplaced_summary.tsv").read())
