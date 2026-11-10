#!/usr/bin/env python
"""wg19 helper: (1) sequence identity gate, (2) protein clean-up and counts.

identity <out.tsv> <chrs: all|NN[,NN]> <hap> <helixer_fa> <asm_fa> [<hap> <helixer_fa> <asm_fa> ...]
    Chromosome number = the two digits after '_chr' at the end of the FASTA name (Qwis_hap1_chr01,
    Qwislizeni_h2_chr01, Qwis_hap2_chr01; unplaced sequences are ignored).
    md5_fwd_equal = md5(upper(helixer chr)) == md5(upper(asm chr))
    md5_rc_equal  = md5(upper(helixer chr)) == md5(revcomp(upper(asm chr)))   (A<->T, C<->G, N->N, others kept)
    Gate: every expected row (12 per hap, or the requested chrs) present AND md5_fwd_equal == True, else exit 1
    (the TSV is still written so the failure can be inspected).

proteins <counts.tsv> <smoke 0|1> <hap> <gff3> <raw.fa> <out.fa> [<hap> <gff3> <raw.fa> <out.fa> ...]
    raw = gffread -y output (one record per mRNA). Clean-up: remove ONE terminal '.' or '*' (stop codon);
    a protein with any remaining '.' or '*' has an internal stop and is excluded (count logged).
    Checks: n_raw == number of mRNA features in the GFF3; full mode (smoke 0) also n_mRNA == 36682 (hap1) / 36582 (hap2).
    counts.tsv: hap n_mrna_gff n_raw n_internal_stop n_written
"""
import sys, hashlib

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
EXPECTED_MRNA = {"hap1": 36682, "hap2": 36582}   # cg05 / wg01 GFF3 counts


def fasta(path):
    name, buf = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def chrnum(name):
    i = name.rfind("_chr")
    if i < 0:
        return None
    s = name[i + 4:]
    return s if len(s) == 2 and s.isdigit() else None


def identity(args):
    out, chrs = args[0], args[1]
    want = None if chrs == "all" else set(chrs.split(","))
    rows, ok = [], True
    trip = args[2:]
    for k in range(0, len(trip), 3):
        hap, hfa, afa = trip[k:k + 3]
        hel = {}
        for n, s in fasta(hfa):
            c = chrnum(n)
            if c and (want is None or c in want):
                s = s.upper(); hel[c] = (len(s), hashlib.md5(s.encode()).hexdigest())
        asm = {}
        for n, s in fasta(afa):
            c = chrnum(n)
            if c and (want is None or c in want):
                s = s.upper()
                asm[c] = (len(s), hashlib.md5(s.encode()).hexdigest(),
                          hashlib.md5(s.translate(COMP)[::-1].encode()).hexdigest())
                del s
        exp = sorted(want) if want else ["%02d" % i for i in range(1, 13)]
        for c in exp:
            h = hel.get(c); a = asm.get(c)
            fwd = bool(h and a and h[1] == a[1]); rc = bool(h and a and h[1] == a[2])
            rows.append((hap, "chr" + c, h[0] if h else "NA", a[0] if a else "NA", fwd, rc))
            ok &= fwd
        extra = (set(hel) | set(asm)) - set(exp)
        if extra:
            print("ERROR: %s has unexpected chromosome numbers %s" % (hap, sorted(extra))); ok = False
    with open(out, "w") as o:
        o.write("hap\tchr\tlen_helixer\tlen_asm\tmd5_fwd_equal\tmd5_rc_equal\n")
        for r in rows:
            o.write("\t".join(map(str, r)) + "\n")
    n_fwd = sum(r[4] for r in rows)
    print("identity: %d of %d rows forward-equal" % (n_fwd, len(rows)))
    if not ok:
        print("GATE FAILED: not all chromosomes forward-identical; stopping")
        sys.exit(1)


def proteins(args):
    out, smoke = args[0], args[1] == "1"
    quad = args[2:]
    lines, ok = [], True
    for k in range(0, len(quad), 4):
        hap, gff, raw, fa = quad[k:k + 4]
        n_mrna = 0
        with open(gff) as f:
            for l in f:
                if not l.startswith("#"):
                    t = l.split("\t", 3)
                    if len(t) > 2 and t[2] == "mRNA":
                        n_mrna += 1
        n_raw = n_int = n_out = 0
        with open(fa, "w") as o:
            for n, s in fasta(raw):
                n_raw += 1
                if s and s[-1] in ".*":
                    s = s[:-1]
                if "." in s or "*" in s:
                    n_int += 1; continue
                o.write(">%s\n%s\n" % (n, s)); n_out += 1
        lines.append((hap, n_mrna, n_raw, n_int, n_out))
        print("%s: mRNA %d, raw proteins %d, internal stop %d, written %d" % (hap, n_mrna, n_raw, n_int, n_out))
        if n_raw != n_mrna:
            print("ERROR: %s raw protein count != mRNA count" % hap); ok = False
        if not smoke and n_mrna != EXPECTED_MRNA[hap]:
            print("ERROR: %s mRNA count %d != expected %d" % (hap, n_mrna, EXPECTED_MRNA[hap])); ok = False
    with open(out, "w") as o:
        o.write("hap\tn_mrna_gff\tn_raw\tn_internal_stop\tn_written\n")
        for r in lines:
            o.write("\t".join(map(str, r)) + "\n")
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    {"identity": identity, "proteins": proteins}[sys.argv[1]](sys.argv[2:])
