#!/usr/bin/env python
"""wg21b: summarise the wg21 InterProScan 5.73-104.0 TSVs per gene.

Usage: wg21b_interproscan_summary.py <out_dir> <smoke 0|1> <gene_te_overlap.tsv>
           hap1 <gff3> <prot.fa> <tsv_dir> <n_chunks> [hap2 <gff3> <prot.fa> <tsv_dir> <n_chunks>]
InterProScan TSV columns (1-based): 1 protein, 4 analysis (member DB), 5 signature, 6 signature description,
  9 score (e-value for Pfam/PANTHER), 12 InterPro accession, 13 InterPro description, 14 GO terms ('|'-separated).
protein -> gene: GFF3 mRNA Parent. Sets: hap1_all, hap1_nonTE (gene_te_overlap.tsv frac_all <= 0.5; TE-derived =
  frac_all > 0.5), hap2_all.
Denominator per set: all genes of the set in the GFF3 (genes whose protein was excluded for an internal stop count as
  'no hit'); in smoke: the genes whose protein is in the chunk .ids file.
  pct_interpro = 100 * genes with >= 1 row with col12 != '-' / genes
  pct_pfam     = 100 * genes with >= 1 Pfam row / genes
  pct_go       = 100 * genes with >= 1 row with a GO term in col14 / genes
  pct_any_db   = 100 * genes with >= 1 row from any of the 6 member DBs / genes
Products (interproscan_gene_products_hap<h>.tsv): per gene, the Pfam or PANTHER row with the lowest e-value
  (ties: Pfam before PANTHER, then signature accession); description = col6, or col13 if col6 is '-'/empty;
  go_terms = union of col14 over all rows of the gene's protein(s). Genes without a Pfam/PANTHER row are listed with '-'.
Checks (exit 1): n_chunks TSVs + .ids present per hap, each TSV non-empty; distinct proteins in the TSVs <= input
  proteins; every TSV protein is in the input FASTA.
"""
import sys, os, glob, math

OUTD, SMOKE, TEO = sys.argv[1:4]
smoke = SMOKE == "1"
rest = sys.argv[4:]
haps = [rest[i:i + 5] for i in range(0, len(rest), 5)]
ok = True

te = {}
with open(TEO) as f:
    h = f.readline().rstrip("\n").split("\t"); gi, fi = h.index("gene"), h.index("frac_all")
    for l in f:
        t = l.rstrip("\n").split("\t"); te[t[gi]] = int(t[h.index("te_bp_all")]) / int(t[h.index("cds_bp")]) > 0.5  # exact, not rounded frac_all (fix 2026-11-05)


def attr(a, key):
    for kv in a.strip().split(";"):
        if kv.startswith(key + "="):
            return kv[len(key) + 1:]
    return None


def evalue(x):
    try:
        return float(x)
    except ValueError:
        return math.inf


summary_rows = []
for hap, gff, prot, tdir, nch in haps:
    nch = int(nch)
    genes, m2g = [], {}
    with open(gff) as f:
        for l in f:
            if l.startswith("#"): continue
            t = l.rstrip("\n").split("\t")
            if len(t) < 9: continue
            if t[2] == "gene": genes.append(attr(t[8], "ID"))
            elif t[2] == "mRNA": m2g[attr(t[8], "ID")] = attr(t[8], "Parent")
    inprot = set(l[1:].split()[0] for l in open(prot) if l.startswith(">"))
    tsvs = sorted(glob.glob(os.path.join(tdir, "chunk_*.tsv")))
    ids = sorted(glob.glob(os.path.join(tdir, "chunk_*.ids")))
    if len(tsvs) != nch or len(ids) != nch:
        print("ERROR %s: %d TSVs / %d ids files, expected %d" % (hap, len(tsvs), len(ids), nch)); ok = False
    chunk_prot = set()
    for p in ids:
        chunk_prot |= set(x.strip() for x in open(p) if x.strip())
    hit = {"interpro": set(), "pfam": set(), "go": set(), "any": set()}
    best, gos, seenp = {}, {}, set()
    for p in tsvs:
        if os.path.getsize(p) == 0:
            print("ERROR %s: empty %s" % (hap, p)); ok = False
        for l in open(p):
            t = l.rstrip("\n").split("\t")
            if len(t) < 13: continue
            pr = t[0]; seenp.add(pr)
            g = m2g.get(pr)
            if g is None:
                print("ERROR %s: protein %s not an mRNA of the GFF3" % (hap, pr)); ok = False; continue
            db = t[3]; go = t[13] if len(t) > 13 else "-"
            hit["any"].add(g)
            if t[11] not in ("-", ""): hit["interpro"].add(g)
            if db == "Pfam": hit["pfam"].add(g)
            if go not in ("-", ""):
                hit["go"].add(g); gos.setdefault(g, set()).update(go.split("|"))
            if db in ("Pfam", "PANTHER"):
                desc = t[5] if t[5] not in ("-", "") else (t[12] if len(t) > 12 else "-")
                key = (evalue(t[8]), 0 if db == "Pfam" else 1, t[4])
                if g not in best or key < best[g][0]:
                    best[g] = (key, pr, db, t[4], desc, t[8], t[11], t[12])
    bad = seenp - inprot
    if bad:
        print("ERROR %s: %d TSV proteins not in input FASTA" % (hap, len(bad))); ok = False
    if len(seenp) > len(inprot):
        print("ERROR %s: distinct TSV proteins %d > input %d" % (hap, len(seenp), len(inprot))); ok = False
    print("%s: input proteins %d, chunk proteins %d, proteins with >=1 hit %d" % (hap, len(inprot), len(chunk_prot), len(seenp)))
    if smoke:
        genes = sorted({m2g[p] for p in chunk_prot if p in m2g})
    sets = [(hap + "_all", genes)]
    if hap == "hap1":
        sets.append(("hap1_nonTE", [g for g in genes if g in te and not te[g]]))
    for name, gs in sets:
        n = len(gs); gs = set(gs)
        row = [name, n] + ["%.2f" % (100.0 * len(gs & hit[k]) / n if n else float("nan"))
                           for k in ("interpro", "pfam", "go", "any")]
        summary_rows.append(row)
    with open(os.path.join(OUTD, "interproscan_gene_products_%s.tsv" % hap), "w") as o:
        o.write("gene\tmrna\tsource_db\tsignature\tdescription\tevalue\tinterpro_acc\tinterpro_desc\tgo_terms\tTE_derived\n")
        for g in genes:
            ted = ("yes" if te.get(g) else "no") if hap == "hap1" else "NA"
            go = "|".join(sorted(gos.get(g, []))) or "-"
            if g in best:
                _, pr, db, sig, desc, ev, ipr, iprd = best[g]
                o.write("\t".join([g, pr, db, sig, desc, ev, ipr, iprd, go, ted]) + "\n")
            else:
                o.write("\t".join([g, "-", "-", "-", "-", "-", "-", "-", go, ted]) + "\n")

out = os.path.join(OUTD, "interproscan_summary.tsv" + ("" if ok else ".FAILED"))
with open(out, "w") as o:
    o.write("set\tn_genes\tpct_interpro\tpct_pfam\tpct_go\tpct_any_memberdb\n")
    for r in summary_rows:
        o.write("\t".join(map(str, r)) + "\n")
print(open(out).read())
if not ok:
    sys.exit(1)
