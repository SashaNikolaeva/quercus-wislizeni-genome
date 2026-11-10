#!/usr/bin/env python3
"""wg22b: tRNA and rRNA summary tables for the wislizeni haplotypes.

Called by wg22b_ncrna_summary.sh, which prepares in --work, per haplotype h:
  hap{h}.trnascan.out           concatenated per-chromosome tRNAscan-SE 2.0.12 --detail output (12 chromosomes)
  hc_hap{h}/hap{h}.out          EukHighConfidenceFilter (tag mode) run once on the concatenated haplotype output
  hap{h}.hc_perchr.out          the per-chromosome filter outputs from wg22, concatenated (check column only)
  hap{h}_NC_062091.paf, hap{h}_CP129458.paf   minimap2 -cx asm20 -N 1000 -p 0.1, organelle query -> 12 chromosomes
  hap{h}.barrnap.gff3           barrnap 0.9 --kingdom euk on the whole assembly (chromosomes + unplaced)

Definitions (tRNA, per haplotype, 12 chromosomes only):
  trnascan_total    = data rows of the raw tRNAscan-SE output (every prediction, incl. pseudogenes).
  pseudogenes       = raw rows whose Note contains "pseudo".
  SeC / suppressor / undetermined = raw rows with Type SeC / Sup / Undet.
  hc_total          = rows of the haplotype-level filtered output whose Note contains "high confidence set".
  hc_perchr_filter  = same, summed over the 12 per-chromosome filter runs (check only; the filter's anticodon-count
                      cutoff of 40 and isotype score cutoffs depend on its input, so the haplotype-level set is reported).
  organelle blocks  = PAF rows with col11 (alignment block length) >= 500; target intervals [col8, col9) merged per
                      chromosome, separately for plastid (NC_062091.1) and mito (CP129458.1).
  plastid-like      = overlap(tRNA, merged plastid blocks) >= 0.5 x tRNA length; tRNA interval = [min(Begin,End)-1,
                      max(Begin,End)) in 0-based half-open coordinates. mito-like: the same with mito blocks.
  organelle-like    = plastid-like OR mito-like;  hc_excl_organelle = hc_total - hc_organelle_like.
  standard isotypes = Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val (20); Type iMet and
                      fMet count as Met, Ile2 as Ile. n_std_isotypes_hc = how many of the 20 occur in the hc set.
rRNA (barrnap 0.9): unit from Name (18S_rRNA, 5_8S_rRNA, 28S_rRNA, 5S_rRNA); partial = attributes contain "(partial)"
  (barrnap writes it into product, for hits < 0.8 x expected length); otherwise full. scope = chromosomes (seqid ends in
  _chrNN), unplaced (contains _unplaced_), all = both. k-mer estimate = haploid column of rdna_kmer_copy_number_summary.tsv
  (18S -> 18S, 28S -> 25S (plant LSU), 5S -> 5S, 5.8S -> NA, no k-mer query for it).
Success checks: hc_total per haplotype in the hundreds to low thousands,
  read as 100 <= hc_total <= 5000; all 20 standard isotypes in the hc set; >= 1 full 18S, 28S and 5S in scope "all".
  Results are written first; the script then exits 2 if any check fails (not in smoke mode).
"""
import argparse
import os
import re
import sys
from collections import Counter, defaultdict

STD = ["Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
       "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"]
ALIAS = {"iMet": "Met", "fMet": "Met", "Ile2": "Ile"}
ORGS = {"plastid": "NC_062091", "mito": "CP129458"}
RRNA = {"18S_rRNA": "18S", "5_8S_rRNA": "5.8S", "28S_rRNA": "28S", "5S_rRNA": "5S"}
KMER_UNIT = {"18S": "18S", "28S": "25S", "5S": "5S", "5.8S": None}


def read_trnascan(path):
    """Return list of dicts (seq, num, s0, e, type, anti, note) from a (possibly concatenated) tRNAscan-SE .out."""
    rows, note_i = [], None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            cols = [c.strip() for c in line.split("\t")]
            if cols[0].startswith("Sequence") or cols[0].startswith("-----"):
                continue
            if cols[0].startswith("Name"):
                note_i = cols.index("Note") if "Note" in cols else None
                continue
            b, e = int(cols[2]), int(cols[3])
            note = cols[note_i] if note_i is not None and note_i < len(cols) else ""
            rows.append(dict(seq=cols[0], num=cols[1], s0=min(b, e) - 1, e=max(b, e),
                             type=cols[4], anti=cols[5], note=note))
    return rows


def read_blocks(paf, minlen=500):
    iv = defaultdict(list)
    with open(paf) as f:
        for line in f:
            c = line.split("\t")
            if len(c) < 12 or int(c[10]) < minlen:
                continue
            iv[c[5]].append((int(c[7]), int(c[8])))
    merged = {}
    for t, lst in iv.items():
        lst.sort()
        m = [list(lst[0])]
        for s, e in lst[1:]:
            if s <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        merged[t] = m
    return merged


def overlap(blocks, seq, s, e):
    tot = 0
    for bs, be in blocks.get(seq, []):
        if bs >= e:
            break
        if be > s:
            tot += min(be, e) - max(bs, s)
    return tot


def std_iso(t):
    t = ALIAS.get(t, t)
    return t if t in STD else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    ap.add_argument("--haps", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--rdna", required=True)
    ap.add_argument("--fai-dir", required=True)
    a = ap.parse_args()
    smoke = os.path.basename(os.path.normpath(a.work)) == "smoke"
    W = a.work

    kmer = {}
    with open(a.rdna) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            d = dict(zip(hdr, line.rstrip("\n").split("\t")))
            kmer[d["unit"]] = d["haploid"]

    trna_rows, rrna_rows, checks = [], [], []
    for h in a.haps:
        raw = read_trnascan(f"{W}/hap{h}.trnascan.out")
        hc = [r for r in read_trnascan(f"{W}/hc_hap{h}/hap{h}.out") if "high confidence set" in r["note"]]
        hc_pc = [r for r in read_trnascan(f"{W}/hap{h}.hc_perchr.out") if "high confidence set" in r["note"]]
        with open(f"{W}/hap{h}.chrlist") as f:
            nchr = sum(1 for x in f if x.strip())
        blocks = {k: read_blocks(f"{W}/hap{h}_{v}.paf") for k, v in ORGS.items()}
        hc_keys = {(r["seq"], r["num"]) for r in hc}
        flag = {}
        with open(f"{W}/hap{h}.trna_flags.tsv", "w") as fo:
            fo.write("seq\ttrna_num\tstart0\tend\ttype\tanticodon\thigh_conf\tplastid_frac\tmito_frac\tplastid_like\tmito_like\tnote\n")
            for r in raw:
                L = r["e"] - r["s0"]
                fr = {k: overlap(b, r["seq"], r["s0"], r["e"]) / L for k, b in blocks.items()}
                fl = {k: fr[k] >= 0.5 for k in fr}
                flag[(r["seq"], r["num"])] = fl
                fo.write(f"{r['seq']}\t{r['num']}\t{r['s0']}\t{r['e']}\t{r['type']}\t{r['anti']}\t"
                         f"{int((r['seq'], r['num']) in hc_keys)}\t{fr['plastid']:.3f}\t{fr['mito']:.3f}\t"
                         f"{int(fl['plastid'])}\t{int(fl['mito'])}\t{r['note']}\n")
        types = Counter(r["type"] for r in raw)
        pl = sum(flag[(r["seq"], r["num"])]["plastid"] for r in hc)
        mi = sum(flag[(r["seq"], r["num"])]["mito"] for r in hc)
        org = sum(flag[(r["seq"], r["num"])]["plastid"] or flag[(r["seq"], r["num"])]["mito"] for r in hc)
        iso_hc = {std_iso(r["type"]) for r in hc} - {None}
        missing = [x for x in STD if x not in iso_hc]
        bp = {k: sum(e - s for v in b.values() for s, e in v) for k, b in blocks.items()}
        trna_rows.append([f"hap{h}", nchr, len(raw), len(hc), len(hc_pc),
                          sum("pseudo" in r["note"] for r in raw), types.get("SeC", 0), types.get("Sup", 0),
                          types.get("Undet", 0), pl, mi, org, len(hc) - org, len(iso_hc),
                          ",".join(missing) if missing else "none", bp["plastid"], bp["mito"]])
        checks.append([f"hap{h}", "hc_total_100_to_5000", len(hc), "PASS" if 100 <= len(hc) <= 5000 else "FAIL"])
        checks.append([f"hap{h}", "all_20_std_isotypes_in_hc", len(iso_hc), "PASS" if not missing else "FAIL"])

        # rRNA
        cnt = Counter()
        unknown = Counter()
        with open(f"{W}/hap{h}.barrnap.gff3") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                c = line.rstrip("\n").split("\t")
                attrs = dict(x.split("=", 1) for x in c[8].split(";") if "=" in x)
                unit = RRNA.get(attrs.get("Name", ""))
                if unit is None:
                    unknown[attrs.get("Name", "")] += 1
                    continue
                if re.search(r"_chr\d\d$", c[0]):
                    scope = "chromosomes"
                elif "_unplaced_" in c[0]:
                    scope = "unplaced"
                else:
                    sys.exit(f"ERROR: seqid {c[0]} is neither chromosome nor unplaced")
                status = "partial" if "(partial)" in c[8] else "full"
                cnt[(scope, unit, status)] += 1
                cnt[("all", unit, status)] += 1
        if unknown:
            print(f"hap{h}: barrnap Names not in the euk set (ignored): {dict(unknown)}")
        for scope in ("chromosomes", "unplaced", "all"):
            for unit in ("18S", "5.8S", "28S", "5S"):
                for status in ("full", "partial"):
                    ku = KMER_UNIT[unit]
                    rrna_rows.append([f"hap{h}", scope, unit, status, cnt[(scope, unit, status)],
                                      kmer.get(ku, "NA") if ku else "NA", ku or "NA",
                                      "assembled copies are collapsed arrays, not copy number"])
        for unit in ("18S", "28S", "5S"):
            n = cnt[("all", unit, "full")]
            checks.append([f"hap{h}", f"full_{unit}_ge1_assembly", n, "PASS" if n >= 1 else "FAIL"])

    os.makedirs(a.out, exist_ok=True)
    with open(f"{a.out}/trna_summary.tsv", "w") as fo:
        fo.write("# tRNAscan-SE 2.0.12 -E (12 chromosomes only); hc = EukHighConfidenceFilter run per haplotype; "
                 "organelle-like = >= 50% of tRNA length in minimap2 asm20 blocks >= 500 bp of NC_062091.1 "
                 "(plastid) or CP129458.1 (mito); definitions in scripts/wg22b_ncrna_summary.py\n")
        fo.write("hap\tn_chromosomes\ttrnascan_total\thc_total\thc_perchr_filter_sum\tpseudogenes\tSeC\tsuppressor\t"
                 "undetermined\thc_plastid_like\thc_mito_like\thc_organelle_like\thc_excl_organelle\t"
                 "n_std_isotypes_hc\tmissing_std_isotypes_hc\tplastid_block_bp\tmito_block_bp\n")
        for r in trna_rows:
            fo.write("\t".join(map(str, r)) + "\n")
    with open(f"{a.out}/rrna_summary.tsv", "w") as fo:
        fo.write("# barrnap 0.9 --kingdom euk, whole assembly; partial = < 0.8 x expected length; "
                 "kmer_haploid_estimate from results/rdna/rdna_kmer_copy_number_summary.tsv (28S <- 25S k-mer query)\n")
        fo.write("hap\tscope\trRNA\tstatus\tcount\tkmer_haploid_estimate\tkmer_unit\tnote\n")
        for r in rrna_rows:
            fo.write("\t".join(map(str, r)) + "\n")
    nfail = sum(c[3] == "FAIL" for c in checks)
    with open(f"{a.out}/ncrna_checks.tsv", "w") as fo:
        fo.write("hap\tcheck\tvalue\tverdict\n")
        for c in checks:
            fo.write("\t".join(map(str, c)) + "\n")
    for fn in ("trna_summary.tsv", "rrna_summary.tsv", "ncrna_checks.tsv"):
        print(f"== {fn} ==")
        with open(f"{a.out}/{fn}") as f:
            sys.stdout.write(f.read())
    print(f"checks failed: {nfail}" + (" (smoke: checks are informational)" if smoke else ""))
    if nfail and not smoke:
        sys.exit(2)


if __name__ == "__main__":
    main()
