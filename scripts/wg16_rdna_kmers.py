#!/usr/bin/env python
"""wg16: rDNA copy number from HiFi k-mers (see wg16_rdna_kmers.sh header for the method).
select <w>: write q_<unit>_<i>.fa query sequences (3 copies each of 18S, 25S, 5S from hap1; 5 single-copy controls).
count <w> <res>: parse the .wig counts, median per query, copy number = median / KMERCOV (diploid) and / 2 (haploid).
Absent k-mers (count 0, not reported by meryl-lookup) are counted as 0 in the median.
"""
import sys, os, re, subprocess, glob, statistics
H1 = "/home/sashanik/projects/def-santang3/sashanik/california_oaks/genomes/Qwislizeni_reasm_hap1.fa"
C8 = "/home/sashanik/scratch/california_oaks/pattern_hunt/close/out/c08"
FEAT = "/home/sashanik/scratch/california_oaks/pattern_hunt/broad/out/ph03_features_100k.tsv"
KMERCOV, K = 22.79, 21

def fa(reg, name, out):
    s = subprocess.run(f"samtools faidx {H1} {reg}", shell=True, check=True, capture_output=True, text=True).stdout
    open(out, "w").write(re.sub(r"^>.*", f">{name}", s, count=1))

def select(w):
    hits = {"18S": [], "25S": []}
    for l in open(f"{C8}/rdna45S_hap1.tsv"):
        t = l.split("\t")
        q, sid, ln, qs, qe, ss, se = t[0], t[1], int(t[3]), int(t[6 - 2]), int(t[5]), int(t[6]), int(t[7])
        if q.startswith("X16077") and qs <= 80 and qe >= 1899: hits["18S"].append((sid, min(ss, se), max(ss, se)))
        if q.startswith("X52320") and qs <= 850 and qe >= 4228: hits["25S"].append((sid, min(ss, se), max(ss, se)))
    five = []
    for l in open(f"{C8}/rdna5S_hap1.tbl"):
        if l[0] == "#": continue
        t = l.split()
        if int(t[5]) <= 5 and int(t[6]) >= 114: five.append((t[0], min(int(t[7]), int(t[8])), max(int(t[7]), int(t[8]))))
    hits["5S"] = five
    for u, hh in hits.items():
        # prefer chromosome-placed copies, spread along the list
        hh = sorted(hh, key=lambda x: ("_chr" not in x[0], x))
        pick = [hh[i] for i in (0, len(hh) // 2, len(hh) - 1)] if len(hh) >= 3 else hh
        print(u, len(hh), "hits; picked", pick)
        for i, (sid, s, e) in enumerate(pick):
            fa(f"{sid}:{s}-{e}", f"{u}_{i}_{sid}:{s}-{e}", f"{w}/q_{u}_{i}.fa")
    # single-copy controls: 100-kb windows with syn_hap2 >= 0.99, het 4-10 SNP/kb, genes >= 60/Mb, low TE; 2 kb from the middle
    import pandas as pd
    F = pd.read_csv(FEAT, sep="\t", usecols=["chrom", "start", "syn_hap2", "d_hap2", "genes_perMb", "TE_TEall"])
    F = F[(F.syn_hap2 >= 0.99) & F.d_hap2.between(4, 10) & (F.genes_perMb >= 60) & (F.TE_TEall < 0.3)].sample(5, random_state=16)
    for i, (c, s) in enumerate(zip(F.chrom, F.start)):
        fa(f"Qwis_hap1_{c}:{s + 49001}-{s + 51000}", f"control_{i}_{c}:{s + 49001}", f"{w}/q_control_{i}.fa")

def count(w, res):
    os.makedirs(res, exist_ok=True)
    rows = []
    for f in sorted(glob.glob(f"{w}/q_*.fa")):
        seq = "".join(l.strip() for l in open(f) if l[0] != ">"); name = open(f).readline()[1:].strip()
        nk = len(seq) - K + 1
        c = {}
        for l in open(f[:-3] + ".wig"):
            t = l.split()
            if len(t) == 2 and t[0].isdigit(): c[int(t[0])] = int(t[1])
        vals = [c.get(i, 0) for i in range(1, nk + 1)] if c and min(c) >= 1 else [c.get(i, 0) for i in range(nk)]
        med = statistics.median(vals); unit = os.path.basename(f)[2:-3].rsplit("_", 1)[0]
        rows.append((unit, name, len(seq), nk, sum(v > 0 for v in vals), med, round(med / KMERCOV, 1), round(med / KMERCOV / 2, 1)))
    import pandas as pd
    T = pd.DataFrame(rows, columns=["unit", "query", "bp", "kmers", "kmers_found", "median_kmer_count", "copies_diploid", "copies_per_haploid"])
    T.to_csv(f"{res}/rdna_kmer_copy_number.tsv", sep="\t", index=False)
    S = T.groupby("unit").agg(median_count=("median_kmer_count", "median"), diploid=("copies_diploid", "median"), haploid=("copies_per_haploid", "median"))
    S.to_csv(f"{res}/rdna_kmer_copy_number_summary.tsv", sep="\t")
    # assembled floor (c08): full-length 18S copies in hap1 / hap2
    def nfull(p):
        try: return sum(1 for l in open(p) if l.startswith("X16077") and int(l.split("\t")[3]) >= 1800)
        except FileNotFoundError: return "NA"
    a1, a2 = nfull(f"{C8}/rdna45S_hap1.tsv"), nfull(f"{C8}/rdna45S_hap2.tsv")
    open(f"{res}/README.md", "w").write(
        "# rDNA copy number from HiFi k-mers\n\n"
        f"meryl-lookup -wig-count of each query's 21-mers in hifi_k21.meryl; copies (diploid) = median k-mer count / {KMERCOV}\n"
        "(GenomeScope2 per-haploid-copy k-mer coverage); per haploid = /2. Single-copy controls should give ~2 diploid copies.\n"
        f"Assembled full-length 18S (alignment >= 1800 bp, c08): hap1 {a1}, hap2 {a2} (a floor).\n"
        "k-mer totals are diploid: they cannot say which haplotype is under-assembled.\n\n" + S.to_string() + "\n")
    print(T.to_string()); print(S.to_string()); print("assembled 18S hap1/hap2:", a1, a2)

if __name__ == "__main__":
    select(sys.argv[2]) if sys.argv[1] == "select" else count(sys.argv[2], sys.argv[3])
