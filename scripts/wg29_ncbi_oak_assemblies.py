#!/usr/bin/env python
"""wg29: published oak (Quercus) genome assemblies on NCBI, for the the genome report Discussion.
Input: JSON lines from `datasets summary genome taxon Quercus --assembly-source all --as-json-lines` (run 2026-11-10).
Output (stdout): TSV, one row per assembly (GenBank only; RefSeq copies dropped). Fields: organism, accession, level,
assembly name, release date, assembly type, total length, contig N50, scaffold N50, n scaffolds.
Then (stderr) a summary: species with >= 1 chromosome-level assembly; for each species the chromosome-level assembly
with the largest contig N50.
Usage: wg29_ncbi_oak_assemblies.py ncbi_quercus.jsonl > ncbi_quercus.tsv 2> summary.txt
"""
import sys, json, re
rows = []
for l in open(sys.argv[1]):
    d = json.loads(l)
    acc = d["accession"]
    if acc.startswith("GCF_"): continue
    ai, st = d.get("assembly_info", {}), d.get("assembly_stats", {})
    org = d["organism"]["organism_name"]
    rows.append((org, acc, ai.get("assembly_level", ""), ai.get("assembly_name", ""), ai.get("release_date", ""),
                 ai.get("assembly_type", ""), int(st.get("total_sequence_length", 0)), int(st.get("contig_n50", 0)),
                 int(st.get("scaffold_n50", 0)), int(st.get("number_of_scaffolds", 0))))
rows.sort()
print("organism\taccession\tlevel\tassembly_name\trelease\ttype\ttotal_len\tcontig_N50\tscaffold_N50\tn_scaffolds")
for r in rows: print("\t".join(map(str, r)))
sp = lambda o: re.sub(r" (subsp\.|var\.|x ).*", "", o)
allsp = sorted({sp(r[0]) for r in rows}); chrom = {}
for r in rows:
    if r[2] == "Chromosome" and (sp(r[0]) not in chrom or r[7] > chrom[sp(r[0])][7]): chrom[sp(r[0])] = r
print(f"{len(rows)} GenBank assemblies, {len(allsp)} species; {len(chrom)} species with a chromosome-level assembly", file=sys.stderr)
for s, r in sorted(chrom.items(), key=lambda x: -x[1][7]):
    print(f"  {s}\t{r[1]}\t{r[3]}\t{r[4]}\ttotal {r[6]/1e6:.1f} Mb\tcontig N50 {r[7]/1e6:.2f} Mb\tscaffold N50 {r[8]/1e6:.2f} Mb", file=sys.stderr)
