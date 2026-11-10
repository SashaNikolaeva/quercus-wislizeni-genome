# Quercus wislizeni genome assembly: analysis scripts

Scripts for the haplotype-resolved, chromosome-scale genome assembly of interior live oak (*Quercus wislizeni*,
isolate QUWIS_CCGP_1), reported in the genome report by Nikolaeva et al. (in preparation).

- Raw data: NCBI BioProject PRJNA808371 (PacBio HiFi SRX15304001; Omni-C SRX15304002, SRX15304003)
- Assemblies: GenBank GCA_023055345.2 (haplotype 1) and GCA_023055485.2 (haplotype 2)

All steps were run on the Digital Research Alliance of Canada cluster Narval as SLURM shell scripts and Python
scripts. Script prefixes: `rp` = assembly, `wg` = genome report analyses,
`cg` = annotation steps shared with a comparative project.

## Scripts by paper section

| Paper section | Scripts |
|---|---|
| Sequencing data, read QC, genome size | `rp07_download_wislizeni_raw.sh`, `rp08_wislizeni_read_qc.sh` |
| Assembly (hifiasm Hi-C mode) | `rp09_wislizeni_hifiasm.sh`, `rp10_wislizeni_contig_qc.sh` |
| Omni-C scaffolding and contact maps | `rp11_wislizeni_yahs.sh`, `rp12_wislizeni_contact_maps.sh`, `rp12b_contact_png.py` |
| Chromosome naming, synteny, telomeres, unplaced sequence | `rp13_wislizeni_validate.sh`, `rp13_synteny.py`, `rp13c_hap_sizes.*`, `rp13d_big_unplaced.*`, `rp13e_blast_unplaced.sh`, `rp13f_missing_telomeres.sh`, `rp13f_telomere_scan.py` |
| Assembly statistics, Table 2, Figure 2 | `wg02_assembly_stats.sh`, `wg03_fig2.py`, `wg24_compleasm_odb12.sh` |
| Contamination and organelle screening, cleaned assembly, Figure 5 | `wg04a_fcs_download.sh`, `wg04b_fcs_adaptor.sh`, `wg04c_fcs_gx.sh`, `wg07_hifi_junctions.*` (HiFi mapping to the diploid assembly, for read depth), `wg09_organelle_screen.sh`, `wg25_organelle_classify.*`, `wg28_clean_fasta.sh`, `wg27_blobtoolkit.*` |
| Organelle genomes | `wg26_oatk.sh`, `wg26b_organelle_check.*` |
| Gene annotation, Table 4 | `upstream/cg05_helixer.sh`, `wg01_helixer_haps.sh`, `wg10_proteome_te_genes.*`, `wg19_annotation_prep.*`, `wg20_annotation_stats.*`, `wg21_interproscan.sh`, `wg21b_interproscan_summary.*`, `wg24b_busco_dup_check.*` |
| Non-coding RNA genes, rDNA copy number | `wg22_trnascan.sh`, `wg23_barrnap.sh`, `wg22b_ncrna_summary.*`, `wg16_rdna_kmers.*` |
| Transposable elements, Table 3, Figures 3-4 | `upstream/cg04_edta.sh`, `wg05_te_summary.*` |
| Discussion (published oak assemblies on NCBI) | `wg29_ncbi_oak_assemblies.py` |
