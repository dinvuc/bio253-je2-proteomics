# Poster references (Bio253, group mval)

The numbers match the numbers in the poster text. The sources of each figure are listed by figure further down.

## Literature

1. M. Débarbouillé *et al.*, "Characterization of a Serine/Threonine Kinase Involved in Virulence of
   *Staphylococcus aureus*," *J. Bacteriol.*, vol. 191, no. 13, pp. 4070–4081, Jul. 2009,
   doi: [10.1128/JB.01813-08](https://doi.org/10.1128/JB.01813-08).
2. F. Sun *et al.*, "Protein cysteine phosphorylation of SarA/MgrA family transcriptional regulators
   mediates bacterial virulence and antibiotic resistance," *Proceedings of the National Academy of Sciences*,
   vol. 109, no. 38, pp. 15461–15466, Sep. 2012,
   doi: [10.1073/pnas.1205952109](https://doi.org/10.1073/pnas.1205952109).
3. J. Makhlin *et al.*, "*Staphylococcus aureus* ArcR Controls Expression of the Arginine Deiminase Operon,"
   *J. Bacteriol.*, vol. 189, no. 16, pp. 5976–5986, Aug. 2007,
   doi: [10.1128/JB.00592-07](https://doi.org/10.1128/JB.00592-07).
4. N. Prust, S. van der Laarse, H. W. P. van den Toorn, N. M. van Sorge, and S. Lemeer, "In-Depth
   Characterization of the *Staphylococcus aureus* Phosphoproteome Reveals New Targets of Stk1,"
   *Molecular & Cellular Proteomics*, vol. 20, p. 100034, 2021,
   doi: [10.1074/mcp.RA120.002232](https://doi.org/10.1074/mcp.RA120.002232).
5. P. D. Fey *et al.*, *mBio*, vol. 4, e00537-12, 2013,
   doi: [10.1128/mBio.00537-12](https://doi.org/10.1128/mBio.00537-12).
   Cited for: the Nebraska Transposon Mutant Library, the source of our two mutants.

## Figure sources

We made all figures from our own analysis of the Bio253 course data (unpublished, not in this repository):
*S. aureus* JE2 anc, *stk* and *stp*, total proteome, n = 3 per group, abundances normalised (vsn) in the
supplied prolfquapp output. Below is what was needed to make each figure.

### Gene map (introduction) – the *stp1–stk1* region

- **Data:** whether each protein was detected in the total proteome
- **Genome positions and strands:** [NCBI GenBank CP020619.1][genbank]
- **Gene symbols:** [UniProt UP000008816][uniprot]
- **Mutants:** Nebraska Transposon Mutant Library ([Fey et al. 2013][fey])
- **Insertion sites:** [NTML Tn insertion map][ntml]
- **Transcription unit:** AureoWiki [SAOUHSC_01186][aw86] and [SAOUHSC_01187][aw87] (accessed 12 Nov 2026), citing
  [Mäder et al. 2016][mader]. Known for NCTC 8325, inferred for JE2.
- **Software:** R 4.5.0, [ggplot2][ggplot]
- **Scripts:** [04][s04], [04d][s04d]
- **Figure in the repository:** [gene map](../results/step04_polar_effect/fig_gene_context_premium_wide_mval.png)

### Fig. 1A – Peptide evidence for both enzymes

- **Data:** number of peptides quantified per protein and sample
- **External databases:** none
- **Mutants:** Nebraska Transposon Mutant Library ([Fey et al. 2013][fey])
- **Method:** counting only, no statistical test
- **Software:** R 4.5.0, [ggplot2][ggplot]
- **Script:** [03][s03]
- **Figure in the repository:** [peptides per sample](../results/step03_validation/fig_enzyme_peptides_mval.png)
  (the poster version has a different layout and labels)

### Fig. 1B – The neighbouring genes

- **Data:** normalised log2 protein abundance per sample
- **Gene order and locus tags:** [NCBI GenBank CP020619.1][genbank]
- **Gene symbols:** [UniProt UP000008816][uniprot]
- **Mutants:** Nebraska Transposon Mutant Library ([Fey et al. 2013][fey])
- **Method:** values per sample, no statistical test
- **Software:** R 4.5.0, [ggplot2][ggplot]
- **Scripts:** [04][s04], [04b][s04b]
- **Figure in the repository:** [neighbours per strain](../results/step04_polar_effect/fig_neighbours_by_genotype_mval.png)
  (the poster version has a different layout and labels)

### Fig. 2A – PASN response, each mutant vs anc

- **Data:** normalised log2 protein abundance in TSB and PASN
- **Locus tags:** [NCBI GenBank CP020619.1][genbank]
- **Gene symbols:** [UniProt UP000008816][uniprot]
- **Statistics:** [limma][limma] (Ritchie et al. 2015), [treat][treat] test against a 2-fold change
  (McCarthy & Smyth 2009), [Benjamini–Hochberg][bh] FDR < 0.05. Pink = exploratory (FDR 0.05–0.1), not a hit.
- **Software:** R 4.5.0, [ggplot2][ggplot]
- **Scripts:** [05][s05], [06][s06]
- **Figure in the repository:** [response vs anc](../results/step06_responses/fig_response_vs_anc_mval.png)
  (the poster version has a different layout and labels)

### Fig. 2B – PASN response, *stp* vs *stk*

- **Data:** normalised log2 protein abundance in TSB and PASN
- **Locus tags:** [NCBI GenBank CP020619.1][genbank]
- **Gene symbols:** [UniProt UP000008816][uniprot]
- **Statistics:** [limma][limma] (Ritchie et al. 2015), [treat][treat] test against a 2-fold change
  (McCarthy & Smyth 2009), [Benjamini–Hochberg][bh] FDR < 0.05. Pink = exploratory (FDR 0.05–0.1), not a hit.
- **Software:** R 4.5.0, [ggplot2][ggplot]
- **Scripts:** [05][s05], [06][s06]
- **Figure in the repository:** [stp vs stk response](../results/step06_responses/fig_response_stk_vs_stp_mval.png)
  (the poster version has a different layout and labels)

### Fig. 2C – *arcA* and its neighbours

- **Data:** normalised log2 protein abundance per sample in TSB and PASN
- **Genome positions and strands:** [NCBI GenBank CP020619.1][genbank]
- **Gene symbols:** [UniProt UP000008816][uniprot]
- **Statistics:** [limma][limma] (Ritchie et al. 2015), [treat][treat] test against a 2-fold change
  (McCarthy & Smyth 2009), [Benjamini–Hochberg][bh] FDR < 0.05
- **Own prediction:** "same-strand run" = same strand, gaps < 150 bp. No operon database was used.
- **Software:** R 4.5.0, [ggplot2][ggplot], [patchwork][patchwork]
- **Scripts:** [05][s05], [06b][s06b], [06c][s06c]
- **Figure in the repository:** [arcA neighbourhood](../results/step06_responses/hit_neighbourhoods/fig_neighbourhood_arcA_14675_mval.png)
  (the poster version has a different layout and labels)

### Fig. 3A – *stp* effects vs *stk* effects

- **Data:** log2 effects of each mutant vs anc, in TSB and on the PASN response
- **Protein names:** [NCBI GenBank CP020619.1][genbank]
- **Gene symbols:** [UniProt UP000008816][uniprot]
- **Statistics:** [limma][limma] (Ritchie et al. 2015), [treat][treat] test against a 2-fold change
  (McCarthy & Smyth 2009), [Benjamini–Hochberg][bh] FDR < 0.05
- **Own method:** r corrected for the noise of the shared anc samples, with a bootstrap 95 % CI
- **Software:** R 4.5.0, [ggplot2][ggplot], [patchwork][patchwork]
- **Scripts:** [05][s05], [07][s07], [12][s12]
- **Figure in the repository:** [Fig. 3A](../results/step12_poster/fig3_h3_shared_distinct_compact_mval.png)

### Fig. 3B – Genome scan

- **Data:** log2 differences of each mutant vs anc in TSB, in PASN and for the PASN response
- **Genome positions and product names:** [NCBI GenBank CP020619.1][genbank]
- **Gene symbols:** [UniProt UP000008816][uniprot]
- **Statistics:** [limma][limma] (Ritchie et al. 2015) for the log2 differences
- **Own method:** mean of 10 neighbouring proteins, 5 % threshold from 1000 permutations per panel.
  The names of regions R1–R8 are ours.
- **Software:** R 4.5.0, [ggplot2][ggplot], [patchwork][patchwork]
- **Scripts:** [05][s05], [11][s11], [11b][s11b]
- **Figure in the repository:** [Fig. 3B](../results/step11_robustness/fig_genome_scan_labelled_compact_mval.png)

### Fig. 3C – Close-up of R6 and R7

- **Data:** log2 differences of each mutant vs anc in TSB
- **Genome positions, strands and product names:** [NCBI GenBank CP020619.1][genbank]
- **Gene symbols:** [UniProt UP000008816][uniprot] (the match for *hlb* is weak)
- **Statistics:** [limma][limma] (Ritchie et al. 2015) for the log2 differences
- **Own prediction:** "putative operon" = same strand, gaps < 150 bp. No operon database was used.
- **Software:** R 4.5.0, [ggplot2][ggplot], [patchwork][patchwork]
- **Scripts:** [05][s05], [11b][s11b], [11c][s11c]
- **Figure in the repository:** [Fig. 3C](../results/step11_robustness/fig_R6_R7_genome_context_compact_mval.png)

## Full details of the figure sources

- **AureoWiki:** entries [SAOUHSC_01186][aw86] and [SAOUHSC_01187][aw87], section Operon.
- **Benjamini Y, Hochberg Y (1995).** Controlling the false discovery rate: a practical and powerful approach to
  multiple testing. *Journal of the Royal Statistical Society Series B* 57, 289–300.
  doi: [10.1111/j.2517-6161.1995.tb02031.x][bh]
- **GenBank:** *S. aureus* JE2, [CP020619.1][genbank], assembly GCA_002085525.1 (ASM208552v1), as supplied with the
  course.
- **Mäder U et al. (2016).** *PLoS Genetics*. PMID [27035918][mader]
- **McCarthy DJ, Smyth GK (2009).** *Bioinformatics* 25, 765–771. doi: [10.1093/bioinformatics/btp053][treat]
- **Pedersen T (2024).** *patchwork: The Composer of Plots*. R package version 1.3.0.
  doi: [10.32614/CRAN.package.patchwork][patchwork]
- **R Core Team (2025).** *R: A Language and Environment for Statistical Computing*. R Foundation for Statistical
  Computing, Vienna, Austria. <https://www.R-project.org/>
- **Ritchie ME et al. (2015).** limma powers differential expression analyses for RNA-sequencing and microarray
  studies. *Nucleic Acids Research* 43(7), e47. doi: [10.1093/nar/gkv007][limma]
- **UniProt:** *S. aureus* NCTC 8325 reference proteome [UP000008816][uniprot], via the course annotation table
  (download date not recorded). Gene symbols are transferred from the closest NCTC 8325 protein (BLAST), not
  determined in JE2.
- **Wickham H (2016).** *ggplot2: Elegant Graphics for Data Analysis*. Springer-Verlag New York.
  ISBN 978-3-319-24277-4. <https://ggplot2.tidyverse.org>

No figure from 1A to 3C uses BioCyc, STRING, KEGG or phosphoproteome data. The figure sources were compiled from our
scripts and R's built-in citations with Claude Code (AI), which did not open the linked websites or papers.

## Our analysis

- [Code and reports](https://github.com/dinvuc/bio253-je2-proteomics) (no raw data)

[genbank]: https://www.ncbi.nlm.nih.gov/nuccore/CP020619.1
[uniprot]: https://www.uniprot.org/proteomes/UP000008816
[fey]: https://doi.org/10.1128/mBio.00537-12
[ntml]: https://ntml.unmc.edu/tninsertionmap
[aw86]: https://aureowiki.med.uni-greifswald.de/SAOUHSC_01186
[aw87]: https://aureowiki.med.uni-greifswald.de/SAOUHSC_01187
[mader]: https://pubmed.ncbi.nlm.nih.gov/27035918/
[limma]: https://doi.org/10.1093/nar/gkv007
[treat]: https://doi.org/10.1093/bioinformatics/btp053
[bh]: https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
[ggplot]: https://ggplot2.tidyverse.org
[patchwork]: https://doi.org/10.32614/CRAN.package.patchwork
[s03]: ../../scripts/03_mutant_validation_mval.R
[s04]: ../../scripts/04_polar_effect_mval.R
[s04b]: ../../scripts/04b_neighbours_by_genotype_mval.R
[s04d]: ../../scripts/04d_gene_context_premium_mval.R
[s05]: ../../scripts/05_formal_test_mval.R
[s06]: ../../scripts/06_altered_responses_mval.R
[s06b]: ../../scripts/06b_hit_genome_context_mval.R
[s06c]: ../../scripts/06c_hit_neighbourhoods_mval.R
[s07]: ../../scripts/07_shared_distinct_mval.R
[s11]: ../../scripts/11_robustness_mval.R
[s11b]: ../../scripts/11b_genome_blocks_mval.R
[s11c]: ../../scripts/11c_R6_R7_context_mval.R
[s12]: ../../scripts/12_poster_figures_mval.R
