# Poster references (Bio253, group mval)

The numbers match the numbers in the poster text.

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

We drew all figures from our own analysis of the course data. They use these external sources and tools.

### Databases and annotation

6. NCBI GenBank, *S. aureus* JE2 genome [CP020619.1](https://www.ncbi.nlm.nih.gov/nuccore/CP020619.1)
   (assembly GCA_002085525.1, ASM208552v1), as supplied with the course.
   Used for: locus tags, gene positions, strands, product names.
7. UniProt, *S. aureus* NCTC 8325 reference proteome
   [UP000008816](https://www.uniprot.org/proteomes/UP000008816), via the course annotation table
   (download date not recorded).
   Used for: gene symbols. They are transferred from the closest NCTC 8325 protein (BLAST), not determined in JE2.
8. Nebraska Transposon Mutant Library, [Tn insertion map](https://ntml.unmc.edu/tninsertionmap).
   Used for: the insertion sites of NE217 and NE1919.
9. AureoWiki, entries [SAOUHSC_01186](https://aureowiki.med.uni-greifswald.de/SAOUHSC_01186) and
   [SAOUHSC_01187](https://aureowiki.med.uni-greifswald.de/SAOUHSC_01187), section Operon (accessed 12 Nov 2026).
   Used for: the *stp1–stk1* transcription unit in NCTC 8325. Its position in JE2 is inferred.
10. U. Mäder *et al.*, *PLoS Genet.*, 2016, PMID [27035918](https://pubmed.ncbi.nlm.nih.gov/27035918/).
    Used for: the same transcription unit, as cited by AureoWiki.

### Statistics and software

11. R Core Team, *R: A Language and Environment for Statistical Computing* (version 4.5.0),
    R Foundation for Statistical Computing, Vienna, Austria, 2025. <https://www.R-project.org/>
12. M. E. Ritchie *et al.*, "limma powers differential expression analyses for RNA-sequencing and microarray
    studies," *Nucleic Acids Research*, vol. 43, no. 7, e47, 2015,
    doi: [10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007).
13. D. J. McCarthy and G. K. Smyth, *Bioinformatics*, vol. 25, pp. 765–771, 2009,
    doi: [10.1093/bioinformatics/btp053](https://doi.org/10.1093/bioinformatics/btp053).
    Used for: the test against a 2-fold change (`treat`).
14. Y. Benjamini and Y. Hochberg, "Controlling the false discovery rate: a practical and powerful approach to
    multiple testing," *Journal of the Royal Statistical Society Series B*, vol. 57, pp. 289–300, 1995,
    doi: [10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x).
15. H. Wickham, *ggplot2: Elegant Graphics for Data Analysis*, Springer-Verlag New York, 2016,
    ISBN 978-3-319-24277-4. <https://ggplot2.tidyverse.org>
16. T. Pedersen, *patchwork: The Composer of Plots*, R package version 1.3.0, 2024,
    doi: [10.32614/CRAN.package.patchwork](https://doi.org/10.32614/CRAN.package.patchwork).

### Data

17. Bio253 course data set (unpublished): total proteome of *S. aureus* JE2 anc, stk and stp, n = 3 per group,
    with abundances normalised (vsn) in the supplied prolfquapp output. The raw data are not in this repository.

## Sources by figure

| Figure | What it shows | External content | Method and software | Data | Our script |
|---|---|---|---|---|---|
| Gene map (introduction) | *stp1–stk1* region and insertions | positions [6], gene symbols [7], mutants [5], insertion sites [8], transcription unit [9, 10] | [11, 15] | detection only [17] | [04](../../scripts/04_polar_effect_mval.R), [04d](../../scripts/04d_gene_context_premium_mval.R) |
| 1A | peptides per sample | mutants [5] | [11, 15] | [17] | [03](../../scripts/03_mutant_validation_mval.R) |
| 1B | proteins next to *stp1/stk1* | gene order [6], gene symbols [7], mutants [5] | [11, 15] | [17] | [04](../../scripts/04_polar_effect_mval.R), [04b](../../scripts/04b_neighbours_by_genotype_mval.R) |
| 2A | PASN response, each mutant vs anc | locus tags [6], gene symbols [7] | [11–15] | [17] | [05](../../scripts/05_formal_test_mval.R), [06](../../scripts/06_altered_responses_mval.R) |
| 2B | PASN response, stp vs stk | locus tags [6], gene symbols [7] | [11–15] | [17] | [05](../../scripts/05_formal_test_mval.R), [06](../../scripts/06_altered_responses_mval.R) |
| 2C | *arcA* and its neighbours | positions and strands [6], gene symbols [7] | [11–16] | [17] | [05](../../scripts/05_formal_test_mval.R), [06c](../../scripts/06c_hit_neighbourhoods_mval.R) |
| 3A | stp effects vs stk effects | product names [6], gene symbols [7] | [11–16], own noise correction | [17] | [07](../../scripts/07_shared_distinct_mval.R), [12](../../scripts/12_poster_figures_mval.R) |
| 3B | genome scan | positions and product names [6], gene symbols [7], *stp1/stk1* locus [6] | [11, 12, 15, 16], own permutation threshold | [17] | [11](../../scripts/11_robustness_mval.R), [11b](../../scripts/11b_genome_blocks_mval.R) |
| 3C | close-up of R6 and R7 | positions, strands and product names [6], gene symbols [7] | [11, 12, 15, 16] | [17] | [11c](../../scripts/11c_R6_R7_context_mval.R) |

**Notes**

- The scripts produce the repository versions of the figures. The poster versions of 1A, 1B, 2A, 2B and 2C have
  a different layout and labels.
- "Same-strand run" (2C) and "putative operon" (3C) are our own prediction (same strand, gaps < 150 bp).
  They do not come from an operon database.
- The region names R1–R8 and their short descriptions are ours, based on [6] and [7].
- No figure from 1A to 3C uses BioCyc, STRING, KEGG or phosphoproteome data.
- [6]–[17] were compiled from our scripts and R's built-in citations with Claude Code (AI). The AI did not open
  the linked websites or papers.

## Our analysis

- [Code and reports](https://github.com/dinvuc/bio253-je2-proteomics) (no raw data)
