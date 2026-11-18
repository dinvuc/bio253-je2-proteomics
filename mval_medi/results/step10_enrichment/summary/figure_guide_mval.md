# Step 10 outputs: what each figure and table compares

Made by [`scripts/10_enrichment_mval.R`](../../../../scripts/10_enrichment_mval.R) (one folder
per contrast, next to this one) and
[`scripts/10b_step10_summary_mval.R`](../../../../scripts/10b_step10_summary_mval.R) (this folder).
Question: is Stk1 "broadly involved" (presentation feedback): **how many functional
categories** shift in `stk` vs `anc`, and do the same categories shift in `stp`
(H3 at pathway level)?

## Methods in short

- **Background** = the proteins tested in that contrast (≥ 2 of 3 in every group =
  the Step 5 eligible set), not the genome. Terms: GO, KEGG, UniProt keywords,
  Reactome, COMPARTMENTS, STRING local clusters (`read_protein_terms()`); size 5–500;
  BH within each annotation category.
- **ORA** (Fisher, one-sided) on the Step 5 formal hits, up and down separately.
  With 1–20 hits per contrast many "significant" ORA terms rest on **1–2 proteins**:
  they restate the hit proteins, not a pathway. They are counted separately.
- **Rank test** (Wilcoxon, all tested proteins ranked by the moderated t). It treats
  proteins as independent, which is anti-conservative when pathway members are
  correlated (Wu & Smyth 2012, *Nucleic Acids Res* 40:e133,
  [doi:10.1093/nar/gks461](https://doi.org/10.1093/nar/gks461)). We therefore also run
  **limma `cameraPR`** (inter-gene correlation 0.01) on the same statistic.
  **Robust shift = FDR < 0.05 in both, same direction.** Rules set 17 Nov before this
  summary was run; the group still has to confirm them.
- **Direction:** `auc` > 0.5 = term members rank towards the top = **higher in the
  mutant** (baseline contrasts) or **SN response lower in the mutant than in `anc`**
  (interaction contrasts). Every row of `robust_terms_mval.csv` says this in words.

## Tables

### `term_counts_mval.csv`: how many terms shift?

| mutant | layer | formal protein hits | ORA ≥ 3 hits | ORA 1–2 hits | Wilcoxon | cameraPR | **robust** |
|---|---|---|---|---|---|---|---|
| `stk` | baseline TSB | 20 | 13 | 17 | 51 | 75 | **35** |
| `stk` | baseline SN | 12 | 5 | 18 | 316 | 174 | **143** |
| `stk` | SN response | 5 | 1 | 13 | 97 | 62 | **44** |
| `stp` | baseline TSB | 14 | 1 | 24 | 284 | 202 | **192** |
| `stp` | baseline SN | 4 | 11 | 38 | 294 | 179 | **168** |
| `stp` | SN response | 1 | 0 | 21 | 41 | 14 | **8** |

- The counts are **not** numbers of independent pathways: GO terms are nested and STRING
  clusters overlap.
- **"Stk1 broadly involved"?** Many terms shift in `stk`, but `stp` has as many or more
  in both media, with fewer formal hits. Broad pathway shifts are therefore not specific
  to the kinase mutant. What is specific is the **kind** of terms (see `robust_terms`):
  - `stk`: small, specific modules (median term size 8 in TSB, 11.5 in the response).
  - `stp`: large, general terms (median size 32.5 in TSB, 40 in SN), e.g. nucleic-acid
    metabolism higher and surface proteins lower.
- Only in the **SN response** does `stk` clearly have more robust shifts than `stp`
  (44 vs 8).

### `robust_terms_mval.csv`
All robust terms with AUC, both FDRs and the direction in words. Top terms (cameraPR FDR):
- `stk_vs_anc_at_TSB` / `_at_SN`, higher in `stk`: the STRING cluster CL:1356
  ("…Protein serine/threonine phosphatase activity, and Sigma factor activity");
  vancomycin resistance (map01502); peptidoglycan synthesis (CL:932); ATP synthase F1
  (GO:0033178, highest in SN).
  - **Caveat: CL:1356 is not Stp1/Stk1.** Its members are the SigB operon (`rsbV`,
    `rsbW`, SigB, and the RsbU-like phosphatase `B7H15_11485`) plus `mazEF` and
    neighbours, locus tags B7H15_11465–11525, all +0.7 to +1.75 log2 in `stk`.
    `alr1`, `murF` and `ddl` (vancomycin/peptidoglycan terms) and `kdpD` (Step 6) lie in
    the same stretch (B7H15_11465–11570). Part of the `stk` "pathway" signal may be
    **one genomic region**. Checked in Step 11.
- `interaction_stk_vs_anc_at_TSB_vs_SN`, SN response lower in `stk`: de novo IMP/purine
  biosynthesis (CL:699, AUC 0.98), branched-chain amino-acid and diaminopimelate biosynthesis.
- `stp_vs_anc_at_TSB`, higher in `anc` than `stp`: signal peptide, cell wall,
  peptidoglycan anchor, secreted, cell adhesion. Higher in `stp`: nucleic-acid-related terms.
- `interaction_stp_vs_anc_at_TSB_vs_SN` (8 terms, all `auc` < 0.5, i.e. SN response
  higher in `stp` than in `anc`): cell wall, peptidoglycan anchor, secreted, signal,
  extracellular region, cell adhesion.

### `auc_correlation_stk_stp_mval.csv`: pathway-level similarity of the two mutants
Pearson r of the term AUCs (`stk` vs `stp`) over all tested terms: **0.23** baseline TSB,
**0.53** baseline SN, **0.30** SN response. This is descriptive: terms overlap, and both
contrasts share the `anc` reference (Step 7), so even pure noise gives r > 0.

### `candidate_annotation_mval.csv`: source-linked candidate table
All formal hits (Step 5 rule) in the H1–H3 contrasts (`stk`/`stp`/`stp_vs_stk` at TSB
and the three interactions), plus the genotype phosphosite hits (Step 9): 71 rows,
49 proteins, 43 with a good UniProt ortholog. Columns: estimate, SE, FDR, **direction in
words**, `gene_symbol_best` (only if `ortholog_quality == "good"`), UniProt name,
accession and **UniProt link**, EC, KEGG, keywords. Names come from
`read_protein_metadata()`, not from the AI.

### `phospho_ora_mval.csv`: phospho enrichment, own background
Protein level: 25 proteins with ≥ 1 site hit (Step 8) against the 133 proteins with ≥ 1
supported site, each protein counted once. 12 terms with FDR < 0.05, all with ≥ 3 hit
proteins: carbohydrate metabolism, glycolysis and the pentose phosphate pathway
(STRING clusters, GO). Caveat: 55 of the 62 site hits are TSB-vs-SN contrasts, so this
describes the **medium** response of the phosphoproteome, not a genotype effect.

### STRING files
`string_upload_candidates_mval.txt` (49 locus tags) for "Multiple proteins", organism
`STRG0B09NBK`. **Set `string_background_tested_mval.txt` as the background** in STRING;
the website's default (whole genome) inflates easy-to-measure categories.

## Figures

### `fig_pathway_counts_mval.png`: Fig 4a
- **Compares:** the number of robust terms per contrast. Left `stk` vs `anc`, right `stp`
  vs `anc`. Above 0 = higher in the mutant / SN response lower in the mutant; below 0 =
  the opposite. Colour = annotation source.
- **Read with:** the counts table above (overlapping terms; `stp` ≥ `stk` except in the
  response).
- **Key points:** robust terms `stk` 35 / 143 / 44 and `stp` 192 / 168 / 8 (baseline TSB /
  baseline SN / SN response). `stk` terms are small, specific modules; `stp` terms are
  large, general ones.
- **Conclusion:** both mutants shift many functional categories, so "Stk1 broadly involved"
  is not specific to the kinase mutant. Only the change of the PASN response is broader
  in `stk` (44 vs 8 terms).

### `fig_pathway_shifts_mval.png`: Fig 4b
- **Compares:** the 16 top robust `stk` terms (size ≤ 200; a term nested ≥ 80 % in one
  already shown is dropped) in all three `stk` contrasts (left) and the same terms in
  `stp` (right).
- **How to read:** colour = AUC − 0.5 (orange = higher in the mutant / SN response lower
  in the mutant). Circle = robust, square = one test only, triangle = not significant.
  Size = cameraPR significance.
- **Result:** `stp` shares the peptidoglycan-synthesis shift (higher) and the
  cell-killing / haemolysis clusters (lower). It does **not** share the SigB/`mazEF`
  cluster, ATP synthase F1 or vancomycin resistance, and not the purine-response
  difference.
- **Conclusion (H3 at pathway level):** `stp` shares part of the `stk` pattern (cell-wall
  synthesis up, toxins/cytolysis down), consistent with its reduced Stk1. The most specific
  `stk` signals are not shared, and the SigB/*mazEF* one coincides with a genomic block of
  the `stk` strain (Step 11).

## Interpretation (marked as interpretation; group to check)
- Pathway level supports Step 7: `stp` shares **part** of the `stk` pattern (cell wall
  synthesis, cytolysis-related terms) and has its own, broader changes.
- The cytolysis / cell-killing terms are lower in `stk` (and `stp`) than in `anc`.
  Burnside et al. 2010 (*PLOS ONE* 5:e11071, in
  [our evidence review](../../../docs/stk_stp_rescue_evidence_mval.md)) report altered
  alpha-toxin expression and haemolysis in Newman *stk1*/*stp1* deletions. This is the
  same **topic** in a different strain and assay. The direction must be checked in the
  paper before any claim.
- Nothing here shows direct Stk1 substrates or rescue. Pathway shifts are measured
  protein-abundance changes.
