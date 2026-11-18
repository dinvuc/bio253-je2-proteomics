# Step 6d outputs: what each figure and table compares

Made by [`scripts/06d_response_hits_string_mval.R`](../../../../scripts/06d_response_hits_string_mval.R)
(~1 min; queries STRING once, then uses the cached `string_*_mval.tsv`).
Proteins: the 7 coloured points of
[`../fig_response_vs_anc_mval.png`](../fig_response_vs_anc_mval.png) (6 formal interaction
hits + B7H15_06000, exploratory). Biology and sources: [`hits_biology_mval.md`](hits_biology_mval.md).

**Reading key.** response = log2(SN) − log2(TSB). Node colour / delta = response in the
mutant − response in `anc` (= −interaction diff). **Red = rises more (falls less) in SN in
the mutant; blue = rises less / falls more.** Our limma `treat(lfc = 1)`, BH-FDR < 0.05,
observed estimates only, n = 3. `stp` also lacks most Stk1.

## Figures

### `fig_string_network_mval.png`: STRING neighbourhood coloured by our data
- **Compares:** for every protein in the network, the PASN response of `stk` (left) or
  `stp` (right) against the PASN response of `anc`.
- **Data:** STRING v12.5, *S. aureus* NCTC 8325 (taxon 93061), edges with combined score
  ≥ 0.7. Network = the 7 proteins + up to 6 best partners each (30 nodes, 64 edges).
  Each STRING node is mapped back to a JE2 locus tag via the UniProt/KEGG ortholog
  (`up_kegg`, `ortholog_quality == "good"`).
- **How to read:** square = Fig 2 protein, circle = STRING partner, thick outline =
  formal interaction hit in that panel, grey = not testable (< 2/3 replicates in a group),
  line width = STRING score. Label = gene (or short product) + locus number. Layout: one
  circle per connected module (hits in the middle), positions are arbitrary.
- **Result:** 5 separate modules + 1 isolated protein (B7H15_06000). Only arcA–argF are
  linked to each other. Partners that move with their hit: betA with B7H15_14515 (`stk`),
  ilvD / B7H15_11390 with ilvA (`stk`). Partners that do not: rbsK (rbsD), arcC1 / argG /
  argH (arcA/argF). The kdp partners are almost all grey (not detected in SN).
- **Caveat:** STRING edges include text-mining and co-expression in other organisms, so
  an edge means "functionally associated", not a physical contact. Mapping choice: JE2
  has a second copy of arcA / argF / arcC1 (B7H15_00330–00350); when two JE2 proteins
  map to one STRING node, the Fig 2 protein is kept, otherwise the best ortholog, then the
  one detected in more samples (listed in `checks_mval.txt`; argF B7H15_00335 vs 06085
  was a tie).
- **Conclusion:** the Fig 2 proteins are not one pathway. Only the arc operon, and weakly
  betaine and ilv–leu, show partners that move together.

### `fig_string_official_mval.png`: STRING's own picture of the same 30 nodes
- **Compares:** nothing from our data; it is the STRING website rendering (evidence view,
  coloured lines = evidence types) for comparison and for the poster.
- **How to read:** STRING names, e.g. `argF-2` = our argF hit B7H15_14670, `ABD31928.1` =
  B7H15_14515, `ABD30224.1` = B7H15_06000. The background is transparent (looks black in
  dark viewers).

### `fig_kdp_operon_raw_mval.png`: kdp proteins per sample, 9 groups
- **Compares:** kdpC, kdpB, kdpA, kdpD, kdpE (KdpF = B7H15_11545 not detected) in
  TSB / SNpulse / SN for `anc`, `stk`, `stp` (`plot_protein()`).
- **Result:** all kdp proteins are higher in `stk` in TSB and SNpulse (the genomic block),
  not in `stp`. No change after the 15-min pulse. In SN: many missing values, KdpE down in
  all strains, KdpD 2/3 · 2/3 · 1/3.
- **Caveat:** SN values above ~18 for KdpA/KdpB come from 1–2 peptides in single samples.

### `fig_kdpD_peptides_mval.png`: is the kdpD response real at peptide level?
- **Compares:** (A) peptides detected per sample for kdpB, kdpD, kdpE, TSB vs SN;
  (B) each kdpD / kdpE peptide seen in SN against the **same peptide** in TSB of the same
  strain (log2 difference). × = protein-level response used in Fig 2.
- **Result:** KdpD drops from 12–27 peptides per sample in TSB to 0–3 in SN in every
  strain; 5 % of peptide observations remain (rank 9 of 1,121 proteins). The `anc`
  protein value (+1.03) comes from 2 rising peptides (one only seen in one TSB sample)
  next to 2 falling ones. All `stk` and all KdpE peptides fall.
- **Conclusion:** the kdp system is switched off by PASN in **all** strains. The
  `stk`-specific kdpD interaction is not supported at peptide level.

## Tables

| file | content |
|---|---|
| `response_hits_function_mval.csv` | the 7 proteins: gene, product, UniProt, SAOUHSC, UniProt name / EC / keywords, response per genotype, where they are hits (delta, FDR) |
| `string_partners_mval.tsv` | STRING `interaction_partners` answer (score ≥ 0.7, 6 per protein), raw |
| `string_edges_mval.tsv` | STRING `network` answer among all 30 nodes (score ≥ 0.7), raw, with evidence channels |
| `string_version_mval.txt` | STRING version, species and query date |
| `string_nodes_mval.csv` | each node: SAOUHSC, JE2 locus tag, gene, product, samples detected, responses, delta and FDR for `stk` and `stp` vs `anc` |
| `kdp_contrasts_mval.csv` | kdp proteins: SN responses, TSB baselines and interactions from our limma fit (diff, SE, FDR, tested) |
| `kdpD_peptide_matched_mval.csv` | every kdpD / kdpE peptide seen in SN with the mean of the same peptide in TSB (same genotype) and the difference |
| `checks_mval.txt` | CHECK lines (counts, mapping, peptide fractions, depth) |
