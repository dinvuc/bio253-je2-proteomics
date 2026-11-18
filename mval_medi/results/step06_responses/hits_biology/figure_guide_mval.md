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
- **Data:** STRING v12.5, the course's own JE2 organism `STRG0B09NBK` ("JE2_Bio253_HS26",
  the same genome and protein ids as our FASTA). Links with combined score ≥ 0.7. Network =
  the 7 proteins + up to 6 best partners each (37 proteins, 70 links). STRING's protein id
  is our `protein_key`, so every node joins directly to our data (no ortholog mapping).
- **How to read:** square = Fig 2 protein, circle = STRING partner, thick outline =
  formal interaction hit in that panel, grey = not testable (< 2/3 replicates in a group),
  line width = STRING score. Label = gene (or short product) + locus number. Layout: one
  circle per **hit module** (Fig 2 proteins linked to each other, here only arcA + argF,
  in the middle; each partner sits with the hit it scores highest with). Solid lines =
  within a module, dashed = between modules (6 of 70). Positions are arbitrary.
- **Result:** 6 modules; no direct link between Fig 2 proteins except arcA–argF. The
  dashed links run through shared partners (aspartate aminotransferase 14470,
  B7H15_11190, lacD 12200 / lactose transporter 12195). Partners that move with their hit
  in `stk`: betA (B7H15_14515), acetolactate synthase 11390 (ilvA), Y1819 (kdpD). Partners
  that do not: the arc / arginine partners, rbsK and the LacI regulators. The kdp partners
  are grey (too few SN values).
- **Caveat:** STRING scores for an uploaded genome are transferred from related organisms
  (text-mining, co-expression, databases), so a link means "functionally associated", not
  a physical contact. The partner list differs slightly from STRING's NCTC 8325 entry.
- **Conclusion:** the Fig 2 proteins are not one pathway. Only the arc operon, and weakly
  betaine and ilv–leu, show partners that move together.

### `fig_string_official_mval.png`: STRING's own picture of the same 37 proteins
- **Compares:** nothing from our data; it is the STRING website rendering for the course
  JE2 organism (evidence view, coloured lines = evidence types), for comparison and for
  the poster.
- **How to read:** nodes are labelled with JE2 locus tags (B7H15_…). The background is
  transparent (looks black in dark viewers).

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
| `response_hits_function_mval.csv` | the 7 proteins: protein_key, gene, product, UniProt, in STRING, UniProt name / EC / keywords, response per genotype, where they are hits (delta, FDR) |
| `string_partners_mval.tsv` | STRING `interaction_partners` answer, JE2 organism (score ≥ 0.7, 6 per protein), raw |
| `string_edges_mval.tsv` | STRING `network` answer among all 37 proteins (score ≥ 0.7), raw, with evidence channels |
| `string_version_mval.txt` | STRING version, organism and query date |
| `string_nodes_mval.csv` | each node: protein_key (= STRING id), JE2 locus tag, gene, product, samples detected, responses, delta and FDR for `stk` and `stp` vs `anc` |
| `kdp_contrasts_mval.csv` | kdp proteins: SN responses, TSB baselines and interactions from our limma fit (diff, SE, FDR, tested) |
| `kdpD_peptide_matched_mval.csv` | every kdpD / kdpE peptide seen in SN with the mean of the same peptide in TSB (same genotype) and the difference |
| `checks_mval.txt` | CHECK lines (counts, mapping, peptide fractions, depth) |
