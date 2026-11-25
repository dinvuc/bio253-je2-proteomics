# Step 12 poster figures: what each figure compares

Made by [`C_mutant_validation_DinMete/scripts/12_poster_figures_mval.R`](../../scripts/12_poster_figures_mval.R) from saved
results of Steps 3–11 (no new statistics), 300 dpi. All per-sample panels show the 18
primary samples (TSB, SN), n = 3 per group. Full write-up: [`reports/project_c_report_mval.html`](../../../reports/project_c_report_mval.html).

### `fig1_h1_enzymes_mval.png`: H1, what the mutants lack
- **a** Stk1 and Stp1 normalised abundance (log2) per sample. **×** at the bottom = not
  detected (no value; never drawn as 0 in log space).
- **b** Peptides quantified per sample (0 = protein not quantified).
- **c** ErmB spectral counts (transposon marker, separate search incl. transposon proteins, Step 4b).
- **Read:** Stk1 absent in `stk`, 1–2 peptides (vs 29–32) in `stp_TSB`, absent in `stp_SN`;
  Stp1 0–2 peptides in `stp` (vs 9–15). Remaining Stp1 values in `stp` are erratic (1–2 peptides).
  ErmB: 0 spectra in every `anc` sample, present in every mutant sample.
- **Conclusion (H1 supported):** `stk` lacks Stk1. `stp` lacks Stp1 and most of Stk1, so it
  is a phosphatase mutant with strongly reduced kinase. The transposon marker confirms
  both insertions.

### `fig2_h2_responses_mval.png`: H2, PASN response relative to `anc`
- **Compares:** per protein, the `anc` response (x, log2 TSB − SN; > 0 = lower in SN) with the
  mutant response (y). Only proteins tested in the interaction contrast (1,882 `stk`, 1,931 `stp`).
- **Colours:** formal interaction hits (treat lfc = 1, BH-FDR < 0.05), coloured by Step 11
  robustness: green = stable in 18/18 leave-one-out refits, orange = fragile, pink = lies in a
  genomic block that shifts as a unit (kdpD, rbsD).
- **Read:** points on the diagonal = same response as `anc` (r = 0.85 / 0.86).
- **Key points:** 5 altered responses in `stk` (argF and B7H15_14515 stable; arcA fragile;
  kdpD and rbsD in shifted blocks) and 1 in `stp` (ilvA, fragile).
- **Caption caveat:** kdpD also fails the peptide check of Step 6d (*kdp* is switched off by
  PASN in all three strains, 0–3 peptides left per SN sample), so do not caption it as an
  `stk`-specific response.
- **Conclusion (H2 supported, but small):** both mutants keep the PASN response of the
  parent almost completely. Only a handful of proteins respond differently, and only two
  of them are fully robust.

### `fig3_h3_shared_distinct_mval.png`: H3, `stp` effects vs `stk` effects
- **a** Baseline in TSB (mutant − `anc`), **b** PASN response (interaction mutant vs `anc`).
  Proteins eligible in both mutants (1,988 / 1,853).
- **Colours:** blue = hit in both mutants, same direction; orange = direct `stp_vs_stk` hit.
- **Read:** the box gives the noise-corrected r (0.37 / 0.38) next to the r expected from the
  shared `anc` noise alone (0.50). Do not judge the cloud by eye: the shared reference tilts it.
- **Key points:** 3 proteins are altered in both mutants in the same direction (coagulase,
  `B7H15_02275`, fibronectin-binding protein; all lower than `anc` in TSB). 20 (baseline) and
  4 (response) proteins differ directly between `stp` and `stk`.
- **Conclusion (H3 partly supported):** `stp` shares a moderate part of the `stk` pattern,
  consistent with its reduced Stk1, but most of its pattern is its own. `stp` is not
  simply a kinase mutant.

### `fig3_h3_shared_distinct_compact_mval.png`: Fig 3, narrow version for the poster
- **Compares** the same as Fig 3 (`stp` vs `anc` against `stk` vs `anc`; a baseline in TSB, b PASN
  response), with the same proteins, colours and r (7 x 6.15 in instead of 12 x 7.5 in).
- **Axes are swapped:** `stp` is on x and `stk` on y, so above the diagonal = more positive in `stk`.
  a: > 0 = higher in the mutant. b: > 0 = response lower in the mutant.
- **Named:** coloured proteins with a gene symbol, the 3 shared ones and the 4 response hits. The
  noise-only r (0.50) is left out, because next to 0.37 it was misread (notebook 18 Nov).
- `fig3_h3_shared_distinct_compact_transparent_mval.png` **compares** exactly the same, without
  background (for the poster). Gene names keep a thin white outline so they stay readable.

### `fig4_pathways_genome_mval.png`: pathway level and genome position
- **a** Robust shifted functional terms (Wilcoxon rank test and cameraPR both FDR < 0.05, Step 10)
  per mutant and layer; above 0 = higher in the mutant / SN response lower in the mutant.
  Terms overlap, so bars are not counts of independent pathways.
- **b** Mean effect of 10 neighbouring proteins along the genome (`stk` and `stp` vs `anc`, TSB);
  dashed = genome-wide 5 % threshold from 1000 permutations (Step 11); pink = *stp1/stk1* locus.
- **Read:** two blocks change as a unit: SigB/*mazEF*/*kdp* in `stk` (+1.5) and a prophage
  in `stp` (−2.1). Co-regulation or a strain difference (copy number, prophage); DNA data needed.
- **Key points:** robust shifted terms `stk` 35 / 143 / 44, `stp` 192 / 168 / 8 (TSB / SN /
  response). The `stk` block is the strongest window in the genome (family-wise p < 0.001).
- **Conclusion:** at pathway level both mutants change broadly, `stp` at least as much as
  `stk`, so "Stk1 broadly involved" is not specific to the kinase mutant. Part of each
  mutant's signal comes from one genomic region, which needs DNA-level verification before
  biological interpretation.

### Not on the poster as a separate figure
Gene context (Step 4, `step04_polar_effect/fig_gene_context_mval.png`, 150 dpi; replaced on the
poster by the two-gene sketch `poster/p0_operon_mval`) and the phospho panels (Steps 8–9).

## `poster/`: the A0 poster versions (Step 12b, 18 Nov)

Made by [`C_mutant_validation_DinMete/scripts/12b_poster_a0_figures_mval.R`](../../scripts/12b_poster_a0_figures_mval.R) . Same data and
numbers as `fig1`–`fig4` above, redrawn for an A0 poster:
drawn at the printed column width (243 mm), Arial ≥ 22–24 pt, no gridlines, direct labels instead
of legends, no figure titles (the poster heading carries the message), vector PDF + 300 dpi PNG.

| file | poster section | compares | differences from the report figure |
|---|---|---|---|
| `p0_operon_mval` | 1 | *stp1* and *stk1* drawn to scale (247 / 664 aa) with the two insertion sites (codon 189 in *stp1*, codon 46 in *stk1*), from `TRANSPOSON_INSERTIONS` | new; a two-gene sketch instead of the full Step 4 gene context |
| `p1_h1_enzymes_mval` | 2 | as `fig1`: a abundance, b peptides, c ErmB, 18 samples | protein in columns, medium in rows; strains on the x axis (no legend); "n.d." row for not detected; the Stk1 fold change `stp` vs `anc` in TSB (log2 −5.42, 43× lower; computed in the script) written on the panel |
| `p2_h2_responses_mval` | 3 | as `fig2`: mutant vs `anc` PASN response, per protein | panels stacked; hits filled = robust (stable in 18/18 refits and not in a shifted block), open = fragile or in a block; key written in the panel |
| `p3_h3_shared_distinct_mval` | 4 | as `fig3`: `stp` vs `stk` effects, baseline and response | panels stacked, each with a text column on the right: noise-corrected r **with its bootstrap 95 % CI and the words "corrected for the shared anc noise"** (18 Nov: "r = 0.37" alone invited a comparison with the noise-only r = 0.50; the corrected r must be compared with 0), key with counts, and the 3 "altered in both" proteins numbered 1–3 in the panel and named in the key; direction of each axis written under the panel title; axes ±5.6 so no protein is clipped; 20 baseline direct hits counted, not named |
| `p4_pathways_genome_mval` | 5 | a total robust terms per mutant and layer; b genome scan | a: totals only (no up/down split), horizontal bars; b: filled area instead of a thin line, labels on the panels |
| `qr_repo_mval.svg` | fine print | link to the repository | new |
| `qr_references_mval.{svg,png}` (Step 12d, `C_mutant_validation_DinMete/scripts/12d_qr_references_mval.R`) | fine print | link to [`references_mval.md`](../../docs/references_mval.md) on `main` (works once merged) | new |
