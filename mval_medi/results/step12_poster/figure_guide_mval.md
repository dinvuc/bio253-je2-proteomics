# Step 12 poster figures: what each figure compares

Made by [`scripts/12_poster_figures_mval.R`](../../../scripts/12_poster_figures_mval.R) from saved
results of Steps 3–11 (no new statistics), 300 dpi. All per-sample panels show the 18
primary samples (TSB, SN), n = 3 per group. Text for the poster: [`poster_text_mval.md`](poster_text_mval.md);
abstract: [`abstract_mval.md`](abstract_mval.md). Full write-up: [`reports/project_c_report_mval.html`](../../../reports/project_c_report_mval.html).

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
Gene context (Step 4, `step04_polar_effect/fig_gene_context_mval.png`, 150 dpi; re-run Step 4
with `dpi = 300` if it is used on the A0 poster) and the phospho panels (Steps 8–9).
