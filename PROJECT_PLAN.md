# Revised topic C: Does the stp mutant behave like a kinase mutant during supernatant exposure?

**Current research goal, updated 5 November 2026.** This focused subsection of course topic C uses total-proteome protein abundance in TSB and SN. Phosphorylation responses, phosphosite analysis and protein adjustment of phospho signals are outside the agreed scope. The hypotheses below are predictions to test, not results from our dataset.

## Research question

Does disruption of `stp1` also reduce Stk1 protein abundance, consistent with a possible polar effect, and do the `stp` and `stk` mutants show similar protein-abundance responses to *Pseudomonas aeruginosa* supernatant?

## Background

Stp1 removes phosphate groups from proteins, whereas Stk1 adds them. The neighboring `stp1` and `stk1` genes are described as an operon in the course material and literature reviewed in [the Stk1/Stp1 evidence review](docs/stk_stp_rescue_evidence.md). Disrupting `stp1` may also reduce expression of the neighboring `stk1` gene: a possible polar effect.

The `stp` mutant may therefore involve reduced function of both enzymes. Its behavior cannot automatically be explained by loss of the phosphatase alone. Total-proteome measurements can establish protein-abundance patterns consistent with this explanation; they cannot prove transcriptional polarity or enzyme activity.

## Hypothesis

The `stp` mutant has reduced amounts of both Stp1 and Stk1. Some of its protein-abundance responses to supernatant therefore resemble those of the `stk` mutant.

- **H1 - mutant profiles:** `stp` has reduced Stp1 and Stk1 relative to the parental strain. `stk` is expected to have reduced Stk1; Stp1 is checked rather than assumed unchanged.
- **H2 - altered exposure response:** some protein-abundance changes from TSB to SN differ between a mutant and the parental strain.
- **H3 - shared response component:** some proteins show same-direction deviations from the parental response in both mutants, consistent with a contribution from reduced Stk1 in `stp`. Other responses may distinguish the mutants.

## Part 1: Establish what is disrupted

Use the total-proteome dataset to examine Stp1 (`B7H15_06350`) and Stk1 (`B7H15_06355`) abundance in `anc`, `stk` and `stp`. Plot individual replicates and assess consistency across TSB and SN conditions. Inspect missingness and outliers; examine neighboring-gene products only as supporting context.

Not detected is not proof of absence. Flag `estimate_type == "lod_imputed"`: such estimates support detection-pattern descriptions rather than precise fold-change magnitudes. Keep missing values as NA and investigate unexpected mutant signals rather than silently excluding samples.

## Part 2: Test whether the mutants respond similarly

Use genotype-by-condition interaction contrasts to compare the change from TSB to SN:

- Does the `stk` response differ from the parental response? `interaction_stk_vs_anc_at_TSB_vs_SN`
- Does the `stp` response differ from the parental response? `interaction_stp_vs_anc_at_TSB_vs_SN`
- Do the mutant responses differ from each other? `interaction_stp_vs_stk_at_TSB_vs_SN`

First examine baseline differences with `stk_vs_anc_at_TSB` and `stp_vs_anc_at_TSB`. These distinguish existing mutant differences from altered exposure responses.

Compare interaction effect sizes and uncertainty rather than only lists of significant results. Plot individual protein measurements for selected shared or distinct responses. Prioritize reliably observed measurements and report missing replicates and imputed estimates. Any formal screen must account for multiple testing; the course default is FDR < 0.1 and absolute log2 effect > 1, reported alongside estimates and uncertainty.

The first two interactions share the parental reference, so their correlation is descriptive and cannot independently establish a shared mechanism. Use the direct `stp`-versus-`stk` interaction to investigate differences. A nonsignificant difference does not establish equivalence, especially with three replicates per group.

### Contrast direction

The supplied contrasts use `A_vs_B = log2(A) - log2(B)`. For example:

`interaction_stp_vs_anc_at_TSB_vs_SN = (stp_TSB - stp_SN) - (anc_TSB - anc_SN)`

A positive value means the TSB-minus-SN difference is more positive in `stp`. Equivalently, its SN-minus-TSB response is less positive or more negative than the parental response. Use `explain_contrast()` and state which group is higher when reporting abundance differences.

## Scope and design

- Organism: *S. aureus* JE2; parental `anc`, `stk1::Tn` and `stp1::Tn`.
- Conditions: TSB versus SN, three replicates in each of six genotype-condition groups (18 samples).
- SN: approximately 15 h growth in 70% TSB plus 30% cell-free, sterile-filtered *P. aeruginosa* PAO1 spent medium. Secreted products and nutrient depletion are not separated by this design.
- Dataset: total-proteome abundance and its differential-analysis tables, accessed with the repository's R helpers.
- SNpulse and the phosphoproteome remain available in the course dataset but are outside this focused project.

## Anticipated figures

1. Stp1 and Stk1 protein abundance across individual TSB and SN samples, grouped by genotype.
2. Comparison of the two mutants' protein-abundance interaction effects, with uncertainty and observed/imputed status made clear.
3. Response plots for selected proteins showing shared or distinct TSB-to-SN responses.

## Significance and interpretation limits

The project tests how specifically these mutants represent kinase or phosphatase disruption and how that affects interpretation of their protein-abundance responses to supernatant. Similar mutant responses would support a contribution from reduced Stk1 in `stp`, but would not prove that a polar effect causes those responses.

Published suppression studies provide context in [the evidence review](docs/stk_stp_rescue_evidence.md). This dataset lacks a clean Stp1-only deletion retaining Stk1, a verified double deletion and complemented strains. It cannot demonstrate causal genetic rescue, restored growth, survival or antibiotic tolerance. Alternative explanations include altered protein stability or regulation, additional mutations and sample problems; contradictory observations must be retained.

## Outputs and next step

Write R analysis and a rendered report; keep dated entries in `notebook/`. Publish code, documentation and reports only. Raw and derived course measurements stay in ignored `data/`.

Begin with the per-sample Stp1/Stk1 total-proteome plots in TSB and SN, then evaluate interaction effects. This scope update does not add biological results.
