# Selected BIO253 project: mutant validation, polar effects and exposure responses

Status: chosen project, agreed 5 November 2026. This replaces the earlier provisional topic-H plan. Protein-abundance adjustment remains a supporting method for selected phosphosites.

## Working title

Does the stp mutant behave like a kinase mutant during Pseudomonas-supernatant exposure?

## Main question

Does the stp1 transposon mutant have reduced Stk1 as well as Stp1, and do its protein-abundance and phosphorylation responses to supernatant resemble those of the stk1 mutant?

We will establish the measured mutant profiles, then use genotype-by-condition interactions to distinguish baseline differences from altered exposure responses. A polar effect is a proposed explanation for reduced Stk1, not a mechanism established by proteomics alone.

## Predictions recorded before analysis

1. The stp mutant shows reduced Stp1 and accompanying reduction of Stk1 relative to the parental strain. The consistency and magnitude must be checked in individual samples.
2. Some exposure responses in stp and stk are similar, consistent with shared consequences of reduced kinase function. Other responses may distinguish them.
3. Similarity and differences may depend on treatment. A failure to reject a difference will not be called equivalence.

Alternative explanations include altered protein stability or regulation, additional mutations, sample problems and changes driven by protein abundance rather than phosphorylation. We will retain contradictory observations.

## Scope and design

- Organism: S. aureus JE2; reference anc, stk1::Tn and stp1::Tn.
- Primary comparison: TSB versus SN, with three samples per genotype-condition combination.
- SN: approximately 15 h growth with 30% cell-free, sterile-filtered P. aeruginosa PAO1 spent medium. This is exposure to spent medium, not co-culture; secreted products and nutrient depletion are not separated by the design.
- Optional extension: SNpulse, TSB-grown cells followed by 15 min exposure. This is a different growth history, not a longitudinal recovery experiment. Two stk pulse phospho measurements were re-injected; investigate batch effects before interpreting them.
- Layers: total proteome and phosphoproteome. STY and STYH are alternative searches, not biological replicates. Record the search and statistical model used.

## Analysis sequence

### 1. Validate the mutant profiles

Use the repository's R helpers and the existing mutant-validation notebook as the starting point. Plot Stp1 (B7H15_06350) and Stk1 (B7H15_06355) in every sample. Check missingness, outliers and consistency across treatments. Examine nearby genes as supporting context, without assuming every neighboring gene is cotranscribed.

Not detected is not proof of absence. Estimates marked lod_imputed support detection-pattern descriptions rather than precise loss magnitudes. Investigate unexpected high mutant signals with the lecturer; do not silently exclude replicates.

### 2. Separate baseline differences from exposure-response differences

First examine stk_vs_anc_at_TSB and stp_vs_anc_at_TSB. Then evaluate:

- interaction_stk_vs_anc_at_TSB_vs_SN
- interaction_stp_vs_anc_at_TSB_vs_SN
- interaction_stp_vs_stk_at_TSB_vs_SN

The first two assess each mutant's response relative to the parental response; the third directly compares the two mutant responses. Compare effect sizes and uncertainty, not just overlaps between significant-hit lists. Because the first two comparisons share the parental reference, their scatterplot correlation is descriptive and not independent proof of a shared mechanism.

For any feature, the stp-versus-anc interaction is:

`(stp_TSB - stp_SN) - (anc_TSB - anc_SN)` on the log2 scale.

Positive means the TSB-minus-SN difference is more positive in stp. Equivalently, its SN-minus-TSB response is less positive or more negative than the parental response. Use explain_contrast() and translate selected results into plain words.

### 3. Examine supported phosphorylation candidates

Start with reproducibly measured features, initially requiring detection in at least two of three samples in each group needed for a quantitative comparison. Inspect localization confidence and individual values. Keep consistently one-sided detections in a separate exploratory analysis rather than treating missing values as zeros.

For selected sites, inspect the corresponding total protein under the same contrast. Adjust only supported observed measurements; do not subtract imputed effects. Propagate uncertainty with stated assumptions if formal adjusted tests are used. Adjustment does not measure absolute phosphorylation occupancy.

### 4. Interpret shared and distinct responses

Group supported candidates by function and interpret against literature and measured-protein coverage. Use direct stp-versus-stk contrasts to investigate differences. Do not label stp-specific effects as direct phosphatase substrates. If enrichment is performed, use the eligible measured background and control multiple testing.

Course default screens use FDR < 0.1 and absolute log2 effect > 1; report estimates and uncertainty alongside the screen. With three replicates per group, absence of significance is weak evidence of no effect.

## Rescue and suppression: literature-guided interpretation

See [the rescue evidence review](docs/stk_stp_rescue_evidence.md). Published double-mutant studies make suppression of some phosphatase-loss effects by kinase loss plausible. Our stp transposon strain is not a verified clean double deletion, and the dataset lacks a clean stp-only deletion and complemented rescue strains.

Therefore, this project can identify molecular patterns consistent with combined disruption and possible suppression. It cannot demonstrate genetic rescue, prove a transcriptional polar effect, or establish restored growth, survival or antibiotic tolerance. Confirmation would require additional strains and measurements.

## Anticipated figures

1. Gene arrangement and per-sample Stp1/Stk1 protein abundances.
2. Comparison of mutant interaction effects, including uncertainty and observed/imputed status.
3. Heatmap and per-sample plots of selected shared and distinct responses.
4. Matched protein/phosphosite plots for a small, supported biological module, if coverage permits.

## Outputs and first next step

Write R analysis and a rendered report, keep dated notebook entries, and publish only code and reports. Raw and derived course measurements stay in ignored data/.

First establish a working R environment, then validate Stp1/Stk1 measurements before evaluating downstream hypotheses. No biological results have yet been generated for this project.
