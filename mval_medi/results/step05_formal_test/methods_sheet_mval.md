# Step 5 methods sheet: contrasts, direction, eligibility, model, thresholds

Written 13 Nov 2026 from the outputs of
[`scripts/05_formal_test_mval.R`](../../../scripts/05_formal_test_mval.R) (the test) and
[`scripts/05b_contrast_sheet_mval.R`](../../../scripts/05b_contrast_sheet_mval.R)
(definitions and checks). Every number below is in `contrast_sheet_mval.csv`
or `contrast_sheet_checks_mval.txt`. Thresholds were fixed by the group on
12 Nov 2026, before results were seen.

## 1. Data and model

| item | value |
|---|---|
| features | 2,232 JE2 proteins (`keep_je2()`: has a locus tag, no `rev_` decoy) |
| values | `normalized_abundance` from `read_abundances("total")`, already log2 and normalised; not changed. NA = not detected, never 0 |
| samples | all 27 (3 genotypes × 3 media × n = 3), joined to the sample sheet by name |
| model | limma, `~ 0 + group` (one mean per genotype × medium, no intercept) |
| contrasts | all 33 course contrasts written out as `A - B` (`mval_contrasts()`) |
| test | `limma::treat(lfc = 1)`: H0 is \|log2FC\| ≤ 1 |
| hit | BH-FDR < 0.05, **within each contrast**, over the eligible proteins of that contrast |
| batch | none: one acquisition day; run order blocked by group (confounded, not modellable) |

## 2. Contrast families (primary = the 18 TSB/SN samples)

| family | primary contrasts | meaning, sign |
|---|---|---|
| baseline: genotype in TSB | `stk_vs_anc_at_TSB`, `stp_vs_anc_at_TSB`, `stp_vs_stk_at_TSB` | diff > 0 = higher in the first-named strain |
| exposure within a genotype | `TSB_vs_SN_at_anc`, `TSB_vs_SN_at_stk`, `TSB_vs_SN_at_stp` | diff > 0 = higher in TSB = **lower in SN** |
| interaction | `interaction_stk_vs_anc_at_TSB_vs_SN`, `interaction_stp_vs_anc_at_TSB_vs_SN`, `interaction_stp_vs_stk_at_TSB_vs_SN` | (X_TSB − X_SN) − (Y_TSB − Y_SN). Positive = the SN response (SN minus TSB) is **less positive / more negative in X** than in Y |

The other 24 contrasts (SNpulse, averaged over media or genotypes) are listed in
`contrast_sheet_mval.csv` with their family; they are not used for H2/H3.

**Reminder for every `stp` contrast:** the `stp` strain also lacks most of Stk1
(Step 3: median 0–1 Stk1 peptides vs 28–32 in `anc`). `stp_vs_anc` is therefore
"phosphatase lost + kinase mostly lost", not phosphatase alone.

## 3. Direction check

- For every one of the 64,876 tested protein × contrast rows, the limma estimate
  equals the contrast formula evaluated on plain group means of the observed values
  (max deviation 2e-14). So the sign is as the name says, in all 33 contrasts.
- Worked example, Stk1 (`B7H15_06355`), `stp_vs_anc_at_TSB`: mean stp_TSB = 10.98,
  mean anc_TSB = 16.40, difference = −5.42 = limma diff. **Stk1 is higher in anc,
  i.e. lower in stp.**

## 4. Eligibility (comparison-specific, fixed before selecting hits)

A protein is **eligible** for a contrast if it has ≥ 2 of 3 observed values in
**every** group that contrast uses (= `tested` in the formal test). Excluded
proteins are kept, in two classes:

- **on/off (one-sided):** ≥ 2/3 in one group, 0/3 in another: presence/absence
  evidence, listed in `one_sided_signals_mval.csv` (never given a quantitative fold change).
- **sparse:** anything else (some group 1/3, or no group with ≥ 2/3).

| contrast | eligible | on/off | sparse | formal hits (up / down) | supplied screen hits |
|---|---|---|---|---|---|
| `stk_vs_anc_at_TSB` | 2002 | 43 | 187 | 7 / 13 | 134 |
| `stp_vs_anc_at_TSB` | 2065 | 33 | 134 | 2 / 12 | 161 |
| `stp_vs_stk_at_TSB` | 2004 | 55 | 173 | 4 / 16 | 197 |
| `TSB_vs_SN_at_anc` | 2001 | 58 | 173 | 136 / 68 | 445 |
| `TSB_vs_SN_at_stk` | 1918 | 75 | 239 | 158 / 72 | 480 |
| `TSB_vs_SN_at_stp` | 1996 | 69 | 167 | 131 / 82 | 447 |
| `interaction_stk_vs_anc_at_TSB_vs_SN` | 1882 | 153 | 197 | 2 / 3 | 96 |
| `interaction_stp_vs_anc_at_TSB_vs_SN` | 1931 | 142 | 159 | 1 / 0 | 90 |
| `interaction_stp_vs_stk_at_TSB_vs_SN` | 1876 | 155 | 201 | 3 / 1 | 120 |

"up" = diff > 0 (higher in the first-named group, or for interactions a positive
difference of differences). Interactions use four groups, so more proteins fall
out of them (on/off 142–155 vs 33–75).

## 5. Supplied DEA: what its labels and FDR mean

- **Supplied FDR scope:** BH within each contrast over all 2,324 supplied rows,
  including contaminants and decoys (recomputed exactly). Model `lm_impute`.
- **Finding: `estimate_type` is a per-protein label.** If any of the 9 groups of a
  protein is 0/3, prolfquapp labels *all 33* contrasts of that protein
  `lod_imputed` (240 proteins), also contrasts whose own groups are fully
  observed. E.g. 99 of the 2,065 eligible proteins of `stp_vs_anc_at_TSB` carry the
  supplied `lod_imputed` label. For those rows our observed-only estimate is the
  one to use. "Observed rather than `lod_imputed`" in the workflow is therefore
  applied through our own test (which never imputes), not through the supplied label.
- The supplied screen (FDR < 0.1, then |diff| > 1) is kept as a **labelled
  comparison only**. Its p-values/FDR are never reused for a new quantity.

## 6. Testing families: what the FDR does and does not cover

- One family = one contrast. BH controls the false discovery rate *within* that
  contrast's hit list. There is no correction *across* the 9 primary contrasts.
- Combined statements in Steps 6–7 ("hit in both interactions", "hit in `stk` but
  not `stp`") are **descriptive classifications**, not tested claims; no FDR
  applies to them. A hit in one contrast and not in another does not show a
  difference: only the direct contrast (e.g. `interaction_stp_vs_stk_…`) tests that.
- Any new test in a later step gets its own BH correction across all features it tests.
- n = 3 per group: non-significance is weak evidence for "no effect".
