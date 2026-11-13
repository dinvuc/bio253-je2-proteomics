# Step 7 outputs: what each figure and table compares

Made by [`scripts/07_shared_distinct_mval.R`](../../../scripts/07_shared_distinct_mval.R).
Question (H3): does `stp` look like `stk`? The polar-effect idea (Steps 3–4: the
`stp` strain also lacks most Stk1) predicts that the `stp`-vs-`anc` effects lie near
the `stk`-vs-`anc` effects. Proteins far from that are candidate phosphatase-specific
effects. Rules: [step 5 methods sheet](../step05_formal_test/methods_sheet_mval.md).
No new test is run here. All hits are the Step 5 formal hits (treat lfc = 1, BH-FDR < 0.05).

## Two traps when reading a stp-vs-stk comparison

1. **Shared reference.** `stk_vs_anc` and `stp_vs_anc` subtract the *same* `anc`
   mean. Its noise moves both estimates the same way, so pure noise already gives a
   positive r. We compute this **noise-only r from the data**: each group mean of n
   observed values has noise variance s²/n, with s² the protein's pooled within-group
   variance (the variance the `~ 0 + group` model uses). A simulation with the same
   s² and n reproduces it (difference ≤ 0.004). It is 0.50 in both layers, but the
   simulated 95% range is wide (0.39–0.61 baseline, 0.37–0.63 response), because a few
   high-variance proteins dominate Pearson r. **So the observed r alone cannot be
   compared with "noise".** We use instead the **noise-corrected r of the true
   effects**: observed covariance minus the mean noise covariance, divided by the
   noise-free variances (method of moments; bootstrap 95% CI over proteins). It assumes
   equal variance across groups within a protein, as limma does.
   Sensitivity check: the *split-reference r* compares `stk` with one `anc` replicate
   and `stp` with the other two (independent noise, null r = 0). It is much noisier, so
   it is pulled towards 0: a lower bound, not an estimate.
2. **"Hit in one, not the other" is not a difference.** Only the direct contrasts
   `stp_vs_stk_at_TSB` / `interaction_stp_vs_stk_at_TSB_vs_SN` test whether the two
   mutants differ. A non-significant direct contrast does not show that they are equal (n = 3).

## Figures

### `fig_mutant_effects_scatter_mval.png`: Fig 3a, stp effect vs stk effect
- **Compares:** for every protein eligible in both mutants, the `stp` effect (y) against
  the `stk` effect (x), both relative to `anc`. Left: baseline in TSB
  (`stp_vs_anc_at_TSB` vs `stk_vs_anc_at_TSB`, > 0 = higher in the mutant). Right:
  PASN response (the two interactions, > 0 = SN response lower in the mutant than in `anc`).
- **How to read:** diagonal = same effect in both mutants. Colours = Step 5 hit
  classes; orange = direct `stp_vs_stk` hit (the mutants differ). The box gives n,
  observed Pearson r with its 95% CI, the noise-only r, the noise-corrected r of the
  true effects with its bootstrap 95% CI, and the split-reference r range (3 ways to
  pick the single `anc` replicate). The total-least-squares slope is in
  `similarity_stats_mval.csv` only (see caveat).
- **Result:**

  | layer | observed r | noise-only r (sim. 95% range) | **noise-corrected r (95% CI)** | split-reference r |
  |---|---|---|---|---|
  | baseline (TSB) | 0.41 | 0.50 (0.39–0.61) | **0.37 (0.29–0.45)** | 0.11–0.33 |
  | response (interaction) | 0.43 | 0.50 (0.37–0.63) | **0.38 (0.24–0.52)** | 0.03–0.15 |

  The true `stk` and `stp` effects are **moderately correlated**: clearly above 0 and
  far below 1. Part of the `stp` pattern is shared with `stk`, but much of it is not.
  Only 3 proteins are formal hits in both mutants (same direction, all baseline).
  20 proteins differ directly between `stp` and `stk` in TSB, and 4 in the response.
  The low count of shared hits also reflects the strict test (treat lfc = 1, n = 3).
- **Caveat:** the TLS slope (baseline 1.18, response 0.72) mostly reflects the relative
  spread of the two axes when r is moderate. Do not read it as "stp has x % of the stk
  effect". The correction assumes noise independent of the true effects and equal
  variance across groups within a protein. It is dominated by high-variance proteins,
  hence the wide CI.

### `fig_heatmap_mval.png`: Fig 3b, all 46 genotype-hit proteins
- **Compares:** observed group means of the 6 TSB/SN groups for every protein that is
  a formal hit in any baseline or interaction contrast of the two mutants vs `anc`, or
  in the direct `stp_vs_stk` contrasts.
- **How to read:** each row is centred on its own mean over the 6 groups (display
  only, no scaling); orange = above the protein's mean, blue = below. Dark grey = 0/3
  detected. Rows are clustered for display (NA set to 0 only for the ordering).
- **Check:** Stk1 (`B7H15_06355`) is dark grey in `stk_TSB`, `stk_SN`, `stp_SN` and
  low in `stp_TSB`, as Step 3 found.

### `fig_shared_distinct_raw_mval.png`: per-sample values, shared vs distinct
- **Compares:** all 27 samples for (a) the 3 proteins that are hits in both mutants in
  the same direction (`B7H15_01230` coagulase, `B7H15_02275`, `B7H15_13900`; all lower in
  both mutants than in `anc` in TSB) and (b) the 4 strongest direct `stp_vs_stk` hits
  (`B7H15_11645` SceD, `B7H15_11900`, `B7H15_06000`, `B7H15_14865` Lipase 1).
- **How to read:** one point = one sample; `n/3` = detected replicates.
- **Inspected 13 Nov:** all shown groups 3/3, replicates tight; the effects are visible
  in the raw values.

## Tables

| file | content |
|---|---|
| `data/derived/effect_pairs_step07_mval.tsv.gz` (local only, git-ignored) | per protein × layer (baseline / response): `stk`, `stp` and direct `stp_vs_stk` estimate, SE, treat FDR, eligibility, hit class |
| `class_counts_mval.csv` | counts per class: hit in both same/opposite direction, one mutant only, no hit; eligible in one mutant only; direct hits |
| `similarity_stats_mval.csv` | observed r + CI, TLS slope + bootstrap CI, noise share of the variance on each axis, noise-only r (formula and simulation with 95% range), noise-corrected r + bootstrap CI, split-reference r range |
| `direct_contrast_candidates_mval.csv` | every direct `stp_vs_stk` hit with both mutant effects and direction in words |
| `shared_distinct_checks_mval.txt` | PASS/FAIL checks |
