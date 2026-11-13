# Step 6 outputs: what each figure and table compares

Made by [`scripts/06_altered_responses_mval.R`](../../../scripts/06_altered_responses_mval.R).
Question (H2): does a mutation change the response to PASN (SN vs TSB),
beyond a difference already present in TSB? Every response is shown relative to
the ancestor `anc`. Rules (model, eligibility, testing families):
[step 5 methods sheet](../step05_formal_test/methods_sheet_mval.md).

**Reading key.** A *response* is written as SN minus TSB:
`response_X = log2(X_SN) − log2(X_TSB)`, positive = higher in SN.
`interaction_X_vs_anc_at_TSB_vs_SN = response_anc − response_X`
(checked numerically for all 3,813 eligible rows), so a **positive interaction =
the SN response is lower (less up / more down) in the mutant than in `anc`**.
Every `stp` result: that strain also lacks most of Stk1 (Step 3).

## Figures

### `fig_response_vs_anc_mval.png`: Fig 2, mutant response vs ancestor response
- **Compares:** for each eligible protein, the SN response in `stk` (left) or `stp`
  (right) against the SN response in `anc`.
- **Data:** our limma fit (Step 5), observed values only, 27 samples, n = 3.
- **How to read:** one point = one protein. On the solid diagonal = same response
  as `anc`; dashed = diagonal ± 1 log2. Below the diagonal = response lower in the
  mutant (interaction > 0). Orange = formal interaction hit (treat lfc = 1,
  BH-FDR < 0.05). **Pink = exploratory only** (0.05 ≤ FDR < 0.1, see below). The box
  gives n eligible, Pearson r, the hit count and the extra exploratory count.
- **Result:** r = 0.85 (`stk`) and 0.86 (`stp`): most of the PASN response is
  shared with `anc`. Formal hits: 5 in `stk`, 1 in `stp`. Exploratory: +0 in `stk`,
  +1 in `stp` (`B7H15_06000`, fibrinogen-binding protein).
- **Caveat:** r is descriptive (no CI here); both axes contain noise, so r is lower
  than the true similarity. A point outside the dashed lines but grey is NOT evidence
  of a difference > 1 (the interaction SE is large, median 0.17).
- **About the pink layer (added 13 Nov on request):** our hit rule stays BH-FDR < 0.05,
  fixed on 12 Nov before results (rules section 5). FDR < 0.1 was chosen *after*
  seeing the results, so the pink proteins are **not hits**. They are not counted in
  any hit number and should be called "exploratory" on the poster. For `B7H15_06000`
  the raw values (`fig_candidates_raw_mval.png`) show that `stp` is low in TSB and
  normal in SN, i.e. mainly a TSB baseline difference.

### `fig_response_stk_vs_stp_mval.png`: stp response vs stk response (direct comparison)
- **Compares:** for each protein eligible in the direct comparison (1,876), the SN
  response in `stp` (y) against the SN response in `stk` (x). This is
  `interaction_stp_vs_stk_at_TSB_vs_SN` = response_stk − response_stp (checked for
  all 1,876 rows).
- **How to read:** same as Fig 2, but now the reference is `stk`, not `anc`. On the
  diagonal = both mutants react the same. Below = response lower in `stp` than in
  `stk` (interaction > 0); above = higher in `stp`. Orange = formal hit
  (FDR < 0.05); pink = exploratory (0.05 ≤ FDR < 0.1, not a hit).
- **Result:** r = 0.84: the two mutants mostly react the same way. 4 formal hits:
  `B7H15_06000` (above: rises more in `stp`), lipA, arcA, argF (below: lower response
  in `stp`). Exploratory +1: `B7H15_14515`.
- **Caveat:** no `anc` here, so this r has no shared-reference problem, but it still
  says only how similar the *responses* are; whether each mutant differs from `anc`
  is in Fig 2. Several hits are mainly TSB baseline differences (see the raw plots).
  The `stp` strain also lacks most Stk1. This figure belongs thematically to H3 (Step 7).

### `fig_effect_distributions_mval.png`: effect sizes, reference vs interactions
- **Compares:** the distribution of `TSB_vs_SN_at_anc` (the ancestor's own response)
  with the two mutant-vs-`anc` interactions, eligible proteins only.
- **How to read:** dashed lines at |log2| = 1, the treat threshold. Note the sign:
  `TSB_vs_SN_at_anc > 0` = higher in TSB = lower in SN.
- **Result:** median |diff| 0.45 for the `anc` response vs 0.21 (`stk`) and 0.22
  (`stp`) for the interactions. Responses that differ between strains are rare and small.

### `fig_candidates_raw_mval.png`: every formal interaction hit, per sample
- **Compares:** normalised log2 abundance of the 8 distinct proteins behind the 10
  interaction hits (incl. `stp_vs_stk`), all 27 samples (SNpulse for context only).
- **How to read:** one point = one sample; `n/3` above = detected replicates;
  missing samples are not drawn.
- **What the raw values show (checked by eye, 13 Nov):**
  - argF (`B7H15_14670`), `B7H15_14515`: `stk` differs from `anc` **in TSB only**
    (stk − anc = −1.63 and +1.79; in SN 0.74 and −0.28). The interaction is mainly a
    baseline difference that disappears in SN.
  - arcA (`B7H15_14675`), kdpD (`B7H15_11550`): `stk` differs from `anc` in **both** media,
    with opposite sign (arcA: −3.74 in TSB, +1.23 in SN; kdpD: +1.35 in TSB, −3.78 in SN).
  - kdpD: `stk_SN` lower than `anc_SN` (2/3 vs 2/3 detected); `stp_SN` 1/3, so not testable
    in `stp` (supplied estimate 0.60, FDR 0.60: no evidence there).
  - rbsD (`B7H15_01460`): high in `stk_SN` (2/3); `stp_SN` replicates spread from about 12.5 to 16.
  - ilvA (`B7H15_11425`, the only `stp` hit): `stk_SN` drops too, but only 1/3 detected,
    so our test could not score `stk`. The supplied model (which uses that one value)
    gives `interaction_stk_vs_anc` = 3.01, FDR 0.001, the same direction as `stp` (2.36).
    The lower SN response is therefore probably **shared by both mutants**; in `stk` it
    rests on one replicate.
  - arcA in `stk_SNpulse`: one replicate far above the other two.

## Tables

| file | content |
|---|---|
| `response_summary_mval.csv` | per primary contrast (baseline / exposure / interaction): eligible, hits up/down, median and 95th percentile \|diff\|, median SE, number of estimates with \|diff\| > 1 (descriptive, NOT a hit rule) |
| `interaction_candidates_mval.csv` | each formal interaction hit: diff, SE, treat FDR, model, `estimate_type` (always observed), responses in anc/stk/stp, strain difference in TSB and in SN, `where_difference` (descriptive, on estimates: mainly TSB baseline / mainly SN response / both), direction in words, `detected_TSB_SN` (replicates per group), `other_mutant` (is the same interaction a hit / testable in the other mutant; if not testable, the supplied `lm_impute` estimate as presence/absence evidence) |
| `one_sided_responses_mval.csv` | on/off proteins of the three interactions (exploratory), with `pattern_class`: *on/off in one strain only* (a possible altered response), *shared* (both strains empty in the same medium), *absent from a strain in both media* (not a response), *mixed* |
| `responses_checks_mval.txt` | PASS/FAIL checks |

The on/off list is presence/absence evidence only. It gets no fold change and no FDR,
and it is not ranked.
