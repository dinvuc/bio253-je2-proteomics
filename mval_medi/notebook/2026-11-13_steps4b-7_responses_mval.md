# 2026-11-13 – Steps 4b, 5 (rest), 6, 7: transposon control, contrast rules, H2, H3 (group mval)

## Question for today
- 4b: is the transposon marker ErmB present in `stk`/`stp` and absent from `anc`?
- 5: fix contrast meanings, direction, eligibility and testing families before
  reading results.
- 6 (H2): does a mutation change the response to PASN, relative to `anc`?
- 7 (H3): does `stp` look like `stk` (polar loss of Stk1)?

## What we did
- branch `steps6-7-responses-mval`; one commit per step
- `scripts/04b_transposon_products_mval.R`: Jonas's new FragPipe 24 search incl. the
  pBursa/pFA545 proteins (`data/raw/FP24_transposonDB/`, local, not in the checksummed bundle)
- `scripts/05b_contrast_sheet_mval.R` → `results/step05_formal_test/methods_sheet_mval.md`
- `scripts/06_altered_responses_mval.R` → `results/step06_responses/`
- `scripts/07_shared_distinct_mval.R` → `results/step07_shared_distinct/`
- every results folder has a `figure_guide_mval.md`; all checks PASS

## Results (key numbers)
- **4b:** ErmB 0 spectra in all 9 `anc` samples, 11–26 in `stk`, 19–30 in `stp`
  (7–10 peptides, 68 % coverage). Prediction met. GFP was in the database but not
  detected (probability 0).
- **5:** direction verified for all 64,876 tested rows (limma diff = group-mean
  formula). Eligible per primary contrast: 1,876–2,065 of 2,232. Supplied FDR = BH
  per contrast over all 2,324 rows. **Finding:** the supplied `lod_imputed` is a
  per-protein label (240 proteins with any group 0/3 get it in all 33 contrasts).
- **6 (H2):** formal interaction hits `stk` vs `anc` = 5, `stp` vs `anc` = 1
  (supplied screen: 96 / 90). Mutant and `anc` SN responses correlate (r = 0.85 / 0.86).
  - `interaction_stk_vs_anc_at_TSB_vs_SN`, kdpD (`B7H15_11550`): diff = 5.13,
    SE 0.18, FDR 9.5e-11 → SN response lower in `stk` than in `anc`
  - `interaction_stp_vs_anc_at_TSB_vs_SN`, ilvA (`B7H15_11425`): diff = 2.36,
    SE 0.25, FDR 0.039 → SN response lower in `stp`. `stk_SN` only 1/3 detected, so
    our test could not score `stk`; the supplied `lm_impute` estimate for `stk` is
    3.01 (FDR 0.001), same direction → probably shared by both mutants (1 `stk` replicate)
  - argF, `B7H15_14515`: the `stk` interaction is mainly a TSB baseline difference;
    arcA, kdpD: `stk` differs from `anc` in both media with opposite sign; rbsD: mainly SN
- **7 (H3):** observed r of `stp` vs `stk` effects = 0.41 (baseline) / 0.43 (response).
  Shared `anc` noise alone gives r = 0.50 (simulated range 0.39–0.61 / 0.37–0.63), so
  the observed r cannot be judged alone. **Noise-corrected r of the true effects =
  0.37 (95% CI 0.29–0.45) / 0.38 (0.24–0.52): a moderate shared component.**
  Split-reference r (lower bound) 0.11–0.33 / 0.03–0.15. 3 proteins are hits in both
  mutants (same direction, all lower than `anc` in TSB); 20 + 4 direct `stp_vs_stk` hits.
- all estimates in 6–7 are observed values (our limma refit), not `lod_imputed`.
  On/off proteins are listed separately as presence/absence only.

- **NTML check (open since 12 Nov), done by hand by the group:** NE217 = `stk`
  (stk1::Tn), NE1919 = `stp` (stp1::Tn). Matches `R/constants.R`, CLAUDE.md and
  notebook 07. The insertion codons (46 / 189) were not re-read by hand.

## Problems / surprises
- **Corrected the same day:** the first Step 7 version compared the observed r with a
  fixed 0.5 ("noise alone") and called the similarity "weak". The simulated noise-only
  r turned out to vary widely (a few high-variance proteins dominate Pearson r), so
  that comparison was not valid. Solved with a noise-corrected r (method of moments
  from per-protein variances, bootstrap CI): 0.37 / 0.38, a moderate shared component.
- `stp` is only partly like `stk`, although Stk1 is nearly gone in both. Possible
  readings (not tested): residual Stk1 in `stp`, a real Stp1-specific effect, or other
  differences between the strains.
- a first draft called arcA "mostly a TSB baseline difference"; the computed split
  shows a difference in both media with opposite sign. Fixed in workflow and guide.
- the supplied `lod_imputed` label is per protein, not per contrast (see 5)
- the first 4b check wrongly assumed every searched protein appears in the prot.xml.
  Proteins without any match show only as decoys. Fixed.
- the full Step 7 per-protein table is kept local (`data/derived/`), like Step 5

## AI use
Claude Code wrote and ran the four scripts, figure guides, methods sheet and this
draft. It found the FP24 transposon folder in `data/raw/`. Gene symbols come from
`read_protein_metadata()` (ortholog_quality good), not from memory. No biological
interpretation was made. **A group member still has to check every number above.**

## Next steps
- group review of the branch / PR; agree on the H3 wording for the poster
  ("`stp` shares a moderate part of the `stk` pattern, r ≈ 0.37–0.38 after noise
  correction, plus 24 direct differences")
- ask Jonas: is FP24 the final transposon search, and may we cite it on the poster?
- Step 10 enrichment for the `stk` contrasts; Step 12 poster figures (Fig 2 = step 6
  scatter, Fig 3 = step 7 scatter + heatmap, Fig 1e = ErmB)
