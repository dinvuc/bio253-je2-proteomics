# 2026-11-13 – Steps 8–9: phosphosite coverage and site-vs-protein matching (group mval)

## Question for today
- 8: is the phosphosite layer measured often enough and localised well enough to
  support a follow-up? Which search (STY or STYH) do we use?
- 9: for the supported site changes, does the site follow its protein (same contrast),
  or is a change left over once protein amount is taken out?

## What we did
- branch `steps8-9-phospho-mval`; one commit per step
- **Rules written into the script headers BEFORE any site was selected** (AI proposal,
  marked as checked by the group on 13 Nov 2026, rules section 10):
  - eligible = 2 of 3 per group (Step 5 rule); localised = median per-sample
    `loc_prob` ≥ 0.75 (class I cutoff, Olsen et al. 2006 *Cell*, PMID 17081983,
    doi:10.1016/j.cell.2006.09.026);
    pH exploratory only; primary search = more supported sites
  - Step 9 categories: consistent (|site − protein| < 1), residual (≥ 1 and ≥ 2 SE),
    unresolved; SE assumes independent layers; no p-value/FDR for the adjusted effect
- `C_mutant_validation_DinMete/scripts/08_phospho_coverage_mval.R` → `results/step08_phospho_coverage/`
  (audit, depth, localisation, coverage funnel, our site refit with treat lfc = 1,
  shortlist, STY vs STYH and reinjection sensitivity)
- `C_mutant_validation_DinMete/scripts/09_site_vs_protein_mval.R` → `results/step09_site_vs_protein/`
  (match by exact `protein_key`, adjusted effect, categories, offset sensitivity,
  paired per-sample plots)
- new helper `add_x_primary()` / `primary_x_axis()` in `R/helpers_mval.R` + test;
  `Rscript tests/run_tests.R` passes
- every results folder has a `figure_guide_mval.md`; all checks PASS

## Results (key numbers)
- **Coverage (STY, primary):** 3,774 JE2 sites; per primary contrast only 119–296
  eligible and **91–184 supported**. Shortlist: 236 S/T/Y sites on 133 proteins.
- **Search choice:** STY 1,159 vs STYH 1,150 supported site × contrast pairs → STY
  (near tie). Same site, same contrast: r = 0.96; 51 hits shared, 9 STY-only, 2 STYH-only.
- **Localisation:** on the sites detected in the primary samples, the best-of-27 score
  would call ~1.7× as many sites "localised" as our per-sample median (STY S/T/Y
  2,109 vs 1,272). 401 sites have a median of 0, 221 of them pass on the best score alone.
- **Site hits (our test, new BH family):** 62 site × contrast hits = 37 sites on 25 proteins.
  55 are within-genotype TSB vs SN, 6 genotype differences in TSB, 1 interaction
  (`stp_vs_stk`). **0 hits in `interaction_stk_vs_anc` and `interaction_stp_vs_anc`.**
- **Step 9 categories:** 17 consistent with protein amount, 40 residual, 5 unresolved.
  - `stk_vs_anc_at_TSB`, deoB T85 (`B7H15_00775`): site diff = −3.81 (lower in `stk`),
    protein +0.05 → adjusted −3.86 ± 0.68. `stk_TSB` 2/3 detected. T87 shares its peptides (near-identical values): not separate evidence.
  - `stp_vs_anc_at_TSB`, pyk T537 (`B7H15_09250`): site +1.99 (higher in `stp`),
    protein +0.14 → adjusted +1.84 ± 0.25.
  - `stp_vs_anc_at_TSB`, Stk1 T164 (`B7H15_06355`): site −4.22, protein −5.42 (both
    lower in `stp`) → adjusted +1.21 ± 0.61, unresolved: the site loss mostly follows
    the protein loss (fits Step 3, polar effect).
- **Layer offset:** median (site − protein) in TSB_vs_SN = −0.65 (`anc`), −0.43 (`stk`),
  −0.17 (`stp`): sites lean towards "higher in SN". Removing it: 5 residual calls out, 5 in.

## Problems / surprises
- **Bug caught by a figure:** `MEDIA` is ordered TSB, SNpulse, SN, so on a TSB/SN-only
  axis all SN points were silently dropped (only a ggplot warning). Fixed with
  `add_x_primary()`, plus a test.
- **Bug caught by reading the output:** in the offset sensitivity, a `summarise()` column
  reused the name of the per-site flag, so lost/gained lists were wrong (fda S25 listed as
  "gained", impossible). Renamed; the lists now match a hand check.
- The STY abundance table has no NA rows (observed values only): must be completed to all
  27 samples before counting "2 of 3".
- The supplied STY DEA calls every row `observed`, although most rows come from
  `Imputed_Mean_moderated`. Its `estimate_type` cannot be used to filter.
- In the total table `B7H15_13880` has a 6-frame pseudogene entry as well → we match by
  exact `protein_key`, not by locus tag.
- `B7H15_04940` T269: +10.8 log2 on two `stp_SN` values near the detection floor. We
  added a QC flag (within-group range > 3 log2, 3 hits) after seeing this. It changes no
  category or hit rule.
- In `anc` and `stk`, SN phospho samples detect ~2× more sites than TSB (`stp` only
  ~1.1×): missing-not-at-random, and a likely contributor to the layer offset (which is
  smallest in `stp`).

## Recheck (16 Nov)
An independent recalculation (separate code: plain group means and a direct limma call)
reproduced every count, estimate, category and the offset sensitivity exactly. It found
three wrong statements, now fixed in the scripts and documents:
1. "Best score roughly doubles the localised sites": compared different site sets.
   Correct: ~1.7× (2,109 vs 1,272 on the same sites).
2. "SN detects ~2× more sites": true for `anc` and `stk`, not for `stp` (673 vs 623).
3. Audit row "534 / 492 peptides carrying ≥ 2 sites" compared whole peptide-list
   strings. Correct count: 2,780 (STYH) / 2,445 (STY) sites share a peptide sequence
   with another site, i.e. about two thirds are not independent.
Also: the Olsen et al. 2006 reference was checked in PubMed (doi:10.1016/j.cell.2006.09.026);
the 0.75 "class I" definition itself is in the full text, not the abstract (marked as checked by the group on 13 Nov 2026).
The Step 9 figure's label placement is now seeded, so the figure is reproducible.
- The project plan (5 Nov) had put phospho out of scope; the workflow keeps it as a
  supporting layer. Result: it supports H1 (Stk1 site), adds nothing to H2.

## AI use
Claude Code (Opus 5.5) proposed the Step 8/9 rules, wrote and ran both scripts, made the
figures and guides, and found and fixed the two bugs above. The localisation rule,
the primary-search rule, the Step 9 categories and the numbers above:
**marked as checked by the group on 13 Nov 2026.** Gene symbols come from `read_protein_metadata()`
(`ortholog_quality == "good"`), not from the AI. No biological interpretation of the
site changes was made.

## Next steps
- Step 8/9 rules, unchanged: marked as checked by the group on 13 Nov 2026.
- Step 10 enrichment for the `stk` contrasts; Step 11 robustness (deoB T85 rests on
  2/3 `stk_TSB` values; T269 flagged).
- Poster: Fig 4 only as a small descriptive panel (e.g. Stk1 T164 + deoB T85 paired plots),
  if the group wants it.
