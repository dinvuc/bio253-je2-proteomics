# 2026-11-12 – formal statistical test, step 5 (group mval)

## Question for today
Make our hit definition a real test instead of "FDR < 0.1, then filter
|log2FC| > 1", following the lecture's guardrail slide (explicit thresholds,
explicit model formula, guardrail checks, script + Quarto report).

## Decision (made BEFORE seeing results)
- Hit = `limma::treat(lfc = 1)` (H0: |log2FC| ≤ 1), BH-FDR < 0.05 per contrast.
- This is now the primary hit list. The supplied prolfquapp screen
  (FDR < 0.1 & |diff| > 1) is kept as a labelled comparison only.
- Rules section 5, plan and workflow step 5 updated accordingly.

## What we did
- branch `stats-formal-treat-mval`
- new helpers in `R/helpers_mval.R`: `abundance_matrix()` (join by sample name),
  `mval_contrasts()` (all 33 contrasts written out as A - B), `fit_limma_treat()`;
  constants `MVAL_FDR = 0.05`, `MVAL_LFC = 1`, `MVAL_MIN_OBS = 2`
- model `~ 0 + group`, 27 samples, n = 3; a protein is tested in a contrast only
  if every group of it has ≥ 2 observed values; no imputation
- no batch term: all 27 total-proteome runs on 20260612, `reinjected` all FALSE;
  run order is blocked by group, so drift is confounded (limitation)
- ran `Rscript scripts/05_formal_test_mval.R` -> all 7 checks PASS;
  `quarto render analyses/05_formal_test_mval.qmd` (16 s)
- step 10 enrichment now uses the formal hits; re-ran it for `stk_vs_anc`
- 2 new unit tests (toy data: sign, treat threshold, min_obs); all tests pass

## Results (key numbers)
- 2232 JE2 proteins; 1853–2103 tested per contrast
- sign check: cor(our diff, supplied diff) ≥ 0.971 in all 33 contrasts (observed rows)
- formal hits vs supplied screen, e.g. `stk_vs_anc_at_TSB` 20 (7 higher in stk,
  13 higher in anc) vs 134; `stp_vs_anc_at_TSB` 14 vs 161; `TSB_vs_SN_at_anc`
  204 vs 445; `interaction_stp_vs_anc_at_TSB_vs_SN` 0 vs 90
- Stk1, `stp_vs_anc_at_TSB`: diff = -5.42, FDR = 0.0013 -> lower in stp
  (polar effect, stp strain lacks most of the kinase)
- Stk1 not observed at all in stk, Stp1 observed in only 1/3 stp TSB samples:
  not testable (presence/absence evidence only, see supplied `lod_imputed`)
- step 10 `stk_vs_anc`: 17 formal hits; cytolysis / cell killing / extracellular
  region still enriched among proteins higher in anc (= lower in stk)

## Problems / surprises
- far fewer hits than the screen: expected, the screen never tested |log2FC| > 1
- "Partial NA coefficients for 240 probe(s)" warning from limma = proteins
  absent from at least one whole group; they are marked `tested = FALSE`
- lod_imputed (presence/absence) proteins drop out of the formal hit list

## AI use
Claude Code: wrote helpers, script, tests and report, checked the treat
reference in PubMed (PMID 19176553). The group chose the thresholds.
A group member still has to check every number above.

## Next steps
- group review of the PR; decide if steps 6/7 report the presence/absence list
  next to the formal hits
- re-run step 10 for the primary 18-sample contrasts
