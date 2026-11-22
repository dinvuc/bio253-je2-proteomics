# 2026-11-17 – Steps 10–12: pathways, robustness, poster and report (group mval)

## Question for today
- 10: is Stk1 "broadly involved" (pathway level), and does `stp` share the `stk` pathways?
- 11: which headline conclusions survive leave-one-out, other thresholds and a check
  for genomic regions?
- 12: poster figures, abstract, report and project notebooks.

## What we did
- branch `steps10-12-final-mval`; one commit per step
- `scripts/10_enrichment_mval.R` on 6 contrasts (`stk`/`stp` vs `anc` in TSB, SN, response);
  `scripts/10b_step10_summary_mval.R` (cameraPR check, counts, candidate table, phospho ORA)
- `scripts/11_robustness_mval.R` (leave-one-sample-out ×18, thresholds, genome scan)
- `scripts/12_poster_figures_mval.R` (4 figures, 300 dpi); abstract + poster text;
  `reports/project_c_report_mval.qmd` (rendered, 24 s); three project notebooks in
  `analyses/` (each < 15 s)
- References checked in PubMed: Olsen 2006, Wu & Smyth 2012, Ritchie 2015, McCarthy & Smyth 2009
- all checks PASS; `Rscript tests/run_tests.R` passes

## Results (key numbers)
- **Pathways (robust = rank test and cameraPR both FDR < 0.05):** `stk` 35 / 143 / 44,
  `stp` 192 / 168 / 8 terms (TSB / SN / response). "Broad" shifts are **not** specific to
  `stk`. `stp` shares peptidoglycan synthesis (higher) and cytolysis (lower) with `stk`.
- The top `stk` "Ser/Thr phosphatase" STRING cluster is **not** Stp1/Stk1: it is the SigB
  operon + `mazEF` (B7H15_11465–11525), all +0.7 to +1.75 log2 in `stk`.
- **Genome scan:** `stk` block B7H15_11470–11590 (SigB, `mazEF`, `kdpCBA-DE`, `murF`, `ddl`)
  +1.51 log2, top window in the genome (family-wise p < 0.001); `stp` prophage block
  B7H15_10910–11150 −2.08 log2 (p < 0.001). 14 of 64 H1–H3 hits lie in such blocks
  (incl. kdpD and rbsD, two of the 5 `stk` response hits).
- **Leave-one-out:** 34 of 64 hits stable in 18/18 refits. Stk1 drop in `stp` 18/18;
  argF, B7H15_14515 18/18; ilvA (`stp`) 9/18.
- **Thresholds:** 1.5-fold triples the lists (e.g. `stk` response 5 → 15, `stp` 1 → 5),
  same shape.
- **Hand check:** Stk1 `stp_vs_anc_at_TSB` = −5.42 from group means = model (SE 0.58,
  FDR 0.0013), higher in `anc`, ~43-fold.

## Problems / surprises
- **Wrong numbers reported mid-session, corrected:** I read the term-count table in the
  wrong row order (it was sorted alphabetically). Correct: `stp` has most shifts in TSB,
  `stk` only leads in the response.
- **Fig 4b label bug:** term size in the label split one term into several rows (and NA).
  Fixed by labelling with the term ID; nested terms removed by overlap coefficient.
- **Broken null model:** the first genome scan used circular shifts, which keep every
  neighbourhood intact; nothing could ever be flagged (p = 1 everywhere). Replaced by
  random permutations.
- Interaction contrasts were printed as "higher in interaction_stk" by the enrichment
  tool: fixed.
- My Step 11 table first called the Stp1 reduction "unresolved"; Step 3 had called it a
  supported reduction (0–2 vs 9–15 peptides). Fixed: reduced = stable, size unresolved.
- Some H3 numbers in my first draft of the Step 11 table were wrong (3 instead of 6
  direct hits in blocks): fixed against the CSV.

## AI use
Claude Code wrote and ran all scripts, proposed the Step 10/11 rules (marked as checked by the group on 22 Nov 2026), drafted the abstract, poster text and report, and found the errors above by
checking figures and CSVs. **The numbers in the abstract and poster text were marked as checked by the group on 22 Nov 2026.** Biological interpretation is limited
to annotation-table facts and our evidence review; the "copy number vs regulation" reading
of the genome blocks is a hypothesis.

## Next steps
- Rules: marked as checked by the group on 22 Nov 2026. Group: poster layout + QR code; presentation; submission requirements.
- Ask Jonas: genome data for NE217/NE1919? (`stk` block, prophage); "stk1 broadly involved".
- Optional: STRING pictures with the written lists and background.
