# 2026-11-18 – kdpD added to the refined N-to-C plots (group mval)

## Question for today
What do the protein and phosphosite layers show for kdpD (`B7H15_11550`), the Step 6
interaction hit (`interaction_stk_vs_anc_at_TSB_vs_SN`, `fig_response_vs_anc_mval.png`),
when it is drawn in the same N-to-C style as the refined-hit proteins of Step 9c?

## What we did
- added kdpD as an extra protein (`EXTRA`) to `C_mutant_validation_DinMete/scripts/09c_refined_hits_n_to_c_mval.R`
  and re-ran it (~1 min); product now taken from the STY table, so proteins without a
  site hit also get one
- new figure `C_mutant_validation_DinMete/results/step09_site_vs_protein/n_to_c_refined/fig_n_to_c_refined_B7H15_11550_kdpD_mval.png`
  (also a page in `n_to_c_refined_all_mval.pdf`); figure guide updated
- checks: all PASS (62 Step 8 hits reproduced, 18 LOSO refits, 5 figures written);
  `refined_hits_mval.csv` and `refined_summary_mval.csv` unchanged in content (20/62 refined)

## Results (key numbers)
- kdpD protein, `TSB_vs_SN_at_stk`: diff = +4.09, FDR = 1e-12 -> lower in SN (PASN) in `stk`
- kdpD protein, `TSB_vs_SN_at_anc`: diff = −1.03, FDR = 1.0 -> no change in `anc`
- kdpD protein, `TSB_vs_SN_at_stp`: not testable (`stp_SN` 1/3 detected)
- kdpD protein in TSB: `stk` − `anc` +1.35 (FDR 0.28), no genotype difference in TSB
- phosphosites: 10 in the STY table, only S830 testable, only in `stk_vs_anc_at_TSB`:
  +2.16, SE 0.94, FDR 1.0 (not a hit). S830 is never detected in SN samples.
- all estimates observed (our limma treat model, no imputation)

## Problems / surprises
- the phospho layer has almost no kdpD data, so the plot mainly shows protein rectangles
- the `stk` drop rests on 2/3 replicates in `anc_SN` and `stk_SN`; Step 11: 14/18 LOSO,
  inside a genomic block that shifts in `stk` -> tentative

## AI use
Claude Code: edited the script, re-ran it, checked the kdpD site/protein numbers and wrote
the figure-guide section. The numbers and the figure were marked as checked by the group on 18 Nov 2026.

## Part 2: deoB, pyk, Y1819 (non-refined site hits) added to the same folder
- question: what do the site hits that failed the refinement look like, and why did they fail?
- added deoB (`B7H15_00775`), pyk (`B7H15_09250`), Y1819 (`B7H15_09315`) to `EXTRA`, labelled
  "Step 8 site hit(s), NOT refined"; site hits outside the 6-panel grid now get an extra
  4th-row panel (only Y1819 T157, `interaction_stp_vs_stk_at_TSB_vs_SN`); subtitle now
  "Step 8 site hits: n, refined: m"; fixed an empty caption line. 8 figures, all checks PASS.
- deoB T85/T87: all 5 hits come from `stk_TSB` being low (2/3 detected, 16.2/16.3 vs ~20–21
  elsewhere); e.g. `stp_vs_stk_at_TSB` T85 +5.08 (higher in `stp`), FDR 6e-4; protein flat.
  LOSO 13–16/18 -> not robust.
- pyk T537, `stp_vs_anc_at_TSB`: +1.99 (higher in `stp`), SE 0.24, FDR 0.013, protein +0.14;
  STYH hit, LOSO 17/18 (lost only without `JE2anc_TSB_1`) -> near-robust candidate.
- Y1819 T157, `interaction_stp_vs_stk_at_TSB_vs_SN`: −2.28, FDR 0.045 (lower in SN in `stk`,
  higher in SN in `stp`); protein interaction −1.21; LOSO 5/18 -> fragile.
- all estimates observed (no imputation); n = 3 per group.

## Next steps
- decide whether kdpD belongs on the poster as a protein-level `stk` result
- decide whether pyk T537 is worth mentioning as a near-robust `stp` candidate
