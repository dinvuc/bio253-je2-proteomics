# 2026-11-18 – kdpD added to the refined N-to-C plots (group mval)

## Question for today
What do the protein and phosphosite layers show for kdpD (`B7H15_11550`), the Step 6
interaction hit (`interaction_stk_vs_anc_at_TSB_vs_SN`, `fig_response_vs_anc_mval.png`),
when it is drawn in the same N-to-C style as the refined-hit proteins of Step 9c?

## What we did
- added kdpD as an extra protein (`EXTRA`) to `scripts/09c_refined_hits_n_to_c_mval.R`
  and re-ran it (~1 min); product now taken from the STY table, so proteins without a
  site hit also get one
- new figure `mval_medi/results/step09_site_vs_protein/n_to_c_refined/fig_n_to_c_refined_B7H15_11550_kdpD_mval.png`
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
the figure-guide section. A group member still has to check the numbers and the figure.

## Next steps
- decide whether kdpD belongs on the poster as a protein-level `stk` result
