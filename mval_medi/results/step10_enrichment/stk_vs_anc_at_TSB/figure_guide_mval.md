# Step 10: enrichment of `stk_vs_anc_at_TSB` (baseline: `stk` vs `anc` in TSB)

Made by [`scripts/10_enrichment_mval.R`](../../../../scripts/10_enrichment_mval.R). Methods, robust-term rule
(Wilcoxon and cameraPR both FDR < 0.05) and the cross-contrast figures: [summary guide](../summary/figure_guide_mval.md).
Direction: up / top = **higher in `stk`** than `anc`; down / bottom = **higher in `anc`**.

### `fig_top_terms_mval.png`
- **Description:** the most significant terms of the three tests for this contrast: ORA up
  (top), ORA down (middle) and the rank test (bottom), as −log10 FDR bars. Dashed = FDR 0.05.
  Background = the 2,002 proteins tested in this contrast.
- **How to read:** a longer bar = a stronger enrichment. ORA uses only the 20 formal hits
  (7 up, 13 down); the rank test uses all 2,002 proteins. Check `ora_*_mval.csv` for how many
  hit proteins a term rests on (`hits_in_term`).
- **Key points:**
  - ORA: 30 terms FDR < 0.05 (10 up, 20 down), but only 13 rest on ≥ 3 hit proteins (12 of
    them down). The down terms are cytolysis / cell killing / haemolysis, driven by three
    haemolysins: alpha-haemolysin `B7H15_06050` and the gamma-haemolysin subunits
    `B7H15_11180` and `B7H15_13435`.
  - Rank test: 51 terms; 35 robust (28 higher in `stk`, 7 higher in `anc`). Higher in
    `stk`: a STRING cluster that is the SigB/*mazEF* region (AUC 0.95), vancomycin
    resistance, peptidoglycan synthesis, ATP synthase F1. Higher in `anc`: cell killing,
    cytolysis, penicillin binding (AUC 0.04), a haemolysis cluster.
- **Conclusion:** in TSB, `stk` has lower toxin/cytolysis and surface proteins and higher
  cell-wall synthesis, ATP synthase and SigB-region proteins than `anc`. The SigB-region
  signal coincides with a genomic block that is higher in `stk` as a whole (Step 11), so part
  of it may be regional rather than a pathway effect.
