# Step 10: enrichment of `stp_vs_anc_at_SN` (baseline: `stp` vs `anc` in PASN)

Made by [`C_mutant_validation_DinMete/scripts/10_enrichment_mval.R`](../../../scripts/10_enrichment_mval.R); methods in the [summary guide](../summary/figure_guide_mval.md).
Direction: up / top = **higher in `stp`** than `anc` (in SN); down / bottom = **higher in `anc`**.

### `fig_top_terms_mval.png`
- **Description:** top terms of ORA (up has no hits, so its panel is empty) and the rank
  test. Background = 1,979 proteins tested in this contrast.
- **How to read:** ORA rests on only 4 formal hits (all down), so it mostly restates them.
- **Key points:**
  - ORA down: 49 terms FDR < 0.05, 11 with ≥ 3 hit proteins (from the 4 hits, incl. ilvA).
  - Rank test: 294 terms; 168 robust (150 higher in `stp`, 18 higher in `anc`). Top: RNA,
    ncRNA and tRNA metabolism (AUC 0.72–0.76). Bottom: ornithine metabolism, peroxidase /
    antioxidant activity, signal peptide.
- **Conclusion:** in PASN, `stp` keeps the broad nucleic-acid shift seen in TSB and has
  lower antioxidant/peroxidase proteins than `anc`. These are many small shifts, not a few
  strong protein changes (only 4 formal hits).
