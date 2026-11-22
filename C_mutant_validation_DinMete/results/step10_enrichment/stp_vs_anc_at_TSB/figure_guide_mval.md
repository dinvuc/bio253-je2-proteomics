# Step 10: enrichment of `stp_vs_anc_at_TSB` (baseline: `stp` vs `anc` in TSB)

Made by [`C_mutant_validation_DinMete/scripts/10_enrichment_mval.R`](../../../scripts/10_enrichment_mval.R); methods in the [summary guide](../summary/figure_guide_mval.md).
Direction: up / top = **higher in `stp`** than `anc`; down / bottom = **higher in `anc`**.
Remember: `stp` also lacks most Stk1 (Steps 3–4).

### `fig_top_terms_mval.png`
- **Description:** top terms of ORA up, ORA down and the rank test. Background = 2,065
  proteins tested in this contrast.
- **How to read:** ORA rests on 14 formal hits (2 up, 12 down). The 21 "up" ORA terms all
  rest on 1–2 proteins (e.g. heme biosynthesis from a single hit), so they restate the hits.
- **Key points:**
  - ORA: 25 terms FDR < 0.05, only 1 with ≥ 3 hit proteins.
  - Rank test: 284 terms; 192 robust (148 higher in `stp`, 44 higher in `anc`). Top: broad
    nucleic-acid terms (DNA/tRNA-acting catalytic activity, ligases; AUC 0.71–0.81). Bottom:
    signal peptide, cell wall, peptidoglycan anchor (AUC 0.07), secreted, extracellular region.
  - Hits include Stk1 (lower in `stp`) and three prophage proteins (B7H15_10925/10935/11085).
- **Conclusion:** in TSB, `stp` has broadly lower surface and secreted proteins and
  slightly higher nucleic-acid metabolism than `anc`. The surface-protein decrease is shared
  with `stk` in direction (Step 10 summary, Fig 4b). Part of the `stp` signal sits in a
  prophage block that is lower in `stp` as a whole (Step 11).
