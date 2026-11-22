# Step 10: enrichment of `stk_vs_anc_at_SN` (baseline: `stk` vs `anc` in PASN)

Made by [`C_mutant_validation_DinMete/scripts/10_enrichment_mval.R`](../../../scripts/10_enrichment_mval.R); methods in the [summary guide](../summary/figure_guide_mval.md).
Direction: up / top = **higher in `stk`** than `anc` (in SN); down / bottom = **higher in `anc`**.

### `fig_top_terms_mval.png`
- **Description:** top terms of ORA up, ORA down and the rank test (−log10 FDR; dashed =
  FDR 0.05). Background = 1,964 proteins tested in this contrast.
- **How to read:** as in the other contrast folders. ORA rests on only 12 formal hits
  (4 up, 8 down), so check `hits_in_term`.
- **Key points:**
  - ORA: 23 terms FDR < 0.05, of which 5 rest on ≥ 3 hit proteins. The down terms (higher in
    `anc`) are glycine cleavage / lipoylation, driven by 2–3 proteins.
  - Rank test: 316 terms FDR < 0.05; 143 robust (96 higher in `stk`, 47 higher in `anc`). Top
    (higher in `stk`): ATP synthase F1 (AUC 0.99), plus many broad biosynthesis terms. Bottom
    (higher in `anc`): de novo IMP / purine biosynthesis (AUC 0.08), transition-metal transport.
- **Conclusion:** in PASN, `stk` shows broad, small shifts across many processes, with
  ATP synthase higher and purine synthesis lower than in `anc`. Many terms overlap (GO
  hierarchy), so "143 terms" means a few biological themes, not 143 independent pathways.
