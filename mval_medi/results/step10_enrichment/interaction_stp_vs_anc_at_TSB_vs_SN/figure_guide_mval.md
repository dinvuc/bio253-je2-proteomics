# Step 10: enrichment of `interaction_stp_vs_anc_at_TSB_vs_SN` (H2: altered PASN response in `stp`)

Made by [`scripts/10_enrichment_mval.R`](../../../../scripts/10_enrichment_mval.R); methods in the [summary guide](../summary/figure_guide_mval.md).
Direction: up / top = **SN response lower in `stp`** than in `anc`; down / bottom = **SN
response higher in `stp`**. (The figure title was fixed on 17 Nov.)

### `fig_top_terms_mval.png`
- **Description:** top terms of ORA up, ORA down (no hits, so empty) and the rank test.
  Background = 1,931 proteins tested in the interaction.
- **How to read:** ORA rests on a single formal hit (ilvA), so all 21 ORA terms are simply
  the annotation terms of that one protein. Use the rank test.
- **Key points:** rank test 41 terms FDR < 0.05; only **8 robust, all bottom** (SN response
  higher in `stp`): cell wall, peptidoglycan anchor (AUC 0.06), secreted, signal,
  extracellular region, cell adhesion.
- **Conclusion:** the only robust pathway-level change in the `stp` response concerns
  surface proteins. They are low in `stp` in TSB and catch up in SN, so the "higher SN
  response" mostly reflects a lower TSB baseline. Otherwise the `stp` response is like
  `anc` (H2: small).
