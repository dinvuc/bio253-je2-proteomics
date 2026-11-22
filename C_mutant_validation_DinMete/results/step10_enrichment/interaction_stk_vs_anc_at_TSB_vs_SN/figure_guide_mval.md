# Step 10: enrichment of `interaction_stk_vs_anc_at_TSB_vs_SN` (H2: altered PASN response in `stk`)

Made by [`C_mutant_validation_DinMete/10_enrichment_mval.R`](../../../10_enrichment_mval.R); methods in the [summary guide](../summary/figure_guide_mval.md).
Direction: up / top = **SN response lower in `stk`** than in `anc` (interaction > 0); down /
bottom = **SN response higher in `stk`**. (The figure title was fixed on 17 Nov; it used to say
"higher in interaction_stk".)

### `fig_top_terms_mval.png`
- **Description:** top terms of ORA up, ORA down and the rank test. Background = 1,882
  proteins tested in the interaction.
- **How to read:** ORA rests on only 5 formal hits (2 up, 3 down), so most ORA terms are the
  hit proteins themselves. The rank test is the informative panel.
- **Key points:**
  - ORA: 14 terms, all "down". The only one with ≥ 3 hit proteins is the broad keyword
    "Cytoplasm" (3/255). The arginine and ornithine metabolism terms rest on 2 proteins
    (argF, arcA).
  - Rank test: 97 terms; 44 robust (18 top, 26 bottom). Top (SN response lower in `stk`): de
    novo IMP / purine biosynthesis (AUC 0.98). Bottom (SN response higher in `stk`): pyrimidine
    nucleotide biosynthesis, glycerol metabolism.
- **Conclusion:** at pathway level, `stk` changes its PASN response in nucleotide
  synthesis: purine synthesis falls more in `stk` than in `anc`, pyrimidine synthesis rises
  more. The arginine/ornithine signal comes from the argF–arcA pair (Step 6).
