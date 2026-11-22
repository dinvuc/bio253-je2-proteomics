# Step 10 plausibility check: enrichment of `stk_vs_anc` (averaged over all media)

Made on 11 Nov by [`C_mutant_validation_DinMete/scripts/10_enrichment_mval.R`](../../../scripts/10_enrichment_mval.R) (default contrast)
to check the method against the lecture's example (lecture 06, slides 36–38). This is
**not** one of our primary contrasts: it averages TSB, SNpulse and SN. The primary analysis
is in the per-medium folders and the [summary](../summary/figure_guide_mval.md).
Direction: up / top = **higher in `stk`**; down / bottom = **higher in `anc`**.

### `fig_top_terms_mval.png`
- **Description:** top terms of ORA up (5 hits), ORA down (12 hits) and the rank test.
- **Key points:**
  - ORA: 32 terms FDR < 0.05, 16 with ≥ 3 hit proteins.
  - Rank test: 201 terms. Top includes broad biosynthesis / translation terms and the
    SigB-region STRING cluster (AUC 0.97).
  - Lecture terms reproduced: cytolysis / killing of cells and copper ion transport lower
    in `stk`, vancomycin resistance higher in `stk`. The "Virulence" keyword did not reproduce
    with our measured background.
- **Conclusion:** our local enrichment reproduces the lecture's example, so the method is
  sound. The result depends on the background (measured proteins vs genome), which is why we
  always use the tested proteins as background.
