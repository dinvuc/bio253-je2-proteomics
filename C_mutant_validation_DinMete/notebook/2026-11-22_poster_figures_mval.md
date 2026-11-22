# 2026-11-22 – poster versions of the H1/H2 figures (group mval)

## Question for today
Make the figures used on the poster readable at poster size, with transparent backgrounds.

## What we did
- Redrew six figures with the same data and thresholds, larger text (base 17 pt) and no titles
  (`C_mutant_validation_DinMete/scripts/12c_poster_*_mval.R`, shared helper `C_mutant_validation_DinMete/scripts/12c_poster_transparent_helper_mval.R`):
  Step 3 enzyme peptides, Step 4b neighbourhood by strain, Step 6 arcA neighbourhood, the two Step 6
  response scatter plots and one combined 3-panel response figure.
- The Step 5 table `data/derived/limma_treat_total_mval.tsv.gz` was missing locally; it was recreated
  with the fit part of `C_mutant_validation_DinMete/scripts/05_formal_test_mval.R` (git-ignored, not committed). Check: the
  interaction estimates equal the stored Step 6 values (58 genes, max difference 0).

## Results (key numbers)
- Unchanged from Step 6: r = 0.85 / 0.86 / 0.84; hits 5 (stk vs anc), 1 (stp vs anc), 4 (stp vs stk)
  of 1,882 / 1,931 / 1,876 eligible proteins.
- The figures rendered from the repo are byte-identical to the versions used on the poster.

## Problems / surprises
- Equal axis units with different value ranges made the combined figure's panels unequal; fixed by
  square value ranges per plot and fixed 3.6 in panels.

## AI use
Claude Code wrote and ran the scripts and checked each figure by eye and the numbers against Step 6.
The figures were marked as checked by the group on 22 Nov 2026.

## Next steps
- Use the transparent PNGs on the poster background; captions in the poster text.
