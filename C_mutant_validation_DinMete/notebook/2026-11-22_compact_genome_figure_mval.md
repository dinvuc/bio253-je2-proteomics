# 2026-11-22 – compact R6/R7 figure for the poster (group mval: Din Vuckic & Mete Kabranlar)

## Question for today
Can the R6/R7 genome figure (poster Fig. 3C) be narrower with the same message?

## What we did
- Plan: aim = narrow version of panels A + B; input = the Step 11c objects; analysis = none new;
  output = `fig_R6_R7_genome_context_compact_mval.png`; check = Step 11c checks still pass, labels readable.
- Added section 7 to `C_mutant_validation_DinMete/scripts/11c_R6_R7_context_mval.R` and ran it.
- One-off check, not saved as a result: Step 7 noise-corrected r without the block proteins.

## Results (key numbers)
- Figure is 7 x 6.4 in (ratio 1.09).
- Baseline r: 0.373 with all proteins, 0.395 without R6/R7, 0.408 without R1-R8.
- PASN-response r: 0.378, 0.364, 0.366 (same three sets).

## Problems / surprises
- The blocks explain little of the limited stk/stp similarity.
- Poster wording to check: R7 in `stk` (higher than `anc`) is +1.12 TSB, +0.70 SN, +0.96 SNpulse, so not
  "the same" in every medium; ATP synthase is in R8, not R7; R6 is also +0.49 in `stk` in SN.

## AI use
Claude Code changed the script, ran the check and drafted this entry. A group member must check it.

## Next steps
- Replace the figure on the poster and decide on the wording.

## Addendum: compact genome scan
- **Question:** can the labelled genome scan be narrower without crowded labels?
- **Run:** new section 4b in `C_mutant_validation_DinMete/scripts/11b_genome_blocks_mval.R`; no new statistics.
- **Output:** `fig_genome_scan_labelled_compact_mval.png` (7 x 7.4 in; mutants = columns, media = rows).
- **To check on the poster:** R6 is flagged in the PASN response of BOTH mutants, so "each block
  appears only in one mutant" holds at baseline only. The 5% threshold is per panel (six panels,
  not corrected), so blocks just over the line are weak evidence.
- **AI use:** Claude Code changed the script and drafted this addendum.

## Addendum: compact Fig 3 (H3)
- **Question:** can the stp-vs-stk scatter be narrower and less crowded?
- **Run:** new compact section in `C_mutant_validation_DinMete/scripts/12_poster_figures_mval.R`; no new statistics.
- **Output:** `step12_poster/fig3_h3_shared_distinct_compact_mval.png` (7 x 6.15 in; axes swapped).
- **To check on the poster:** "the largest differences sit in strain-specific genome blocks" is too
  strong. Of the direct stp-vs-stk hits, 6 of 20 at baseline (3 in R6/R7) and 0 of 4 in the PASN
  response lie in R1-R8. The largest gap overall is arcA in the response (5.36 log2), outside all blocks.
- **AI use:** Claude Code changed the script, ran the count and drafted this addendum.
