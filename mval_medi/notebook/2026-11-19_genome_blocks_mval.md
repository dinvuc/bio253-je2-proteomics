# 2026-11-19 – Step 11b: what is inside the flagged genome blocks (group mval)

## Question for today
The Step 11 genome scan shows WHERE neighbouring proteins change together, but not WHAT
they are. Which genes lie in the flagged blocks, and how much of the Step 10 pathway
shifts (Fig 4b) is one genomic region? Asked while choosing the H3 figures.

## What we did
- ran `Rscript scripts/11b_genome_blocks_mval.R` (~20 s; no new test). It merges the
  Step 11 blocks into 8 regions R1–R8, gives each a label from the annotation, shows the
  region means in the Fig 4b layout and counts the Fig 4b term members inside the regions.
- checks (`genome_blocks_checks_mval.txt`): 8/8 PASS, incl. "every gene named in a label lies
  in that region", "recomputed windows reproduce the Step 11 flags", "Fig 4b rebuilt with
  the same 16 terms".
- new outputs in `mval_medi/results/step11_robustness/`: `fig_genome_blocks_mval.png`,
  `fig_pathway_terms_in_blocks_mval.png`, `fig_genome_scan_labelled_mval.png` + 3 tables;
  figure guide section 3b.

## Results (key numbers)
- `stk` regions, all higher in `stk` than `anc`: R7 (*sigB*, *mazEF*, *alr1*, *kdp*, *murF*,
  *ddl*) +1.1 in TSB; R8 (*sceD*, *fabZ*, *murA*, ATP synthase) +1.1 in SN; R2 (*rbs*, *lytM*,
  *esxA*) +1.0 in SN. The same regions are −0.2 to +0.2 in `stp`.
- `stp` regions, all lower in `stp` than `anc`: R5 (*spl* proteases, lantibiotic genes) −1.1
  and R6 (prophage) −1.1 in TSB; R3 −1.0 and R4 (prophage) −0.8 in SN.
- R6 is flagged in the PASN response of both mutants (−0.6 / −0.8: SN response higher in the
  mutant than in `anc`). It may be an `anc` feature seen through the shared reference.
- Fig 4b terms inside `stk` blocks: ATP synthase 5/5 (R8), SigB/*mazEF* cluster CL:1356 9/10
  (R7), vancomycin resistance 4/12, peptidoglycan 4/13 and 3/15. The `stk` response terms
  (purine, branched-chain amino acids, diaminopimelate) and the cell-killing terms: 0.
- Inside `stp` blocks: cell killing 3/12, CL:4241 3/30 (*lukD*, *hlgC*).

## Problems / surprises
- The ATP synthase shift in `stk` (Fig 4b) is entirely one block (R8, next to R7). We had
  read it as a pathway result. It is a regional signal; its cause is open.
- Region means (tested proteins of the region) differ from the scan's window means; both are
  stated in the figure guide.

## AI use
Claude Code (Opus 5.5) wrote the 11b script, the figures and figure-guide section 3b after we
asked what changes inside the blocks. Region labels were read off the member table, and the
script checks them. To check by the group: the labels R1–R8, especially
"prophage" for R4/R6 (based on product names only).

## Next steps
- Decide whether the H3 report shows `fig_genome_blocks_mval.png` next to Fig 4b.
- Optional sensitivity check: re-run the Fig 4b term tests without the block proteins.
- DNA-level check of R7/R8 (`stk`) and R5/R6 (`stp`), as in the poster's "Next" line.
