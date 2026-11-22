# 2026-11-19 – Step 11b: what is inside the flagged genome blocks (group mval)

## Question for today
The Step 11 genome scan shows WHERE neighbouring proteins change together, but not WHAT
they are. Which genes lie in the flagged blocks, and how much of the Step 10 pathway
shifts (Fig 4b) is one genomic region? Asked while choosing the H3 figures.

## What we did
- ran `Rscript C_mutant_validation_DinMete/11b_genome_blocks_mval.R` (~20 s; no new test). It merges the
  Step 11 blocks into 8 regions R1–R8, gives each a label from the annotation, shows the
  region means in the Fig 4b layout and counts the Fig 4b term members inside the regions.
- checks (`genome_blocks_checks_mval.txt`): 8/8 PASS, incl. "every gene named in a label lies
  in that region", "recomputed windows reproduce the Step 11 flags", "Fig 4b rebuilt with
  the same 16 terms".
- new outputs in `C_mutant_validation_DinMete/results/step11_robustness/`: `fig_genome_blocks_mval.png`,
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
script checks them. The labels R1–R8, including
"prophage" for R4/R6 (based on product names only), were marked as checked by the group on 19 Nov 2026.

## Next steps
- Decide whether the H3 report shows `fig_genome_blocks_mval.png` next to Fig 4b.
- Optional sensitivity check: re-run the Fig 4b term tests without the block proteins.
- DNA-level check of R7/R8 (`stk`) and R5/R6 (`stp`), as in the poster's "Next" line.

---

# 2026-11-19 (later) – Step 11c: genome context of R7 and R6 at baseline TSB

## Question
Do the changed proteins in R7 (`stk`) and R6 (`stp`) sit in one operon, many operons or a
whole stretch of DNA? Do `stk` effects follow operons genome-wide? Could a polar effect of
the transposon explain R6/R7?

## What we did
- ran `Rscript C_mutant_validation_DinMete/11c_R6_R7_context_mval.R` (~20 s; no new test). Putative TU = same
  strand, gap < 150 bp (a hint from the sequence, not a measured operon).
- checks (`R6_R7_context_checks_mval.txt`): 9/9 PASS, incl. "R7 mean reproduces Step 11b
  (1.12)", "R6 prophage sits inside hlb", "stk TSB: same-TU pairs more alike than divergent pairs".
- new: `fig_R6_R7_genome_context_mval.png` + 3 tables; figure guide section 3c.

## Results (key numbers)
- R7 (~23 kb): in `stk` TSB 7 of 8 tested TUs higher by > 0.5 log2 on both strands; *kdpABC*
  +2.3 (higher in `stk`). All 3 `stk` replicates ≈ +1.1; still up in SN (+0.6) and SNpulse (+1.0).
  `stp` ≈ 0.
- SigB test: *sigB* operon +0.7 on average in `stk`, target Asp23 −0.01 → no sign of more SigB
  activity (cf. Tamber et al. 2010, doi:10.1128/IAI.00296-10: lower SigB activity in a USA300
  *pknB* mutant).
- R6: phiSa3-like prophage in *hlb* (*scn*, *chp*, staphylokinase, integrase). Most phage
  proteins about −2 in `stp` TSB (lower in `stp`), ≈ 0 in `stk`. Per-sample −0.6 to −1.3.
- Genome-wide (without R6/R7), r of adjacent-pair diffs in `stk` TSB: same TU 0.31,
  divergent −0.05. `stp`: 0.19 vs 0.08.
- Polar effect excluded: R6/R7 lie 860/953 kb from the insertions (401 putative TUs between).
- No LOD-imputed values in R6/R7.

## Problems / interpretation
- R7 looks like a regional, strain-specific effect (copy-number change? local regulator?).
  It is not the Stk1-loss signature, because `stp` also lacks most Stk1 and R7 is flat there.
- Proteomics only: "follows TUs" suggests transcription-level effects but is indirect.

## AI use
Claude Code (Opus 5.5) wrote the 11c script, the figure and guide section 3c; it looked up
Tamber et al. 2010 on PubMed. The TU rule, the choice of SigB reporters (Asp23, CrtM/N)
and the "phiSa3-like" label (based on gene content only) were marked as checked by the group on 19 Nov 2026.

## Next steps
- WGS / genomic qPCR across B7H15_11460–11600 in NE217 vs JE2 (copy number), and the R6
  prophage in all three strains.
- Optional: AureoWiki / published transcriptome operon maps instead of the 150 bp rule.

## Addendum: are the blocks "an effect of stk"? (discussion with AI)
- Our professor's hypothesis was that the R blocks come from `stk`. Our data: the blocks belong to
  one strain only (R7 +1.12 in `stk`, −0.08 in `stp`; R6 −1.09 in `stp`, −0.12 in `stk`),
  although both strains lack most Stk1. So the data support an effect of the `stk` **strain**,
  but not of Stk1 loss. Genome-wide, `stk` effects still follow TUs.
- Definition of "strain background" (any genetic difference besides the insertion),
  together with a table of direct / polar / trans / background explanations, added to the figure
  guide (section 3, under the conclusion).
- NTML reference (Fey et al. 2013, doi:10.1128/mBio.00537-12): marked as checked by the group on 19 Nov 2026.
- Next: propose complementation + WGS / genomic qPCR to the professor.

## Addendum: STRING network of R6 / R7 (Step 11d)
- **Question:** which functional groups do the R7 and R6 block proteins form, and does each
  group move with the block?
- **Run:** `Rscript C_mutant_validation_DinMete/11d_R6_R7_string_mval.R` (STRING API, course genome STRG0B09NBK,
  score ≥ 0.4; Louvain clusters; label = most specific shared STRING term + member products).
  Figure `fig_R6_R7_string_network_mval.png`, guide section 3d.
- **Results:** 123 links (R7 95, R6 27, R6–R7 1; 56 neighbourhood-only). R7: 4 clusters, all
  higher in `stk` at TSB (sigB/mazEF +0.91, kdp +1.89, alr1/murF/ddl +0.73, CsoR +1.13), flat
  in `stp`. R6: 35/49 proteins unlinked. The phage clusters (capsid/portal −2.48 in `stp`) rest
  on ≤ 2 tested proteins each.
- **Problems:** the course file `string_physical_links.txt` has only 4 links in R6, so we used
  the API. Weak-ortholog gene symbols put bacterial names on phage proteins, so we dropped them
  as labels. The KEGG term "Vancomycin resistance" is a map name (map01502), not a phenotype.
- **AI use:** Claude Code wrote and ran the script and drafted this entry and guide section 3d.
  The numbers were marked as checked by the group on 19 Nov 2026.
