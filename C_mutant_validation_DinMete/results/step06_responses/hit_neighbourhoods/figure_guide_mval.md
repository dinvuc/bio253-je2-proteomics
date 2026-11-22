# Step 6c: one neighbourhood figure per interaction hit

Made by [`C_mutant_validation_DinMete/06c_hit_neighbourhoods_mval.R`](../../../06c_hit_neighbourhoods_mval.R).
Full method in the [Step 6 figure guide](../figure_guide_mval.md#hit_neighbourhoodsfig_neighbourhood_gene_locus_mvalpng-one-figure-per-hit-with-its-neighbours);
numbers in `neighbourhood_genes_mval.csv` and `neighbourhood_summary_mval.csv`.

**Question:** does the strain effect stay on the hit, or do the genes around it change too?

**Every figure has the same three parts:**
- **(a) Gene map:** the hit (orange) with the 5 nearest protein-coding genes on each side,
  drawn to scale with arrows for the strand. A blue box marks a same-strand run, which is a
  hint of co-transcription only.
- **(b) Per-sample values:** TSB (top) and SN (bottom), the three strains side by side.
  Red `n/3` = fewer than 3 replicates detected.
- **(c) The three interactions per gene:** estimate ± 2 SE. Filled = formal hit; `n.t.` = not
  testable. Dashed = ± 1 log2, the test threshold. > 0 = the SN response is lower in the
  first-named strain.

**Overall conclusion:** apart from argF and arcA, which move together, every altered response is
specific to the hit gene. Weaker same-direction trends in single neighbours (betA, leuA)
are not supported by the test (n = 3).

| figure | hit (contrast) | key points | conclusion |
|---|---|---|---|
| `fig_neighbourhood_argF_14670_mval.png` | argF (stk vs anc; stp vs stk) | neighbour arcA is also a hit, same direction (stk vs anc −2.37 / −4.98; stp vs stk +1.88 / +5.36); both in one 5-gene same-strand run with arcR, arcC2 and an antiporter | the arc/arg cluster behaves as a unit in `stk`. argF is stable in 18/18 refits (Step 11) |
| `fig_neighbourhood_arcA_14675_mval.png` | arcA (stk vs anc; stp vs stk) | same cluster as argF; arcC2 no hit; arcR not testable; one `stk_SNpulse` replicate far off | the same unit as argF. arcA holds in 16/18 refits, so slightly less robust |
| `fig_neighbourhood_14515_mval.png` | `B7H15_14515` (stk vs anc) | neighbour betA shows the same direction in stk vs anc (+1.76 vs +2.07) but is not a hit (FDR 0.8) | gene-specific by the test; mainly a TSB baseline difference. Stable 18/18 |
| `fig_neighbourhood_kdpD_11550_mval.png` | kdpD (stk vs anc) | kdpE not testable; kdpA/B/C have few SN values; kdpB stp vs stk −6.62 with a very wide SE (FDR 0.41) | the neighbourhood cannot be judged well. Step 11 shows the whole region B7H15_11470–11590 is higher in `stk` in TSB, so this hit may be part of a regional effect |
| `fig_neighbourhood_rbsD_01460_mval.png` | rbsD (stk vs anc) | measured neighbours show no effect above 1; `stp_SN` replicates spread from ~12.5 to 16 | gene-specific by the test, but rbsD lies in a block flagged in `stk` (SN) in Step 11; tentative |
| `fig_neighbourhood_ilvA_11425_mval.png` | ilvA (stp vs anc) | last gene of a 9-gene run (ilvD … leuA–D, ilvA); leuA stk vs anc +2.53 (FDR 0.73, wide SE); leuC/leuD not testable | a shared SN drop of the branched-chain amino acid operon is possible but not supported. ilvA is fragile (9/18 refits) |
| `fig_neighbourhood_06000_mval.png` | `B7H15_06000` (stp vs stk) | neighbours no effect above 1 (except `B7H15_06020` stp vs stk −1.04, FDR 1) | gene-specific; 67 kb from the insertions, so not a local transposon effect |
| `fig_neighbourhood_lipA_14865_mval.png` | lipA (stp vs stk) | most neighbours (ica genes) never detected | the neighbourhood is mostly unmeasured, so little can be said about co-regulation |

## Poster version (Step 12c)

### `fig_neighbourhood_arcA_poster_mval.png` — arcA at poster size, transparent background
- **Made by:** `C_mutant_validation_DinMete/12c_poster_arcA_mval.R` (needs the Step 5 limma table in `data/derived/`).
- **Same three panels and data as `fig_neighbourhood_arcA_14675_mval.png`**: 10 × 11 in, text ≈ 1.7×
  larger; title/subtitle replaced by a two-line caption; legends above/below the panels; gene names
  only; panel names TSB / PASN; red digit = samples detected of 3. Backgrounds transparent.
