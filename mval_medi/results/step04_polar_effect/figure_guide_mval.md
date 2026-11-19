# Step 4 outputs: what each figure and table compares

Made by [`scripts/04_polar_effect_mval.R`](../../../scripts/04_polar_effect_mval.R)
(needs the Step 3 output). Question: is the loss of Stk1 in `stp` consistent
with a polar transposon insertion, and are other nearby proteins changed?
Results and numbers: [notebook 2026-11-12](../../notebook/2026-11-12_steps3-4_mutant_validation_mval.md).

All outputs are descriptive (no new statistical test). A pattern that fits a
polar effect is consistent with it, not proof of it.

## Figures

### `fig_gene_context_mval.png` — Fig 1d: where the insertions sit
- **Shows (no measurements):** the JE2 genome around *stp1*–*stk1*
  (CP020619.1, ~1.227–1.238 Mb). Arrows = genes and their direction; red
  triangles = the two bursa aurealis insertions (NE1919 in *stp1* at codon 189,
  NE217 in *stk1* at codon 46).
- **Data:** gene coordinates and strand from `JE2_genome.gff`. Insertion
  positions from `TRANSPOSON_INSERTIONS` (`R/constants.R`, NTML map; not yet
  re-checked on the NTML website). Gene symbols from `read_protein_metadata()`.
- **How to read:** pink = *stp1*/*stk1*; dark grey = protein measured in the
  total proteome; light grey = not measured. Grey box = the transcription unit
  listed by AureoWiki for strain **NCTC 8325** (SAOUHSC_01182–01187;
  MicrobesOnline and Mäder et al. 2016, PMID 27035918). We placed it on JE2
  genes by gene order, which is an assumption.
- **Key points:** *stk1* starts 4 bp before the end of *stp1* (overlapping genes, same
  strand) and is the **last gene** of the transcription unit. The NE1919 insertion sits at
  codon 189 of 247 in *stp1*, i.e. upstream of the whole *stk1* gene.
- **Conclusion:** the genome layout makes a polar effect plausible: a transposon in
  *stp1* can stop transcription or translation of *stk1*. This is a mechanism to test,
  not a measurement.

### `fig_gene_context_premium_mval.png` / `.pdf` — the same map, presentation version
- **Made by** [`scripts/04d_gene_context_premium_mval.R`](../../../scripts/04d_gene_context_premium_mval.R)
  from `gene_context_mval.csv` and `TRANSPOSON_INSERTIONS` only (no `data/` needed; ~5 s).
- **Shows (no measurements):** the same content as `fig_gene_context_mval.png`, in two tracks.
  Top: the ~11 kb window with gene arrows, the transcription unit as a bent arrow (promoter
  position schematic) and a light "lens" to the zoom. Bottom: *stp1*–*stk1* enlarged, with
  both insertions (codon, % into the gene by nucleotide), protein lengths and the 4 bp overlap.
- **How to read:** insertion colours = genotype colours of every other figure (`stp` green,
  `stk` orange); *stp1*/*stk1* arrows use the same colours. Blue-grey = same transcription
  unit (NCTC 8325, mapped to JE2 by gene order), light grey = flanking gene, dashed outline =
  protein never detected in our 27 total-proteome samples (B7H15_06325).
- **Conclusion:** unchanged: the `stp` insertion (codon 189 of 247) lies upstream of the
  whole *stk1* gene, so a polar effect is plausible, not shown.
- **Wide variant** `fig_gene_context_premium_wide_mval.{png,pdf}` (14 × 3.4 in, same script):
  one row at true genome scale, no zoom; insertions drawn below the genes, legend on the right.
  Same content and sources; for slides or a wide poster slot.

### `fig_neighbours_mval.png` — are nearby proteins changed?
- **Compares:** per-sample log2 abundance of the 9 measured proteins around
  *stp1*/*stk1* (`B7H15_06320`–`06375`), each mutant against `anc` in the same
  medium. The genes upstream in the same unit (06330–06345) are the ones a
  local disturbance would hit; 06360–06375 are outside the unit.
- **Data:** `read_abundances("total")`, normalised log2.
- **How to read:** one point = one sample, colour = genotype. All panels share
  one y axis, so small shifts look small. All 9 proteins are detected in 27/27 samples.
- **Key points:** all 9 neighbours are detected in 27/27 samples. The largest mutant-vs-`anc`
  difference is |log2| = 0.74 (`neighbour_contrasts_mval.csv`), with no formal hit.
- **Conclusion:** only *stp1*/*stk1* themselves change. The insertions do not disturb
  the expression of the surrounding genes, so the Stk1 loss in `stp` is specific to the
  gene downstream of the insertion, as a polar effect predicts.

### `fig_neighbours_by_genotype_mval.png` — the same neighbourhood, one strain per panel
- **Made by:** `scripts/04b_neighbours_by_genotype_mval.R`.
- **Compares:** the 11 measured proteins from `B7H15_06320` to `06375`, **including
  *stp1* and *stk1***, in genome order along the x axis. There is one panel per strain
  (`anc`, `stk`, `stp`), and all panels share the same y axis, so you can compare
  strains by reading straight down.
- **Data:** `read_abundances("total")`, normalised log2; gene order and labels come
  from `gene_context_mval.csv`.
- **How to read:** one point = one sample; shape = medium (TSB / SNpulse / SN).
  The grey band marks *stp1* and *stk1*. A red `n/3` means fewer than 3 replicates
  were detected. Missing samples are NA and are not drawn (never 0).
- **What to look at:** the neighbours are flat in all three strains. Only the grey
  band changes. In `stk`, Stk1 is detected in 1 of 9 samples. In `stp`, both Stp1
  and Stk1 are low and patchy, and the replicates scatter widely: Stp1 SN_2 = 18.4
  and pulseSN_3 = 18.1, while the other `stp` samples are around 12–14. Check those
  two samples before you trust a group mean.
- **Key points:** neighbours flat in every strain; in `stk` only Stk1 drops (1/9 detected);
  in `stp` both Stp1 and Stk1 drop, with two outlying Stp1 values.
- **Conclusion:** the same picture as `fig_neighbours_mval.png`, now with the two enzymes
  in place. The defect is confined to the two genes, and `stp` loses both.

### `fig_neighbours_by_genotype_medium_mval.png` — 6 panels: strain × medium
- **Made by:** `scripts/04b_neighbours_by_genotype_mval.R` (figure 2).
- **Compares:** the same 11 proteins, split into rows by medium (TSB on top, SN below)
  and into columns by strain (`anc`, `stk`, `stp`). Reading along a row compares
  strains within one medium. Reading down a column compares TSB with SN for one strain.
  **SNpulse is not shown** (18 of 27 samples).
- **How to read:** one point = one sample, and each panel holds the 3 replicates of
  one group. The grey band marks *stp1* and *stk1*. A red `n/3` means fewer than 3
  replicates were detected (NA, not drawn).
- **What to look at:** Stk1 is 0/3 in `stk` in both media and 0/3 in `stp` SN.
  In `stp` TSB, Stp1 is 1/3 (~11.9) and Stk1 is 3/3 at ~10.5–11.4, about 5 log2
  below `anc`. In `stp` SN, Stp1 splits into one sample at 18.4 (`JE2stp_SN_2`)
  and two at ~13.
- **Key points:** the pattern is the same in TSB and SN. PASN does not rescue or deepen
  the loss of either enzyme.
- **Conclusion:** the enzyme defects are stable across the two media of the primary
  analysis, so H2/H3 compare mutants with a constant enzyme state.

### `fig_stk1_peptides_by_strain_mval.png` — polar effect at peptide level
- **Compares:** Stk1 peptides quantified per sample in `anc` vs `stk` vs `stp`,
  within each medium.
- **Prediction tested by eye:** if the *stp1* insertion is polar, `stp` should
  drop to (near) the `stk` level, far below `anc`.
- **Data:** Step 3 table `enzyme_per_sample_mval.csv` (`nr_children_protein_Id`).
- **How to read:** filled dot = detected; `x` = not detected (0 peptides, count
  scale, never logged).
- **Key points:** Stk1 peptides: `anc` 27–32, `stk` 0, `stp` 1–2 in TSB and 0 in SN.
- **Conclusion:** `stp` falls to almost the `stk` level, far below `anc`, which is exactly
  what a polar insertion predicts. It is consistent with polarity but does not prove the
  mechanism: reduced stability, regulation or a secondary mutation would look the same in
  protein data.

## Tables

| file | one row per | compares / contains |
|---|---|---|
| `gene_context_mval.csv` | gene in the window (12) | start, end, strand, gap to next gene (negative = overlap), in AureoWiki unit?, UniProt ortholog, samples detected |
| `neighbour_contrasts_mval.csv` | neighbour × contrast (54) | `stk_vs_anc` / `stp_vs_anc` per medium: supplied diff, FDR, `estimate_type`, our treat FDR, `formal_hit` |
| `stk1_peptides_by_strain_mval.csv` | genotype × medium (9) | Stk1 detections, peptides per replicate, median peptides (NA counted as 0) |
| `stk1_site_drop_stp_mval.csv` | Stk1 site × medium seen in both anc and stp | n and mean log2 intensity in anc and stp (observed only, STYH), `diff_stp_minus_anc` (negative = lower in stp) |
| `polar_checks_mval.txt` | check | PASS/FLAG lines printed by the script |
