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
- **Point of the figure:** *stk1* is directly downstream of *stp1* (4 bp overlap)
  and is the last gene of the unit, so an insertion in *stp1* can block *stk1*.

### `fig_neighbours_mval.png` — are nearby proteins changed?
- **Compares:** per-sample log2 abundance of the 9 measured proteins around
  *stp1*/*stk1* (`B7H15_06320`–`06375`), each mutant against `anc` in the same
  medium. The genes upstream in the same unit (06330–06345) are the ones a
  local disturbance would hit; 06360–06375 are outside the unit.
- **Data:** `read_abundances("total")`, normalised log2.
- **How to read:** one point = one sample, colour = genotype. All panels share
  one y axis, so small shifts look small. All 9 proteins are detected in 27/27 samples.
- **Matching numbers:** `neighbour_contrasts_mval.csv` (largest |diff| = 0.74, no formal hit).

### `fig_stk1_peptides_by_strain_mval.png` — polar effect at peptide level
- **Compares:** Stk1 peptides quantified per sample in `anc` vs `stk` vs `stp`,
  within each medium.
- **Prediction tested by eye:** if the *stp1* insertion is polar, `stp` should
  drop to (near) the `stk` level, far below `anc`.
- **Data:** Step 3 table `enzyme_per_sample_mval.csv` (`nr_children_protein_Id`).
- **How to read:** filled dot = detected; `x` = not detected (0 peptides, count
  scale, never logged).

## Tables

| file | one row per | compares / contains |
|---|---|---|
| `gene_context_mval.csv` | gene in the window (12) | start, end, strand, gap to next gene (negative = overlap), in AureoWiki unit?, UniProt ortholog, samples detected |
| `neighbour_contrasts_mval.csv` | neighbour × contrast (54) | `stk_vs_anc` / `stp_vs_anc` per medium: supplied diff, FDR, `estimate_type`, our treat FDR, `formal_hit` |
| `stk1_peptides_by_strain_mval.csv` | genotype × medium (9) | Stk1 detections, peptides per replicate, median peptides (NA counted as 0) |
| `stk1_site_drop_stp_mval.csv` | Stk1 site × medium seen in both anc and stp | n and mean log2 intensity in anc and stp (observed only, STYH), `diff_stp_minus_anc` (negative = lower in stp) |
| `polar_checks_mval.txt` | check | PASS/FLAG lines printed by the script |
