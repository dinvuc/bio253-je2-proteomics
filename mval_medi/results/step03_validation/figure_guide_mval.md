# Step 3 outputs: what each figure and table compares

Made by [`scripts/03_mutant_validation_mval.R`](../../../scripts/03_mutant_validation_mval.R).
Question (H1): do the mutants lack the enzyme they should lack, and how much
peptide evidence is behind each "detected" value? Results and numbers:
[notebook 2026-11-12](../../notebook/2026-11-12_steps3-4_mutant_validation_mval.md).

None of the figures is a statistical test. They show the raw per-sample values
that every fold change for Stp1/Stk1 rests on (CLAUDE.md rule 4).

## Figures

### `fig_enzyme_abundance_mval.png` — Fig 1a: enzyme amount per sample
- **Compares:** Stp1 (`B7H15_06350`) and Stk1 (`B7H15_06355`) in all 27 samples:
  3 genotypes (`anc`, `stk`, `stp`) × 3 media (TSB, SNpulse, SN) × 3 replicates.
  The main comparison is each mutant against `anc` in the same medium.
- **Data:** normalised log2 protein abundance, `read_abundances("total")`
  (already log2 and normalised, not changed). The number above each point is
  the peptides quantified in that sample (`nr_children_protein_Id`, sheet
  `normalized_abundances` of `DE_total.xlsx`).
- **How to read:** one point = one sample. Filled = ≥ 2 peptides, open = 1 peptide.
  Red `n/3` under a group = only n of 3 samples detected; missing samples are
  NOT drawn (NA, never 0). The panel title gives the `anc` peptide range as the reference.
- **Caveat:** a value from 1 peptide is not comparable with one from 30. The
  "wild-type-like" values in `stp` (Stp1) and `stk_pulseSN_1` (Stk1) are single-
  or two-peptide values.

### `fig_enzyme_peptides_mval.png` — Fig 1b: peptide evidence per sample
- **Compares:** the number of peptides quantified per sample for Stp1 and Stk1
  across the 9 groups. The maximum seen in the whole data set is 20 (Stp1) and
  37 (Stk1).
- **Data:** `nr_children_protein_Id` from `DE_total.xlsx`, joined to the sample
  sheet by sample name.
- **How to read:** filled dot = detected, `x` = not detected, drawn at 0 peptides.
  This is a count scale; the 0 is never used in a log2 calculation. Y axes differ
  between the two panels.

### `fig_stk1_phosphosite_map_mval.png` — Fig 1c: which parts of Stk1 are seen
- **Compares:** Stk1 phosphosite detections per sample along the Stk1 sequence
  (x = amino-acid position, 1–664), between strains, and between the two phospho
  searches (STY left, STYH right; these are alternative searches, not replicates).
- **Data:** `read_sites("STY")` and `read_sites("STYH")`, `log2_intensity`
  (FragPipe site intensity, not normalised). Stp1 has no phosphosites in either
  search, so there is no Stp1 map.
- **How to read:** a coloured tile = site detected in that sample (colour =
  log2 intensity); blank = not detected. All 27 samples are listed, so an
  empty row is a sample with no Stk1 site at all. The red dashed line is the `stk`
  insertion (codon 46). The `stp` insertion is in Stp1 and is not on this axis.
- **Caveat:** localisation probability is not shown (see
  `hiq_peptide_check_mval.csv`). `JE2stk_pulseSN_1` was re-injected on another day.

## Tables

| file | one row per | compares / contains |
|---|---|---|
| `enzyme_per_sample_mval.csv` | protein × sample (54) | log2 abundance, peptides in the sample, peptides in the data set |
| `enzyme_per_group_mval.csv` | protein × group (18) | detections out of 3, samples with ≥ 2 peptides, peptides per replicate, mean of detected values |
| `enzyme_contrasts_mval.csv` | protein × contrast (18) | `stk_vs_anc`, `stp_vs_anc`, `stp_vs_stk` within each medium: supplied DEA (diff, SE, FDR, `estimate_type`, model) next to our treat test (diff, SE, FDR, tested?). `direction` says which group is higher |
| `stk1_sites_mval.csv` | search × site (27) | detections out of 9 per genotype, `best_loc_prob` |
| `hiq_peptide_check_mval.csv` | detected T231/S240/S242 value × sample | intensity, per-sample `loc_prob`, `reinjected`: is the `HIQDSVPNVTtDVR` peptide real evidence? |
| `evidence_table_mval.csv` | enzyme × mutant (4) | rule-based call: "supported reduction", "no reduction" or "inconsistent pattern" (rule written in the script, section 7, before results were read) |
| `validation_checks_mval.txt` | check | PASS/FLAG/INFO lines printed by the script |

Diff convention: `A_vs_B` = log2(A) − log2(B); negative = higher in B.
`lod_imputed` = at least one group was never measured: presence/absence
evidence only, and the size of the fold change is not meaningful.
