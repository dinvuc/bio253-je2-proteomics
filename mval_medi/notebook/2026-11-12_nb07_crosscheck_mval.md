# 2026-11-12 – notebook 07 peptide layer vs steps 3–4, plus neighbour figures (group mval)

## Question for today
1. Show the proteins around stp1/stk1 one strain at a time (and strain × medium).
2. Does the starter-pack notebook `analyses/07_transposon_mutants.qmd` (peptide layer)
   contradict our Step 3 (mutant validation) and Step 4 (polar effect) results?

## What we did
- `Rscript scripts/04b_neighbours_by_genotype_mval.R` → `fig_neighbours_by_genotype_mval.png`
  (3 panels: anc / stk / stp, genes in genome order incl. stp1/stk1) and
  `fig_neighbours_by_genotype_medium_mval.png` (6 panels: TSB and SN × 3 strains; SNpulse left out).
  Merged to `main` (c061bb0).
- Notebook 07 could not render: `data/derived/peptides_total.tsv.gz` was missing from our data zip.
  We extracted the course `DEAs.zip` into `data/raw/DEAs/` (git-ignored, 1.7 GB) and built only the
  peptide table from `DEA_20261029_WUtotal_JE2_Oct26_vsn/.../msstats.csv` with `parse_msstats_peptides()`.
  The full `00_fetch_data.R` was **not** run (it also needs `phospho_combined_sites/` and
  `more_resources/`, which are not in that zip). Tests pass; notebook 07 renders.
- `Rscript scripts/04c_peptide_crosscheck_mval.R` → `peptide_crosscheck_mval.csv`,
  `peptide_crosscheck_checks_mval.txt`; written comparison in
  `mval_medi/results/step04_polar_effect/nb07_comparison_mval.md`.

## Results (key numbers)
- **No contradiction** between notebook 07 and steps 3–4. Peptide counts per sample equal our
  Step 3 `n_peptides` in all 54 Stp1/Stk1 × sample rows (PASS).
- Neighbours (06320–06375): flat in all strains and media in the new figures; only stp1/stk1 change.
- Stk1 in stk: 1/9 samples (`stk_pulseSN_1`, 1 peptide, residues 332–340), same as Step 3.
- Stk1 in stp: peptides per sample TSB 1/1/2, SNpulse 2/1/–, SN 0, same as Step 4.
- Stp1, peptides after the insertion (codon 189): anc 24, stk 23, **stp 1** detection.
  Notebook 07 reads this as a truncated Stp1 (N-terminal part only). Our "supported reduction" call stands.
- The wild-type-like Stp1 value `JE2stp_SN_2` (18.4 log2) rests on residues 51–57 and **233–246**.
  The second peptide lies after the insertion, so it should not exist in stp1::Tn → most likely an
  artefact. `JE2stp_pulseSN_3` (18.1) rests on one 8-residue peptide (174–181).

## Problems / surprises
- Our rebuilt `peptides_total.tsv.gz` has a different md5 from the lecturer's
  (`data/check_data.R` will flag it). Expected for a re-compressed file; ask the lecturer for theirs.
- First version of the cross-check FLAGged mismatches: it counted unique sequences, while Step 3
  counts peptide forms (e.g. oxidised M separately). Fixed; now all rows match.
- Notebook 07 says "for Stk1, stk and stp are simply empty". This refers to the N-terminal part only.
  stp still has some Stk1 (TSB/SNpulse), so we keep "strong but not complete knock-down".
- Still unexplained: Stk1 signals in stk after codon 46 (1 total-proteome detection, phospho T231 in
  3/9 samples). A truncated product cannot explain them (carry-over? mis-assignment?).
- Notebook 07 calls the polar effect "exactly" what the proteome shows; we keep "consistent with,
  not proof of".

## AI use
An AI coding assistant (Claude Code) wrote and ran the two scripts and the figures, set up the peptide
table, and wrote the comparison. It found and fixed two of its own bugs: n/3 labels under the wrong
medium, and peptide forms vs. sequences. The group still has to check the numbers against the CSVs.

## Next steps
- Look up the run order around `JE2stp_SN_2` and the stk samples with Stk1 signals (carry-over?).
- Ask the lecturer for the original `peptides_total.tsv.gz`.
- Group: check `nb07_comparison_mval.md` and the two new figures.
