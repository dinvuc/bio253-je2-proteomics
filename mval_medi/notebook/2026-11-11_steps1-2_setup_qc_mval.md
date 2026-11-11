# 2026-11-11 – Workflow steps 1 & 2: setup and sample QC (group mval)

## Question for today
Is our working copy complete and reproducible (Step 1), and are the 27 samples
correctly labelled and of usable quality before we interpret any biology (Step 2)?

## What we did
- Created branch `steps1-2-setup-qc-mval` from `main` at revision `93f72a0`.
- Ran `Rscript scripts/01_setup_check_mval.R`. It records git provenance, checksums,
  software versions and git exclusions, runs `tests/run_tests.R` and renders
  `analyses/01_explore_data.qmd`.
  Record: `mval_medi/results/step01_setup/setup_record_mval.txt`.
- Ran `Rscript scripts/02_sample_qc_mval.R`. It checks the design, sample-name repair,
  sample mapping across all 5 data layers, identifier uniqueness and contaminants,
  detection per sample, distributions, correlation, PCA, the reinjected phospho
  samples and DEA table metadata.
  Results: `mval_medi/results/step02_qc/` (`qc_checks_mval.txt`, CSV summaries, `fig_*_mval.png`).
- How missing values were handled: NA stays NA. Detection counts and boxplots use
  observed values only. Correlation and PCA use only proteins observed in *all*
  samples of the subset (complete-case), which biases them towards abundant proteins.

## Results (key numbers)
**Step 1**
- Checksums: 22 OK, 1 missing (`derived/peptides_total.tsv.gz`), 0 corrupt.
  The peptide table needs `msstats.csv` from the course bundle and is only needed for Step 4.
- R 4.5.0, Quarto 1.6.42, limma 3.66.0, ggplot2 3.5.2, dplyr 1.1.4.
- Helper tests pass (1 skipped: needs the peptide table). Notebook 01 renders in 11 s.
- `data/` and `analyses/_output/` are git-ignored.

**Step 2**
- 27 unique samples, 3 per genotype × medium group; primary subset (TSB + SN) = 18. PASS.
- 9 truncated names (`2anc_pulseSN_1` …) are repaired, and every sample-sheet entry
  matches its raw file name. Note: sample names keep the raw spelling
  `JE2xxx_pulseSN_n`, while the `medium` column says `SNpulse`.
- The same 27 samples with the same genotype/medium/replicate appear in total,
  phospho_STYH, phospho_STY, sites_STYH and sites_STY. PASS.
- Reinjected phospho samples: `JE2stk_pulseSN_1` and `_2` (as expected).
- Total proteome: **2,232 JE2 proteins**; 1,742 are observed in all 27 samples,
  1,758 in all of the primary 18.
- Detection per sample (total): 1,981–2,147 proteins. Group means: SN ≈ 2,025–2,044,
  TSB/SNpulse ≈ 2,025–2,138. Fewer proteins are detected in SN.
- Median log2 abundance is 16.02–16.10 in every sample, so normalisation looks fine.
- Replicate correlation (Pearson, complete proteins): within group r = 0.964–0.991,
  vs other groups r = 0.817–0.896. No sample flagged as an outlier.
- PCA (27 samples): PC1 = 62.0 % of variance, explained by medium (R² = 0.98).
  SN separates from TSB; SNpulse sits next to TSB. PC2 = 10.7 %, explained by
  genotype (R² = 0.88): `stp` separates from `anc` and `stk`, which lie close together.
  The 18-sample PCA looks the same (PC1 66.4 %, PC2 10.1 %).
  This is a QC description only, not yet a tested result.
- DEA tables: all 33 contrasts. total and phospho_STYH = `lm_impute`; `lod_imputed`
  rows = 9,042 of 76,692 (total) and 134,145 of 139,755 (STYH, i.e. 96 % imputed!).
  phospho_STY = older moderated models, no imputation. All tables carry
  `std.error`, `conf.low/high`, `df`, `p.value`, `FDR`.

## Decision on sample inclusion
**Keep all 27 samples** (all 18 for the primary analysis). No sample fails any
check. The weakest is `JE2stk_TSB_3` (lowest detection, 1,981 proteins,
z = −2.1; lowest within-group r = 0.964). It is still within normal
range, so it stays in. We will keep an eye on it when we plot `stk` TSB values.

## Problems / surprises
- **Contaminant filter is incomplete.** `read_abundances("total")` with the default
  `drop_contaminants = TRUE` still contains 68 contaminant proteins (keratins etc.,
  not flagged by `is_contaminant`) and 24 reversed decoys (`rev_…`). The total DEA
  also contains the 68 contaminants. Our scripts keep only rows with a `locus_tag`
  that are not `rev_`. Every later script must do the same → tell the lecturer.
- **Phospho detection is very uneven between replicates.** Example: anc TSB rep 1
  has ≈ 1,000 sites, reps 2–3 ≈ 450. The reinjected `stk_pulseSN_1/2` see more sites
  (≈ 400) than rep 3 (≈ 285), which fits a batch effect. Phospho is outside our
  focused scope, but this matters if we ever use it.
- `R/setup.R` does not load purrr/tibble; our script loads them itself.
- Stray `*_unzipped*` copies in `data/derived/` and `reference/` (about 130 MB). The 3 in
  `reference/` are **not** git-ignored → do not `git add .` until they are removed.

## AI use
Claude Code (Opus 5.5): wrote and ran `scripts/01_setup_check_mval.R` and
`scripts/02_sample_qc_mval.R` and drafted this entry. Corrections during the
session: the reinjected-sample check first expected `SNpulse` in the sample names
(they are `pulseSN`), and the script needed purrr loaded explicitly. To check
the numbers yourselves, look at `qc_checks_mval.txt` and the PCA/detection figures.

## Addendum (same day): lecture 06 "web resources" applied to our plan

**Question:** Does the lecture `06_web_resources_2026.pdf` help us, or solve one of our problems?
It does not solve the missing peptide table, the contaminant flag or the phospho
detection problem. It does help with Step 4 (sources) and Step 10 (enrichment).

**What we did**
- Added `R/helpers_mval.R` (loaded automatically by `R/setup.R`):
  - `keep_je2()`: our shared JE2-protein filter.
  - `hit_list()`: up/down hits; imputed hits labelled.
  - `run_ora()`: Fisher test against the measured background.
  - `run_rank_enrichment()`: GSEA-like Wilcoxon rank test.
  - `write_string_upload()`: locus-tag lists for STRING.
- Added `tests/testthat/test-helpers_mval.R` (toy data). `Rscript tests/run_tests.R`: all pass
  (data, helpers_mval 11 checks, parsing; 1 skip = peptide table).
- Added `scripts/10_enrichment_mval.R` and ran it on `stk_vs_anc` as a method check
  (outputs in `mval_medi/results/step10_enrichment/stk_vs_anc/`). Runtime ≈ 7 s.
- `scripts/02_sample_qc_mval.R` now uses `keep_je2()`. The rerun gives identical numbers.
- Updated `PROJECT_WORKFLOW_mval.md`:
  - JE2-filter rule (no `grep 'ARG'`);
  - STRING organism `STRG0B09NBK` and AureoWiki/NTML sources;
  - enrichment method in Step 10;
  - status section.

**Results (method check, `stk_vs_anc`; diff > 0 = higher in stk, diff < 0 = higher in anc)**
- Hits (FDR < 0.1, |log2FC| > 1): 36 higher in stk (4 `lod_imputed`) and 51 higher in anc (7 `lod_imputed`),
  out of 2,232 tested JE2 proteins. The rank test used 1,992 observed estimates (240 imputed left out).
- ORA, higher in anc (= lower in stk): extracellular region (GO:0005576, 13 hits, FDR 1e-6),
  cytolysis (GO:0019835, 6/18, FDR 2.5e-4), killing of cells of another organism (GO:0031640, 5/13, FDR 2.5e-4).
- Rank test: vancomycin resistance (map01502, 12 proteins, AUC 0.90, FDR 1.4e-4) is shifted towards
  higher in stk. Translation-related STRING clusters lean slightly towards higher in stk (AUC ≈ 0.62).
- **Agreement with the lecture example (slides 36–38):** cytolysis/killing (lower in stk) and
  vancomycin resistance (higher in stk) reproduce. Copper ion transport (GO:0006825) reproduces only
  if 3-protein terms are allowed (rank test, FDR 0.049). The two copper ATPases B7H15_00415 and
  B7H15_14230 and CopZ (B7H15_14235) are all lower in stk (diff −1.1 to −1.5, observed).
  "Virulence" (KW-0843) does **not** reproduce (1/22 hits). The lecture probably used a whole-genome
  background and/or a different cut-off.
- This is a test of the tool, not our Step 10 result: `stk_vs_anc` averages over all three media,
  and we have not yet applied the Step 5 eligibility rules.

**Problems / notes**
- The local STRING files use organism id `STRG0A35SHG`; the lecture uses `STRG0B09NBK`. Match on locus tags.
- The BioCyc slide shows a shared course login. We do not copy it into the repo (public on GitHub).

## Next steps
- Remove the stray `*_unzipped*` files; get `msstats.csv` → peptide table (for Step 4).
- Ask the lecturer about the unflagged contaminants/decoys.
- Step 3: validate Stp1 and Stk1 levels (H1) with `plot_protein()` and the DEA,
  using JE2 proteins only.
