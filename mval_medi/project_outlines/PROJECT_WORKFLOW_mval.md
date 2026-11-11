# Project C workflow: from mutant validation to a defensible conclusion

Created 10 November 2026. This is the practical companion to [PROJECT_PLAN_mval.md](PROJECT_PLAN_mval.md), which defines our agreed research questions and hypotheses. Check off a step only after its output has been inspected and recorded.

## What we want to produce

Answer three questions with evidence:

1. **H1: What protein profiles do the mutants actually have?** Does `stp` show reduced Stp1 and Stk1, and what happens to both proteins in `stk`?
2. **H2: Do mutations alter the response to supernatant?** Separate differences already present in TSB from differences in the TSB-to-SN response.
3. **H3: Which response changes are shared or distinct?** Assess whether `stp` has a kinase-mutant-like component, while examining differences directly.

The wanted result is a reliable answer to each question, including uncertainty or contradictory evidence. We are testing these hypotheses, not trying to force a predicted outcome. Phosphorylation is a supporting layer; causal rescue and transcriptional polarity cannot be established by this dataset alone.

## Scope and working rules

- Primary analysis: `anc`, `stk`, `stp` in TSB and SN, three samples per group: **18 samples**. Inspect all 27 samples for mutant validation and sample quality; reserve the nine SNpulse samples for an optional extension.
- SN is approximately 15 h growth with 30% cell-free, sterile-filtered PAO1 spent medium. This is not co-culture; secreted products and nutrient depletion are not separated.
- Use R and the existing helpers, starting with `source("R/setup.R")` at the repository root or `source("../R/setup.R")` inside `analyses/`.
- Normalized abundance is already on the log2 scale. Do not normalize or log-transform it again. Missing values stay `NA`, never zero.
- Keep course inputs and full machine-readable measurement extracts local and ignored. Share analysis code, selected figures, reports and concise notebook findings; never stage course measurements under `data/`.
- Record decisions and AI assistance in a dated `mval_medi/notebook/` entry after each session. Save methods even when a result is negative.
- **Real JE2 proteins only:** pipe every abundance or DEA table through `keep_je2()` (`R/helpers_mval.R`). It keeps rows with a locus tag and drops reversed decoys. The default `drop_contaminants = TRUE` misses 68 contaminants and 24 `rev_` decoys (Step 2 finding, 11 Nov). Do **not** build hit lists with the lecture's `grep 'ARG'` shortcut: it silently drops 15 real JE2 proteins whose key has no "ARG".

## What we need, and what each item is for

| Requirement | Existing source or location | Purpose |
|---|---|---|
| Agreed hypotheses and design | `PROJECT_PLAN_mval.md` | Keep the analysis focused on Topic C and distinguish predictions from findings. |
| R, packages and Quarto | `install_packages.R`, `R/setup.R`, `tests/` | Run reproducible analysis and render readable reports. Record versions with `sessionInfo()`. |
| Authentic course bundle | `data/checksums.md5`, `data/check_data.R` | Verify inputs before interpreting unexpected measurements. A previous check passed all 22 expected files on 5 November; recheck this working copy. |
| Sample identities and exposure metadata | `read_samples()`, `data/derived/samples.tsv` | Assign genotype, medium, replicate and reinjection information correctly; prevent sample-join errors. |
| Individual total-protein measurements | `read_abundances("total")` | Validate Stp1/Stk1, inspect replicate consistency and check selected statistical effects. |
| Total-protein statistical results | `read_dea("total")` | Obtain baseline and interaction estimates, standard errors, FDR and observed/imputed labels. |
| Individual phosphosite measurements and localization | `read_sites("STYH")` / `read_sites("STY")` | Decide which sites have enough replicate support and credible residue assignment. |
| Matching phospho abundances and statistics | `read_abundances("phospho_STYH")`, `read_dea("phospho_STYH")`, and corresponding STY calls | Compare supported phosphorylation responses using the search/model actually used for inference. |
| JE2 protein and genome annotation | `read_annotation()`, `data/raw/resources/JE2_proteins.fasta`, `JE2_genome.gff` | Resolve proteins, site coordinates, neighboring genes and gene orientation. |
| Bacterial functional mappings | `read_protein_terms()` (GO, KEGG, UniProt keywords, STRING clusters; 2,144 of 2,232 measured JE2 proteins have ≥ 1 term); raw `string_aliases.txt` / `string_enrichment_terms.txt` in `data/raw/resources/`; STRING website with the course organism **`STRG0B09NBK` (JE2_Bio253_HS26)**, which accepts `B7H15_…` locus tags | Interpret candidates and test enrichment locally against a measured background (`scripts/10_enrichment_mval.R`); use the STRING website for network pictures and as a cross-check. The local raw files carry a different STRING organism id (`STRG0A35SHG`), so match on locus tags, not on STRING ids. |
| Gene context and mutant sources | AureoWiki (aureowiki.med.uni-greifswald.de, searchable by NCTC 8325 locus tags / gene names), NTML methods (ntml.unmc.edu/methods), NCBI genome GCF_002085525.1, BLAST for doubtful orthologs | Citable sources for the *stp1*–*stk1* operon structure, the NE217/NE1919 insertions and ortholog checks (Step 4). Lecture `06_web_resources_2026.pdf`. Do not copy the shared BioCyc login into the repo. |
| Verified biological evidence | `mval_medi/docs/stk_stp_rescue_evidence_mval.md`, course papers; checked UniProt/primary papers as needed | Support protein-function and suppression claims, while checking species, strain and experimental context. |
| Lecturer clarification where necessary | Sample/batch records, insertion details, peptide evidence and original QC reports if available | Resolve anomalies that processed protein tables alone cannot explain. These are conditional requests, not reasons to postpone all work. |

The bundle contains processed course results, not instrument-level spectra. Spectrum reprocessing is outside the main workflow. The previously recorded environment lacked R; Rscript and Quarto are now discoverable on this PC, but package loading and rendering still need verification.

## Step 1 — Establish a reproducible starting point

**Do**

- [x] Confirm the working branch and record the starting revision and input provenance.
- [x] Run `Rscript data/check_data.R`; resolve missing or mismatched inputs before analysis.
- [x] Check R/package availability; use `Rscript install_packages.R` if needed. Verify `source("R/setup.R")`, run `Rscript tests/run_tests.R`, and render a starter notebook.
- [x] Record R, package and Quarto versions. Confirm the local data/output locations and Git exclusions.

**Why:** Later results must be traceable to intact inputs and a working analysis environment.

**Output / completion check:** A dated setup record, checksum outcome, helper-test outcome and one successfully rendered notebook. Setup checks establish technical readiness, not biological validity.

## Step 2 — Audit sample identities and measurement quality

**Do**

- [x] Load the sample sheet, annotation, abundance and DEA tables with the helpers. Inspect actual columns and identifier uniqueness before joining.
- [x] Verify 27 unique samples and three samples in each genotype-by-medium group; flag any discrepancy. Confirm the primary subset contains 18 samples.
- [x] Verify repaired `pulseSN`/`SNpulse` names and consistent sample mapping across layers. A repeated replicate number alone does not establish cross-layer pairing.
- [x] Summarize detected proteins/sites and missingness per sample and group. Inspect abundance distributions, sample correlations and PCA using `01_explore_data.qmd` as a starting point; document how missing values are handled for exploratory plots.
- [x] Identify outliers, potential batch structure and reinjected pulse phospho samples. Retain samples by default; document any proposed exclusion, its independent evidence and a sensitivity comparison.
- [x] Record the search, `modelName`, `estimate_type` and available uncertainty fields for each statistical table.

**Why:** Poor sample quality, mislabeled samples and unequal detection can mimic a genotype effect.

**Output / completion check:** Sample-design summary, QC plots and a written decision on sample inclusion. Investigate serious anomalies with the lecturer before making claims that depend on those samples.

## Step 3 — Validate Stp1 and Stk1 first (H1)

**Do**

- [ ] Use `03_mutant_validation.qmd` as the starting point. Write the predicted enzyme profiles before examining results.
- [ ] Plot individual normalized protein abundances for Stp1 (`B7H15_06350`) and Stk1 (`B7H15_06355`) across all nine groups. Start interpretation with TSB, then check consistency in SN and SNpulse.
- [ ] Count detections out of three samples per group and examine unexpected high mutant signals. Record missingness and outliers explicitly.
- [ ] Inspect corresponding mutant-versus-parental contrasts within each medium. Keep `lod_imputed` estimates separate: they describe detection patterns, not precise loss magnitudes.
- [ ] If unexpected signal persists, ask about insertion position, peptide uniqueness/coverage, possible protein fragments, carry-over or sample identity. Use peptide-level evidence only if it is available; protein summaries alone cannot prove full-length functional protein.

**Why:** Downstream `stp_vs_anc` effects cannot be interpreted as phosphatase-only effects if Stk1 is also reduced.

**Output / completion check:** Figure 1 enzyme-abundance panels and a short evidence table for each enzyme/strain: supported reduction, inconsistent pattern, or insufficient information. Do not describe nondetection as proven absence or assume a clean knockout.

## Step 4 — Examine the possible polar effect

**Do**

- [ ] Inspect JE2 genome coordinates, strand/orientation and neighboring genes around `stp1` and `stk1`; draw a small annotated gene diagram.
- [ ] Cite the operon/transcription-unit evidence for *stp1*–*stk1* from AureoWiki (look up the NCTC 8325 orthologs via `staa8_uniprot` / `read_protein_metadata()`), and the NE217/NE1919 insertion positions from NTML (ntml.unmc.edu) and `TRANSPOSON_INSERTIONS` in `R/constants.R`. Record which strain the AureoWiki evidence refers to (NCTC 8325, not JE2).
- [ ] When the peptide table (`peptides_total.tsv.gz`, needs `msstats.csv`) arrives, use `plot_peptide_coverage()` / notebook 07 to see which part of Stp1/Stk1 is still detected.
- [ ] Plot measured neighboring proteins and assess whether any accompanying changes are consistent across samples and conditions.
- [ ] Compare the pattern with the verified insertion/polar-effect evidence in `mval_medi/docs/stk_stp_rescue_evidence_mval.md`. Check the original source before using a specific biological claim.
- [ ] List alternatives: protein stability/regulation, secondary mutations, peptide assignment and sample problems. Retain observations that contradict polarity.

**Why:** Gene proximity and reduced Stk1 can support a possible polar explanation, but do not demonstrate shared transcription or the causal mechanism.

**Output / completion check:** Gene-context panel and a paragraph separating measured protein changes, the proposed explanation and unresolved alternatives. If Stk1 reduction is unsupported, revise the interpretation and continue the response comparisons without assuming combined disruption.

## Step 5 — Fix comparison definitions and analysis rules

**Do**

- [ ] Separate baseline contrasts (`stk_vs_anc_at_TSB`, `stp_vs_anc_at_TSB`) from within-genotype exposure contrasts and the three interaction contrasts below.
- [ ] Use `explain_contrast()` and check one example against group means to verify the direction. Say which group is higher whenever reporting a simple fold change.
- [ ] Begin with supplied DEA results. Document the tested feature universe, filtering, missingness and original FDR adjustment scope. Do not silently refit or renormalize.
- [ ] Use the course screen **FDR < 0.1 and absolute log2 effect > 1**, alongside effect estimates and uncertainty. It is a screen, not a formal test that the true effect exceeds one log2 unit.
- [ ] Define comparison-specific quantitative eligibility before selecting hits: initially require detection in at least two of three samples in every group needed for that comparison, and observed rather than `lod_imputed` estimates. Report eligible and excluded counts; retain sparse/one-sided signals separately.
- [ ] Keep supplied FDR values labeled as supplied. Declare any additional testing families if new tests are performed; do not imply the original FDR controls every later selection or combined claim.

| Contrast | What it answers |
|---|---|
| `interaction_stk_vs_anc_at_TSB_vs_SN` | Does the stk response differ from the parental response? |
| `interaction_stp_vs_anc_at_TSB_vs_SN` | Does the stp response differ from the parental response? |
| `interaction_stp_vs_stk_at_TSB_vs_SN` | Do the two mutant responses differ directly? |

For example, the stp-versus-parent interaction is `(stp_TSB - stp_SN) - (anc_TSB - anc_SN)` on the log2 scale. Positive means the **SN-minus-TSB response is less positive or more negative in stp** than in anc. All four genotype-condition groups are needed to assess this interaction.

**Why:** Baseline differences, exposure effects and altered exposure responses answer different questions. Clear rules limit selection based on attractive results.

**Output / completion check:** A methods/contrast sheet that states direction, eligible features, model and statistical thresholds. Use standard errors for uncertainty; construct intervals only with the appropriate model information and label any approximation.

## Step 6 — Identify altered total-protein responses (H2)

**Do**

- [ ] Summarize baseline mutant differences separately from TSB-to-SN changes within each genotype.
- [ ] Examine the two mutant-versus-parent interaction contrasts across eligible proteins. Report effect distributions, eligible counts and screened-hit counts.
- [ ] Retain effect estimates, standard errors, FDR, model and observed/imputed status for each selected candidate.
- [ ] Plot per-sample values for representative candidates before interpreting them. Use consistent axes and explain the response difference in words.
- [ ] Keep a separate exploratory summary of one-sided detection patterns; do not rank large imputed fold changes as precise quantitative effects.

**Why:** This establishes whether a mutation changes the exposure response beyond any pre-existing TSB difference.

**Output / completion check:** A reproducible interaction summary and a small set of inspected candidate plots. If few or no proteins pass the screen, report that outcome with effect sizes, coverage and uncertainty; do not relax cutoffs just to generate hits.

## Step 7 — Compare shared and distinct mutant responses (H3)

**Do**

- [ ] Compare the two mutant-versus-parent interaction estimates for the same eligible proteins, with direction and uncertainty visible. Record proteins eligible in one comparison but not the other.
- [ ] Describe same-direction changes and opposite-direction changes. Label significance in both comparisons as a descriptive classification with its testing limits.
- [ ] Inspect the direct `stp_vs_stk` interaction for candidates and report its uncertainty. A hit in one mutant but not the other is not by itself evidence of a difference between mutants.
- [ ] Use a heatmap and individual-value plots to show representative shared and distinct patterns. State any centering/scaling used only for display.
- [ ] Explain that correlations share the parental reference, and that a nonsignificant direct contrast does not establish equivalence.

**Why:** Similar response components could be consistent with reduced kinase function, while direct comparisons reveal supported differences.

**Output / completion check:** Figure 2 interaction-effect comparison and Figure 3 selected response patterns, with a conclusion about similarity limited to the measured features. Similarity does not demonstrate rescue or a return to normal.

## Step 8 — Assess whether phosphosite coverage supports a follow-up

**Do**

- [ ] Use `04_phospho_two_roadmaps.qmd` as the starting point. Keep STY and STYH separate: they are alternative searches with different supplied models, not extra biological replicates.
- [ ] Audit site identifiers, JE2 protein mapping, residue, position and peptide. Verify join cardinality; retain ambiguous mappings instead of silently duplicating rows.
- [ ] Count detections per group and apply the two-of-three rule in every group required for the chosen comparison. Inspect `plot_sites()` values and the matching normalized phospho abundances used for DEA.
- [ ] Choose and record a localization rule before site selection. Assess per-sample `loc_prob` where available; `best_loc_prob` alone is only the best score across samples. A conventional STY cutoff such as 0.75 is not a measured false-localization rate or a validated histidine cutoff.
- [ ] Choose the primary search for follow-up based on usable coverage and localization, with the reason documented. Treat the other search as a sensitivity check and record model differences.
- [ ] Keep histidine phosphorylation optional and exploratory, with chemical/localization limitations. Do not expand this into a separate main project.

**Why:** An apparent phosphosite effect is useful only if the site assignment and measurements support the comparison.

**Output / completion check:** Coverage and mapping audit plus a supported-site shortlist. If reliable coverage is inadequate, retain a descriptive phospho section and finish the main total-protein project without a quantitative phosphorylation claim.

## Step 9 — Match selected phosphosites to total-protein changes

**Do**

- [ ] Match each supported site to its total protein using the JE2 identifiers and **exactly the same contrast**. Check sample pairing before any sample-level adjustment; retain unmatched/ambiguous sites separately.
- [ ] Plot both measurement layers. Ask whether the phosphosite change could track a change in protein amount.
- [ ] For supported observed effects only, calculate the descriptive adjusted effect `diff_site - diff_protein`. Do not subtract an imputed protein or site estimate.
- [ ] Label adjustment as relative protein-adjusted phosphorylation signal, not absolute occupancy. Do not reuse the original site p-value/FDR for the new adjusted effect.
- [ ] If a formal adjusted test is justified, account for covariance: `Var(site - protein) = Var(site) + Var(protein) - 2*Cov(site, protein)`. Document any independence assumption; otherwise keep the adjustment descriptive.

**Why:** More phosphosite signal can reflect more protein rather than altered phosphorylation regulation.

**Output / completion check:** Figure 4 paired protein/site plots for a small supported module, if available, and clear categories: consistent with protein amount, residual adjusted change, or unresolved. These categories do not identify direct kinase/phosphatase substrates.

## Step 10 — Interpret functions and verify the biological story

**Do**

- [ ] Annotate reliable shared and distinct candidates using JE2 annotations, bacterial STRING mappings and verified protein/literature evidence. Record mapping ambiguity and accession provenance.
- [ ] Use `mval_medi/docs/stk_stp_rescue_evidence_mval.md` to distinguish complementation, genetic suppression and compensation. Verify cited claims in the source and note differences from our strain, exposure and measurements.
- [ ] If candidate numbers and mapping coverage permit enrichment, define the foreground and the **eligible measured/tested background under the same filters**. Analyze directions separately where meaningful; deduplicate proteins and term memberships.
- [ ] Use a separate phosphoproteome background for phospho candidates and count a protein only once despite multiple sites. Record mapping losses, tested terms, term-size filters and multiple-testing correction.
- [ ] If coverage is poor, provide an annotated candidate summary instead of an uninformative enrichment test.
- [ ] Run enrichment with `Rscript scripts/10_enrichment_mval.R <contrast> ...` (helpers in `R/helpers_mval.R`, tested in `tests/testthat/test-helpers_mval.R`). It does both methods from lecture 06 (slide 25):
  - **ORA:** hits FDR < 0.1 and |log2FC| > 1, up and down separately, one-sided Fisher test. The background is the JE2 proteins tested in that contrast that carry ≥ 1 term of the category. BH-FDR per category, term size 5–500. `lod_imputed` hits are included and counted (`hits_imputed`).
  - **Rank-based (GSEA-like):** observed estimates only, ranked by the moderated t `statistic`, Wilcoxon test per term. `auc` > 0.5 = shifted towards higher in the first group of the contrast.
  - For the real Step 10 run, restrict hits *and* background to the Step 5 eligible proteins, and use the Step 6–7 candidate sets as foreground.
- [ ] Paste the written `string_upload_*_mval.txt` lists into STRING ("Multiple proteins", organism `STRG0B09NBK`) for network pictures. Note that the website's default background is the whole genome; set a custom background where possible, or report the difference.
- [ ] Plausibility check, done 11 Nov on `stk_vs_anc`: our local method reproduces the lecture's example (slides 36–38). Cytolysis / killing of cells and copper ion transport come out lower in `stk`; vancomycin resistance comes out higher in `stk`. The "Virulence" keyword (KW-0843) did not reproduce with our measured background. Treat it as background-dependent.

**Why:** Functional themes can help explain results, but annotation coverage and the background determine what an enrichment analysis can support.

**Output / completion check:** A source-linked interpretation table and, if justified, an enrichment summary with overlaps and denominators. Distinguish measured changes from inferred mechanisms; do not infer growth rescue, survival, antibiotic tolerance or direct substrates from these omics alone.

## Step 11 — Check robustness and resolve contradictions

**Do**

- [ ] Reinspect the individual samples behind each headline conclusion, including contradictory features and unexpected enzyme signals.
- [ ] Check whether conclusions rely on one replicate, imputed values, ambiguous mapping, localization choices or marginal thresholds. Where appropriate, compare reasonable prespecified alternatives and document why they were chosen.
- [ ] Keep any justified sample-exclusion or alternative-search analysis labeled as sensitivity analysis. Refit only if necessary, with explicit design, filtering, batch and missingness assumptions.
- [ ] Treat SNpulse as optional after the primary story is complete. Investigate the two reinjected stk pulse phospho measurements and batch/genotype confounding first; do not claim a longitudinal recovery trajectory.

**Why:** A conclusion should survive reasonable checks and acknowledge where the design cannot resolve an alternative explanation.

**Output / completion check:** A short robustness/limitations section stating which conclusions are stable, tentative or unresolved. With three samples per group, nonsignificance is weak evidence for no effect.

## Step 12 — Assemble and verify the final deliverables

- [ ] Create short project notebooks in `analyses/` for QC/validation, response comparisons and phospho follow-up, reusing existing helpers and templates. Keep each notebook independently renderable in under one minute.
- [ ] Write `reports/project_c_report.qmd` and render a final HTML report with explicit setup paths/output settings. Include the question, design, methods, figures, H1-H3 answers, limitations, citations and software versions.
- [ ] Deliver Figure 1 (gene context and enzyme profiles), Figure 2 (interaction comparisons), Figure 3 (selected shared/distinct responses), and Figure 4 only if protein/site coverage supports it.
- [ ] Give each hypothesis one conclusion: supported within stated limits, unsupported by the measurements, or unresolved. Separate expected results from actual findings.
- [ ] Include a short follow-up proposal: verify insertion/genotype, measure `stk1`/`stp1` transcripts, and compare clean single/double deletions plus complemented strains to test polarity/suppression causally. These are future experiments, not prerequisites for completing this course analysis.
- [ ] Check one numerical example by hand, render all new/changed reports, inspect every figure and confirm contrast directions. Run `Rscript tests/run_tests.R` if `R/` helpers change.
- [ ] Update the lab notebook. Review the Git diff and staged-file list before any commit or publication; include code/reports and exclude course measurement files and unintended files.
- [ ] Prepare the final presentation from verified findings, including negative results and limitations. Confirm course submission requirements with the group/lecturer.

**Why:** The group must be able to explain, reproduce and defend what it submits.

**Completion standard:** A reproducible report and presentation answer H1-H3 as far as the data allow, every headline finding has an inspected sample plot and traceable method, and unsupported causal claims are excluded.

## Suggested order for the group

| Work block | Steps | What must be ready before moving on |
|---|---|---|
| Foundation | 1-2 | Working environment, intact inputs, understood sample design and QC concerns. |
| Mutant identity | 3-4 | Enzyme profiles assessed; polar-effect explanation clearly separated from observations. |
| Main response analysis | 5-7 | Defined contrasts and eligibility; shared/distinct effects checked against individual samples. |
| Supporting biology | 8-10 | Phospho coverage decision, supported protein/site matches and sourced interpretation. |
| Completion | 11-12 | Robustness checks, readable figures, rendered report and presentation. |

If time is short, finish Steps 1-7 and the relevant interpretation/reporting in Steps 10-12. Add quantitative phospho adjustment, enrichment and SNpulse only when coverage and time justify them. A complete, well-supported total-protein answer is sufficient to address the core mutant-validation question.

## Current status and next action

**11 November 2026:** Steps 1–2 are done (`scripts/01_setup_check_mval.R`, `scripts/02_sample_qc_mval.R`; notebook `2026-11-11_steps1-2_setup_qc_mval.md`). All 27 samples are kept. The peptide table is still missing (needed for Step 4). The enrichment tooling for Step 10 (`scripts/10_enrichment_mval.R`) is built, and it was checked against the lecture's `stk_vs_anc` example. This is a method check only; the Step 10 interpretation still comes after Steps 5–7.

**Next action: Step 3, inspect Stp1/Stk1 individual measurements (`plot_protein()`, `03_mutant_validation.qmd`).**
