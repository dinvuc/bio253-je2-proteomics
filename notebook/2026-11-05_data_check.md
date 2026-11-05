# 2026-11-05 - Course data completeness check

## Question for today
Does gh_data.zip contain the inputs needed for topic H (phosphosite changes beyond protein abundance), with optional Stk/Stp comparisons?

## What we did
- Compared archive contents with data/README.md and data/checksums.md5.
- Verified all expected files against the repository's MD5 checksums before extraction.
- Extracted only the 22 expected files into the local, git-ignored data folder, preserving any existing matching files.
- Inspected derived-table headers and the sample metadata. No hypothesis tests or biological analyses were performed.

## Results (key numbers)
- 22 expected files present, 22 matching checksums, zero missing or mismatched files.
- 27 unique samples: anc, stk and stp in TSB, SN and SNpulse, three samples per combination.
- Total-proteome and both STY/STYH phospho differential-analysis and abundance tables are available.
- Differential-analysis tables include contrasts, effect estimates, standard errors, adjusted significance values and observed/imputed estimate labels.
- Site tables include residues, positions, localization probabilities and per-sample intensities. Protein identifiers and annotations support matching the layers.
- The bundle is complete for the proposed analysis of processed course data. Availability of usable observed site-protein pairs remains to be assessed.

## Problems / surprises
- Windows Rscript was not found on PATH or in the standard Program Files/R location during this check. Base PowerShell/.NET was used for archive integrity and metadata inspection only; biological analysis will use R and the repository helpers.
- The bundle contains original course result tables, not instrument-level mass-spectrometry raw files. Reprocessing spectra is outside the proposed project's scope.
- STY and STYH are alternative searches, not additional biological replicates.

## AI use
Codex inspected repository documentation, archive contents, checksums and metadata, and recorded this entry. No raw-data tables were copied into the notebook or published.

## Next steps
- Locate or configure R and required packages.
- Start with TSB_vs_SN_at_anc, check per-sample detection and observed versus imputed estimates, then quantify reliable site-protein matches.
- Correct only supported observed measurements and retain unmatched sites separately.
- If coverage permits, extend selected candidates to genotype-by-condition comparisons, accounting for documented mutant and reinjection caveats.

## Follow-up: alternative project questions
- Considered five extensions of the lecturer's topics: distinct responses of different sites within one protein; mutant baseline profiles resembling the parental supernatant response; opposite responses to short versus prolonged exposure; coordinated defense and virulence-associated protein changes; and alternative cystine-importer abundance patterns.
- These are candidate hypotheses, not findings. Feature coverage and observed replicate measurements have not yet been assessed for these questions.
- Prioritized within-protein site differences for a focused extension of topic H, and baseline-versus-response comparisons for an alternative biological question.
- Important interpretation limits: stp also has reduced Stk1; short and prolonged exposures are different growth histories, not a longitudinal time course; two stk pulse phospho runs were re-injected. No claims of direct substrates, survival, transporter flux or resource allocation follow from these omics data alone.
- Codex used the course topic slides, repository design notes and the Prust/Schwyter papers to develop these suggestions. No additional data analysis was run.

## Project decision and rescue literature review
- The user selected revised topic C: validate Stp1/Stk1 protein profiles, investigate a possible polar effect, and compare genotype-by-condition responses. This supersedes the earlier provisional topic-H project; protein correction is a supporting analysis.
- Saved PROJECT_PLAN.md and linked it from README.md. Saved primary-source evidence and strain/assay limitations in docs/stk_stp_rescue_evidence.md.
- Huemer et al. 2023 provides double-deletion evidence relevant to suppression under acid stress. Sun et al. 2012 reports a polar stp1 insertion affecting both genes. Burnside et al. 2010 documents gene-add-back complementation; Beltramini et al. 2009 shows double loss is not universally restorative. See the evidence review for full citations.
- Explicitly distinguished complementation, genetic suppression and compensation. Our course dataset lacks clean stp-only and complemented rescue comparisons, so it cannot demonstrate causal genetic rescue.
- Codex searched primary literature and checked the course Huemer PDF. No hypothesis tests or new biological results were generated. No raw data were staged or uploaded.

## Cross-PC continuity
- User requested the same current chat on both PCs and reports a dataset copy on the other PC.
- Confirmed the current chat is local Codex, named Assess BIO253 datasets and projects. Only the current local host is exposed to the app tools.
- Documented supported remote access, its host-online requirement, separate Git file updates, and the distinction from a static share link. Pairing and second-PC access remain unverified and require user setup in the app.
- Published the selected project documents and connection guide on GitHub main so the second checkout can obtain the same saved decisions. Raw data remain excluded. This does not synchronize chat history; device pairing remains unverified.

## Updated project questions and hypotheses
- At the user's request, clarified PROJECT_PLAN.md as the current project reference for both PCs.
- RQ1/H1: assess the actual Stp1/Stk1 protein profiles, including accompanying Stk1 reduction in stp.
- RQ2/H2: test genotype-dependent TSB-to-SN responses through interaction contrasts.
- RQ3/H3: evaluate shared kinase-mutant-like and distinct stp responses, including the direct stp-versus-stk interaction.
- Kept the proposed rescue/suppression interpretation explicitly separate from what this three-genotype observational comparison can demonstrate.
- Codex edited documentation and checked Git formatting and publication. No biological analysis or new results; course measurements remain excluded from GitHub.
