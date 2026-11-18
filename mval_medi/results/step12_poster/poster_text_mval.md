# Poster text: long version (for questions at the poster)

> **Since 18 Nov the poster text lives in [`poster_a0_mval.typ`](poster_a0_mval.typ)**
> (rendered: [`poster_a0_mval.pdf`](poster_a0_mval.pdf); see [`poster_guide_mval.md`](poster_guide_mval.md)).
> Following the poster lecture, most of the text below was cut from the poster ("edit
> ruthlessly"). It stays here as the long version: details to explain in the 5-minute tour and
> to answer questions. Edit wording on the poster in the `.typ` file, not here.

Report figures: `fig1_h1_enzymes_mval.png` … `fig4_pathways_genome_mval.png` in this folder
(300 dpi, made by [`scripts/12_poster_figures_mval.R`](../../../scripts/12_poster_figures_mval.R)).
Poster figures: [`poster/`](poster/) (made by [`scripts/12b_poster_a0_figures_mval.R`](../../../scripts/12b_poster_a0_figures_mval.R)).
Every number below is in a figure guide of Steps 3–11; check it there before printing.

## 1. Question and design
- **Question:** Does the *stp1*::Tn mutant also lack the kinase Stk1 (polar effect), and do
  the two mutants respond like the parent to *P. aeruginosa* spent medium (PASN)?
- **H1** `stp` has reduced Stp1 **and** Stk1. **H2** Some PASN responses differ from `anc`.
  **H3** Part of the `stp` pattern resembles `stk`.
- **Design:** *S. aureus* JE2 (`anc`), NE217 *stk1*::Tn (`stk`), NE1919 *stp1*::Tn (`stp`);
  TSB vs ~15 h in 70 % TSB + 30 % cell-free, sterile-filtered PAO1 spent medium (SN);
  n = 3; total proteome (2,232 JE2 proteins) + Fe-IMAC phosphoproteome. PASN contains only
  secreted molecules; *S. aureus* never met living *P. aeruginosa*.

## 2. Methods box
- Normalised log2 abundances (prolfquapp); missing = not detected, never 0.
- limma, model `~ 0 + group` (9 means), contrasts written as A − B.
- **Hit = `treat` test of |log2FC| > 1, BH-FDR < 0.05**; tested only if ≥ 2 of 3 values in
  every group of the contrast. Fixed before results were seen (12 Nov).
- Interaction = (mutant TSB − SN) − (anc TSB − SN).
- H3 similarity corrected for the shared `anc` reference (noise-only r = 0.50).
- Phospho: ≥ 2 of 3 per group and median per-sample localisation ≥ 0.75; site minus
  protein effect in the same contrast.
- Pathways: rank test + limma cameraPR; robust = both FDR < 0.05.
- Robustness: leave-one-sample-out, alternative thresholds, genome-position scan.

## 3. Results (one line per figure)
- **Fig 1 (H1):** Stk1 is absent in `stk` and ~43-fold lower in `stp` (1–2 of ~30
  peptides left). Stp1 is down to 0–2 peptides in `stp`. The transposon marker ErmB is in
  all mutant samples and in no `anc` sample.
- **Fig 2 (H2):** Both mutants keep most of the PASN response (r = 0.85 / 0.86). Only
  5 (`stk`) and 1 (`stp`) of ~1,900 proteins respond differently. 2 of the 5 `stk` hits
  survive every leave-one-out refit and lie outside shifted genomic blocks (argF, B7H15_14515).
- **Fig 3 (H3):** Beyond shared noise, `stp` and `stk` effects correlate moderately:
  r = 0.37 (TSB baseline) and 0.38 (response). 20 + 4 proteins differ directly between the mutants.
- **Fig 4:** Many pathway terms shift in both mutants (`stp` as many as `stk`). Two genomic
  blocks change as a unit: SigB/*mazEF*/*kdp* (+1.5 log2 in `stk`) and a prophage
  (−2.1 log2 in `stp`).
- **Phospho (text only):** 236 usable sites; none respond differently to PASN in a mutant.
  The Stk1 phosphosite drops with the Stk1 protein.

## 4. Answers
- **H1 supported:** `stp` lacks most Stk1 as well as Stp1. This fits a polar insertion,
  but protein data cannot prove the mechanism.
- **H2 supported, but small:** few altered responses; the PASN response is largely intact.
- **H3 partly supported:** a moderate shared component, plus clear differences between the mutants.

## 5. Limitations
n = 3; no complemented or clean deletion strains; phosphosite coverage is sparse; blocks of
neighbouring proteins (copy number? prophage?) are unexplained without DNA data; PASN mixes
secreted molecules with nutrient depletion.

## 6. Next experiments
Genome sequencing or qPCR of the *stk* block and the prophage region; *stk1*/*stp1*
transcripts; clean Δ*stp1*, Δ*stk1* and double deletions with complemented strains.

## 7. Code
All figures: github.com/dinvuc/bio253-je2-proteomics (code and reports only, no raw data).
QR code: `poster/qr_repo_mval.svg` (R package `qrcode`, made by Step 12b).
