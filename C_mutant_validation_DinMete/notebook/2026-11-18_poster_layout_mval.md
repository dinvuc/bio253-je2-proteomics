# 2026-11-18 – A0 poster layout after the poster lecture (group mval: Din Vuckic & Mete Kabranlar)

## Question for today
Apply the lecture "Academic posters" (BIO253, J. L. Ruiz, Nov 2026) to our poster: turn the
Step 12 figures and text blocks into a readable A0 portrait poster.

## What we did
- New branch `poster-guidelines-mval`.
- Plan (coding rules, section 4), written at the top of the new script:
  aim = figures readable at A0; input = the same saved step results as Step 12;
  analysis = no new statistics; output = `step12_poster/poster/*` + `poster_a0_mval.pdf`;
  check = 9 PASS lines (sample counts, Stk1 fold change = Step 5, hit counts = Step 12,
  term counts = Step 10, genome-scan top window = Step 11).
- `C_mutant_validation_DinMete/scripts/12b_poster_a0_figures_mval.R`: figures redrawn at their printed size (one column =
  243 mm), Arial ≥ 22–24 pt, no gridlines, direct labels, vector PDF + 300 dpi PNG, QR code
  (installed the R package `qrcode` from CRAN). The Step 12 report figures are unchanged.
- `C_mutant_validation_DinMete/results/step12_poster/poster_a0_mval.typ` (Typst, bundled with Quarto 1.4+):
  layout and all poster text → `poster_a0_mval.pdf` (1 page, A0 portrait, fully vector).
- Lecture changes, in short: message title + subtitle instead of a question; numbered
  message headings; 3-column grid read top→bottom, left→right; take-home box bottom right;
  fine print last; text cut to 1–2 lines per figure; one font, consistent strain colours;
  filled genome-scan area instead of a thin line. Full checklist: `poster_guide_mval.md`.
- Arm's length test on a 1/4-scale preview: title, headings, take-home box and figure text
  readable; fine print small by design.
- Checked the transposon-library reference in PubMed: Fey et al. 2013, mBio 4:e00537-12,
  PMID 23404398, doi:10.1128/mBio.00537-12.

## Results (key numbers)
No new results. Every number on the poster comes from Steps 3–11 and is recomputed or
checked in the script: Stk1 `stp_vs_anc_at_TSB` log2 −5.42 (43× lower in `stp` than `anc`,
observed values); altered PASN responses 5 (`stk`, 1,882 tested) and 1 (`stp`, 1,931);
direct `stp_vs_stk` hits 20 (baseline) + 4 (response); robust terms `stk` 35/143/44,
`stp` 192/168/8.

## Problems / surprises
- Report figures at A0 would have had ~10 pt text and legends far from the data; redrawing at
  print size was simpler than scaling.
- Several labels collided or were clipped at 24 pt (facet strips, "43× lower", block labels in
  the genome scan); fixed by moving labels and stacking panels. Column 2 overflowed to a second
  page once; fixed by giving Fig 3 the empty height of Fig 2.
- In Typst a `~` is a non-breaking space ("~43×" printed as " 43×"); it must be written `\~`.

## AI use
Claude Code (AI agent) read the lecture PDF, wrote the figure script and the Typst layout,
checked the reference in PubMed, and drafted the poster guide and this entry. The numbers and
wording were marked as checked by the group on 18 Nov 2026.

## Next steps
- Fill in author names and e-mail (red on the poster).
- Numbers and QR code: marked as checked by the group on 18 Nov 2026. Print an A4 test copy.
- Rehearse the 10-s pitch and the 1/2/5-min tours (`poster_guide_mval.md`).
- Ask Jonas about printing and the deadline. Open a pull request into `main` for a group review.

## Addendum (18 Nov, later): Fig 3 (H3) label fix and redesign
- **Question:** could a visitor misread the "r = 0.37" on Fig 3?
- **Problem:** yes. The poster text gave the noise-only r (0.50) next to it, which invites
  "0.37 < 0.50, so no similarity". But 0.37 is the *noise-corrected* r of the true effects and
  must be compared with 0, not with 0.50.
- **What was run:** `Rscript C_mutant_validation_DinMete/scripts/12b_poster_a0_figures_mval.R`, then
  `quarto typst compile C_mutant_validation_DinMete/results/step12_poster/poster_a0_mval.typ`.
- **Changes to `p3_h3_shared_distinct_mval`:** r now shown with its bootstrap 95 % CI
  (baseline 0.37, 0.29–0.45; response 0.38, 0.24–0.52; read from
  `step07_shared_distinct/similarity_stats_mval.csv`) and the words "corrected for the shared
  anc noise"; r, CI and key moved to a text column next to each panel (inside the panel they
  covered points); key with counts (20 / 4 differ directly, 3 altered in both); the 3 "both"
  proteins numbered 1–3 and named in the key (fibronectin-binding, B7H15_02275, coagulase);
  direction of the axes under each panel title; axes ±5.6 (±5 clipped one protein at −5.36).
  New PASS checks in the script for the r values and the axis range.
- **Poster text:** the "how to read" line of section 4 now explains r in words
  (0 = unrelated, 1 = identical) instead of quoting the noise-only r; the dot key moved
  into the figure. Poster still one page.
- **AI use:** Claude Code redrew the figure, edited the poster line and this entry.
