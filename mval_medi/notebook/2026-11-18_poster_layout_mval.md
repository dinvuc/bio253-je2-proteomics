# 2026-11-18 – A0 poster layout after the poster lecture (group mval: [names])

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
- `scripts/12b_poster_a0_figures_mval.R`: figures redrawn at their printed size (one column =
  243 mm), Arial ≥ 22–24 pt, no gridlines, direct labels, vector PDF + 300 dpi PNG, QR code
  (installed the R package `qrcode` from CRAN). The Step 12 report figures are unchanged.
- `mval_medi/results/step12_poster/poster_a0_mval.typ` (Typst, bundled with Quarto 1.4+):
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
checked the reference in PubMed, and drafted the poster guide and this entry. We still have to
check every number and all wording before printing.

## Next steps
- Fill in author names and e-mail (red on the poster).
- Check every number against the figure guides; scan the QR code; print an A4 test copy.
- Rehearse the 10-s pitch and the 1/2/5-min tours (`poster_guide_mval.md`).
- Ask Jonas about printing and the deadline. Open a pull request into `main` for a group review.
