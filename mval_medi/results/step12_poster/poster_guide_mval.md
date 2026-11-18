# A0 poster: how it is built, what the lecture asks, how to present it

**Poster:** [`poster_a0_mval.pdf`](poster_a0_mval.pdf) (A0 portrait, 841 × 1189 mm, vector PDF).
**Source:** [`poster_a0_mval.typ`](poster_a0_mval.typ) (Typst, all text and layout) +
[`scripts/12b_poster_a0_figures_mval.R`](../../../scripts/12b_poster_a0_figures_mval.R) (figures in
[`poster/`](poster/), QR code). The Step 12 report figures (`fig1`…`fig4`) stay as they are; the
final report uses them.

```bash
Rscript scripts/12b_poster_a0_figures_mval.R                              # ~25 s, 9 checks must say PASS
quarto typst compile mval_medi/results/step12_poster/poster_a0_mval.typ   # -> poster_a0_mval.pdf (1 page)
```

To change a word, edit the `.typ` file and compile again. Typst is a small markup
language: `*bold*`, `_italic_`, `#stp` prints a coloured strain name. A `~` must be written
`\~`, otherwise Typst prints a space.

## Before printing (group)

- [ ] Fill in the author names and the contact e-mail (printed in **red** until done).
- [ ] Check every number on the poster against the figure guides of Steps 3–12
      (fold change 43×, 5 / 1 altered, 1,882 / 1,931 tested, r = 0.85 / 0.86, r = 0.37 / 0.38,
      noise-only r = 0.50, 20 + 4 direct hits, 236 sites, term counts, 2,232 proteins).
- [ ] Scan the QR code with a phone: it must open the repository.
- [ ] Print the PDF on A4 and do the arm's length test (below). Have someone outside the group
      read it: are the question and the main message obvious?
- [ ] Ask Jonas: print deadline and whether the course prints it or we do.

## Lecture checklist ("Academic posters", BIO253, Nov 2026) → what we did

| Lecture recommendation | On our poster |
|---|---|
| Focus on one message; the title carries most of it | The title states the main result (*stp1* mutant also loses Stk1); the subtitle states the second (PASN response intact). The old title was a question. |
| Headings summarise the work; a hurried reader gets the main points from headings alone | Every section heading is a result: "One operon, two mutants", "The stp mutant also loses Stk1", "Both mutants keep their response to PASN", "stp resembles stk only in part", "Two gene blocks shift as a unit". |
| Grid, three columns, reader gravity (top→bottom, left→right), organisation cues | 3 equal columns (243 mm, 20 mm gutters, 35 mm side margins); numbered sections 1–5 in reading order; H1–H3 in section 1 point to the sections that answer them (→ 2, 3, 4). |
| Hot zone top left, take-home messages bottom right, fine print last | Question and design top left; take-home box (coloured frame) bottom right; methods, limits, next steps, reference and QR code below it. |
| Hierarchy: the more important, the larger the type | Title 104 pt → subtitle 56 pt → headings 46 pt → take-home 34 pt → body 28 pt → "how to read" lines 26 pt → fine print 25 pt. All ≥ 24 pt, as the lecture asks. |
| Graphics dominate; text sparingly; phrases, active voice, left-aligned | About 80 % of the area is figures. Text is 1–2 lines per figure ("Each dot = one protein …"). The long text blocks (old `poster_text_mval.md`) are cut; the details stay in that file for questions. No justified or centred text. |
| Show, don't tell (visual grammar) | Gene sketch with the two insertion sites (from `TRANSPOSON_INSERTIONS`) instead of a sentence; design shown as a strains × media grid of dots (● = one culture). |
| Text in figures must be large, too | Figures are drawn at their printed size (243 mm wide), so ggplot text sizes are real sizes: ≥ 22–24 pt, Arial (the poster font). Previously the figures had ~10 pt text and would have been scaled arbitrarily. |
| Simple, clean graphs; explanations on the figure; no generic titles, no unnecessary gridlines | No gridlines (`theme_classic`). Direct labels instead of legends: strain names on the x axis, "43× lower" on the Stk1 panel, "n.d." on the y axis for not detected, "same response as anc" on the diagonal, the r value, the number of hits and the block names in the panels. Generic figure titles are removed; the poster heading carries the message. |
| Few complex graphs better than many simple ones; label only what matters | Fig 3 labels only the 3 "altered in both" proteins and the 4 response-layer hits; the 20 baseline hits are counted, not named. Fig 4a shows totals as horizontal bars (readable labels). |
| Lines are hard to see from a distance; filled shapes are better | The genome scan is a filled area instead of a 0.5 pt line; all lines ≥ 0.8 mm; points ≥ 5 mm. |
| Consistent fonts, colours, widths | Arial only. Strain colours are the same everywhere, in text and figures (anc blue, stk orange, stp green; colour-blind safe Okabe–Ito). One accent colour (orange-red) for the take-home box and "differ directly". |
| Resolution / vector formats (PDF/SVG over JPG/PNG), 150 dpi prints well | Figures embedded as vector PDF, QR code as SVG, so the poster is fully vector and sharp at any size. 300 dpi PNG previews are written next to them. |
| Leave white space; balance text and graphics | 20 mm gutters, 35/30 mm margins, space between sections; columns end at about the same height. |
| Know your audience: avoid jargon, explain abbreviations | PASN is defined in section 1; "Tn", "n.d." are explained where they first appear; "interaction" is explained in words ("respond differently"); statistics jargon is in the fine print. |
| Be bold: strongest statement the data allow | Title states the Stk1 loss as a result (it is measured). The *mechanism* (polar effect) stays a question ("is Stk1 lost too?"), because protein data cannot prove it. |
| Include e-mail, foster contact; link to more | QR code + URL to the repository (code and reports only, no raw data); contact e-mail field. |
| Cite sources | Mutant library: Fey et al. 2013, *mBio* 4:e00537-12 (PMID 23404398, doi:10.1128/mBio.00537-12). |
| A0 portrait (BIO253) | `paper: "a0"`, portrait, one page. |

## Arm's length test (done 18 Nov on the screen, repeat on paper)

At A4 size (1/4 of A0) the title, headings, take-home box, axis labels and figure annotations
are readable. The fine print is small by design. Print on A4 and check by eye before ordering the A0.

## Presenting the poster

**Use the poster as a visual aid: point at the figures, do not read the text.**
Present in the order of the numbers. Everyone in the group should be able to give each tour.

### 10-second pitch (someone walks by)
"We knocked out the phosphatase Stp1 of *S. aureus* with a transposon. The mutant turns out to
lose the kinase Stk1 as well. Still, both mutants respond almost normally to *Pseudomonas*
spent medium."

### 1-minute tour
1. **Section 1:** Stk1 (kinase) and Stp1 (phosphatase) sit in one operon. The *stp1* insertion lies
   upstream of *stk1*. Question: does the *stp* mutant lose Stk1 too?
2. **Section 2:** Yes. Point at "43× lower" and at the peptide panel (1–2 peptides instead of ~30).
   ErmB confirms both insertions.
3. **Section 3:** Point at the diagonal: almost every protein responds to PASN as in the wild
   type. Only 5 and 1 proteins differ.
4. **Take-home box.**

### 2-minute tour
The 1-minute tour, plus:
- **Section 4:** *stp* and *stk* share part of their pattern (r ≈ 0.37 after removing the
  correlation that the shared wild-type reference creates by itself), but 24 proteins differ directly.
- **Section 5:** Many pathway terms shift in both mutants. Part of the signal comes from two
  blocks of neighbouring genes; that can be a strain difference, which we would check by sequencing.

### 5-minute tour (interested visitor)
The 2-minute tour, plus the details from [`poster_text_mval.md`](poster_text_mval.md) and the
figure guide:
- Hit rule: limma `treat` tests |log2FC| > 1 directly (not FDR first, then a fold-change cut-off),
  BH-FDR < 0.05, fixed before we saw results.
- Robustness: refits leaving out each sample; filled vs open circles in section 3.
- Why "not detected" is not zero (n.d. row in section 2).
- Limits: n = 3, no complemented or clean deletion strains, sparse phosphosites, PASN = secreted
  molecules + nutrient depletion (not co-culture).
- Next experiments (fine print).

### Likely questions
- *"Is it really a polar effect?"* Consistent with it, not proven: we measured proteins, not
  transcripts. Next: *stk1* transcripts and a clean Δ*stp1* with complementation.
- *"Why are some Stp1 values in stp wild-type-like?"* They rest on 1–2 peptides (section 2b), so
  they are unreliable.
- *"What is PASN?"* Cell-free, sterile-filtered spent medium of *P. aeruginosa* PAO1 (30 %,
  ~15 h). *S. aureus* never met living *Pseudomonas*.
- *"Phosphorylation?"* 236 usable sites; none responds differently to PASN in a mutant; the
  Stk1 phosphosite drops together with the Stk1 protein.
