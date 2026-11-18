// BIO253 poster, group mval: A0 portrait (841 x 1189 mm), the course format.
//
// Build from the repository root:
//   Rscript scripts/12b_poster_a0_figures_mval.R                                # figures, QR code
//   quarto typst compile mval_medi/results/step12_poster/poster_a0_mval.typ     # -> poster_a0_mval.pdf
//
// Layout follows the lecture "Academic posters" (BIO253, Nov 2026); see poster_guide_mval.md:
//   - the title states the main result; the headings state the section results
//   - 3-column grid, numbered sections 1-5, read top to bottom, then left to right
//   - take-home messages in the bottom-right box; methods and other fine print last
//   - one font (Arial), body >= 24 pt, headings 46 pt, title 104 pt
// The figures are made at the printed column width (243 mm), so their text prints at >= 22 pt.
// If you change the margins or gutters below, change COL_W_MM in the R script too.
// Every number on the poster is in a figure guide of Steps 3-12. Check it there before printing.
// TODO (group): author names and contact e-mail (printed in red until filled in).

#let C-anc = rgb("#0072B2")
#let C-stk = rgb("#E69F00")
#let C-stp = rgb("#009E73")
#let navy = rgb("#0b3c5d")
#let red = rgb("#D55E00")
#let grey = rgb("#333333")
#let todo(body) = text(fill: rgb("#cc0000"), weight: "bold", body)

#set page(paper: "a0", margin: (x: 35mm, y: 30mm))
#set text(font: "Arial", size: 28pt, lang: "en", hyphenate: false)
#set par(justify: false, leading: 0.45em, spacing: 0.9em)

// strain names always in their figure colour
#let anc = text(fill: C-anc, weight: "bold")[anc]
#let stk = text(fill: C-stk, weight: "bold")[stk]
#let stp = text(fill: C-stp, weight: "bold")[stp]

// numbered section heading: the number guides the reading order, the text is the result
#let sec(n, title) = block(below: 7mm, grid(
  columns: (24mm, 1fr), column-gutter: 7mm, align: (center + horizon, left + horizon),
  box(width: 24mm, height: 24mm, radius: 12mm, fill: navy,
      align(center + horizon, text(fill: white, weight: "bold", size: 44pt)[#n])),
  text(size: 46pt, weight: "bold", fill: navy, title)))
// one or two lines on how to read the figure
#let how(body) = block(below: 6mm, text(size: 26pt, fill: grey, body))
#let fig(path) = block(below: 6mm, image(path, width: 100%))

// ---- Title band ("billboard zone") ------------------------------------------------------------
#block(width: 100%, below: 14mm, inset: (bottom: 10mm), stroke: (bottom: 4pt + navy))[
  #set par(leading: 0.3em)
  #text(size: 104pt, weight: "bold", fill: navy)[The _stp1_ transposon mutant of _S. aureus_ also loses the kinase Stk1]
  #v(6mm)
  #text(size: 56pt)[… yet both mutants keep most of their response to _P. aeruginosa_ spent medium]
  #v(8mm)
  #text(size: 32pt)[#todo[Author 1, Author 2, Author 3] · group mval · BIO253 Research Cycle in Omics, UZH / ETH Zürich (FGCZ) · November 2026]
]

// ---- Column 1: question, design, H1 ----------------------------------------------------------------
#let dots(c) = text(fill: c, size: 34pt)[● ● ●]
#let col1 = [
  #sec(1)[One operon, two mutants]
  #fig("poster/p0_operon_mval.pdf")
  _S. aureus_ JE2 (#anc) and two transposon (Tn) mutants. _stk1_ lies downstream of the #stp insertion: is Stk1 lost too?
  #v(2mm)
  #table(
    columns: (auto, 1fr, 1fr), stroke: none, inset: (x: 3mm, y: 2.5mm), align: (left + horizon, center, center),
    [], text(weight: "bold")[TSB], text(weight: "bold")[PASN],
    [#anc (wild type)], dots(C-anc), dots(C-anc),
    [#stk (_stk1_::Tn)], dots(C-stk), dots(C-stk),
    [#stp (_stp1_::Tn)], dots(C-stp), dots(C-stp),
  )
  #how[● = one culture. PASN: \~15 h in 30 % cell-free spent medium of _P. aeruginosa_ PAO1 (secreted molecules only). Total proteome (2,232 proteins) + phosphoproteome.]
  #v(2mm)
  #grid(columns: (auto, 1fr), column-gutter: 5mm, row-gutter: 4mm,
    text(weight: "bold", fill: navy)[H1], [#stp also lacks Stk1 #h(1fr) → *2*],
    text(weight: "bold", fill: navy)[H2], [mutants respond differently to PASN #h(1fr) → *3*],
    text(weight: "bold", fill: navy)[H3], [#stp resembles #stk #h(1fr) → *4*],
  )
  #v(12mm)
  #sec(2)[The #stp mutant also loses Stk1]
  #how[Each dot = one sample. Stk1 rests on 29–32 peptides in #anc, but only 1–2 in #stp.]
  #fig("poster/p1_h1_enzymes_mval.pdf")
]

// ---- Column 2: H2, H3 ---------------------------------------------------------------------------------
#let col2 = [
  #sec(3)[Both mutants keep their response to PASN]
  #how[Each dot = one protein. Diagonal = same response as #anc. Labelled = altered response (limma treat, |log2FC| > 1, FDR < 0.05).]
  #fig("poster/p2_h2_responses_mval.pdf")
  #how[Phosphoproteome: 236 usable sites; none responds differently to PASN in a mutant.]
  #v(10mm)
  #sec(4)[#stp resembles #stk only in part]
  #how[r corrected for the shared #anc reference (noise alone would give r = 0.50). #text(fill: black)[●] altered in both mutants, same direction; #text(fill: red)[●] #stp and #stk differ directly.]
  #fig("poster/p3_h3_shared_distinct_mval.pdf")
]

// ---- Column 3: pathways + genome, take-home, fine print ------------------------------------------
#let col3 = [
  #sec(5)[Two gene blocks shift as a unit]
  #how[Many functional terms shift in both mutants, in #stp as many as in #stk. Part of the signal maps to two blocks of neighbouring genes.]
  #fig("poster/p4_pathways_genome_mval.pdf")
  #v(6mm)
  #block(width: 100%, fill: rgb("#fdf0e6"), stroke: 3pt + red, radius: 4mm, inset: 9mm, below: 10mm)[
    #text(size: 46pt, weight: "bold", fill: red)[Take-home messages]
    #v(1mm)
    #set text(size: 34pt)
    #enum(spacing: 7mm,
      [*#stp is not a phosphatase-only mutant:* Stk1 is \~43× lower too.],
      [*The PASN response is largely intact:* only 5 (#stk) and 1 (#stp) of \~1,900 proteins respond differently.],
      [*#stp shares only part of the #stk pattern;* two gene blocks need DNA-level checks.],
    )
  ]
  #set text(size: 25pt)
  #set par(leading: 0.4em, spacing: 0.75em)
  *Methods* prolfquapp-normalised log2 abundances; not detected = missing, never 0 · limma \~0 + group, contrasts A − B · hit = treat test |log2FC| > 1, BH-FDR < 0.05, fixed before results · tested if ≥ 2 of 3 values per group · interaction = (mutant TSB − PASN) − (#anc TSB − PASN) · pathways: rank test + cameraPR, both FDR < 0.05 · genome scan: means of 10 neighbouring proteins, permutation threshold · robustness: leave-one-sample-out refits.

  *Limits* n = 3 · no clean deletion or complemented strains · sparse phosphosites · PASN mixes secreted molecules with nutrient depletion.

  *Next* sequence both gene blocks · _stp1_ / _stk1_ transcripts · clean Δ#emph[stp1], Δ#emph[stk1] and complemented strains.

  *Mutants* Nebraska Transposon Mutant Library: Fey et al. 2013, _mBio_ 4:e00537-12.
  #v(3mm)
  #grid(columns: (50mm, 1fr), column-gutter: 8mm, align: horizon,
    image("poster/qr_repo_mval.svg", width: 50mm),
    [*Code and reports* (no raw data)\ github.com/dinvuc/bio253-je2-proteomics\ *Contact* #todo[e-mail]],
  )
]

#grid(columns: (1fr, 1fr, 1fr), column-gutter: 20mm, col1, col2, col3)
