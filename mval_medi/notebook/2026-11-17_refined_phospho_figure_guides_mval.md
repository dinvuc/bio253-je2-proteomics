# 2026-11-17 – Refined phosphosite hits, paper-style N-to-C plots, figure guides (group mval)

## Question for today
- Which phosphosite hits are robust, and can the N-to-C plots look like Fig. 1 of the
  Wolski/Grossmann chapter (ch. 8, ISBN 978-1-0716-5519-1)?
- Every figure of every step needs a guide entry with description, key points and conclusion.

## What we did
- `scripts/09b_n_to_c_plots_mval.R`: N-to-C plots for all 25 proteins with a site hit (`n_to_c/`)
- `scripts/09c_refined_hits_n_to_c_mval.R`: refined hits plus N-to-C plots in the chapter's
  style (`n_to_c_refined/`), restricted on request to proteins that are important (in
  `fig_paired_layers_mval.png`) **and** significant (refined hit). The paired-layer figure is
  copied into the folder.
- Installed `pdftools` and `remotes` (CRAN) to render the chapter's figure page and to read the
  `prophosqua` source (github.com/prolfqua/prophosqua, MIT). `prophosqua` itself could not be
  installed (its dependencies need Rtools), so its plotting style was rebuilt in ggplot.
- Figure guides: new guides for `step02_qc`, `step06_responses/hit_neighbourhoods`,
  `step09_site_vs_protein/n_to_c` and the 7 per-contrast `step10_enrichment` folders.
  Every existing guide got key points and a conclusion per figure. A check confirmed that
  every PNG appears in a guide.

## Results (key numbers)
- Refined site hit = leave-one-sample-out 18/18 + also a hit in STYH + no QC flag:
  **20 of 62 hits (14 sites, 11 proteins)**. Of the genotype site hits only **Stk1
  T164/T166** (lower in `stp`, following the Stk1 protein, −4.2 vs −5.4) survives. deoB T85
  fails (13–16/18 refits, not a hit in STYH).
- The robust "usage" changes (site beyond its protein) are all PASN responses: pflB, fda,
  qoxB, traP, pgi, glyA, eno, B7H15_04690, B7H15_04940.

## Problems / surprises
- Repeated the `summarise()` name-reuse bug (refined counts per category); caught by checking
  that the counts add up, then fixed.
- Guide drafting: three claims written from memory were wrong and were corrected against the
  data before committing: the cytolysis terms are driven by three haemolysins (not coagulase /
  FnBPB); the only interaction ORA term with ≥ 3 proteins is "Cytoplasm"; the reinjected-sample
  remark in the Step 2 guide had not been checked against the figure.
- First paper-style version drew the yellow rectangle over the bars (tinted colours) and
  listed absent residues in the legends; fixed after comparing with the rendered book page.

## AI use
Claude Code wrote the scripts, guides and this entry; it rendered the book page and read the
prophosqua source to match the figure. The conclusions in the guides, including those marked as
interpretation, were marked as checked by the group on 17 Nov 2026.

## Next steps
- Guides (conclusions per figure) and refinement rule (proposed by the AI): marked as checked by the group on 17 Nov 2026.
- Optional: install Rtools to use `prophosqua` directly.
