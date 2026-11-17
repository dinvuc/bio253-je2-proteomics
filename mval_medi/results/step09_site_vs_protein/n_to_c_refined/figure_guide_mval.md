# Refined phosphosite hits and N-to-C plots in the style of the chapter's Fig. 1

Made by [`scripts/09c_refined_hits_n_to_c_mval.R`](../../../../scripts/09c_refined_hits_n_to_c_mval.R) (~1 min).
Reference: Wolski, Dittmann, Panse, Kunz & Grossmann, "Integrated Analysis of
Posttranslational Modifications and Total Proteome: Methods for Distinguishing Abundance
from Usage Changes", Chapter 8 (ISBN 978-1-0716-5519-1), Fig. 1, made there with
`prophosqua::n_to_c_expression_multicontrast`. The first version of these plots is in
[`../n_to_c/`](../n_to_c/).

## Refined hits (`refined_hits_mval.csv`, `refined_summary_mval.csv`)

A Step 8 site hit (treat |log2FC| > 1, BH-FDR < 0.05, STY search) counts as **refined**
only when it passes all three robustness checks. These are sensitivity labels: the
Step 8 hits are not changed.

| check | rule |
|---|---|
| leave-one-sample-out | still a hit in all 18 refits that each drop one primary sample (same model and BH family) |
| search | also a hit in the STYH search (same site, same contrast) |
| QC | no group of the contrast with a within-group range > 3 log2 (Step 9 flag) |

**Result: 20 of 62 site hits are refined (14 sites on 11 proteins).** Most failures come
from leave-one-out: dropping a single sample removes the hit. Per contrast:

| contrast | site hits | refined | residual / consistent / unresolved among refined |
|---|---|---|---|
| TSB_vs_SN_at_anc | 10 | 6 | 5 / 1 / 0 |
| TSB_vs_SN_at_stk | 26 | 9 | 4 / 5 / 0 |
| TSB_vs_SN_at_stp | 19 | 3 | 3 / 0 / 0 |
| stp_vs_anc_at_TSB | 3 | 2 | 0 / 0 / 2 (Stk1 T164, T166) |
| stk_vs_anc_at_TSB, stp_vs_stk_at_TSB, interaction_stp_vs_stk | 1, 2, 1 | 0 | — |

- **Of the genotype site hits, only Stk1 T164/T166 (lower in `stp`, following the Stk1
  protein) survives.** deoB T85/T87 (`stk` vs `anc` / `stp` vs `stk`) fail: 13–16 of 18
  refits and not hits in STYH. pyk T537 fails narrowly (17/18).
- Refined **residual** changes (site beyond its protein) are all in the PASN response:
  pflB S742, fda S238, qoxB T552, traP S150, pgi T143 (`anc`); glyA Y51, fda S238, eno S373,
  pgi T143 (`stk`); B7H15_04690 T22, fda S238, B7H15_04940 S264 (`stp`).
- The table also says, for each hit, why it failed (`failed`) and whether a residual call
  survives the layer-offset sensitivity of Step 9 (`residual_survives_offset`).

## Figures (`fig_n_to_c_refined_<locus>_<gene>_mval.png`, `n_to_c_refined_all_mval.pdf`)

**Which proteins** (`plotted_proteins_mval.csv`, 4 proteins): only proteins that are **important
and significant**, i.e. shown in [`fig_paired_layers_mval.png`](fig_paired_layers_mval.png) (copied
into this folder) **and** carrying ≥ 1 refined hit. That gives glyA (`B7H15_11735`), traP
(`B7H15_10065`), `B7H15_04690` and Stk1 (`B7H15_06355`). deoB and pgcA are in the paired-layer
figure but have no refined hit, so they are not plotted (their hits are still listed with the
reason in `refined_hits_mval.csv`). The other 7 proteins with refined hits (pflB, eno,
`B7H15_04940`, pgi, qoxB, Y2093/`B7H15_10475`, fda) are not in the paired-layer figure; their
plots are in the first version in [`../n_to_c/`](../n_to_c/).

**Note on the copied `fig_paired_layers_mval.png`:** it is the Step 9 figure unchanged, so it
still shows deoB T85/T87 and pgcA T191, which did **not** pass the refinement.

**Layout:** the same 3 × 2 grid for every protein, so proteins can be compared. Left
column: genotype difference in TSB (`stk_vs_anc`, `stp_vs_anc`, `stp_vs_stk`). Right
column: PASN response within each strain (`TSB_vs_SN` at `anc`, `stk`, `stp`; > 0 = lower in SN).

**Style, taken from the printed figure and from the prophosqua source code**
(`github.com/prolfqua/prophosqua`, `R/plotNtoC.R`, MIT licence; the package could not be
installed because its dependencies must be compiled with Rtools):
- bold overall title `Protein: <id> (length: L)`; each panel titled
  `Prot : <id>; length: L; # sites:n; # not localized sites:m`; y `diff : <contrast>`,
  x `posInProtein`; `theme_minimal`; a legend next to every panel (`modAA`,
  `imputation_status`, `Rectangle`), each listing only what that panel contains
- **yellow rectangle** (alpha 0.3) from 0 to the protein log2FC ("diff of protein"), with
  `N`/`C` at its height; a black zero line from 0 to the protein length
- **bars** from 0 to each site log2FC, coloured by `modAA` with prophosqua's colours (S blue,
  T green, Y orange). **Not-localised sites are shown as modAA `NotLoc` (pink)**, which is how
  prophosqua marks them.
- `imputation_status`: our model never imputes, so every bar is "observed" (solid), the
  same as every site in the chapter's example
- red asterisks at the bar tip: `**` FDR < 0.01, `*` FDR < 0.05 (chapter caption)
- protein without an estimate in that contrast: red dashed box "No estimate for diff of
  protein" (as in prophosqua). "no site tested" = no site of this protein had ≥ 2 of 3
  values in every group of the contrast.
- each panel has its own y-scale (as in the chapter). Compare the bar with its own
  rectangle, not bar heights between panels.

**Our additions:** bright red asterisks = refined hit; faded asterisks = hit that failed a
robustness check. Sites too sparse to test are not drawn (prophosqua draws only estimated sites).

**How to read:** a bar that ends inside the yellow rectangle changes no more than its
protein (abundance-driven). A bar that reaches clearly beyond it changes beyond the protein
(the residual / "usage" changes of `../fig_adjusted_effects_mval.png`). Example: Stk1 in
`stp_vs_anc_at_TSB`: the rectangle reaches −5.4 and the T164 bar (−4.2, `**`, refined) ends
inside it, so the site loss follows the protein loss. The two pink `NotLoc` bars pointing up
come from the sporadic `HIQDSVPNVTtDVR` peptide (Step 3) and are not hits.
