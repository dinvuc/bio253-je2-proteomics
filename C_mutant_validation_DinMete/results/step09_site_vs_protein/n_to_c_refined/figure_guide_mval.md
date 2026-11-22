# Refined phosphosite hits and N-to-C plots in the style of the chapter's Fig. 1

Made by [`C_mutant_validation_DinMete/scripts/09c_refined_hits_n_to_c_mval.R`](../../../scripts/09c_refined_hits_n_to_c_mval.R) (~1 min).
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

**Which proteins** (`plotted_proteins_mval.csv`, 8 proteins; column `why`):
- **paired-layer figure + refined hit** (4): proteins that are **important and significant**,
  i.e. shown in [`fig_paired_layers_mval.png`](fig_paired_layers_mval.png) (copied into this
  folder) **and** carrying ≥ 1 refined hit: glyA (`B7H15_11735`), traP (`B7H15_10065`),
  `B7H15_04690` and Stk1 (`B7H15_06355`).
- **Step 6 interaction hit, group request (18 Nov)** (1): kdpD (`B7H15_11550`), a formal
  protein hit of Step 6 ([`fig_response_vs_anc_mval.png`](../../step06_responses/fig_response_vs_anc_mval.png),
  `interaction_stk_vs_anc_at_TSB_vs_SN`). It has no site hit; it is drawn for its protein layer.
- **Step 8 site hit(s), NOT refined, group request (18 Nov)** (3): deoB (`B7H15_00775`, T85/T87),
  pyk (`B7H15_09250`, T537), Y1819 (`B7H15_09315`, T157). Their site hits **failed** the
  refinement (faded asterisks); they are drawn so the reason can be seen, not because they are
  robust. pgcA (in the paired-layer figure, no refined hit) is still not plotted.

The other 7 proteins with refined hits (pflB, eno,
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
If a plotted protein has a site hit in a contrast **outside** the 6-panel grid (so far only
Y1819 T157 in `interaction_stp_vs_stk_at_TSB_vs_SN`), that contrast is added as an extra panel
in a 4th row, so every Step 8 hit of a plotted protein is visible. The subtitle gives
"Step 8 site hits: n, refined: m" for the protein.

**How to read:** a bar that ends inside the yellow rectangle changes no more than its
protein (abundance-driven). A bar that reaches clearly beyond it changes beyond the protein
(the residual / "usage" changes of `../fig_adjusted_effects_mval.png`). Example: Stk1 in
`stp_vs_anc_at_TSB`: the rectangle reaches −5.4 and the T164 bar (−4.2, `**`, refined) ends
inside it, so the site loss follows the protein loss. The two pink `NotLoc` bars pointing up
come from the sporadic `HIQDSVPNVTtDVR` peptide (Step 3) and are not hits.

## Per figure: key points and conclusion

Signs: "PASN response" = TSB − SN (> 0 = lower in SN, < 0 = higher in SN); genotype panels =
first strain minus second, in TSB.

### `fig_n_to_c_refined_B7H15_11735_glyA_mval.png`: glyA (serine hydroxymethyltransferase)
- **Description:** glyA with its phosphosites in the 6 standard contrasts.
- **Key points:** Y51 in the `stk` PASN response: site +4.6 while the protein changes only
  +0.7 (both lower in SN). Refined hit (`**`). No site hit in the genotype panels.
- **Conclusion:** in `stk`, Y51 phosphorylation drops in PASN far more than the glyA
  protein does: a robust usage change of a one-carbon/amino-acid metabolism enzyme in the
  medium response. It is not a strain difference.

### `fig_n_to_c_refined_B7H15_10065_traP_mval.png`: traP (signal transduction protein TRAP)
- **Key points:** S150 in the `anc` PASN response: site −2.7 (higher in SN) while the
  protein is flat (+0.2). Refined. The same site in `stp` (−2.6 vs −0.2) is a hit but not
  refined (faded asterisk).
- **Conclusion:** PASN increases S150 phosphorylation of traP without changing the
  protein: a usage change in the parent, with the same direction in `stp`. It is not a
  strain difference.

### `fig_n_to_c_refined_B7H15_04690_Y845_mval.png`: B7H15_04690 (CsbD family protein)
- **Key points:** 10 site hits, all in the PASN responses (site lower in SN). In `stk`, T52, T17
  and T22 (+3.2 to +4.4) follow the protein (+3.8): refined, abundance-driven. In `stp`, T22
  (+3.8 vs protein +2.8) reaches beyond the protein (refined); T52/T17/S5 also exceed the
  protein in `stp` but are not refined.
- **Conclusion:** most of this protein's phospho signal simply follows its strong drop in
  PASN. Only in `stp` is there a robust extra drop at T22 beyond the protein.

### `fig_n_to_c_refined_B7H15_06355_nogene_mval.png`: Stk1 (serine/threonine kinase)
- **Key points:** `stp_vs_anc_at_TSB`: protein −5.4 (rectangle), T164 −4.2 and T166 −4.1
  (`**`, refined), both inside the rectangle. No Stk1 site is testable in `stk` (the protein
  is absent), and in `TSB_vs_SN_at_stk` and `stp_vs_stk` there is no estimate at all. The
  pink `NotLoc` bars pointing up are the sporadic `HIQ…` peptide, not hits.
- **Conclusion (H1):** the only robust strain difference in the phosphoproteome is the loss
  of Stk1 sites in `stp`, and it follows the loss of the Stk1 protein itself. This is
  consistent with the polar effect, and it gives no evidence that Stp1 loss changes Stk1
  phosphorylation.

### `fig_n_to_c_refined_B7H15_11550_kdpD_mval.png`: kdpD (sensor histidine kinase KdpD)
- **Description:** kdpD (885 aa; UniProt ortholog Q2FWH7, "Sensor histidine kinase KdpD",
  membrane, from `read_protein_metadata()`) in the same 6 contrasts. Added because it is a
  Step 6 interaction hit (`stk` vs `anc`), not because of a phosphosite.
- **Key points (protein, yellow rectangles; total-proteome treat test):**
  - `TSB_vs_SN_at_anc`: −1.03 (slightly higher in SN; FDR 1.0, no change).
  - `TSB_vs_SN_at_stk`: **+4.09 (much lower in SN; FDR 1e-12)**. The difference between these
    two panels is the Step 6 interaction (+5.13, FDR 1e-10): only `stk` loses kdpD in PASN.
  - `TSB_vs_SN_at_stp`: **no estimate** (`stp_SN` 1/3 detected), so `stp` cannot be judged.
  - Genotype in TSB: `stk` − `anc` +1.35 (higher in `stk`, FDR 0.28), `stp` − `anc` +0.08,
    `stp` − `stk` −1.27 (FDR 0.73): no TSB difference.
  - **Phosphosites:** 10 sites in the STY table, but only S830 (localised) has ≥ 2 of 3 values
    in both groups of any of the 6 contrasts: `stk_vs_anc_at_TSB` +2.16 (SE 0.94, FDR 1.0, not
    a hit; `anc_TSB` 2/3, `stk_TSB` 3/3). S830 is never detected in an SN sample, so there is no
    site information for the PASN response at all.
- **Caveats:** the `stk` drop rests on 2/3 detected replicates in both `anc_SN` and `stk_SN`;
  in Step 11 kdpD holds in 14/18 leave-one-out refits and lies in a genomic block
  (`B7H15_11470`–`11590`) that shifts as a unit in `stk`, so the call is **tentative**.
- **Conclusion:** kdpD is a protein-level, not a phosphorylation-level, result: in `stk` the
  KdpD protein drops strongly in PASN, while in `anc` it does not. The phosphoproteome cannot
  test whether KdpD phosphorylation changes too (too few values).

### `fig_n_to_c_refined_B7H15_00775_deoB_mval.png`: deoB (phosphopentomutase) — NOT refined
- **Description:** deoB (392 aa) with its sites T85 and T87, which sit on the same peptide.
  Group request; none of its 5 Step 8 site hits passed the refinement.
- **Key points:** protein flat everywhere (|diff| ≤ 0.66, FDR 1). The three hit contrasts all
  point at **one group**: the sites are low in `stk` TSB.
  - `stk_vs_anc_at_TSB`: T85 −3.81 (FDR 0.034, hit; lower in `stk`); T87 −3.54 (FDR 0.06)
  - `stp_vs_stk_at_TSB`: T85 +5.08, T87 +4.74 (FDR < 0.01; higher in `stp` = lower in `stk`)
  - `TSB_vs_SN_at_stk`: T85 −5.34, T87 −5.01 (FDR < 0.001; higher in SN, because `stk` TSB is low)
  - raw values: `stk_TSB` **2/3 detected**, 16.2 and 16.3; every other group 3/3 at about
    19–23 (e.g. T85 `anc_TSB` 19.9–20.4, `stp_TSB` 21.3–21.4).
- **Why not refined:** dropping `JE2stk_TSB_1` or `JE2stk_TSB_3` removes the hits (16/18;
  T85 `stk_vs_anc` 13/18); T85 is not a hit in STYH for the genotype contrasts, T87 is not
  supported in STYH.
- **Conclusion:** the deoB "strain differences" and the `stk` PASN response rest on two
  low values in one group, with the third replicate missing. Whether the missing value is
  also low (not detected because it is low) or just missed cannot be told. Not a robust result;
  the protein itself does not change.

### `fig_n_to_c_refined_B7H15_09250_pyk_mval.png`: pyk (pyruvate kinase) — NOT refined
- **Description:** pyk (585 aa); one Step 8 site hit, T537 in `stp_vs_anc_at_TSB`.
- **Key points:** T537 `stp` − `anc` in TSB = **+1.99** (SE 0.24, FDR 0.013; higher in `stp`)
  while the protein is +0.14: a change beyond the protein. All groups 3/3 detected and tight
  (`anc_TSB` 24.2–24.6, `stp_TSB` 26.3–26.4). Same direction but not a hit in `stp_vs_stk`
  (+1.67, FDR 0.12) and none in `stk_vs_anc` (+0.31). In SN, `stp` is also highest (26.6–27.0).
- **Why not refined:** it is a hit in STYH and has no QC flag, but fails leave-one-out
  **narrowly: 17/18**, lost only when `JE2anc_TSB_1` is dropped. With n = 3 the evidence is
  close to the treat threshold (|log2FC| > 1).
- **Conclusion:** the most nearly robust of the three: T537 tends to be more phosphorylated
  in `stp` (which lacks Stp1 and most of Stk1) than in `anc`, in TSB, without a protein change.
  Treat as a candidate, not a result.

### `fig_n_to_c_refined_B7H15_09315_Y1819_mval.png`: Y1819 (universal stress protein UspA) — NOT refined
- **Description:** Y1819 (166 aa); one Step 8 site hit, T157 in
  `interaction_stp_vs_stk_at_TSB_vs_SN`, shown in the **extra bottom panel** (4th row).
- **Key points:** interaction = (TSB − SN in `stp`) − (TSB − SN in `stk`) = **−2.28**
  (SE 0.33, FDR 0.045): the PASN response of T157 differs between the mutants. Read from the
  grid: in `stk` T157 is lower in SN (+0.99), in `stp` higher in SN (−1.28); neither alone is
  a hit (FDR 0.92, 0.52). Part of it follows the protein (protein interaction −1.21, FDR 1),
  so the site goes about 1 log2 beyond its protein. Raw values 3/3 everywhere; the pattern
  comes from a low `stp_TSB` (17.7–17.8) and a low `stk_SN` (18.0–18.2).
- **Why not refined:** a hit in STYH and no QC flag, but leave-one-out only **5/18**: dropping
  almost any of 13 samples removes it. The weakest of the three.
- **Conclusion:** a fragile interaction close to the FDR cut-off; partly explained by the
  protein. Do not interpret on its own.

### `fig_paired_layers_mval.png` (copied from Step 9)
- **Description / how to read:** see the [Step 9 guide](../figure_guide_mval.md#fig_paired_layers_mvalpng-fig-4c-per-sample-values-of-both-layers).
  Top row = site, bottom row = its protein, per sample.
- **Conclusion:** it is the selection basis for this folder. Of its 6 proteins, glyA, traP,
  B7H15_04690 and Stk1 survive the refinement; deoB and pgcA (shown there) do not.

## Folder conclusion
Of 62 phosphosite hits, 20 are robust. The robust "usage" changes are all PASN responses of
metabolic, stress or signalling proteins. The only robust strain difference is the Stk1
site loss in `stp`, explained by the Stk1 protein loss. The phosphoproteome therefore
supports H1 and adds nothing to H2/H3.
