# Step 9 outputs: what each figure and table compares

Made by [`scripts/09_site_vs_protein_mval.R`](../../../scripts/09_site_vs_protein_mval.R).
Question: a phosphosite signal goes up when the site is phosphorylated more **or**
when there is more of the protein. For each supported site hit from
[Step 8](../step08_phospho_coverage/figure_guide_mval.md), does the site change follow
the protein change in **exactly the same contrast**?

## Method in short

- **Inputs:** our Step 8 site refit (STY search, primary; STYH = sensitivity) and our
  Step 5 protein refit. Same model for both: limma `~ 0 + group` on 27 samples,
  observed values only. Matched by the **exact `protein_key`** (same search FASTA),
  with a locus-tag check. All 236 supported sites found their protein.
  (`B7H15_13880` has a second, 6-frame pseudogene entry in the total table; the exact
  key avoids picking it.)
- **Adjusted effect** = `diff_site − diff_protein` (log2, A − B of the same contrast).
  It is a **relative, protein-adjusted phosphorylation signal**, not occupancy or
  stoichiometry. Computed only when the protein was tested (≥ 2 of 3 in each group),
  so nothing imputed is ever subtracted.
- **Uncertainty:** `se_adj = sqrt(se_site² + se_protein²)`. This **assumes the two
  estimates are independent**. The covariance term of
  `Var(site − protein) = Var(site) + Var(protein) − 2·Cov` cannot be estimated: the
  phospho and total runs are separate measurements, and Step 2 did not establish that
  "X_1" is the same culture aliquot in both layers. With a positive covariance `se_adj` is
  too large (conservative). **No p-value, no FDR** for the adjusted effect: the site
  FDR is not reused for this new quantity.
- **Categories (fixed before the matched values were seen; AI proposal, group to confirm),
  for site hits only:**

  | category | rule | meaning |
  |---|---|---|
  | consistent with protein amount | protein tested, \|adjusted\| < 1 | the site change is within 2-fold of the protein change |
  | residual adjusted change | protein tested, \|adjusted\| ≥ 1 **and** ≥ 2·se_adj | more site change than the protein change explains |
  | unresolved | protein not tested, or \|adjusted\| ≥ 1 but < 2·se_adj | cannot tell |

  These categories **do not identify Stk1/Stp1 substrates**. A residual change can also
  come from a different peptide form, a co-modification, or layer-specific normalisation.
- **QC flag (added after the first run, does not change any category):**
  `fragile_spread_gt3` = some group of the contrast has a within-group range > 3 log2.
  Reason: `B7H15_04940` T269 came out at +10.8 log2 in `TSB_vs_SN_at_stp`, resting on
  two `stp_SN` values (10.9 and 14.1) near the detection floor.

## Tables

### `site_hit_categories_mval.csv` + `category_counts_mval.csv`
**62 site × contrast hits (37 sites, 25 proteins): 17 consistent with protein amount,
40 residual adjusted change, 5 unresolved (adjusted effect uncertain).** Every hit's
protein was tested, so there are no "protein not tested" cases.

| contrast (STY) | consistent | residual | unresolved |
|---|---|---|---|
| stk_vs_anc_at_TSB | 0 | 1 | 0 |
| stp_vs_anc_at_TSB | 0 | 1 | 2 |
| stp_vs_stk_at_TSB | 0 | 2 | 0 |
| TSB_vs_SN_at_anc | 4 | 6 | 0 |
| TSB_vs_SN_at_stk | 10 | 16 | 0 |
| TSB_vs_SN_at_stp | 3 | 13 | 3 |
| interaction_stp_vs_stk_at_TSB_vs_SN | 0 | 1 | 0 |

The genotype (baseline) and interaction hits, the ones that matter for our question:

| contrast | site | site diff | protein diff | adjusted ± se_adj | category | in words |
|---|---|---|---|---|---|---|
| stk_vs_anc_at_TSB | deoB T85 (`B7H15_00775`) | −3.81 | +0.05 | −3.86 ± 0.68 | residual | site lower in `stk` than `anc` in TSB; protein unchanged. `stk_TSB` 2/3 detected |
| stp_vs_stk_at_TSB | deoB T85 / T87 | +5.08 / +4.74 | +0.02 | +5.06 / +4.71 | residual | site higher in `stp` than `stk`; T85 and T87 sit on the same peptide → one piece of evidence |
| stp_vs_anc_at_TSB | pyk T537 (`B7H15_09250`) | +1.99 | +0.14 | +1.84 ± 0.25 | residual | site higher in `stp` than `anc`; protein unchanged |
| stp_vs_anc_at_TSB | Stk1 T164 / T166 (`B7H15_06355`) | −4.22 / −4.14 | −5.42 | +1.21 ± 0.61 / +1.28 ± 0.71 | unresolved | site and protein both much lower in `stp`; the site drop is mostly the protein drop (H1) |
| interaction_stp_vs_stk_at_TSB_vs_SN | T157 of `B7H15_09315` | −2.28 | −1.21 | −1.06 ± 0.35 | residual (borderline, just past 1) | the `stp`-vs-`stk` difference at this site is larger in SN than in TSB |

No site hit in `interaction_stk_vs_anc_*` or `interaction_stp_vs_anc_*`.
QC-flagged (range > 3 log2): pgcA S193 (`TSB_vs_SN_at_stk`), `B7H15_04940` T269 and
Y845 T52 (`TSB_vs_SN_at_stp`).

### `site_vs_protein_tracking_mval.csv`: all supported sites, not only hits
- `r_site_vs_protein`: how well site effects follow protein effects. Higher in the
  medium contrasts (0.68–0.72) than in genotype contrasts (0.32–0.63) and interactions (0.32–0.44).
- **`median_signed_adjusted` (layer offset):** −0.65 (`TSB_vs_SN_at_anc`), −0.43 (`stk`),
  −0.17 (`stp`). Across all sites, site signals lean towards "higher in SN" relative to
  their proteins. Likely in part technical: SN phospho samples detect ~2× more sites
  (Step 8), and each layer was normalised on its own.

### `offset_sensitivity_mval.csv`
The categories are recomputed after subtracting each contrast's median offset (a labelled
sensitivity, not a correction). **Residual hits stay at 40, but 5 leave** (pflB S742 in
`anc`; fda S22, rplY S67, qoxB T552 in `stk`; fda S22 in `stp`) **and 5 enter**
(`B7H15_04940` S264 and Y845 T22 in `anc`; Y2093 S51 and `B7H15_04940` S264 in `stk`;
trxA S71 in `stp`). So 35 of 40 residual calls do not depend on the offset. Near
|adjusted| = 1 the category is fragile.

### Search sensitivity (`site_vs_protein_checks_mval.txt`)
Same site, same contrast in STYH: same category for 46 of 62 hits, different for 14
(mostly near the |adjusted| = 1 boundary), site not supported in STYH for 2.

### `data/derived/site_protein_matched_mval.tsv.gz` (local only, full table), `unmatched_sites_mval.csv`
All supported site × primary contrast rows with both estimates, SEs, site FDR (ours),
protein FDR (Step 5), adjusted effect and category. `unmatched_sites` lists sites whose
protein was not tested in that contrast (not a hit-level problem: none of the 62 hits).

## Figures

### `fig_site_vs_protein_mval.png`: Fig 4a, site effect vs protein effect
- **Compares:** per primary contrast, every supported site whose protein was tested:
  protein effect (x, total proteome) vs site effect (y), both log2 A − B.
- **How to read:** solid diagonal = site change fully explained by protein amount;
  dotted = ±1 log2 (the "consistent" band). Grey = not a site hit. Blue/orange/pink =
  category. Points far above or below the band = residual changes (labelled).
- **What it shows:** in the three TSB-vs-SN panels many hits lie in or near the band
  (protein-driven), others well below it (site higher in SN than its protein). Genotype
  and interaction panels have almost no hits.

### `fig_adjusted_effects_mval.png`: Fig 4b, adjusted effects of all site hits
- **Compares:** for each of the 62 hits, the adjusted effect (filled, ± 2·se_adj) and the
  unadjusted site effect (open circle). Dotted = ±1.
- **How to read:** a filled point moved towards 0 compared with the open circle = part of
  the site change was protein amount. **x** = QC flag (fragile spread).
- **Caveat:** error bars assume independence (see above); they are not confidence
  intervals from a test.

### `fig_paired_layers_mval.png`: Fig 4c, per-sample values of both layers
- **Compares:** top row = normalised phosphosite values, bottom row = its total protein,
  18 primary samples, numbers = samples with a value (of 3). Every panel has its own y-scale.
- **Shown:** the 5 largest non-flagged residual sites (deoB T87, glyA Y51, deoB T85,
  pgcA T191, traP S150), the strongest "consistent" example (Y845 T52 in `stk`), and Stk1 T164.
- **How to read:** compare the *pattern* across groups between the rows. Same pattern =
  protein-driven; site moves while protein stays flat = residual (e.g. deoB T85: the site
  drops in `stk_TSB`, the protein does not).
- **Stk1 T164:** protein 0/3 in `stk` and 3/3 at a very low level in `stp_TSB` (0/3 in
  `stp_SN`); the site is absent in `stk` and low in `stp` (2/3 in `stp_SN` even though
  the protein is 0/3 there: phospho-enrichment can see a peptide the total proteome misses).

## Limits
- n = 3 per group; one search choice that was nearly a tie; categories rest on a 2-fold
  band and an approximate SE.
- Only 236 sites (133 proteins) are measurable at all; absence of a phospho change at the
  other ~3,500 sites is **not** evidence of no change.
- No direct Stk1/Stp1 substrates can be claimed: the data cannot separate direct from
  indirect phosphorylation changes.
