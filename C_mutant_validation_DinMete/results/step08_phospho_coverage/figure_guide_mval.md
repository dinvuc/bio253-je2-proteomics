# Step 8 outputs: what each figure and table compares

Made by [`C_mutant_validation_DinMete/08_phospho_coverage_mval.R`](../../08_phospho_coverage_mval.R).
Question: is the phosphosite layer measured often enough, and localised well enough,
to support a follow-up (Step 9)? Data: FragPipe site tables (`read_sites()`), the
normalised site abundances the supplied DEA used (`read_abundances("phospho_*")`,
log2, **not** re-normalised) and the supplied DEAs (comparison only).
Primary analysis = 18 samples (`anc`/`stk`/`stp` × TSB/SN) and the 9 primary
contrasts of the [Step 5 contrast sheet](../step05_formal_test/contrast_sheet_mval.csv).

## Rules (fixed before any site was selected; AI proposal, marked as checked by the group on 13 Nov 2026)

| rule | what | why |
|---|---|---|
| JE2 only | `keep_je2()`: locus tag, no `rev_` decoy | same filter as the total proteome |
| eligible | ≥ 2 of 3 observed values in **every** group of the contrast | Step 5 rule; NA stays NA |
| localised | **median** per-sample `loc_prob` over the primary samples where the site was detected ≥ 0.75 | `best_loc_prob` is the best of 27 samples, so one lucky spectrum passes. 0.75 = conventional "class I" cutoff (Olsen et al. 2006, *Cell* 127:635–48, PMID 17081983, [doi:10.1016/j.cell.2006.09.026](https://doi.org/10.1016/j.cell.2006.09.026); reference checked in PubMed, the 0.75 definition is in the full text, not the abstract: marked as checked by the group on 13 Nov 2026). It is **not** a measured false-localisation rate |
| histidine | pH sites counted and listed, never in the shortlist | no validated pH localisation cutoff; pH is acid-labile |
| supported | eligible + localised + residue S/T/Y | |
| primary search | more supported S/T/Y site × contrast pairs; tie → STYH | |
| site test | our refit, limma `~ 0 + group` (27 samples), `treat(lfc = 1)`, **BH within contrast over the supported sites only**, hit = FDR < 0.05 | **new test family, our FDR, not the supplied one**; observed values only, so it can be matched with the protein estimates in Step 9 |

**Result of rule 6: STY is the primary search, by 1,159 vs 1,150 supported pairs, i.e.
practically a tie.** STYH is the sensitivity check; both give nearly the same answer (below).

## Audit findings (tables)

### `mapping_audit_mval.csv`: features, decoys, joins per search

| | STYH | STY |
|---|---|---|
| site features in the abundance table | 4,278 | 3,863 |
| reversed decoys (`rev_`) | 43 | 0 |
| no locus tag (contaminants, 6-frame entries) | 131 | 89 |
| **JE2 sites used** | **4,104** | **3,774** |
| features per `site_id` | 1 (one peptide feature per site) | 1 |
| in abundance table, not in site table | 33 | 0 |
| in site table, not in abundance table | 0 | 32 |
| sites with ≥ 2 peptide forms (e.g. singly + doubly phosphorylated) | 473 | 414 |
| **sites that share a peptide sequence with another site** | **2,780** | **2,445** |

- **The STY abundance table stores observed values only** (16,209 rows, no NA rows).
  The script fills in every site × sample with NA before counting; otherwise "2 of 3"
  would be counted against the wrong denominator.
- Detection agrees 100 % between the site table (`intensity`) and the abundance table
  (`normalized_abundance`) for every shared site × sample cell, in both searches.
- **Sites that share a peptide are not independent evidence.** About two thirds of the
  JE2 sites (STY 2,445 of 3,774) share an unmodified peptide sequence with another
  site. Example: deoB T85 and T87 are measured partly from the same peptides (incl. the
  doubly phosphorylated form) and have almost identical values (`stk_TSB` 16.2/16.3 for
  both). Count proteins/peptides, not sites, when you summarise.
  *(Corrected 16 Nov: the first version reported "534 / 492 peptides carrying ≥ 2 sites";
  that number compared whole peptide-list strings and was not a count of shared peptides.)*

### `supplied_model_audit_mval.csv`: what the supplied phospho DEAs say about themselves
- **STYH** (`lm_impute`): in the 9 primary contrasts 35,523 of 36,936 rows (96 %) are
  labelled `lod_imputed`. As in the total proteome (Step 5), this is a per-site label.
- **STY**: two model names, `Imputed_Mean_moderated` (30,168 rows) and
  `Linear_Model_moderated` (3,798 rows), but **`estimate_type` is "observed" in every
  row**. In STY this label therefore cannot tell imputed from observed estimates. Do not
  use the STY `estimate_type` column to filter. Our refit avoids the problem.

### `detection_per_sample_mval.csv` + `fig_detection_depth_mval.png`
- **Compares:** number of JE2 sites with a value in each of the 27 samples, both searches.
  Triangles = the two reinjected `stk_pulseSN` phospho samples.
- **Result (median per group, STY):** `anc_TSB` 536, `anc_SN` 1,033; `stk_TSB` 424,
  `stk_SN` 874; `stp_TSB` 623, `stp_SN` 673. STYH is similar.
- **Caveat:** in `anc` and `stk`, SN samples see about **twice as many sites** as TSB
  (1,033 vs 536; 874 vs 424). In `stp` the difference is small (673 vs 623). For `anc`
  and `stk`, "only detected in SN" is therefore expected for many sites even without
  biology (missing not at random). It also matters for Step 9 (layer offset, which is
  largest in `anc` and smallest in `stp`, matching this pattern).
  *(Corrected 16 Nov: the first version said "SN ~2× TSB" for all strains.)*
- **Conclusion:** phospho detection depth differs strongly between samples and media,
  much more than in the total proteome (Step 2). Presence/absence of a site between TSB
  and SN is therefore not trustworthy on its own, and site effects must be compared with
  care across media.

### `fig_localisation_mval.png` + `residue_counts_mval.csv`
- **Compares:** for each site detected in ≥ 1 primary sample, the median per-sample
  `loc_prob` (our rule, orange) and `best_loc_prob` (blue), per residue and search. Dashed = 0.75.
- **How to read:** many sites have a **median per-sample loc_prob of 0**: in most samples
  FragPipe placed the phosphate on another residue of the peptide. In STY these are 401
  sites (S 158 / T 161 / Y 82); their best score is still high (median 0.78) and 221 of
  them would pass a `best_loc_prob` ≥ 0.75 filter. A peak at 0.5 = two equally likely positions.
- **Result (STY, sites detected in ≥ 1 primary sample):** localised by our rule
  S 690 / T 412 / Y 170 (1,272), by `best_loc_prob` alone S 1,083 / T 740 / Y 286
  (2,109). **The best score alone would call about 1.7× as many sites "localised".**
  *(Corrected 16 Nov: the first version counted `best_loc_prob` over all sites, not only
  those detected in the primary samples, and said "roughly double".)*
- **Histidine (STYH):** 408 pH sites; 82 pass our localisation rule (170 by best score),
  30 are eligible in some primary contrast, **8 are eligible and localised**. Listed in
  `histidine_sites_mval.csv`, exploratory only.
- **Conclusion:** many sites are placed on a residue only in a minority of samples.
  Using the per-sample median instead of the best score removes about 40 % of the
  "localised" sites and keeps only sites whose position is reproducible. Histidine
  phosphorylation is too sparse to analyse (8 usable sites).

### `coverage_per_contrast_mval.csv` + `fig_coverage_funnel_mval.png`
- **Compares:** per primary contrast, how many JE2 sites exist (grey), are eligible
  (2 of 3 in every group, light blue) and are supported (eligible + localised, S/T/Y, dark blue).
- **Result:** of 3,774 STY sites only **119–296 are eligible and 91–184 supported per
  contrast** (2–5 %). Interactions need 4 groups and have the fewest (91–125). STYH: 93–186.
- **Meaning:** the phospho layer can only speak about a small, abundant subset of sites.
- **Conclusion:** coverage is the limiting factor of the phosphoproteome: 2–5 % of the
  sites per contrast are usable, and the interaction contrasts (H2) are the worst covered.
  Absence of a phospho change is therefore weak evidence.

### `phospho_hit_counts_mval.csv`: our site hits vs the supplied screen
Site hit = treat(lfc = 1) BH-FDR < 0.05 over the supported sites. "up" = higher in the
FIRST group of the contrast name (`TSB_vs_SN` up = higher in TSB, down = higher in SN).

| contrast (STY) | supported | up | down | supplied screen, all sites | supplied screen, supported sites |
|---|---|---|---|---|---|
| stk_vs_anc_at_TSB | 118 | 0 | 1 | 40 | 13 |
| stp_vs_anc_at_TSB | 155 | 1 | 2 | 39 | 10 |
| stp_vs_stk_at_TSB | 114 | 2 | 0 | 64 | 18 |
| TSB_vs_SN_at_anc | 184 | 2 | 8 | 170 | 72 |
| TSB_vs_SN_at_stk | 124 | 11 | 15 | 151 | 54 |
| TSB_vs_SN_at_stp | 148 | 11 | 8 | 147 | 50 |
| interaction_stk_vs_anc_at_TSB_vs_SN | 100 | 0 | 0 | 16 | 3 |
| interaction_stp_vs_anc_at_TSB_vs_SN | 125 | 0 | 0 | 33 | 1 |
| interaction_stp_vs_stk_at_TSB_vs_SN | 91 | 0 | 1 | 51 | 6 |

- **62 site × contrast hits = 37 sites on 25 proteins**, almost all in the within-genotype
  TSB-vs-SN contrasts. **No site hit in the two H2 interactions** (`stk`/`stp` vs `anc`).
- The supplied screen (FDR < 0.1, then |diff| > 1, mostly on imputed values) gives 711
  hits over all sites, 227 on supported sites. It is shown for comparison only.
- Guardrail: our site estimates agree with the supplied ones on supported sites
  (r = 1.000 STY, 0.967 STYH): same sign convention.

### `supported_sites_mval.csv`: the shortlist
236 supported S/T/Y sites (S 164, T 58, Y 14) on 133 proteins, with gene symbol (from
`read_protein_metadata()`, `ortholog_quality == "good"`), median and best loc_prob,
detections per primary group, and the contrasts in which the site is a hit.

## Sensitivity checks (`coverage_checks_mval.txt`)
- **STY vs STYH** (`fig_sty_vs_styh_mval.png`): 1,068 site × contrast pairs supported in
  both searches; estimates r = 0.96, points on the identity line. Hits: 51 in both, 9 STY
  only, 2 STYH only. The choice of search hardly matters.
  - *Figure conclusion:* the points lie on the identity line in every panel, so the
    primary-search choice (a near tie) does not drive any result.
- **Reinjected samples** (`stk_pulseSN_1/2`, not in the primary contrasts but in the
  pooled variance): refit without them → hits 62 → 78 (none lost, 16 new). The main fit
  is the **more conservative** one.

## Figures of the hits (look before trusting)
- `fig_supported_hits_raw_mval.png`: the 9 strongest site hits, normalised values the
  model used, 18 primary samples, numbers = samples with a value (of 3).
- `fig_supported_hits_sites_mval.png`: the same sites from the FragPipe site table
  (`plot_sites()`, all 27 samples, not normalised). Same patterns as the normalised values.
- Things visible there: the Stk1 site (`B7H15_06355` T164) is absent in `stk` and
  ~4 log2 lower in `stp` (compare Step 3: the Stk1 protein is lost in `stk` and nearly
  lost in `stp`); deoB T85 is low in `stk` TSB only (2/3, both near 16).
- **Key points:** most of the strongest site hits are TSB-vs-SN changes in metabolic
  enzymes (e.g. fda, pgi, pflB, glyA). The normalised values and the raw FragPipe intensities
  show the same patterns.
- **Conclusion:** the site hits are visible in the per-sample values and do not come from
  the normalisation. The Stk1 site confirms H1 at the phospho level. deoB T85 rests on two
  low `stk_TSB` values and later fails the refinement checks (Step 9c).

## Decision for Step 9
Coverage is **sufficient for a small, descriptive protein-matching step** on the 236
supported sites, and too thin for a quantitative phosphorylation answer to H2/H3:
with no interaction site hits, the phospho layer adds no evidence about altered PASN
*responses*. Baseline site hits (6) and exposure hits (55) go to Step 9.
