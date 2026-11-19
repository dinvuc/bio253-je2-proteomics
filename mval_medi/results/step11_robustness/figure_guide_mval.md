# Step 11 outputs: robustness of the headline conclusions

Made by [`scripts/11_robustness_mval.R`](../../../scripts/11_robustness_mval.R) (~40 s).
Question: do the headline conclusions of Steps 3–10 survive reasonable checks, or do
they rest on one replicate, a marginal threshold, imputed values or one genomic region?
Every analysis here is a **sensitivity analysis**. The main results and the hit rule
(treat lfc = 1, BH-FDR < 0.05) stay as they are.

## Headline table

Status: **stable** = survives every check; **tentative** = supported but depends on a
check, an assumption or a few samples; **unresolved** = the data cannot decide.

| # | conclusion | evidence | checks | status |
|---|---|---|---|---|
| H1a | Stk1 is lost in `stk` | 1/9 `stk` samples detected, 1 peptide (Step 3); ErmB marker in 18/18 mutants, 0/9 `anc` (Step 4b) | on/off pattern, not a fold change (no formal test possible) | **stable** |
| H1b | Stk1 is strongly reduced in `stp` | `stp_vs_anc_at_TSB` −5.42 (higher in `anc`), FDR 0.0013 | hit in **18/18** LOSO refits, estimate −5.62 to −5.16 | **stable** (as a measurement) |
| H1c | …because of a polar insertion | Stk1 is the last gene of the transcription unit (Step 4), Stk1 site follows the protein (Step 9) | protein data cannot separate polarity from stability/regulation/secondary mutation | **tentative** (consistent, not proven) |
| H1d | Stp1 is reduced in `stp` | 0–2 peptides per sample vs 9–15 in `anc` (TSB, SN); Step 3 call "supported reduction" | the abundance values themselves are erratic (they rest on 1–2 peptides), so the **size** of the loss cannot be estimated | **stable** that it is reduced; size **unresolved** |
| H2a | Few proteins change their PASN response; most of the response is shared with `anc` | 5 `stk` / 1 `stp` interaction hits; response r = 0.85 / 0.86 vs `anc` (Step 6) | shape holds under both alternative rules: `stk` 5 / 15 / 5, `stp` 1 / 5 / 2 (main / 1.5-fold / FDR 0.1) | **stable** |
| H2b | `stk`-specific response hits | argF, B7H15_14515: 18/18; arcA 16/18; **kdpD, rbsD 14/18 and in a flagged genomic block of `stk`** | LOSO, genome scan | argF, B7H15_14515 **stable**; arcA **tentative**; kdpD, rbsD **tentative** (region effect possible) |
| H2c | `stp` response hit ilvA | 2.36 (SN response lower in `stp`), FDR 0.039 | **9/18** LOSO; also lost by 9 different single drops | **tentative / fragile** |
| H3a | `stp` shares a moderate part of the `stk` pattern | noise-corrected r 0.37 (0.29–0.45) baseline, 0.38 (0.24–0.52) response (Step 7); shared pathway shifts in peptidoglycan synthesis and cytolysis (Step 10) | depends on the noise model; split-reference lower bound 0.03–0.33 | **tentative** |
| H3b | `stp` and `stk` differ directly (20 baseline + 4 response proteins) | Step 7 direct contrasts | LOSO stable 11/20 and 2/4; **6 of the 20 baseline hits lie in a flagged block** of one of the mutants: the `stp` prophage (B7H15_10935, B7H15_11085), another `stp` block (azoR, B7H15_01055–01135), or `stk` blocks (mazE in 11470–11590, sceD in 11640–11700, lytM in 01445–01545) | **tentative** |
| P1 | Phosphosites add no evidence about altered PASN responses | 0 site hits in the two H2 interactions (Step 8) | same in STY and STYH; same without reinjected samples | **stable** (as absence of evidence) |
| P2 | Stk1 site drop in `stp` follows the protein loss | adj. +1.2 ± 0.6 (Step 9) | descriptive, SE assumes independent layers | **tentative** |
| E1 | "Stk1 broadly involved": many pathway shifts in `stk` | 35 / 143 / 44 robust terms (Step 10) | `stp` has as many in TSB and SN; **part of the `stk` signal is one genomic block** (below) | **unresolved** |

## 1. Leave-one-sample-out (`loso_hit_stability_mval.csv`, `fig_loso_stability_mval.png`)
- **Compares:** every main formal hit in the six H1–H3 contrasts (64 hits) against
  18 refits, each leaving out one of the 18 primary samples (same model, same rule).
- **How to read:** bar = refits in which the protein is still a hit (dashed = 18). Green
  = hit in all 18 ("stable"). Orange = lost at least once. `lost_when_dropping` names
  the sample(s); `diff_range_loso` = range of the estimate over the refits.
- **Result:** **34 of 64** hits are stable. Per contrast (stable/all): `stk_vs_anc_at_TSB`
  9/20, `stp_vs_anc_at_TSB` 10/14, `stp_vs_stk_at_TSB` 11/20, `stk` interaction 2/5,
  `stp` interaction 0/1, direct interaction 2/4. The estimates move little; most losses are
  FDR crossing 0.05 after a group drops to n = 2. The samples whose removal costs the most
  hits are `stk_TSB_2` and `stk_TSB_3` (16 each).
- **Caveat:** with n = 3, dropping one sample removes a third of a group's data, so some
  losses are expected even for real effects. "Fragile" means "do not lean on it alone",
  not "false".
- **Conclusion:** about half of the formal hits are fully robust. The ones the conclusions
  lean on most hold up: the Stk1 drop in `stp`; the shared surface-protein decreases
  (coagulase and `B7H15_02275` stable in both mutants, the fibronectin-binding protein
  `B7H15_13900` stable in `stk` but 14/18 in `stp`); and the `stk` response hits argF and
  B7H15_14515. ilvA, the single `stp` response hit,
  is fragile.

## 2. Threshold sensitivity (`threshold_sensitivity_mval.csv`)

| contrast | main (lfc 1, FDR 0.05) | lfc log2(1.5), FDR 0.05 | lfc 1, FDR 0.10 |
|---|---|---|---|
| stk_vs_anc_at_TSB | 20 | 67 | 20 |
| stp_vs_anc_at_TSB | 14 | 50 | 18 |
| stp_vs_stk_at_TSB | 20 | 92 | 27 |
| TSB_vs_SN_at_anc / stk / stp | 204 / 230 / 213 | 391 / 440 / 407 | 216 / 245 / 226 |
| interaction stk vs anc | 5 | 15 | 5 |
| interaction stp vs anc | 1 | 5 | 2 |
| interaction stp vs stk | 4 | 19 | 5 |

A lower fold-change threshold roughly triples the hit lists, but the **shape** of every
answer stays the same: few response changes, more in `stk` than in `stp`, and a PASN
response of similar size in all three strains.

## 3. Genome-position scan (`fig_genome_scan_mval.png`, `genome_*_mval.csv`)
- **Method:** the tested proteins in genome order; mean effect of 10 neighbouring
  proteins; null = 1000 random permutations of the effects over the positions; a window
  is flagged when it is more extreme than the 95th percentile of the genome-wide
  maximum. *(Fixed 17 Nov: a first version used circular shifts, which keep every
  neighbourhood intact; its null could never flag anything.)*
- **How to read:** a flagged window means neighbouring proteins change together more
  than random proteins do. Operons and regulons do this naturally. A flag is a reason
  to look, **not** proof of an artefact.
- **Result:** 11 flagged blocks over 6 contrasts. The two that matter
  (`genome_block_members_mval.csv`):
  1. **`stk`, B7H15_11470–11590 (~2.19 Mb), mean +1.51 log2, the top window in the genome
     (family-wise p < 0.001).** Almost every protein from the SigB operon (`rsbV`, `rsbW`,
     SigB, RsbU-like phosphatase) through `mazEF`, `alr1`, `acpS`, `kdpCBA`/`kdpDE`,
     `cshA`, `murF`, `ddl` is +0.7 to +2.6 higher in `stk` in TSB, similarly in SN
     (window +1.0, p = 0.10; a neighbouring block B7H15_11640–11700 with `sceD` is
     flagged in SN). Not seen in `stp` (best window −0.15).
     **Interpretation (group to check):** functionally unrelated neighbours rising together
     by ~2× fit a **regional copy-number change** (e.g. a duplication) in the `stk`
     strain. Regional regulation (SigB regulon, `kdp`) is the alternative. Protein data
     cannot decide; DNA data (genome sequencing or qPCR of this region in NE217 vs JE2)
     can.
  2. **`stp`, B7H15_10910–11150 (~2.1 Mb), mean −2.08 log2 (p < 0.001): a prophage**
     (tail, capsid, portal, terminase, recombinase). Lower in `stp` in TSB; the same block
     drives the most extreme windows of **both** interaction contrasts (−1.18 `stk`,
     −1.33 `stp`). This could be a different prophage state or content of the strains,
     or different induction. It is not a Stp1/Stk1 effect we can claim.
- **Formal hits inside a flagged block of their own mutant: 14 of 64**, e.g. kdpC
  (`stk` TSB), kdpD and rbsD (`stk` response), B7H15_10925/10935/11085 (`stp` TSB, prophage),
  mazE (`stp_vs_stk`). Column `in_own_block` in `loso_hit_stability_mval.csv`.
- **Conclusion:** the strongest positional signals are one block per mutant (SigB/*mazEF*/*kdp*
  region higher in `stk`, prophage lower in `stp`). About a fifth of the hits lie in such
  blocks, so they may reflect strain differences in genome content rather than kinase or
  phosphatase function. DNA-level checks are needed before these hits are interpreted biologically.

## 3b. What is inside the flagged blocks (Step 11b, added 19 Nov)
Made by [`scripts/11b_genome_blocks_mval.R`](../../../scripts/11b_genome_blocks_mval.R).
No new test: it re-reads the Step 11 blocks and the Step 5 fit. Flagged blocks that
overlap in any contrast are merged into **8 regions, R1–R8** (genome order). Region labels
were written by hand from the annotation; the script checks that every gene named in a
label lies in that region (`genome_blocks_checks_mval.txt`, 8/8 PASS).

### `genome_regions_mval.csv`: the 8 regions
| region | locus tags | Mb | content (annotation) | flagged in |
|---|---|---|---|---|
| R1 | 01055–01135 | 0.22 | *ptsG*, *murQ*, *azoR* (sugar uptake) | `stp` response |
| R2 | 01445–01545 | 0.31 | *rbs* operon, *lytM*, *esxA*/*esaA* (type VII secretion) | `stk` SN |
| R3 | 04460–04540 | 0.88 | cold-shock protein, *gcvH* | `stp` SN |
| R4 | 07835–08180 | 1.56 | prophage (10 phage products) incl. *lukD*, *hlgC* annotations | `stp` SN |
| R5 | 09880–09970 | 1.94 | *splA–F* proteases, lantibiotic genes, *lukD* | `stp` TSB |
| R6 | 10910–11150 | 2.09 | prophage (8 phage products) | `stp` TSB, both responses |
| R7 | 11470–11590 | 2.19 | *sigB*/*rsbVW*, *mazEF*, *alr1*, *kdp*, *murF*, *ddl* | `stk` TSB |
| R8 | 11640–11700 | 2.22 | *sceD*, *fabZ*, *murA*, ATP synthase (*atpA/C/D/G/H*) | `stk` SN |

### `fig_genome_blocks_mval.png` + `genome_region_effects_mval.csv`: blocks in the Fig 4b layout
- **Compares:** each region (row) in the same 6 contrasts and columns as the pathway figure
  (Fig 4b): `stk` / `stp` vs `anc` × baseline TSB, baseline SN, SN response (interaction).
- **How to read:** colour and number = **mean log2 diff of the tested proteins in the
  region** (> 0 = higher in the mutant than in `anc`; in the response column > 0 = SN
  response lower in the mutant). This is not the window mean of the scan, so the numbers
  differ a little. Circle = a flagged window of that contrast lies in the region;
  triangle = not flagged (a tendency only).
- **Key points:**
  - `stk` regions: R7 +1.1 in TSB (flagged), +0.7 in SN; R8 +1.1 in SN (flagged), +0.6 in
    TSB; R2 +1.0 in SN. **The same regions are about 0 in `stp`** (−0.2 to +0.2).
  - `stp` regions: R5 −1.1 and R6 −1.1 in TSB (lower in `stp` than `anc`); R3 −1.0 and R4 −0.8
    in SN. In `stk` the baseline values of these regions are close to 0 (R5 −0.8 in SN, not flagged).
  - **R6 is the only region flagged in both mutants**, and only in the response column
    (−0.6 `stk`, −0.8 `stp`: SN response higher in the mutant than in `anc`). Something both
    mutants share relative to `anc` may lie in `anc` (shared reference; interpretation).
- **Conclusion:** every flagged region belongs to one mutant only, apart from the R6 response.
  The `stk` blocks sit together at 2.19–2.22 Mb and are all **higher** in `stk`. The `stp`
  blocks are spread out and **lower** in `stp`.

### `fig_pathway_terms_in_blocks_mval.png` + `block_term_overlap_mval.csv`: direct link to Fig 4b
- **Compares:** for the 16 terms of the pathway figure (rebuilt with the 10b rule; check
  PASS), the share of each term's tested proteins that lie inside a flagged region of
  `stk` (left) or `stp` (right). Same row order as Fig 4b.
- **Result (`stk`):**
  - ATP synthase (GO:0033178): **5/5** members in R8.
  - CL:1356 (SigB/*mazEF* cluster): **9/10** in R7.
  - vancomycin resistance (map01502) 4/12, peptidoglycan CL:932 4/13 and KW-0573 3/15:
    *alr1*, *murF*, *ddl* (R7), *murA* (R8).
  - the `stk` PASN-response terms (de novo IMP, branched-chain amino acids, diaminopimelate),
    the transport terms and the cell-killing / haemolysis terms: **0** members in a block.
- **Result (`stp`):** cell killing (GO:0001906) 3/12 and CL:4241 3/30 (*lukD*, *hlgC* in R4/R5);
  all other terms 0 (pyrimidine 1/16 in R6 for both mutants).
- **How to read it:** a term with all its members in one block cannot be told apart from a
  regional change. The block and the term shift are the **same signal** seen twice. The
  figure cannot say which is the cause (copy number vs. regional regulation; see section 3).
- **Conclusion:** two of the strongest `stk`-only pathway shifts in Fig 4b, ATP synthase and
  the SigB/*mazEF* cluster, are entirely one genomic region each. The vancomycin and
  peptidoglycan shifts are about one third inside these regions. The `stk` PASN-response shifts
  and the cell-killing shifts are not in a block. For H3 this means: part of the "not shared
  by `stp`" pathway signal may be `stk` strain background. The **shared** shifts (cell wall
  up; cell killing down) do not rest on a block in `stk`, and only 3 of 12 cell-killing
  proteins lie in `stp` blocks.

### `fig_genome_scan_labelled_mval.png`
The Step 11 scan with each flagged window labelled R1–R8 and the region contents in the
caption. Window means are recomputed (deterministic); the thresholds are the Step 11 ones
(check: same number of flagged windows). Use this version in the report instead of
`fig_genome_scan_mval.png`.

## Not done (and why)
- **SNpulse** stays optional (workflow): the primary story is TSB vs SN. The two reinjected
  `stk_pulseSN` phospho samples were checked in Step 8: removing them gives *more* site
  hits (62 → 78), so the main fit is conservative.
- No sample was excluded: no check pointed to a single bad sample. `stk_TSB_2/3` cost the
  most hits in LOSO, but their estimates stay within the ranges shown.

## Follow-up this suggests (future experiments, not prerequisites)
Sequence (or qPCR) the `stk` strain over B7H15_11465–11700 and the prophage region in
all three strains; then verify `stk1`/`stp1` transcripts. See the follow-up proposal in
the report.
