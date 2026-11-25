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
`fig_genome_scan_mval.png`. Dotted vertical lines mark every region R1–R8 in all six
panels, so you can see whether a block flagged in one contrast also moves in the others.

### `fig_genome_scan_labelled_compact_mval.png`
The same scan redrawn narrow for the poster (7 x 7.4 in instead of 13 x 7 in). **Compares**
`stk` vs `anc` (left) with `stp` vs `anc` (right) in TSB, PASN (= SN) and the PASN response
(rows). > 0 = **higher in the mutant**; in the response row > 0 = **response lower in the mutant**.
Windows, thresholds and flags are unchanged. R1-R8 are named once on top, with a short key below.
The 5% threshold holds per panel and is not corrected for looking at six panels.
`fig_genome_scan_labelled_compact_transparent_mval.png` **compares** exactly the same, without
background (for the poster). The R labels and the grey row titles stay filled.

**Conclusion (Steps 11, 11b and 11c together).**
1. **A few short blocks carry part of the mutant signal.** Most of the genome sits near 0.
   Eight regions pass the genome-wide 5 % threshold: R7, R8 and R2 are **higher in `stk`
   than in `anc`** (about +1); R3, R4, R5 and R6 are **lower in `stp` than in `anc`** (about −1).
2. **At baseline, every block belongs to one mutant only.** No block is shared by `stk` and
   `stp` in TSB or SN. `stp` also lacks most Stk1 (polar effect inside *stp1-stk1*), so if a
   block came from losing Stk1, both mutants should show it. They don't.
3. **No block is a polar effect.** The *stp1/stk1* locus (pink line, 1.23 Mb) sits in a flat
   stretch of every panel. R6 and R7 are 860 and 953 kb away. A polar effect only reaches genes
   downstream in the same transcript.
4. **R7 is a stretch of DNA, not an operon or a regulon (Step 11c).** About 23 kb and 8 putative
   TUs on both strands go up together in `stk`, in all replicates and all media. The SigB targets
   are not up, so this is not more SigB activity. It is most likely a regional effect of the
   `stk` strain background (e.g. an extra copy of that DNA). R8 lies right next to it (2.22 Mb),
   and `stk` values stay at about +0.5 between the two, so R7–R8 may be one larger region
   (B7H15_11460–11700). Not tested.
5. **R6 is a prophage-level change in `stp`.** The phiSa3-like prophage in *hlb* is lower in
   `stp` at TSB, not in `stk`. In the SN-response (interaction) panels it is flagged in *both*
   mutants. Both contrasts share the same `anc` reference, so that flag may say more about how
   the prophage in `anc` responds to PASN than about either mutant.
6. **Consequence for H3.** Two of the "`stk`-only" pathway shifts in Fig 4b are essentially one
   block each: the SigB/*mazEF* cluster (9/10 in R7) and ATP synthase (5/5 in R8). Treat them as
   positional or strain-background signals, not as Stk1 signalling. The **shared** shifts (cell
   wall up, cell killing down) and the `stk` PASN-response terms do not rest on a block.
7. **Genome-wide, effects follow transcription units** (same-TU neighbours r = 0.31, divergent
   neighbours r ≈ 0). Some co-changing neighbours are therefore expected. The flagged blocks
   stand out because they span many TUs at once.

Caveats: n = 3 per group; the scan tests windows, not genes. The causes in 4 and 5 are
hypotheses. A DNA-level check (WGS or genomic qPCR of B7H15_11460–11700 in `stk`, and of the
R6 prophage in all three strains) can decide them.

*One-line version for the report:* "The genome scan finds a few mutant-specific blocks far
from the transposon. The largest (R7, `stk`) spans many operons on both strands without
activating SigB. Part of the `stk`-only pathway signal (SigB/*mazEF*, ATP synthase) is therefore
likely positional or due to strain background, whereas the shared pathway shifts are not."

**Are the blocks "an effect of stk"?** Our professor's hypothesis was that the blocks are
caused by `stk`. The data support an effect of the **`stk` strain**, but not an effect of
losing **Stk1 kinase activity**:
- The blocks belong to one strain only. `stp` also lacks most Stk1, yet R7 is +1.12 (higher)
  in `stk` and −0.08 (unchanged) in `stp`. The reverse holds for R6: −1.09 (lower) in `stp` and
  −0.12 in `stk`.
- R7 spans many putative TUs on both strands, does not activate SigB, and is the same in every
  medium. A signalling response would usually depend on the condition.
- Genome-wide, `stk` effects do follow TUs (point 7). So Stk1 may well regulate genes; the
  blocks just are not the evidence for it.
- Limits: `stp` has *reduced*, not zero, Stk1 and also lacks the phosphatase, so the comparison
  between the two strains is strong but not airtight. n = 3, protein data only.
- What would decide it: complementation of *stk1* (if R7 stays high, it is not Stk1), WGS or
  genomic qPCR of B7H15_11460–11700, and a second, independent *stk1* mutant.

**What we mean by "strain background".** Any genetic difference between a mutant strain and
our `anc` *other than* the transposon insertion we are studying. Each NTML mutant comes from
a single cell and has its own history of growth and freezing, so it can carry unrelated changes:
SNPs or small indels, **copy-number changes** (a duplicated stretch gives about +1 log2; our
hypothesis for R7), deletions, a lost, rearranged or more active **prophage** (a possibility
for R6), plasmid differences, or mutations in often-mutated regulators (e.g. *agr*, *rsbU*).
NTML library: Fey et al. 2013, *mBio*, doi:10.1128/mBio.00537-12 (still to be checked by the group).

| explanation | cause | expected pattern |
|---|---|---|
| direct effect | Stk1 kinase missing | changes in Stk1 substrates or the pathways downstream of them |
| polar effect | the transposon weakens genes downstream in the same operon | only *stp1-stk1* and the genes right after it (1.23 Mb) |
| indirect (trans) effect | Stk1 loss alters regulators, which change genes far away | also seen in `stp` (which lacks most Stk1); depends on the condition |
| strain background | an unrelated difference in that clone's genome | **one** strain only, often a **block** of neighbours, the same in **every medium** |

R7 fits the last row on all three points. "Strain background" names where the cause
probably lies, not the cause itself. It is not proven, only better supported than an Stk1
effect for R7.

## 3c. Genome context of R7 and R6 at baseline TSB (Step 11c)

Made by [`scripts/11c_R6_R7_context_mval.R`](../../../scripts/11c_R6_R7_context_mval.R) (~20 s;
no new test). Question: do the changed proteins in the two baseline-TSB blocks sit in one
operon, many operons or a whole stretch of DNA, and could a polar effect of the transposon
explain them? **Putative TU** (transcription unit) = neighbouring genes on the same strand
with gaps < 150 bp, as in Steps 6b/6c. This is a hint from the sequence, not a measured operon.

### `fig_R6_R7_genome_context_mval.png`
- **A (R7, 2.19 Mb), stk vs anc and stp vs anc in TSB.** Gene map (arrows = strand, colour =
  `stk_vs_anc_at_TSB`, blue bars = putative TUs) and per-gene log2 diffs (> 0 = **higher in the
  mutant**). In `stk`, 7 of the 8 tested TUs in the ~23 kb block are higher by > 0.5 log2,
  on both strands: *rsbU-V-W-sigB*, *mazEF-alr1-acpS*, *kdpABC* (+2.3, the strongest), *kdpDE*,
  *murF-ddl*. In `stp` the same TUs are about 0. The rise crosses TU boundaries and strands,
  so **it is not one operon**.
- **B (R6, 2.09 Mb), same layout, colour = `stp_vs_anc_at_TSB`.** A prophage inserted in *hlb*
  (*hlb* fragments on both sides; integrase, *scn*, *chp*, staphylokinase at the edges:
  phiSa3-like). Most phage proteins (capsid, portal, replication module; several '−' strand TUs)
  are **lower in `stp` than in `anc`** (about −2). A few phage proteins near the left end are
  higher in `stp`. In `stk` the region is about 0.
- **C (raw per-sample values).** Block score = mean over the block's proteins of
  log2(sample) − mean log2(anc TSB). R7: all three `stk` replicates are about +1.1 in TSB, and
  `stk` stays up in SN (+0.6) and SNpulse (+1.0). R6: all three `stp` TSB replicates are
  negative (−0.6 to −1.3; replicate 3 is the weakest).
- **D (genome-wide).** Pearson r of the log2 diffs of adjacent gene pairs. In `stk` TSB, pairs
  in the same putative TU are alike (r = 0.31 without R6/R7, 95% CI 0.24–0.37). Divergent pairs,
  which have their own promoters, are not (r = −0.05, CI −0.19 to 0.09). `stp` shows the same
  pattern, but weaker (0.19 vs 0.08). So mutant effects follow transcription units, which
  suggests they mostly act at the level of transcription. These are protein data, so this is
  indirect.
- **E (SigB test).** In `stk` the *sigB* operon proteins are higher (+0.2 to +1.0). The
  classic SigB target Asp23 (outside R7) is unchanged (−0.01), and so is CrtN (+0.01; CrtM
  +0.35). More SigB protein, together with more anti-sigma factor RsbW, but **no sign of more
  SigB activity**. Tamber et al. 2010 even report *lower* SigB activity in a USA300
  *pknB*/*stk1* mutant (doi:10.1128/IAI.00296-10).

**Polar effect?** No. The transposons sit in *stp1*/*stk1* at 1.23 Mb, 860 kb (401 putative
TUs) before R6 and 953 kb before R7. A polar effect only reaches genes downstream in the
*same* transcript. Also, `stp` lacks most Stk1 (polar effect inside *stp1-stk1*), yet R7 is
flat in `stp`. So R7 does not simply follow Stk1 loss, and R6 (only in `stp`) does not follow
it either.

**Reading.** R7 is a **regional** effect specific to the `stk` strain. It covers many TUs on
both strands, has sharp-ish edges, holds in every medium and doesn't activate the SigB regulon.
This fits a copy-number change of this stretch of DNA in NE217 (one extra copy ≈ +1 log2), or a
local regulator; *kdp* rises further on top of that. R6 is a prophage-level change in `stp`,
e.g. less spontaneous induction or a different prophage state (hypothesis). Both are
hypotheses that a DNA-level check (WGS or genomic qPCR across B7H15_11460–11600 and the
prophage) can decide. All proteins in R6/R7 are observed in both groups (no LOD imputation).
n = 3 per group. The TU-level numbers are descriptive; the genome-wide test is the Step 11
permutation scan.

Tables: `R6_R7_genes_mval.csv` (gene, strand, gap, TU, diffs, estimate type),
`R6_R7_tu_effects_mval.csv` (TU means per contrast), `operon_coherence_mval.csv` (r per pair
type, contrast, with/without R6–R7), `R6_R7_context_checks_mval.txt` (9 PASS).

### `fig_R6_R7_genome_context_compact_mval.png`
Panels A + B of the figure above, redrawn narrow for the poster (7 x 6.4 in instead of a
16 in wide figure). **Compares** the same: `stk` vs `anc` and `stp` vs `anc` in TSB, > 0 =
**higher in the mutant**, with the same data and genome windows. Arrow colour is `stk` vs `anc`
in A and `stp` vs `anc` in B. Only layout changed: shorter titles, shared legend below,
fewer gene labels, "putative operon" instead of "putative TU". Panels C-E are not in it.
`fig_R6_R7_genome_context_compact_transparent_mval.png` **compares** exactly the same, without
background (for the poster). White arrows are not see-through: white = difference near 0.

## 3d. STRING network of the R7 and R6 proteins (Step 11d, `fig_R6_R7_string_network_mval.png`)

Script: `scripts/11d_R6_R7_string_mval.R`. **Compares:** nothing new is tested. The figure
shows which of the block proteins STRING links to each other, and how they group.
- **Nodes** = the Step 11 block proteins (R7: 25, R6: 49). Node colour = log2 diff in the
  contrast that flagged the block, at TSB (A: `stk_vs_anc`, B: `stp_vs_anc`; > 0 = higher in
  the mutant, i.e. lower in anc). Grey = not tested (not enough observed values). All tested
  values are observed (no LOD imputation).
- **Links** = STRING v12 links for the course genome STRG0B09NBK, combined score ≥ 0.4 (medium
  confidence), from the STRING API (cached in `R6_R7_string_edges_mval.tsv`). Line width =
  score. **Dashed** = the link rests only on gene neighbourhood, so it mostly restates that the
  genes sit next to each other (56 of 123 links). We used the API because the course file
  `string_physical_links.txt` has only 4 links inside R6.
- **Shaded areas** = Louvain clusters (seed fixed). Each label gives (1) the shared STRING term
  with the lowest Fisher p vs the genome (≥ 2 members; terms with > 300 genome proteins and
  "Coiled coil" excluded) and (2) the members' GenBank product names. The p is
  **descriptive only**: the region was already selected.
- Check: kdpABC and sigB–rsbW–rsbV–rsbU each fall into one cluster (known operons).

| cluster | function (STRING term / products) | n (tested) | mean log2 diff, flagging contrast | other mutant |
|---|---|---|---|---|
| R7-1 | Transcription regulator activity: sigB operon, mazEF, acpS, cshA | 11 (11) | +0.91 (stk) | −0.05 (stp) |
| R7-2 | Potassium transport: kdpABC, kdpDE | 5 (5) | +1.89 (stk) | −0.26 |
| R7-3 | KEGG "Vancomycin resistance": alr1, murF, ddl, rodA | 4 (4) | +0.73 (stk) | +0.27 |
| R7-4 | no specific term: CsoR regulator + 2 hypotheticals | 3 (3) | +1.13 (stk) | −0.39 |
| R6-1 | no specific term: head-tail adapter, major capsid, portal | 6 (2) | −2.48 (stp) | −0.05 (stk) |
| R6-2 | Pyrimidine metabolism: terminase, HNH endonuclease, dUTPase, peptidase | 4 (1) | −0.12 (stp) | +0.73 |
| R6-3 | no specific term: tail, tail tape measure | 2 (1) | +1.27 (stp) | +1.24 |
| R6-4 | no specific term: DnaD-like replication protein, ATPase | 2 (1) | −2.36 (stp) | −0.30 |

**How to read it.**
- **R7:** every cluster is higher in `stk` (+0.7 to +1.9), in `stp` they stay near 0. So
  four functionally unrelated groups (sigma B regulation, K⁺ uptake, D-Ala-D-Ala /
  peptidoglycan precursors, a CsoR regulator) move together. That fits the Step 11c
  reading (a regional, DNA-level effect) better than one pathway being switched on.
  "Vancomycin resistance" is only the name of the KEGG map that contains alanine racemase
  and D-Ala-D-Ala ligase (KEGG map01502). It does not mean these strains are resistant.
- **R6:** STRING knows little about phage proteins. 35 of 49 proteins have no link ≥ 0.4,
  and most links are neighbourhood-only. The clusters are phage modules
  (head/capsid/portal, tail, replication). The R6-1 and R6-4 means rest on only 2 and 1
  tested proteins, and each R6 cluster has at most 2 tested proteins. So use the per-gene
  dot track in Step 11c (panel B), not these cluster means.
- Only 1 STRING link joins R6 and R7. They are two separate blocks, functionally as well.
- n = 3 per group. Reminder: `stp` also lacks most Stk1 (polar effect), yet R7 is flat in
  `stp`.

Tables: `R6_R7_string_clusters_mval.csv` (protein → cluster, function, diffs),
`R6_R7_string_cluster_terms_mval.csv` (top 3 terms per cluster with k, term size, Fisher p),
`R6_R7_string_checks_mval.txt`.

**STRING's own images** (`string_original_{both,R7,R6}_mval.png`, string-db.org API
`highres_image`, same proteins, score ≥ 0.4): MCL clustering (inflation 3) with
`show_cluster_names=1`, so each colour is an MCL cluster named after its most typical STRING
terms (cut short by STRING; "Mixed, incl." = no single function). Edge colours = evidence
type. No fold changes are shown. The clusters agree with ours, except that MCL merges R7-1
and R7-3 (sigB/mazEF + murF/ddl) into one cluster.

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
