# The proteins of Fig 2: function, STRING neighbourhood and hypotheses (Step 6d)

Made by [`C_mutant_validation_DinMete/scripts/06d_response_hits_string_mval.R`](../../../scripts/06d_response_hits_string_mval.R),
2026-11-18. Figures: see [`figure_guide_mval.md`](figure_guide_mval.md).
**Status: AI draft (Claude Code). A group member has to check every number and source
before this counts** (rules section 8). The group asked for the biology, so this file
interprets it (rules section 10). Everything under *Hypothesis* is **interpretation**,
not a result.

**Reading key.** response = log2(SN) − log2(TSB), positive = higher in SN (PASN).
delta = response in the mutant − response in `anc` (= −interaction diff). Hit rule: limma
`treat(lfc = 1)`, BH-FDR < 0.05; n = 3 per group; all estimates are observed (no
imputation). Every `stp` statement: the `stp` strain also lacks most Stk1 (Step 3).
Robustness (LOSO, Step 11) is given in brackets.

Names and functions come from `read_protein_metadata()` (UniProt entry of the NCTC 8325
ortholog, `ortholog_quality == "good"` for all 7). Literature was found with PubMed. For
every paper I checked only the abstract, so the full text should be checked before the
poster.

## Overview

| protein (JE2) | function (UniProt, NCTC 8325 ortholog) | what changes | robustness | how much we trust it |
|---|---|---|---|---|
| `B7H15_11550` **kdpD** | Sensor histidine kinase KdpD, K⁺ two-component system (Q2FWH7) | `stk`: response −4.09 vs `anc` +1.03, delta −5.13 (FDR 1e-10) | 14/18, in the flagged `stk` block | **low**: peptide check below |
| `B7H15_14675` **arcA** | Arginine deiminase (Q2FUX7) | `stk`: low in TSB (−3.74 vs anc), rises in SN; delta +4.98 (FDR 0.003) | 16/18 | medium |
| `B7H15_14670` **argF** | Ornithine carbamoyltransferase, same operon as arcA (Q2FUX8) | `stk`: −1.63 in TSB, +0.74 in SN; delta +2.37 (FDR 1e-4) | 18/18 | **high** |
| `B7H15_14515` | Betaine-aldehyde dehydrogenase (FASTA); UniProt "Putative aldehyde dehydrogenase AldA" (Q2FV10) | `stk`: +1.79 in TSB, same as anc in SN; delta −2.07 (FDR 0.002) | 18/18 | **high** (as a TSB difference) |
| `B7H15_01460` **rbsD** | D-ribose pyranase (Q2G1A5) | `stk` SN +5.61 vs anc; delta +5.36 (FDR 0.011) | 14/18, in a flagged `stk` block | low |
| `B7H15_11425` **ilvA** | Threonine deaminase, isoleucine biosynthesis (Q2FWJ9) | `stp`: response −2.85 vs anc −0.49; delta −2.36 (FDR 0.039) | 9/18 | low |
| `B7H15_06000` | Fibrinogen-binding protein, secreted (Q2FZC2, unreviewed) | `stp`: −2.89 in TSB, normal in SN; delta +2.31 (FDR 0.077, **exploratory, not a hit**) | (18/18 in stp vs stk) | low |

**STRING (v12.5, the course's own JE2 organism `STRG0B09NBK` "JE2_Bio253_HS26", score
≥ 0.7, up to 6 partners each; 37 proteins, 70 links).** Among the 7 proteins themselves
there is only one link: arcA–argF (0.918). Each protein has its own neighbourhood:
kdpD → kdpA/B/C/E + Y1819 (B7H15_09315, a Step 9c phospho-site protein);
arcA/argF → the second arc copy (arcC1 00330, argF 00335), argF 06085, argB, rocD, carA,
arginase, nitric oxide synthase, B7H15_11190; ilvA → leuB, acetolactate synthase 11390,
threonine synthase, metE, trpB, glyA; B7H15_14515 → betA, a GbsR/MarR-family regulator
(B7H15_14525, the neighbouring gene), crtP, amidohydrolase, lacD, aspartate aminotransferase
14470; rbsD → rbsK, two LacI-family regulators, a lactose-specific transporter;
B7H15_06000 → none. Only 6 of the 70 links join two modules, and none of them is a direct
link between two Fig 2 proteins. **The hits are not one pathway; they are 6 separate
modules.** That argues against a single Stk1-controlled programme. It fits better with
several local effects (or artefacts) next to each other. In our data, partners that move
with their hit in `stk` are betA (−1.76, not a hit), acetolactate synthase 11390 (−1.16)
and Y1819 (−1.12). The arc partners, rbsK and the LacI regulators stay near 0.

---

## 1. kdpD: stk vs stp (main focus)

### What KdpD is
KdpD is the sensor kinase of the KdpDE two-component system. In *S. aureus*, KdpDE
controls the *kdpFABC* K⁺-transporter operon and its own *kdpDE* genes. The second
messenger c-di-AMP binds the USP domain of KdpD and **inhibits** *kdp* induction
(Moscoso et al. 2015, PMID 26195599, [doi:10.1128/JB.00480-15](https://doi.org/10.1128/JB.00480-15);
structure: Dutta et al. 2021, PMID 33989637, [doi:10.1016/j.jbc.2021.100771](https://doi.org/10.1016/j.jbc.2021.100771)).
KdpDE responds to external K⁺. In *S. aureus* it acts more as a regulator of virulence
genes than as the main K⁺ importer. It is **upregulated by Agr/RNAIII** (Xue et al. 2011,
PMID 21422185, [doi:10.1128/IAI.01180-10](https://doi.org/10.1128/IAI.01180-10)). Under
osmotic stress (2 M NaCl), *kdp* is induced via KdpDE (Price-Whelan et al. 2013,
PMID 23963175, [doi:10.1128/mBio.00407-13](https://doi.org/10.1128/mBio.00407-13)).

### What the data show (numbers: `checks_mval.txt`, `kdp_contrasts_mval.csv`)

**TSB (baseline).** `stk` has more of **every** kdp protein than `anc`: KdpD +1.35
(FDR 0.28, higher in `stk`), and KdpC, KdpB, KdpE also higher. `stp` ≈ `anc` (KdpD +0.08).
The whole region B7H15_11470–11590 (SigB operon … *mazEF* … *kdp* … *cshA*) is
about +1.5 log2 in `stk` (Step 11 genome scan, family-wise p < 0.001).

**SNpulse (15 min).** No change in any strain (fig_kdp_operon_raw_mval.png).

**SN (15 h).** On the protein level (the numbers behind Fig 2):
`anc` +1.03 (2/3 replicates), `stk` −4.09 (2/3), `stp` not testable (1/3).
**On the peptide level the picture is different** (`fig_kdpD_peptides_mval.png`):

- In TSB, KdpD is seen with 12–27 peptides per sample. In SN it is seen with **0–3
  peptides in every strain**. Only 5 % of the TSB peptide observations remain (KdpE 8 %).
  That is the 9th (KdpE: 12th) strongest loss among 1,121 well-covered proteins (median
  protein: 100 %). Sequencing depth is the same in TSB and SN (median 32,439 vs
  30,788 peptides per sample), so this is not a depth artefact.
- Peptide-matched (each SN peptide against the same peptide in TSB, same strain):
  `stk` −1.0 to −4.6 (median −4.0, 4 peptides). In `anc`, two peptides go *down*
  (−0.3, −2.7), and two go *up* (+2.9, +4.3). One of the rising peptides is a Met-oxidised
  form that was seen in only one TSB sample. KdpE falls in **all** strains
  (`anc` −4.1, `stk` −3.1, `stp` −4.9).
- The `stp` value rests on a single peptide in one sample, and that peptide was never
  seen in `stp` TSB.

### Interpretation (my reading; the group should check it)
1. **The dominant effect is shared:** in PASN, the KdpDE system (and KdpB) is largely
   switched off **in all three strains**. This is a strong part of the common PASN
   response. Our protein-level model misses it, because so few peptides remain that the
   protein value becomes unstable.
2. **The `stk`-specific "opposite response" is most likely an artefact of summarising
   1–3 peptides.** The `anc` "rise" comes from two inconsistent peptides. The reported
   SE (0.18) is far too optimistic for an estimate built on 2 peptides × 2 replicates.
   With LOSO 14/18 and the genomic block, kdpD should **not** be presented as an
   `stk`-specific response. On the poster I would say: "kdp shut down by PASN in all
   strains; the stk-specific estimate is not supported at peptide level."
3. **stk vs stp in TSB:** the extra kdp in `stk` comes with its neighbours (SigB, *mazEF*,
   *cshA* …). It is **absent in `stp`**, even though `stp` also lacks most Stk1. If the
   loss of Stk1 kinase activity caused it, `stp` should show at least part of it. This
   supports the Step 11 idea: a **strain-specific regional difference** (e.g. a
   duplication in NE217), not Stk1 signalling. A residual-Stk1 effect in `stp` is the
   alternative we cannot exclude. qPCR / genome sequencing of the region would decide.

### Hypothesis: why does PASN switch kdp off? (interpretation)
- **Agr link (supported by our data, plausible):** KdpDE is Agr/RNAIII-activated
  (Xue 2011, above). In our data, **Hld (δ-toxin, encoded within RNAIII) falls strongly in SN in
  all strains**: `anc` −4.52 (FDR 3e-5), `stk` −3.54 (FDR 0.001), `stp` −2.63
  (FDR 0.04); AgrB −1.4 (not a hit). Lower Agr activity in PASN would explain lower
  KdpDE in every strain. Caveat: the Hld ortholog call is `ortholog_quality == "check"`
  (26-aa peptide, Q2FWM8), so treat the annotation with care. That *hld* lies in
  RNAIII: Novick et al. 1993, PMID 7691599,
  [doi:10.1002/j.1460-2075.1993.tb06074.x](https://doi.org/10.1002/j.1460-2075.1993.tb06074.x)
  (the abstract describes RNAIII as the agr effector; the *hld* location is in the full
  text, **please check**). Reduced expression of global regulators and virulence factors
  in *S. aureus* next to *P. aeruginosa* was also seen in a co-culture RNA-seq study
  (Tognon et al. 2019, PMID 30630428, [doi:10.1186/s12864-018-5398-y](https://doi.org/10.1186/s12864-018-5398-y)).
  That study is direct co-culture, not spent medium.
- **c-di-AMP (speculative):** higher c-di-AMP inhibits *kdp* induction (Moscoso 2015).
  We cannot measure c-di-AMP.
- **K⁺/osmotic status (speculative):** Kdp is induced by K⁺ limitation and high
  osmolarity (Price-Whelan 2013). Its shutdown suggests that the SN cultures were *not*
  K⁺-limited. We do not know the K⁺ content of PASN.

**Possible checks:** Agr targets in our data (PSMs, other RNAIII-dependent proteins);
qPCR of *kdpD*, *kdpA* and RNAIII; KdpD by Western blot instead of 1–3 peptides.

---

## 2. arcA + argF (arginine deiminase pathway), stk

**Function.** Arginine deiminase (arcA) and the catabolic ornithine carbamoyltransferase
(argF, B7H15_14670) are the first two steps of the ADI pathway: arginine → citrulline →
ornithine + carbamoyl-P → ATP. It runs under anaerobic conditions and is activated by ArcR
(arcR is in the same gene run, B7H15_14655) (Makhlin et al. 2007, PMID 17557828,
[doi:10.1128/JB.00592-07](https://doi.org/10.1128/JB.00592-07)).

**Data.** Both proteins are **low in `stk` in TSB** (arcA −3.74, argF −1.63 vs `anc`) and
rise in SN (arcA +4.06 in `stk` vs −0.92 in `anc`). Two enzymes of one operon move together,
so this is not a single-protein accident. `stp` behaves like `anc` (and arcA/argF are also
hits in `stp` vs `stk`). JE2 has a **second arc gene copy** (B7H15_00330–00350; STRING
lists its argF and arcC1 as partners). Its arcA rises in SN in `anc` (+1.42) and less in
`stk` (+0.35), and its argF (00335) responds like in `anc` (delta −0.66). So the `stk` effect is specific to the main operon B7H15_14655–14675.

**Hypothesis (interpretation).** PASN contains HQNO, which pushes *S. aureus* towards
fermentative growth (Orazi & O'Toole 2017, PMID 28720732,
[doi:10.1128/mBio.00873-17](https://doi.org/10.1128/mBio.00873-17); Tognon 2019 above).
Fermentative growth is when the ADI pathway pays off. PknB (= Stk1) influences central
metabolism and purine synthesis (Donat et al. 2009, PMID 19376851,
[doi:10.1128/JB.00117-09](https://doi.org/10.1128/JB.00117-09)). A different energy state
of `stk` could make it repress *arc* in rich TSB and depend more on arginine catabolism
once respiration is blocked. Why `stp` (low Stk1) does not do the same is open. The options
are residual Stk1 in `stp`, a phospho-*balance* effect, or a second mutation in NE217.
This is the most robust `stk`-specific finding (argF 18/18).

## 3. B7H15_14515 (+ betA), stk

**Function.** The FASTA annotates it as betaine-aldehyde dehydrogenase. UniProt (Q2FV10,
unreviewed) only says "Putative aldehyde dehydrogenase AldA". The gene next to it,
B7H15_14510, is betA (choline dehydrogenase). Together they would make glycine betaine
from choline, an osmoprotectant. **The pathway assignment rests on automatic
annotation.** I found no PubMed paper for this gene pair in *S. aureus*.

**Data.** `stk` is +1.79 in TSB, equal to `anc` in SN. betA follows (`stk` delta −1.76,
not a hit). This is a TSB baseline difference, stable 18/18, outside the flagged block.

**Hypothesis (speculative).** `stk` may be under some osmotic or cell-envelope stress in
TSB; PknB loss changes cell-wall metabolism and autolysis (Donat 2009). Note that kdp is
also osmo-related, but its `stk` excess is explained by the genomic block, so I would
not combine the two.

## 4. rbsD, stk

**Function.** D-ribose pyranase (Q2G1A5). Operon with rbsK (ribokinase) and the ribose
transporter RbsU; ribose uptake and use.
**Data.** `anc` falls in SN (−2.56), `stk` rises (+2.79). This rests on 2/3 replicates in
both SN groups, rbsK does not follow (+0.37), the protein lies in a flagged `stk` block
(B7H15_01445–01545), and LOSO is 14/18.
**Hypothesis (weak).** Nutrient-depleted PASN could change sugar use. Given the block and
the missing partner, a regional artefact is at least as likely. No literature claim.

## 5. ilvA, stp (and probably stk)

**Function.** Threonine deaminase, the first step of isoleucine synthesis (Q2FWJ9). It is
in the ilv–leu run (B7H15_11385–11425). The ilv–leu operon is repressed by CodY (Lei et
al. 2012, PMID 23056478, [doi:10.1371/journal.pone.0046836](https://doi.org/10.1371/journal.pone.0046836)).
CodY is more active when GTP and branched-chain amino acids are high (Stenz et al. 2011,
PMID 21539625, [doi:10.1111/j.1574-695X.2011.00812.x](https://doi.org/10.1111/j.1574-695X.2011.00812.x)).

**Data.** In `stp` the SN response is −2.85 (`anc` −0.49). In `stk`, ilvA is detected in only
1/3 SN samples, so our test skips it. The supplied prolfquapp estimate says the
response is also lower in `stk` (3.01, FDR 0.001). ilvD and B7H15_11390 trend the same way
in `stk` (−1.2, not hits). LOSO 9/18: fragile.

**Hypothesis (interpretation).** If both strains lacking Stk1 show it, Stk1 loss is the
shared factor (rules: `stp` lacks most Stk1). Donat 2009 links PknB to purine (GTP)
metabolism. Higher GTP → more active CodY → stronger ilv–leu repression in PASN is a
testable chain, but every link is unmeasured here. Isoleucine-derived branched-chain fatty acids feed back on
SaeRS virulence signalling (Alqahtani et al. 2024, PMID 39382300,
[doi:10.1128/jb.00191-24](https://doi.org/10.1128/jb.00191-24)).

## 6. B7H15_06000 (exploratory only), stp

**Function.** Secreted fibrinogen-binding protein (Q2FZC2, unreviewed; keywords Secreted,
Virulence).
**Data.** Low in `stp` in TSB (−2.89), normal in SN. This is a baseline difference, and
FDR 0.077 makes it **not a hit**.
**Hypothesis (interpretation).** Stp1/Stk1 influence secreted virulence factors: an
*stp1* deletion lowered hemolysin expression (Burnside et al. 2010, PMID 20552019,
[doi:10.1371/journal.pone.0011071](https://doi.org/10.1371/journal.pone.0011071)), and
Stk1/Stp1 phosphorylate SarA/MgrA-family regulators (Sun et al. 2012, PMID 22927394,
[doi:10.1073/pnas.1205952109](https://doi.org/10.1073/pnas.1205952109)). Lower
fibrinogen-binding protein in `stp` would fit, but we did not check that this gene is
controlled by those regulators.

---

## Take-home for the poster (draft)
- The 7 proteins fall into **6 STRING modules** (course JE2 organism) with no direct link
  except arcA–argF, so there is no single
  Stk1-controlled response programme.
- **Most solid:** the arc operon (arcA + argF) in `stk`, i.e. arginine catabolism switched
  off in TSB and switched on in PASN.
- **kdpD:** the PASN effect is real but **shared by all strains** (kdp ~95 % of peptide
  evidence lost, probably via lower Agr activity, shown by Hld). The `stk`-specific
  estimate is not supported at peptide level. The extra kdp in `stk` in TSB is regional
  and absent in `stp`.
