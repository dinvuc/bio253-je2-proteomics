# 2026-11-12 – mutant validation and polar effect, steps 3–4 (group mval)

## Question for today
H1: do the mutants lack the enzyme they should lack? Does `stp` also lose
Stk1, and is that consistent with a polar insertion? The peptide table will not
be provided, so peptide evidence comes from `nr_children_protein_Id`
(peptides quantified per protein per sample) and the Stk1 phosphopeptides.

## Predictions (written before running)
- anc: Stp1 + Stk1 present.
- stk: Stk1 absent (insertion codon 46/664), Stp1 unchanged.
- stp: Stp1 absent (codon 189/247); Stk1 also reduced if the insertion is polar.

## What we did
- branch `steps3-4-mutant-validation-mval`
- `Rscript scripts/03_mutant_validation_mval.R` → `mval_medi/results/step03_validation/`
  (8 checks PASS, 1 INFO line). It joins by sample name, and the Excel abundances equal the derived ones.
- `Rscript scripts/04_polar_effect_mval.R` → `mval_medi/results/step04_polar_effect/`
  (all checks PASS)
- checked the transcription unit on AureoWiki and UniProt (sources below)

## Results (key numbers)
**Stp1 / Stk1 protein (Step 3)**
- anc: both in 9/9 samples; Stp1 9–15 peptides/sample, Stk1 27–32.
- stk, Stk1: detected in **1/9** samples (`stk_pulseSN_1`, **1 peptide**, log2 16.8).
  Supplied `stk_vs_anc` diff −1.7 to −2.8, all **`lod_imputed`**. That is presence/absence
  evidence, and the imputed floor makes the loss look smaller than it is.
- stk, Stp1: unchanged (10–15 peptides; `stk_vs_anc` diff −0.17 to −0.40, FDR ≥ 0.56).
- stp, Stp1: detected 7/9, but on only **1–2 peptides** (anc: 9–15). The two
  wild-type-like values rest on almost nothing: `stp_SN_2` 18.4 (2 peptides),
  `stp_pulseSN_3` 18.1 (1 peptide). Supplied `stp_vs_anc_at_TSB` diff = −4.98, FDR 0.023
  → higher in anc. The `observed` mean here mixes single-peptide values, so it is not a precise size.
- stp, Stk1: TSB 3/3 detected at ~11 log2 on 1–2 peptides, SNpulse 2/3, SN 0/3.
  Our treat test `stp_vs_anc_at_TSB` diff = −5.42, FDR = 0.0013 → higher in anc.
- Rule-based calls (`evidence_table_mval.csv`): Stk1 in stk, Stk1 in stp, Stp1 in stp =
  "supported reduction"; Stp1 in stk = "no reduction". Single-peptide values are kept but flagged.

**Stk1 phosphosites**
- Stp1 has no phosphosites in either search.
- Site cluster S159–Y175 (STYH): T161–T172 in 8–9/9 anc samples (S159 7/9, Y175 1/9),
  **0/9 in stk** for every site, 0–6/9 in stp (T164 6, T166 5, H168 5, T172 3, S159 1, Y175 1, T161/S162 0).
  For T161–T172, where both anc and stp have ≥ 2 values in a medium, stp is lower by a median
  of **−4.4 log2** (6 site × medium rows). (Whether this cluster is the kinase activation loop is
  not checked against a source yet.)
- `HIQDSVPNVTtDVR…` (T231/S240/S242): the peptide is unique to Stk1 in `JE2_proteins.fasta`, but it is
  **sporadic in all strains**: anc 3/9, stk 3/9, stp 4/9 (STYH). In stk it appears in `SN_2`,
  `pulseSN_1` (reinjected) and `pulseSN_3`, and **not in TSB**. Median loc_prob is 0 in anc and
  0.64 in stk. Intensities of 15.7–18.6 log2 are at the low end (for comparison, T164 in anc: median 20.2).
  **Correction of the 12 Nov first look:** it is 1 sample per medium, not "3/3 per medium".
  It is not good evidence for a truncated Stk1 product. Unresolved why it appears in stk at all (carry-over? noise?).

**Polar effect (Step 4)**
- JE2: `B7H15_06330` → `06335 fmt` → `06340` → `06345 rlmN` → `06350 stp1` → `06355 stk1`, all + strand.
  stp1/stk1 overlap by 4 bp. **stk1 is the last gene**; 227 bp gap to `rsgA`.
- Measured neighbours (06320–06375, 9 proteins): 27/27 detected; largest |diff| mutant vs anc
  = 0.74 (rpmB, stp TSB, lower in stp); **0 formal hits** in 54 protein × contrast rows.
  Small consistent shifts in stp: 06320 about −0.6 and rlmN +0.71 in TSB (below the |log2FC| > 1 threshold).
- Stk1 peptides per sample, median per medium (TSB/SNpulse/SN): anc 32/32/28, **stk 0/0/0, stp 1/1/0**.
  So Stk1 in stp drops to almost the stk level, as a polar insertion predicts. This is
  consistent with polarity but does not prove the mechanism.

## Sources
- AureoWiki SAOUHSC_01186 (stp1) and SAOUHSC_01187 (pknB): Operon section lists
  SAOUHSC_01182 > … > 01186 > 01187 (MicrobesOnline; SigA-promoter unit from Mäder et al. 2016,
  PLoS Genet, PMID 27035918). https://aureowiki.med.uni-greifswald.de/SAOUHSC_01186 .
  This is **NCTC 8325, not JE2**. We map 01182–01185 to JE2 genes by gene order (an assumption).
- UniProt Q2FZ65 (= SAOUHSC_01186, 247 aa), Q2FZ64 (= SAOUHSC_01187, 664 aa), rest.uniprot.org.
- Precedent for a polar bursa aurealis insertion in stp1: Sun et al. 2012, PNAS,
  doi:10.1073/pnas.1205952109 (different strain/assay; see `mval_medi/docs/stk_stp_rescue_evidence_mval.md`).
- Insertion coordinates: `TRANSPOSON_INSERTIONS` (`R/constants.R`, NTML map). The NTML website
  could not be read by the AI's fetch tool (the content loads by JavaScript): **to do: check NE217/NE1919 by hand**.

## Problems / surprises
- Every "wild-type-like" value in a mutant rests on 1–2 peptides. Without peptide counts these
  values look like failed knockouts.
- Alternatives to polarity that we cannot exclude: Stk1 unstable without Stp1, feedback
  regulation, a secondary mutation, peptide assignment errors, carry-over between runs
  (DIA), and sample identity. Our data cannot separate transcriptional polarity from these.
- Contradicting observation kept: stp still yields a little Stk1 in TSB/SNpulse (1–2 peptides,
  ~5 log2 lower than anc), so the knockdown is strong but not complete.

## AI use
An AI coding assistant wrote and ran the two scripts and looked up the AureoWiki/UniProt pages.
We still need to check the numbers above against the CSVs and figures.

## Next steps
- Group: check this entry, the figures and the evidence table.
- Check NE217/NE1919 on ntml.unmc.edu by hand.
- Steps 6–7: interactions relative to anc; stp vs stk scatter. Remember that stp lacks most Stk1.
