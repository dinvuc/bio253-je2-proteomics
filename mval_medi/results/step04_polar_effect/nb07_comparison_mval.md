# Notebook 07 vs our steps 3–4: agreement check

**Question:** does the starter-pack notebook `analyses/07_transposon_mutants.qmd` (rendered
12 Nov 2026, peptide layer `read_peptides()`) contradict our Step 3 (mutant validation) and
Step 4 (polar effect) results?

**Short answer: no contradiction.** Every number that both analyses report is the same. Notebook 07
adds one new layer (*where* in the protein the peptides sit), and that layer **strengthens**
two of our cautious statements without overturning any of them.

Reproduce: `Rscript scripts/04c_peptide_crosscheck_mval.R` → `peptide_crosscheck_mval.csv`
(Stp1/Stk1, all 27 samples: our Step 3 peptide count next to the peptide-layer count, plus
which residues each detected peptide covers) and `peptide_crosscheck_checks_mval.txt`.
**Check:** peptide-layer counts equal Step 3 `n_peptides` in all 54 protein × sample rows (PASS).
Both count peptide *forms*, e.g. an oxidised and an unmodified copy of the same sequence count
as two.

## 1. Points where both agree

| point | notebook 07 | our steps 3–4 |
|---|---|---|
| insertion sites | Stk1 codon 46/664 (7 %), Stp1 codon 189/247 (76 %) | same (`validation_checks_mval.txt`, PASS) |
| gene arrangement | stp1/stk1 overlap by 4 bp, one operon | same (`polar_checks_mval.txt`, PASS) |
| Stk1 in stk | 1 detection in 9 samples | 1/9 samples (`stk_pulseSN_1`), 1 peptide |
| Stk1 in stp | strongly reduced (7 detections, 4 peptides) | TSB 1/1/2, SNpulse 2/1/–, SN 0 peptides per sample: identical |
| Stp1 in stk | unchanged (peptides like anc) | "no reduction", diff −0.17 to −0.40 |
| Stp1 in stp | erratic, few peptides | 7/9 samples, only 1–2 peptides each (anc 9–15) |
| why Stk1 is low in stp | polar effect of the stp1 insertion | "consistent with polarity, does not prove the mechanism" |
| warning for every `stp` contrast | stp lacks Stp1 **and** most Stk1 | same, repeated in the step 3–4 notebook entry |

## 2. What notebook 07 adds (refines, does not contradict)

**a) The two "wild-type-like" Stp1 values in stp are explained.** Our Step 3 entry already
flagged them as resting on almost nothing. The peptide positions now show what that "nothing" is:

| sample | log2 abundance | peptides behind it (residues) |
|---|---|---|
| `JE2stp_SN_2` | 18.4 | 51–57 **and 233–246** |
| `JE2stp_pulseSN_3` | 18.1 | 174–181 (one 8-residue peptide) |

Residues 233–246 lie **after** the insertion at codon 189. According to notebook 07, a
stp1::Tn strain cannot translate that part of the protein. It is the only C-terminal Stp1 detection
in all 9 stp samples, whereas anc has 24 and stk 23 such detections. This value is therefore most
likely a technical artefact (carry-over between runs, mis-assignment, or a sample problem) and
not real Stp1. **Consequence:** do not read these two samples as a "failed knock-out", and do not
let them drive a group mean for Stp1 in stp.

**b) Stp1 in stp looks truncated, not simply reduced.** N-terminal peptides (before residue 189)
still appear sporadically, while C-terminal ones are almost never seen. Our Step 3 call "supported
reduction" stays correct; notebook 07 makes it more specific. Its own caveat also applies to us:
the N-terminal detections are sporadic (1–3 samples per peptide). They fit "truncated and then
degraded" or "spill-over between runs" equally well.

## 3. Points to word carefully (no data conflict)

1. **Notebook 07 says "For Stk1, stk and stp are simply empty".** In context, this refers to the
   N-terminal part (before codon 46). Its own table and our Step 3–4 numbers both show **residual
   Stk1 in stp** (TSB 3/3 and SNpulse 2/3 samples, ~5 log2 below anc). Keep our wording: the
   knock-down in stp is strong but not complete. The "before/after codon 46" split means nothing
   for stp, because stk1 itself is intact there.
2. **Stk1 detections in stk that a truncated product cannot explain.** The single stk detection
   (`stk_pulseSN_1`, residues 332–340) and the phospho peptide `HIQDSVPNVTtDVR…` (T231, in 3/9 stk
   samples, Step 3) both lie **after** codon 46. An N-terminal fragment cannot produce them.
   Notebook 07 does not look at the phospho layer, so our open question ("carry-over? noise?")
   stays open. Its "spill-over between runs in DIA" hypothesis is one candidate explanation.
3. **Polar effect: proof vs. consistency.** Notebook 07 writes that the polar effect "is exactly
   what the proteome shows". We keep the more careful formulation: the data are *consistent* with
   polarity, but cannot separate it from Stk1 instability, feedback regulation or a secondary
   mutation. Notebook 07 itself asks for a wet-lab test in its discussion prompts.

## 4. Suggested next step

Look up the run order of `JE2stp_SN_2` (and of the stk samples with Stk1 detections). If an anc
or stk sample with high Stp1/Stk1 was injected just before, carry-over becomes the most likely
explanation.
