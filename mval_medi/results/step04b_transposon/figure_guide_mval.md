# Step 4b outputs: what each figure and table compares

Made by [`scripts/04b_transposon_products_mval.R`](../../../scripts/04b_transposon_products_mval.R).
Question: is the transposon in the `stk` and `stp` strains present and expressed,
and absent from `anc`? This is a positive control for the two insertions, not a
test of Stp1/Stk1 themselves.

**Data source (new, local only):** `data/raw/FP24_transposonDB/`, a FragPipe 24.0
search of the same 27 total-proteome runs by Jonas (12 Nov 2026), against
`p41881_db2_SA_JE2_wTransposon_20261111` = JE2 proteins + the proteins of the
transposon delivery plasmid pBursa (ermB, gfp, cat, repC, orf) and of the
helper plasmid pFA545 (tnp, bla, repC, tetB, xylR). It is NOT in the
checksummed course bundle (`data/checksums.md5`). The course DEA tables were
searched without these proteins, which is why ErmB was "never looked for" before.

The prediction was written in `reference/README.md` before these data came:
ErmB detected in the 18 mutant samples, absent from the 9 `anc` samples.

## Figures

### `fig_ermb_per_sample_mval.png`: Fig 1e, ErmB per sample
- **Compares:** ErmB in all 27 samples, 3 genotypes × 3 media × 3 replicates.
  The main comparison is mutant (`stk`, `stp`) against `anc`.
- **Data:** (a) spectral count per sample (`combined_protein.tsv`, count scale);
  (b) log2 of FragPipe MaxLFQ intensity. These intensities are NOT normalised
  like the course tables, so do not compare them with `read_abundances()` values.
- **How to read:** one point = one sample. In (a) `x` = 0 spectra (not detected).
  In (b) samples without signal are not drawn (NA, never 0 on a log scale), so
  the `anc` panel is empty.
- **Result:** 0 spectra in all 9 `anc` samples; 11–26 spectra in `stk`, 19–30 in
  `stp`. Prediction met.
- **Caveat:** spectral counts are a rough amount measure. `stk` has fewer ErmB
  spectra than `stp` in TSB and SNpulse (11–15 vs 25–30); this is described, not
  tested, and with n = 3 we do not interpret it.

### `fig_ermb_peptide_map_mval.png`: which parts of ErmB were seen
- **Compares:** the ErmB peptides (x = position in the 245-aa protein) per sample.
- **Data:** `combined_peptide.tsv`, per-sample spectral count > 0.
- **How to read:** one bar = one peptide seen in that sample; empty row = no ErmB
  peptide. 12 peptides, 166 of 245 residues covered (68 %).
- **Result:** no ErmB peptide in any `anc` sample; the same core peptides
  (e.g. positions 17–40, 77–111, 168–188) in every `stk` and `stp` sample.
- **Caveat:** `QVTSIELDSHLFNLSSEK` (53–70) is seen in 9/9 `stp` and 0/9 `stk`
  samples. The reason is not known (detection limit, as `stk` has fewer ErmB
  spectra overall, or something else); not interpreted.

## Tables

| file | content |
|---|---|
| `transposon_db_proteins_mval.csv` | all 10 transposon/plasmid proteins in the database: ProteinProphet probability, peptide count, whether they passed into `combined_protein.tsv` |
| `ermb_per_sample_mval.csv` | ErmB per sample: spectra, unique spectra, intensity, MaxLFQ, peptides, detected |
| `ermb_peptides_mval.csv` | each ErmB peptide: start/end, number of samples per genotype in which it was seen |
| `transposon_checks_mval.txt` | the PASS/FAIL checks of the script |

**Other transposon proteins:** only ErmB passed (probability 1.0, 384 peptide
spectra). GFP was in the database but got probability 0 and 0 peptides, so here
"not found" means "searched and not detected" (unlike in the course tables).
Whether *gfp* should be expressed from these insertions was not checked against
a source. Helper-plasmid proteins (pFA545) are not expected in the final strains;
two have weak, sub-threshold matches (tnp p = 0.56, xylR p = 0.37) and are not
counted as detected.

**What this does and does not show:** both mutant strains carry and express the
transposon marker, and the `anc` samples do not, so the sample labels and the
mutant identity agree. ErmB does not say *where* the transposon sits: the
insertion sites (stk1 codon 46, stp1 codon 189) come from the NTML map
(`TRANSPOSON_INSERTIONS` in `R/constants.R`).
