# Step 5 outputs: what each figure and table compares

Made by [`scripts/05_formal_test_mval.R`](../../../scripts/05_formal_test_mval.R) (part 1,
the formal test, 12 Nov) and [`scripts/05b_contrast_sheet_mval.R`](../../../scripts/05b_contrast_sheet_mval.R)
(part 2, definitions and checks, 13 Nov). The rules behind all of them are in
[`methods_sheet_mval.md`](methods_sheet_mval.md).

## Figures

### `fig_volcano_primary_mval.png`: formal test, 9 primary contrasts
- **Compares:** the 9 TSB/SN contrasts (baseline, exposure within genotype, interactions).
- **Data:** our limma `treat(lfc = 1)` refit; x = diff = log2(A) − log2(B), y = −log10(treat BH-FDR).
- **How to read:** coloured = formal hit (orange higher in A, blue higher in B);
  dashed lines at |log2FC| = 1 and FDR = 0.05, the same thresholds as the test.
  Only eligible (tested) proteins are drawn.
- **Caveat:** treat FDR is not the ordinary FDR; a point beyond the lines can still be grey.

### `fig_raw_top_hits_mval.png`: raw values behind the top `stk_vs_anc_at_TSB` hits
- **Compares:** per-sample normalised log2 abundance of the 6 strongest hits, all 27 samples.
- **How to read:** one point = one sample; missing samples are not drawn.

### `fig_eligibility_mval.png`: what each primary contrast can test
- **Compares:** for the 9 primary contrasts, how many of the 2,232 JE2 proteins are
  eligible (≥ 2/3 observed in every group used), on/off (≥ 2/3 in one group, 0/3 in
  another) or sparse.
- **How to read:** the number is the eligible count. Interactions use 4 groups, so
  more proteins drop out of them.
- **Caveat:** on/off proteins are presence/absence evidence, not a measured fold change.

## Tables

| file | content |
|---|---|
| `hit_counts_mval.csv` | per contrast: tested, formal hits up/down, supplied screen hits (part 1) |
| `hits_formal_mval.csv` | every formal hit: diff, SE, treat t, FDR (part 1) |
| `sign_check_mval.csv` | our diff vs supplied diff per contrast (part 1) |
| `contrast_sheet_mval.csv` | all 33 contrasts: family, formula, meaning in words, eligible / on-off / sparse counts, eligible rows that carry the supplied `lod_imputed` label, formal and screen hits |
| `direction_check_mval.csv` | per contrast: max deviation between limma diff and the hand-computed group-mean formula |
| `one_sided_signals_mval.csv` | the 783 on/off rows of the 9 primary contrasts: observation pattern (e.g. `anc_TSB=0/3 anc_SN=3/3`), supplied diff/FDR, whether in the supplied screen, which group is higher (simple contrasts only) |
| `formal_test_checks_mval.txt`, `contrast_sheet_checks_mval.txt` | PASS/FAIL checks |
