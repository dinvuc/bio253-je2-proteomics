# Step 2 outputs: sample identity and measurement quality

Made by [`C_mutant_validation_DinMete/scripts/02_sample_qc_mval.R`](../../scripts/02_sample_qc_mval.R) (11 Nov). Question:
are the 27 samples what the sample sheet says, are they of similar quality, and does any
sample or technical factor look like it could fake a genotype effect? Data: normalised
log2 protein abundances (`read_abundances("total")`, 2,232 JE2 proteins), phospho tables
for detection counts. All checks: [`qc_checks_mval.txt`](qc_checks_mval.txt).

**Step conclusion:** all 27 samples are correctly labelled and comparable. Medium (PC1, 62 %)
and genotype (PC2/PC3) are by far the strongest signals. Replicates agree closely
(r 0.964–0.991). No sample was excluded.

---

### `fig_pca_27_mval.png`: PCA of all 27 samples
- **Description:** principal component analysis of the 1,742 proteins observed in all 27
  samples (no imputation needed); colour = genotype, shape = medium, number = replicate.
- **How to read:** samples close together have similar proteomes. The axis labels give
  the share of the total variation each component explains.
- **Key points:**
  - PC1 (62.0 %) separates SN from TSB/SNpulse. Medium explains 98 % of PC1.
  - PC2 (10.7 %) separates `stp` from `anc`/`stk`. Genotype explains 88 % of PC2, and 97 % of PC3.
  - Replicates of a group lie together; SNpulse lies next to TSB (15 min of PASN changes
    little), and `anc` and `stk` overlap more than either does with `stp`.
- **Conclusion:** the experiment is dominated by the medium, then by genotype, with
  tight replicates. `stp` differs from `anc` more than `stk` does. There is no sign of a
  swapped or outlying sample.

### `fig_pca_18_mval.png`: PCA of the 18 primary samples (TSB, SN)
- **Description:** the same PCA restricted to the samples of the primary analysis
  (1,758 complete proteins).
- **How to read:** as above.
- **Key points:** PC1 66.4 % (medium, R² = 0.99), PC2 10.1 % (genotype, R² = 0.90). The
  structure is the same as with 27 samples.
- **Conclusion:** dropping SNpulse does not change the picture. The 18-sample design
  used for H1–H3 is well behaved.

### `fig_sample_correlation_mval.png`: sample-to-sample correlation
- **Description:** Pearson correlation of every pair of samples over the 1,742 complete
  proteins, in group order (`anc`, `stk`, `stp` × TSB, SNpulse, SN).
- **How to read:** yellow = very similar. Bright 3 × 3 blocks on the diagonal are the
  replicates of one group. Large blocks show which groups resemble each other.
- **Key points:**
  - correlation within groups 0.964–0.991, between groups 0.817–0.896 (`sample_qc_metrics_mval.csv`)
  - SN samples of all strains form one block, TSB + SNpulse another
  - the lowest within-group correlation is `JE2stk_TSB_3` (0.964)
- **Conclusion:** every sample correlates best with its own replicates, so no sample is
  mislabelled. The medium block structure confirms the PCA.

### `fig_detection_per_sample_mval.png`: how many features each sample detects
- **Description:** number of features with a value per sample in each layer: total
  proteome, phospho abundances and site tables of both searches (`detection_per_sample_mval.csv`).
- **How to read:** one point = one sample. Compare samples within a layer, not across layers.
- **Key points:**
  - total proteome: 1,981–2,147 of 2,232 proteins per sample (89–96 %). The lowest is
    `JE2stk_TSB_3` (1,981, z = −2.08), within the |z| ≤ 2.5 rule.
  - phospho layers vary far more: e.g. STY sites 291–933 (TSB), 257–628 (SNpulse) and
    458–1,054 (SN) per sample. The two reinjected phospho samples (`stk_pulseSN_1/2`, run on
    another day) are among the SNpulse values; `qc_checks_mval.txt` confirms they are exactly these two.
- **Conclusion:** the total proteome has even depth, so missing values there are mostly
  low-abundance proteins. The phosphoproteome depth differs strongly between samples and
  media, which Step 8 has to take into account.
  - `stk_TSB_3` is the weakest sample but not an outlier. Step 11 later found that dropping
    `stk_TSB_2` or `stk_TSB_3` costs the most hits.

### `fig_abundance_boxplot_mval.png`: abundance distributions per sample
- **Description:** distribution of normalised log2 abundances per sample (boxplots).
- **How to read:** boxes at the same height = samples are normalised to a common level.
- **Key points:** sample medians are 16.02–16.10 (a spread of 0.08 log2), and no sample is shifted.
- **Conclusion:** the supplied normalisation worked. The values can be compared directly
  and must not be normalised again (coding rules, section 5).

### `fig_missing_vs_abundance_mval.png`: which proteins are missing?
- **Description:** mean observed abundance of each protein, grouped by the number of
  samples (of 27) in which it is missing (0, 1–3, 4–12, 13–27).
- **How to read:** if missing values were random, all boxes would be at the same height.
- **Key points:** proteins with no missing value have a median of ~16.3 log2. Proteins
  missing in 1–3 / 4–12 / 13–27 samples sit lower, at ~14.6 / 14.5 / 14.0.
- **Conclusion:** missingness depends on abundance (missing not at random): a missing
  value usually means "below detection", not "absent at random". This is why we never
  enter 0, treat on/off proteins separately, and test only proteins with ≥ 2 of 3 values
  per group (Step 5).

## Tables
| file | content |
|---|---|
| `sample_qc_metrics_mval.csv` | per sample: within- and between-group correlation, median abundance, proteins detected, detection z-score |
| `detection_per_sample_mval.csv` | per sample × layer: features detected |
| `dea_table_metadata_mval.csv` | supplied DEA tables: model name, estimate types, uncertainty columns |
| `protein_key_types_mval.csv` | protein keys by type (JE2 locus tag, contaminant, decoy): basis of `keep_je2()` |
| `qc_checks_mval.txt` | all PASS/FLAG checks and the PCA drivers |
