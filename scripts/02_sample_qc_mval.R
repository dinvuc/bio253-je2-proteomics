# Step 2 of our workflow (mval_medi/project_outlines/PROJECT_WORKFLOW_mval.md):
# audit sample identities and measurement quality.
#
# Why: mislabelled samples, a bad replicate or unequal detection between groups
# can look exactly like a genotype effect. We check this BEFORE any biology.
#
# Run from the repository root:   Rscript scripts/02_sample_qc_mval.R   (~1-2 min)
#
# Outputs (mval_medi/results/step02_qc/, all small summaries, no raw data):
#   qc_checks_mval.txt                pass/fail list of the design + mapping checks
#   protein_key_types_mval.csv        JE2 proteins vs contaminants vs decoys
#   detection_per_sample_mval.csv     detected features per sample, all layers
#   sample_qc_metrics_mval.csv        per-sample QC numbers (total proteome)
#   dea_table_metadata_mval.csv       model + estimate_type per DEA table
#   fig_*_mval.png                    QC plots
#
# Missing values: stay NA. Boxplots/detection counts use observed values only.
# Correlation and PCA need a complete matrix, so they use only proteins
# observed in ALL samples of the subset (complete-case). That is a bias towards
# abundant proteins, which is fine for "do samples look alike?" questions.

suppressMessages(source("R/setup.R"))
# setup.R loads dplyr/tidyr/readr/stringr/ggplot2 only; we also need map() and column_to_rownames()
suppressMessages({ library(purrr); library(tibble) })

out_dir <- "mval_medi/results/step02_qc"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
save_fig <- function(p, name, w = 8, h = 5) {
  ggsave(file.path(out_dir, name), p, width = w, height = h, dpi = 150)
}

checks <- character()
check <- function(label, ok, detail = "") {
  line <- sprintf("[%s] %s %s", if (isTRUE(ok)) "PASS" else "FLAG", label, detail)
  checks <<- c(checks, line)
  message(line)
}

# ---- 1. Sample sheet: design ----------------------------------------------------
samples <- read_samples()
design  <- count(samples, genotype, medium)
primary <- filter(samples, medium %in% c("TSB", "SN"))   # our 18-sample subset

check("27 unique sample names", n_distinct(samples$sample) == 27 && nrow(samples) == 27,
      paste0("(n = ", n_distinct(samples$sample), ")"))
check("3 replicates in each of the 9 genotype x medium groups",
      nrow(design) == 9 && all(design$n == 3))
check("primary subset (TSB + SN) has 18 samples", nrow(primary) == 18)
check("replicates numbered 1-3 within each group",
      all(count(samples, group, replicate)$n == 1) && all(samples$replicate %in% 1:3))

# Name repair: truncated Excel names like "2anc_pulseSN_1" must map to JE2anc_pulseSN_* (medium column = SNpulse).
repaired <- samples |> filter(sample != sample_name_in_excel)
check("no truncated sample names left in the 'sample' column",
      !any(str_detect(samples$sample, "^[0-9]")))
check("each repaired name matches its raw file (genotype, medium, replicate)",
      all(str_detect(repaired$raw_file_total,
                     paste0("JE2", repaired$genotype, "_pulseSN_", repaired$replicate))),
      paste0("(", nrow(repaired), " names repaired)"))
check("raw file name agrees with the sample sheet for every sample",
      all(str_detect(samples$raw_file_total,
                     paste0("JE2", samples$genotype, "_",
                            if_else(samples$medium == "SNpulse", "pulseSN", as.character(samples$medium)),
                            "_", samples$replicate, "_"))))

# ---- 2. Load all layers ----------------------------------------------------------
layers <- list(
  total        = read_abundances("total",        drop_contaminants = FALSE),
  phospho_STYH = read_abundances("phospho_STYH", drop_contaminants = FALSE),
  phospho_STY  = read_abundances("phospho_STY",  drop_contaminants = FALSE),
  sites_STYH   = read_sites("STYH", drop_contaminants = FALSE) |> rename(normalized_abundance = log2_intensity),
  sites_STY    = read_sites("STY",  drop_contaminants = FALSE) |> rename(normalized_abundance = log2_intensity)
)

# ---- 3. Same samples, same labels in every layer? -------------------------------
for (nm in names(layers)) {
  lab <- layers[[nm]] |> distinct(sample, genotype, medium, replicate)
  joined <- inner_join(lab, samples, by = "sample", suffix = c("", ".sheet"))
  check(sprintf("%-13s has exactly the 27 sample-sheet names", nm),
        setequal(lab$sample, samples$sample))
  check(sprintf("%-13s genotype/medium/replicate agree with sample sheet", nm),
        nrow(lab) == 27 &&
          all(joined$genotype == joined$genotype.sheet, joined$medium == joined$medium.sheet,
              joined$replicate == joined$replicate.sheet))
}

reinj <- layers$phospho_STYH |> distinct(sample, reinjected) |> filter(reinjected)
check("reinjected phospho samples are exactly stk pulseSN replicates 1 and 2",
      setequal(reinj$sample, c("JE2stk_pulseSN_1", "JE2stk_pulseSN_2")),
      paste0("(", paste(reinj$sample, collapse = ", "), ")"))

# ---- 4. Identifier uniqueness and what kind of proteins are in the table ---------
tot_all <- layers$total
check("total: one row per protein_key x sample",
      nrow(distinct(tot_all, protein_key, sample)) == nrow(tot_all))

key_types <- tot_all |>
  distinct(protein_key, locus_tag, is_contaminant) |>
  mutate(type = case_when(str_detect(protein_key, "^rev_")                 ~ "decoy (rev_)",
                          str_detect(protein_key, "contam|^zz|CON_")       ~ "contaminant",
                          !is.na(locus_tag)                                ~ "JE2 protein",
                          TRUE                                             ~ "other")) |>
  count(type, flagged_is_contaminant = is_contaminant, name = "n_protein_keys")
write_csv(key_types, file.path(out_dir, "protein_key_types_mval.csv"))
print(key_types)
n_unflagged <- key_types |> filter(type != "JE2 protein", !flagged_is_contaminant) |> pull(n_protein_keys) |> sum()
check("all contaminants/decoys flagged by is_contaminant", n_unflagged == 0,
      paste0("(", n_unflagged, " non-JE2 keys pass the default filter -> we filter them ourselves)"))

# From here on: JE2 proteins only (has a locus tag, not a reversed decoy).
tot <- keep_je2(tot_all)   # shared filter, R/helpers_mval.R
check("total: locus_tag unique per protein_key",
      nrow(distinct(tot, protein_key, locus_tag)) == n_distinct(tot$protein_key))
message("JE2 proteins in the total proteome: ", n_distinct(tot$protein_key))

# ---- 5. Detection / missingness per sample and layer ----------------------------
detection <- imap(layers, function(x, nm) {
  if (nm == "total") x <- tot          # JE2 proteins only
  x |> filter(!is_contaminant) |>
    group_by(sample, genotype, medium, replicate, reinjected) |>
    summarise(n_detected = sum(!is.na(normalized_abundance)),
              n_features = n(), .groups = "drop") |>
    mutate(layer = nm, pct_detected = round(100 * n_detected / n_features, 1))
}) |> bind_rows()
write_csv(detection, file.path(out_dir, "detection_per_sample_mval.csv"))

det_group <- detection |>
  group_by(layer, genotype, medium) |>
  summarise(mean_detected = round(mean(n_detected)), min = min(n_detected),
            max = max(n_detected), .groups = "drop")
print(det_group |> filter(layer == "total"))

p_det <- ggplot(detection, aes(interaction(replicate, medium, sep = "_"), n_detected,
                               fill = genotype, alpha = reinjected)) +
  geom_col(position = position_dodge(width = 0.9)) +
  facet_wrap(~layer, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = GENOTYPE_COLOURS) +
  scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.4)) +
  labs(x = "replicate_medium", y = "detected features (non-NA)",
       title = "Detection per sample", subtitle = "faded bars = reinjected phospho samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p_det, "fig_detection_per_sample_mval.png", h = 10)

# ---- 6. Abundance distributions (total) ------------------------------------------
p_box <- ggplot(tot, aes(sample, normalized_abundance, fill = genotype)) +
  geom_boxplot(outlier.size = 0.2, na.rm = TRUE) +
  facet_grid(~medium, scales = "free_x") +
  scale_fill_manual(values = GENOTYPE_COLOURS) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
  labs(x = NULL, y = "normalised abundance (log2)", title = "Total proteome: distribution per sample")
save_fig(p_box, "fig_abundance_boxplot_mval.png")

# Missingness is not random: is it concentrated in low-abundance proteins?
p_mis <- tot |>
  group_by(protein_key) |>
  summarise(mean_obs = mean(normalized_abundance, na.rm = TRUE),
            n_missing = sum(is.na(normalized_abundance))) |>
  mutate(missing_class = cut(n_missing, c(-1, 0, 3, 12, 27), labels = c("0", "1-3", "4-12", "13-27"))) |>
  ggplot(aes(missing_class, mean_obs)) + geom_boxplot(na.rm = TRUE) +
  labs(x = "samples (of 27) where protein is missing", y = "mean observed abundance (log2)",
       title = "Missing values are mostly low-abundance proteins?")
save_fig(p_mis, "fig_missing_vs_abundance_mval.png", w = 6, h = 4)

# ---- 7. Correlation + PCA (complete-case matrix) --------------------------------
complete_matrix <- function(x, keep_samples) {
  x |> filter(sample %in% keep_samples) |>
    select(protein_key, sample, normalized_abundance) |>
    pivot_wider(names_from = sample, values_from = normalized_abundance) |>
    drop_na() |> column_to_rownames("protein_key") |> as.matrix()
}

m27 <- complete_matrix(tot, samples$sample)
m18 <- complete_matrix(tot, primary$sample)
message("complete-case proteins: 27 samples = ", nrow(m27), ", 18 samples = ", nrow(m18))

cor27 <- cor(m27)
# Outlier screen: each sample's mean correlation with its own replicates
# vs. with all other samples. A good replicate is closest to its own group.
sample_metrics <- samples |>
  select(sample, group, genotype, medium, replicate) |>
  mutate(r_within_group = map2_dbl(sample, group, \(s, g) {
           mates <- setdiff(samples$sample[samples$group == g], s); mean(cor27[s, mates]) }),
         r_other_groups = map2_dbl(sample, group, \(s, g) {
           mean(cor27[s, samples$sample[samples$group != g]]) }),
         median_abundance = map_dbl(sample, \(s) median(tot$normalized_abundance[tot$sample == s], na.rm = TRUE))) |>
  left_join(detection |> filter(layer == "total") |> select(sample, n_detected), by = "sample") |>
  mutate(detected_z = round((n_detected - mean(n_detected)) / sd(n_detected), 2),
         # flag = sample is less similar to its own replicates than the typical sample
         low_within_r = r_within_group < median(r_within_group) - 3 * mad(r_within_group))
write_csv(sample_metrics |> mutate(across(where(is.double), \(v) round(v, 4))),
          file.path(out_dir, "sample_qc_metrics_mval.csv"))
check("no sample with unusually low correlation to its own replicates",
      !any(sample_metrics$low_within_r),
      paste(sample_metrics$sample[sample_metrics$low_within_r], collapse = ", "))
check("no sample with |detection z| > 2.5",
      !any(abs(sample_metrics$detected_z) > 2.5),
      paste(sample_metrics$sample[abs(sample_metrics$detected_z) > 2.5], collapse = ", "))

p_cor <- as_tibble(cor27, rownames = "a") |>
  pivot_longer(-a, names_to = "b", values_to = "r") |>
  mutate(a = factor(a, samples$sample), b = factor(b, samples$sample)) |>
  ggplot(aes(a, b, fill = r)) + geom_tile() + scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6), axis.text.y = element_text(size = 6)) +
  labs(x = NULL, y = NULL, title = sprintf("Pearson correlation, %d complete proteins", nrow(m27)))
save_fig(p_cor, "fig_sample_correlation_mval.png", w = 8, h = 7)

pca_plot <- function(m, title) {
  pca <- prcomp(t(m), center = TRUE, scale. = FALSE)
  ve  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  sc  <- as_tibble(pca$x[, 1:3], rownames = "sample") |> left_join(samples, by = "sample")
  p <- ggplot(sc, aes(PC1, PC2, colour = genotype, shape = medium)) +
    geom_point(size = 3) +
    geom_text(aes(label = replicate), nudge_x = 0.6, size = 3, show.legend = FALSE) +
    scale_colour_manual(values = GENOTYPE_COLOURS) +
    labs(x = paste0("PC1 (", ve[1], " %)"), y = paste0("PC2 (", ve[2], " %)"),
         title = title, subtitle = sprintf("%d proteins observed in all samples; number = replicate", nrow(m)))
  list(plot = p, var = ve[1:4], scores = sc)
}
pca27 <- pca_plot(m27, "PCA total proteome, all 27 samples")
pca18 <- pca_plot(m18, "PCA total proteome, primary 18 samples (TSB + SN)")
save_fig(pca27$plot, "fig_pca_27_mval.png", w = 7, h = 5)
save_fig(pca18$plot, "fig_pca_18_mval.png", w = 7, h = 5)

# Which factor drives PC1/PC2? R^2 of a one-factor linear model on each PC score.
pc_drivers <- expand_grid(pc = c("PC1", "PC2", "PC3"), factor = c("genotype", "medium")) |>
  mutate(r2_27 = map2_dbl(pc, factor, \(p, f) summary(lm(pca27$scores[[p]] ~ pca27$scores[[f]]))$r.squared),
         r2_18 = map2_dbl(pc, factor, \(p, f) summary(lm(pca18$scores[[p]] ~ pca18$scores[[f]]))$r.squared)) |>
  mutate(across(starts_with("r2"), \(v) round(v, 2)))
print(pc_drivers)

# ---- 8. Reinjected phospho samples: do they stand out? ----------------------------
reinj_check <- detection |>
  filter(layer %in% c("phospho_STYH", "phospho_STY"), genotype == "stk", medium == "SNpulse") |>
  select(layer, sample, reinjected, n_detected)
print(reinj_check)

# ---- 9. DEA tables: which model, how much is imputed ----------------------------
dea_meta <- map(c("total", "phospho_STYH", "phospho_STY"), \(d) {
  read_dea(d, drop_contaminants = FALSE) |>
    count(modelName, estimate_type, name = "n_rows") |>
    mutate(table = d, .before = 1,
           n_contrasts = n_distinct(read_dea(d, drop_contaminants = FALSE)$contrast))
}) |> bind_rows()
write_csv(dea_meta, file.path(out_dir, "dea_table_metadata_mval.csv"))
print(dea_meta)
dea_cols <- names(read_dea("total"))
check("DEA tables carry uncertainty columns (std.error, conf.low/high, df, FDR)",
      all(c("std.error", "conf.low", "conf.high", "df", "p.value", "FDR") %in% dea_cols))

# ---- 10. Write the check list ----------------------------------------------------
writeLines(c("# Step 2 QC checks (group mval)", paste("Run at:", format(Sys.time(), "%Y-%m-%d %H:%M")), "",
             checks, "",
             sprintf("JE2 proteins (total): %d; complete in all 27: %d; complete in primary 18: %d",
                     n_distinct(tot$protein_key), nrow(m27), nrow(m18)),
             sprintf("PCA 27: PC1 %.1f %%, PC2 %.1f %%; PCA 18: PC1 %.1f %%, PC2 %.1f %%",
                     pca27$var[1], pca27$var[2], pca18$var[1], pca18$var[2]), "",
             "PC drivers (R^2 of PC score ~ factor):",
             capture.output(print(as.data.frame(pc_drivers), row.names = FALSE))),
           file.path(out_dir, "qc_checks_mval.txt"))
message("Done. Results in ", out_dir)
