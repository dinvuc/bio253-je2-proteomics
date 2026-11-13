# Step 7 (PROJECT_WORKFLOW_mval.md): shared and distinct mutant effects (H3).
#
# Run from the repository root:   Rscript scripts/07_shared_distinct_mval.R   (~40 s)
# Needs scripts/05_formal_test_mval.R first.
#
# ---- Plan (coding rules, section 4) -----------------------------------------
# Aim:      Does `stp` look like `stk`? If the stp1::Tn insertion has "the same
#           effect as the stk mutation" (polar loss of Stk1, Steps 3-4), the
#           stp-vs-anc effects should lie near the diagonal of the stk-vs-anc
#           effects, both in TSB (baseline) and in the PASN response (interaction).
#           Proteins far off the diagonal = candidate phosphatase-specific effects.
# Input:    data/derived/limma_treat_total_mval.tsv.gz (Step 5 formal test)
#           read_abundances("total") |> keep_je2(), read_samples() (split-reference check)
#           read_protein_metadata() gene symbols (ortholog_quality == "good")
# Analysis: no new hypothesis test, so no new FDR. Descriptive:
#           * Pearson r with 95% CI (cor.test) and total-least-squares slope
#             (errors in both axes) with a bootstrap 95% CI over proteins.
#           * Shared-reference problem: both axes subtract the SAME anc mean, so
#             pure noise already gives r ~ 0.5 (equal variances). Sensitivity:
#             split-reference r, where stk is compared with one anc replicate and
#             stp with the other two (independent noise, null r = 0, but noisier
#             -> lower). Done for all 3 ways to pick the single replicate.
#           * Classes from the Step 5 formal hits (treat lfc = 1, BH-FDR < 0.05):
#             hit in both same sign / both opposite / stk only / stp only / none.
#             "Different between mutants" ONLY from the direct stp_vs_stk contrasts.
# Output:   mval_medi/results/step07_shared_distinct/
#           data/derived/effect_pairs_step07_mval.tsv.gz  (local) per protein: stk/stp/direct estimates, SE, FDR, class
#             similarity_stats_mval.csv      r, CI, TLS slope, null r, split-reference r
#             class_counts_mval.csv          same / opposite / one-only counts
#             direct_contrast_candidates_mval.csv  stp_vs_stk hits (TSB and interaction)
#             fig_mutant_effects_scatter_mval.png   Fig 3a: stp vs stk effects
#             fig_heatmap_mval.png                  Fig 3b: centred group means of candidates
#             fig_shared_distinct_raw_mval.png      per-sample values, shared vs distinct
#             figure_guide_mval.md, shared_distinct_checks_mval.txt
# Check:    samples joined by name; stp_vs_stk diff == stp_vs_anc - stk_vs_anc
#           (same model); class counts add up to the eligible-in-both total.

suppressMessages(source("R/setup.R"))
out <- "mval_medi/results/step07_shared_distinct"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(20261113)
checks <- character()
check <- function(label, ok) {
  line <- sprintf("[%s] %s", if (isTRUE(ok)) "PASS" else "FAIL", label)
  message(line); checks <<- c(checks, line)
  invisible(ok)
}

res  <- read_tsv("data/derived/limma_treat_total_mval.tsv.gz", show_col_types = FALSE)
meta <- read_protein_metadata() |> filter(ortholog_quality == "good") |> distinct(locus_tag, gene = gene_symbol_best)

LAYERS <- tibble(
  layer  = c("baseline (TSB)", "response (TSB vs SN interaction)"),
  stk    = c("stk_vs_anc_at_TSB", "interaction_stk_vs_anc_at_TSB_vs_SN"),
  stp    = c("stp_vs_anc_at_TSB", "interaction_stp_vs_anc_at_TSB_vs_SN"),
  direct = c("stp_vs_stk_at_TSB", "interaction_stp_vs_stk_at_TSB_vs_SN"))

# ---- 1. one row per protein x layer with both mutant effects and the direct contrast
get <- function(ct, prefix) {
  res |> filter(contrast == ct) |>
    select(protein_key, diff, std.error, FDR, tested) |>
    rename_with(\(x) paste0(prefix, "_", x), -protein_key)
}
pairs <- purrr::pmap(LAYERS, \(layer, stk, stp, direct) {
  get(stk, "stk") |> full_join(get(stp, "stp"), by = "protein_key") |>
    full_join(get(direct, "direct"), by = "protein_key") |> mutate(layer = layer, .before = 1)
}) |> bind_rows() |>
  left_join(distinct(res, protein_key, locus_tag, product), by = "protein_key") |>
  left_join(meta, by = "locus_tag")

check("same model: direct stp_vs_stk diff == stp_vs_anc - stk_vs_anc (all rows with both)",
      with(pairs, max(abs(direct_diff - (stp_diff - stk_diff)), na.rm = TRUE)) < 1e-8)

pairs <- pairs |>
  mutate(eligible_both = stk_tested & stp_tested,
         eligibility = case_when(eligible_both ~ "eligible in both",
                                 stk_tested ~ "eligible in stk only", stp_tested ~ "eligible in stp only",
                                 TRUE ~ "eligible in neither"),
         stk_hit = stk_tested & stk_FDR < MVAL_FDR, stp_hit = stp_tested & stp_FDR < MVAL_FDR,
         direct_hit = direct_tested & direct_FDR < MVAL_FDR,
         class = case_when(!eligible_both ~ NA_character_,
                           stk_hit & stp_hit & sign(stk_diff) == sign(stp_diff) ~ "hit in both, same direction",
                           stk_hit & stp_hit ~ "hit in both, opposite direction",
                           stk_hit ~ "hit in stk only", stp_hit ~ "hit in stp only",
                           TRUE ~ "no hit"))
write_tsv(pairs, "data/derived/effect_pairs_step07_mval.tsv.gz")   # full table: local only (git-ignored), like Step 5

class_counts <- pairs |> filter(eligible_both) |>
  count(layer, class) |>
  bind_rows(pairs |> count(layer, class = eligibility) |> filter(class != "eligible in both")) |>
  bind_rows(pairs |> filter(direct_hit) |> count(layer) |> mutate(class = "direct stp_vs_stk hit")) |>
  arrange(layer)
write_csv(class_counts, file.path(out, "class_counts_mval.csv"))
print(class_counts, n = 30)
n_both <- pairs |> filter(eligible_both) |> count(layer)
check("class counts add up to the eligible-in-both totals",
      all(class_counts |> filter(str_detect(class, "^hit|^no hit")) |> group_by(layer) |>
            summarise(n = sum(n)) |> pull(n) == n_both$n))

# ---- 2. similarity: r, TLS slope, and the shared-reference problem --------------
tls_slope <- function(x, y) {      # first principal axis of the centred (x, y) cloud
  v <- eigen(cov(cbind(x, y)))$vectors[, 1]
  v[2] / v[1]
}
sim_full <- pairs |> filter(eligible_both) |> group_by(layer) |>
  summarise(n = n(),
            r = cor(stk_diff, stp_diff),
            r_low = cor.test(stk_diff, stp_diff)$conf.int[1], r_high = cor.test(stk_diff, stp_diff)$conf.int[2],
            tls_slope = tls_slope(stk_diff, stp_diff),
            slope_boot = list(replicate(1000, { i <- sample(n(), replace = TRUE); tls_slope(stk_diff[i], stp_diff[i]) })),
            .groups = "drop") |>
  mutate(slope_low = purrr::map_dbl(slope_boot, quantile, 0.025),
         slope_high = purrr::map_dbl(slope_boot, quantile, 0.975), .keep = "unused") |>
  select(-any_of("slope_boot"))

# Null r from the shared anc reference, using the actual per-protein variances:
# with equal residual variance s2 in all groups and n = 3 per group,
#   baseline:    cov = s2/3 (anc_TSB),            var = 2 s2/3  -> r0 = 0.5
#   interaction: cov = 2 s2/3 (anc_TSB, anc_SN),  var = 4 s2/3  -> r0 = 0.5
# Missing values change this slightly; 0.5 is the reference value.
sim_full$r_expected_from_shared_reference_only <- 0.5

# Split-reference sensitivity: independent anc replicates for the two mutants.
samples <- read_samples()
abund <- read_abundances("total") |> keep_je2() |>
  select(protein_key, sample, normalized_abundance) |>
  inner_join(select(samples, sample, group, replicate), by = "sample")
check("27 samples, 3 per group after the join (split-reference check)",
      n_distinct(abund$sample) == 27 && all(table(distinct(abund, sample, group)$group) == 3))
gm <- abund |> group_by(protein_key, group) |>
  summarise(m = mean(normalized_abundance, na.rm = TRUE), .groups = "drop") |>
  mutate(m = if_else(is.nan(m), NA_real_, m)) |>
  pivot_wider(names_from = group, values_from = m)
anc_split <- function(med, single) {
  a <- abund |> filter(group == paste0("anc_", med)) |>
    mutate(part = if_else(replicate == single, "one", "two")) |>
    group_by(protein_key, part) |>
    summarise(m = mean(normalized_abundance, na.rm = TRUE), .groups = "drop") |>
    mutate(m = if_else(is.nan(m), NA_real_, m)) |>
    pivot_wider(names_from = part, values_from = m, names_prefix = paste0("anc_", med, "_"))
}
split_r <- purrr::map(1:3, \(single) {
  d <- gm |> left_join(anc_split("TSB", single), by = "protein_key") |>
    left_join(anc_split("SN", single), by = "protein_key") |>
    transmute(protein_key,
              base_stk = stk_TSB - anc_TSB_one, base_stp = stp_TSB - anc_TSB_two,
              int_stk = (stk_TSB - stk_SN) - (anc_TSB_one - anc_SN_one),
              int_stp = (stp_TSB - stp_SN) - (anc_TSB_two - anc_SN_two))
  elig <- pairs |> filter(eligible_both) |> select(layer, protein_key)
  bind_rows(
    d |> semi_join(filter(elig, layer == LAYERS$layer[1]), by = "protein_key") |>
      summarise(layer = LAYERS$layer[1], r = cor(base_stk, base_stp, use = "complete.obs")),
    d |> semi_join(filter(elig, layer == LAYERS$layer[2]), by = "protein_key") |>
      summarise(layer = LAYERS$layer[2], r = cor(int_stk, int_stp, use = "complete.obs"))) |>
    mutate(single_anc_replicate = single)
}) |> bind_rows()
sim <- sim_full |>
  left_join(split_r |> group_by(layer) |>
              summarise(r_split_reference_min = min(r), r_split_reference_max = max(r)), by = "layer")
write_csv(sim, file.path(out, "similarity_stats_mval.csv"))
print(sim, width = 200)

# ---- 3. direct stp_vs_stk candidates ----------------------------------------------
direct <- pairs |> filter(direct_hit) |>
  mutate(in_words = case_when(
    layer == LAYERS$layer[1] & direct_diff > 0 ~ "higher in stp than in stk (TSB)",
    layer == LAYERS$layer[1] ~ "higher in stk than in stp (TSB)",
    direct_diff > 0 ~ "SN response lower in stp than in stk",
    TRUE ~ "SN response higher in stp than in stk")) |>
  select(layer, locus_tag, gene, product, direct_diff, direct_std.error, direct_FDR,
         stk_diff, stk_FDR, stp_diff, stp_FDR, class, in_words) |>
  arrange(layer, direct_FDR)
write_csv(direct, file.path(out, "direct_contrast_candidates_mval.csv"))
print(direct, n = 40, width = 200)

# ---- 4. figures -------------------------------------------------------------------
sc <- pairs |> filter(eligible_both) |>
  mutate(show = case_when(direct_hit ~ "differs between mutants (direct stp_vs_stk hit)",
                          class == "hit in both, same direction" ~ "hit in both, same direction",
                          class %in% c("hit in stk only", "hit in stp only") ~ "hit in one mutant only",
                          class == "hit in both, opposite direction" ~ "hit in both, opposite direction",
                          TRUE ~ "no formal hit"),
         show = factor(show, c("no formal hit", "hit in one mutant only", "hit in both, same direction",
                               "hit in both, opposite direction", "differs between mutants (direct stp_vs_stk hit)")),
         label = if_else(direct_hit | class == "hit in both, same direction", coalesce(gene, locus_tag), NA_character_))
txt <- sim |> mutate(text = sprintf(
  "n = %d eligible in both\nr = %.2f (95%% CI %.2f-%.2f)\nnoise alone gives r ~ 0.5 (shared anc)\nsplit-reference r = %.2f-%.2f\nTLS slope = %.2f (%.2f-%.2f)",
  n, r, r_low, r_high, r_split_reference_min, r_split_reference_max, tls_slope, slope_low, slope_high))
p_sc <- ggplot(sc, aes(stk_diff, stp_diff)) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_abline(slope = 1, intercept = 0, colour = "grey40") +
  geom_point(data = filter(sc, show == "no formal hit"), colour = "grey60", alpha = 0.3, size = 0.7) +
  geom_point(data = filter(sc, show != "no formal hit"), aes(colour = show), size = 1.8) +
  geom_text(data = filter(sc, !is.na(label)), aes(label = label, colour = show), size = 2.4, vjust = -0.7,
            show.legend = FALSE, check_overlap = TRUE) +
  geom_text(data = txt, aes(x = -Inf, y = Inf, label = text), hjust = -0.05, vjust = 1.1, size = 2.6) +
  scale_colour_manual(values = c("hit in one mutant only" = "#56B4E9", "hit in both, same direction" = "#009E73",
                                 "hit in both, opposite direction" = "#CC79A7",
                                 "differs between mutants (direct stp_vs_stk hit)" = "#D55E00")) +
  facet_wrap(~layer, scales = "free") +
  labs(x = "stk effect vs anc (log2):  stk_vs_anc_at_TSB  /  interaction_stk_vs_anc_at_TSB_vs_SN",
       y = "stp effect vs anc (log2):  stp_vs_anc / interaction_stp_vs_anc",
       colour = NULL,
       title = "Does stp look like stk? Mutant effects relative to the ancestor",
       subtitle = paste0("Diagonal = same effect in both mutants. Left: > 0 = higher in the mutant than in anc (TSB). ",
                         "Right: > 0 = SN response lower in the mutant than in anc.\n",
                         "Hits: treat lfc = 1, BH-FDR < 0.05 (Step 5). stp also lacks most Stk1 (polar effect, Steps 3-4).")) +
  theme(legend.position = "bottom", plot.subtitle = element_text(size = 7.5)) +
  guides(colour = guide_legend(nrow = 2))
ggsave(file.path(out, "fig_mutant_effects_scatter_mval.png"), p_sc, width = 11, height = 6.5, dpi = 300)

# heatmap of all candidates: proteins that are a formal hit in any of the 6
# genotype contrasts (3 baseline + 3 interaction), observed group means of the
# 6 TSB/SN groups, centred per protein for DISPLAY only (no scaling).
cand_keys <- pairs |> filter(stk_hit | stp_hit | direct_hit) |> distinct(protein_key) |> pull()
hm <- gm |> filter(protein_key %in% cand_keys) |>
  select(protein_key, anc_TSB, stk_TSB, stp_TSB, anc_SN, stk_SN, stp_SN)
hm_mat <- as.matrix(hm[, -1]); rownames(hm_mat) <- hm$protein_key
hm_c <- hm_mat - rowMeans(hm_mat, na.rm = TRUE)
ord <- hclust(dist(replace(hm_c, is.na(hm_c), 0)))$order   # NA -> 0 only to ORDER rows, not shown
lab <- tibble(protein_key = rownames(hm_c)) |>
  left_join(distinct(pairs, protein_key, locus_tag, gene), by = "protein_key") |>
  mutate(label = if_else(is.na(gene), locus_tag, paste(locus_tag, gene)))
hm_long <- as_tibble(hm_c, rownames = "protein_key") |>
  pivot_longer(-protein_key, names_to = "group", values_to = "centred") |>
  left_join(lab, by = "protein_key") |>
  mutate(label = factor(label, lab$label[ord]),
         group = factor(group, c("anc_TSB", "stk_TSB", "stp_TSB", "anc_SN", "stk_SN", "stp_SN")))
lim <- max(abs(hm_long$centred), na.rm = TRUE)
p_hm <- ggplot(hm_long, aes(group, label, fill = centred)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", limits = c(-lim, lim),
                       na.value = "grey30") +
  geom_vline(xintercept = 3.5, colour = "black") +
  labs(x = NULL, y = NULL, fill = "log2,\ncentred\nper protein",
       title = sprintf("All %d proteins with a formal genotype hit (baseline or interaction)", length(cand_keys)),
       subtitle = "Mean of observed values per group, minus the protein's mean over the 6 groups (display only).\nDark grey = group not detected (0/3). Rows clustered for display.") +
  theme(axis.text.y = element_text(size = 5), plot.subtitle = element_text(size = 7.5))
ggsave(file.path(out, "fig_heatmap_mval.png"), p_hm, width = 7, height = max(6, 0.11 * length(cand_keys) + 2), dpi = 300)

# per-sample values: up to 4 shared (hit in both, same direction) + up to 4
# distinct (direct stp_vs_stk hits), strongest first
pick_shared <- pairs |> filter(class == "hit in both, same direction") |>
  arrange(pmax(stk_FDR, stp_FDR)) |> distinct(locus_tag) |> slice_head(n = 4) |> pull()
pick_distinct <- pairs |> filter(direct_hit, !locus_tag %in% pick_shared) |>
  arrange(direct_FDR) |> distinct(locus_tag) |> slice_head(n = 4) |> pull()
p_raw <- plot_protein(c(pick_shared, pick_distinct), abundances = read_abundances("total"), max_proteins = 8) +
  labs(title = "Representative shared and distinct mutant effects: per-sample values",
       subtitle = paste0("Shared (hit in both mutants, same direction): ", paste(pick_shared, collapse = ", "),
                         "\nDistinct (direct stp_vs_stk hit): ", paste(pick_distinct, collapse = ", ")))
ggsave(file.path(out, "fig_shared_distinct_raw_mval.png"), p_raw, width = 11, height = 7.5, dpi = 200)
message("shared shown: ", paste(pick_shared, collapse = ", "), " | distinct shown: ", paste(pick_distinct, collapse = ", "))

writeLines(c(paste("Step 7 shared/distinct, run", Sys.Date()), checks),
           file.path(out, "shared_distinct_checks_mval.txt"))
if (any(grepl("^\\[FAIL", checks))) stop("At least one check FAILED, see ", out)
message("done, written to ", out)
