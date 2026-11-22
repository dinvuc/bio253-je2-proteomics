# Step 8 (PROJECT_WORKFLOW_mval.md): does phosphosite coverage support a follow-up?
#
# Run from the repository root:   Rscript C_mutant_validation_DinMete/scripts/08_phospho_coverage_mval.R   (~1 min)
# Needs C_mutant_validation_DinMete/scripts/05_formal_test_mval.R first (contrast sheet, same model helpers).
#
# ---- Plan (coding rules, section 4) -----------------------------------------
# Aim:      Decide whether the phospho layer is good enough for Step 9 (matching
#           sites to protein changes): which sites are measured often enough and
#           localised well enough in the comparisons we use for H2/H3?
# Input:    read_sites("STYH" / "STY")       per-sample intensity + loc_prob (FragPipe)
#           read_abundances("phospho_STYH" / "phospho_STY")   normalised log2 values
#                                            the supplied DEA used (NOT re-normalised)
#           read_dea("phospho_STYH" / "phospho_STY")   supplied models, comparison only
#           read_samples(); step05 contrast_sheet_mval.csv (the 9 primary contrasts)
# Rules, FIXED BEFORE any site was selected (proposed 13 Nov 2026,
# to be confirmed by the group, rules section 10):
#   1. Real JE2 sites only: keep_je2() (locus tag, no rev_ decoy).
#   2. Eligibility per contrast = the Step 5 rule: >= 2 of 3 observed normalised
#      values in EVERY group the contrast uses. Missing stays NA.
#   3. Localisation: a site is "localised" when the MEDIAN of its per-sample
#      loc_prob, over the 18 primary samples in which it was detected, is >= 0.75.
#      Why median per sample: best_loc_prob is only the best score over all 27
#      samples, so one lucky spectrum would pass. 0.75 is the conventional
#      "class I" cutoff (Olsen et al. 2006, Cell 127:635, PMID 17081983); it is
#      NOT a measured false-localisation rate. best_loc_prob >= 0.75 is reported
#      next to it for comparison only.
#   4. Histidine (pH) sites are counted and listed but are exploratory only and
#      never enter the supported shortlist (no validated pH localisation cutoff;
#      pH is acid-labile).
#   5. Supported site (per contrast) = rules 1 + 2 + 3 and residue S/T/Y.
#   6. Primary search = the one with more supported S/T/Y sites summed over the
#      9 primary contrasts; on a tie STYH (its DEA uses lm_impute like the total
#      proteome). The other search is the sensitivity check.
# Analysis: counts and a coverage funnel; our own refit of the site values with the
#           Step 5 model (limma ~0+group, all 27 samples, treat lfc = 1) so that the
#           site effects are observed-only and comparable with the protein effects
#           in Step 9. NEW TEST FAMILY (rules section 5): BH within each contrast
#           over the SUPPORTED sites only, hit = BH-FDR < 0.05. This FDR is ours,
#           not the supplied one. Sensitivity: refit without the two reinjected
#           stk_pulseSN phospho samples; STY vs STYH estimates on shared sites.
# Output:   C_mutant_validation_DinMete/results/step08_phospho_coverage/
#             mapping_audit_mval.csv          per search: features, decoys, unmapped, joins
#             detection_per_sample_mval.csv   sites detected per sample (both searches)
#             residue_counts_mval.csv         S/T/Y/H: all, eligible, localised
#             coverage_per_contrast_mval.csv  funnel per primary contrast and search
#             supported_sites_mval.csv        the shortlist (primary search), one row per site
#             histidine_sites_mval.csv        pH sites, exploratory
#             phospho_hit_counts_mval.csv     our hits vs supplied screen per contrast
#             supplied_model_audit_mval.csv   modelName / estimate_type of the supplied DEAs
#             fig_detection_depth_mval.png, fig_localisation_mval.png,
#             fig_coverage_funnel_mval.png, fig_supported_hits_raw_mval.png,
#             fig_supported_hits_sites_mval.png, fig_sty_vs_styh_mval.png
#             figure_guide_mval.md, coverage_checks_mval.txt
#           data/derived/limma_treat_phospho_<search>_mval.tsv.gz  (local, git-ignored)
# Check:    27 samples / 3 per group after joining by name; one feature per site_id;
#           sites-table detection agrees with abundance-table detection; our site diff
#           has the same sign as the supplied diff on observed rows; every hit has
#           |diff| > 1.

suppressMessages(source("R/setup.R"))
out <- "C_mutant_validation_DinMete/results/step08_phospho_coverage"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
checks <- character()
check <- function(label, ok) {
  line <- sprintf("[%s] %s", if (isTRUE(ok)) "PASS" else "FAIL", label)
  message(line); checks <<- c(checks, line)
  invisible(ok)
}
note <- function(...) { line <- paste0("[NOTE] ", ...); message(line); checks <<- c(checks, line) }

LOC_CUT  <- 0.75
SEARCHES <- c("STYH", "STY")
# `reinjected` is a property of the phospho runs (stk_pulseSN_1/2), so it comes
# from the phospho table, joined to the sample sheet by sample name.
reinj_tbl <- read_abundances("phospho_STYH") |> distinct(sample, reinjected)
samples  <- read_samples() |> arrange(genotype, medium, replicate) |>
  left_join(reinj_tbl, by = "sample")
check("reinjected flag: one value per sample, known for all 27", !anyDuplicated(reinj_tbl$sample) && !anyNA(samples$reinjected))
check("27 samples, 3 per group", nrow(samples) == 27 && all(table(samples$group) == 3))
PRIMARY_SAMPLES <- samples$sample[samples$medium %in% c("TSB", "SN")]
check("18 primary samples (anc/stk/stp x TSB/SN)", length(PRIMARY_SAMPLES) == 18)

sheet   <- read_csv("C_mutant_validation_DinMete/results/step05_formal_test/contrast_sheet_mval.csv", show_col_types = FALSE)
PRIMARY <- sheet |> filter(primary_18_samples) |> select(contrast, family, formula)
check("9 primary contrasts in the Step 5 sheet", nrow(PRIMARY) == 9)
contr_groups <- PRIMARY |> mutate(groups = purrr::map(formula, contrast_groups)) |>
  select(contrast, groups) |> tidyr::unnest(groups)

meta <- read_protein_metadata() |> filter(ortholog_quality == "good") |>
  distinct(locus_tag, gene = gene_symbol_best)

# ---- 1. load one search, audit mapping and joins --------------------------------
load_search <- function(s) {
  sites <- read_sites(s)
  abund <- read_abundances(paste0("phospho_", s))
  # The STY abundance table stores observed values only (no NA rows). Complete it
  # to every site x sample so that "not detected" is an explicit NA (rule 6).
  feat  <- abund |> distinct(site_id, protein_key, locus_tag, product, residue, position, Peptide)
  abund_full <- feat |>
    tidyr::crossing(sample = samples$sample) |>
    left_join(abund |> select(site_id, sample, normalized_abundance, raw_intensity,
                              Localization_Probability), by = c("site_id", "sample")) |>
    left_join(select(samples, sample, genotype, medium, replicate, group, reinjected), by = "sample")
  audit <- tibble(
    search = s,
    abundance_rows_supplied = nrow(abund),
    abundance_table_has_NA_rows = any(is.na(abund$normalized_abundance)),
    features = n_distinct(abund$site_id),
    features_per_site_id_max = max(count(distinct(abund, site_id, site), site_id)$n),
    decoy_rev = n_distinct(abund$site_id[str_detect(abund$protein_key, "^rev_")]),
    no_locus_tag = n_distinct(abund$site_id[is.na(abund$locus_tag) & !str_detect(abund$protein_key, "^rev_")]),
    je2_features = n_distinct(keep_je2(abund)$site_id),
    site_table_sites = n_distinct(sites$site_id),
    in_abundance_not_in_site_table = length(setdiff(unique(abund$site_id), unique(sites$site_id))),
    in_site_table_not_in_abundance = length(setdiff(unique(sites$site_id), unique(abund$site_id))),
    site_id_starts_with_protein_key = all(str_starts(abund$site_id, fixed(abund$protein_key))),
    # `Peptide` is a ";"-separated list of the peptide forms behind a site. Split
    # it and compare unmodified sequences: sites that share a sequence are NOT
    # independent evidence. (Fixed 16 Nov: the first version compared the whole
    # list string, which is not a count of shared peptides.)
    sites_with_2plus_peptide_forms = keep_je2(abund) |> distinct(site_id, Peptide) |>
      filter(str_detect(Peptide, ";")) |> nrow(),
    sites_sharing_a_peptide_sequence = keep_je2(abund) |> distinct(site_id, Peptide) |>
      tidyr::separate_rows(Peptide, sep = ";") |> mutate(seq = toupper(Peptide)) |>
      distinct(site_id, seq) |> group_by(seq) |> filter(n() > 1) |> ungroup() |>
      pull(site_id) |> n_distinct())
  list(sites = sites, abund = abund_full, audit = audit)
}
dat <- purrr::set_names(SEARCHES) |> purrr::map(load_search)
audit <- purrr::map(dat, "audit") |> bind_rows()
write_csv(audit, file.path(out, "mapping_audit_mval.csv"))
print(audit |> tidyr::pivot_longer(-search, values_transform = as.character) |>
        tidyr::pivot_wider(names_from = search), n = 30)

for (s in SEARCHES) {
  a <- dat[[s]]$abund
  check(paste(s, ": one row per site x sample after completing, 27 samples, 3 per group"),
        !anyDuplicated(a[c("site_id", "sample")]) && n_distinct(a$sample) == 27 &&
          all(table(distinct(a, sample, group)$group) == 3))
  check(paste(s, ": one feature per site_id (no duplicated peptide features)"),
        audit$features_per_site_id_max[audit$search == s] == 1)
  # Does a value in the abundance table correspond to a detection in the site table?
  agree <- a |> keep_je2() |>
    inner_join(dat[[s]]$sites |> select(site_id, sample, intensity, loc_prob), by = c("site_id", "sample")) |>
    summarise(n = n(), same = mean(is.na(normalized_abundance) == is.na(intensity)),
              abund_only = sum(!is.na(normalized_abundance) & is.na(intensity)),
              site_only = sum(is.na(normalized_abundance) & !is.na(intensity)))
  note(sprintf("%s: detection agrees between site table and abundance table in %.2f%% of %d site x sample cells (value only in abundance: %d, only in site table: %d)",
               s, 100 * agree$same, agree$n, agree$abund_only, agree$site_only))
}

# ---- 2. detection depth per sample ----------------------------------------------
depth <- purrr::imap(dat, \(d, s) d$abund |> keep_je2() |>
                       group_by(sample, genotype, medium, group, reinjected) |>
                       summarise(sites_detected = sum(!is.na(normalized_abundance)), .groups = "drop") |>
                       mutate(search = s)) |> bind_rows()
write_csv(depth, file.path(out, "detection_per_sample_mval.csv"))
print(depth |> group_by(search, group) |>
        summarise(min = min(sites_detected), median = median(sites_detected),
                  max = max(sites_detected), .groups = "drop") |>
        tidyr::pivot_wider(names_from = search, values_from = c(min, median, max)), n = 20)

p <- ggplot(depth |> add_x_position(), aes(x_pos, sites_detected, colour = genotype, shape = reinjected)) +
  geom_point(size = 2.3, position = position_jitter(width = 0.04, height = 0, seed = 1)) +
  facet_wrap(~search) + medium_x_axis() +
  scale_colour_manual(values = GENOTYPE_COLOURS) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) +
  expand_limits(y = 0) +
  labs(x = "medium", y = "JE2 phosphosites with a value (per sample)", shape = "reinjected",
       title = "Step 8: phosphosite detection depth per sample",
       subtitle = "Each point = one sample. Primary analysis uses TSB and SN only.")
ggsave(file.path(out, "fig_detection_depth_mval.png"), p, width = 9, height = 4.5, dpi = 150)

# ---- 3. localisation per site (rule 3) and eligibility per contrast (rule 2) ----
site_info <- function(s) {
  d <- dat[[s]]
  loc <- d$sites |> filter(sample %in% PRIMARY_SAMPLES, !is.na(intensity)) |>
    group_by(site_id) |>
    summarise(n_primary_detected = n(), median_loc_prob = median(loc_prob, na.rm = TRUE),
              .groups = "drop")
  best <- d$sites |> distinct(site_id, best_loc_prob)
  n_obs <- d$abund |> keep_je2() |>
    group_by(site_id, group) |> summarise(n = sum(!is.na(normalized_abundance)), .groups = "drop")
  elig <- contr_groups |> inner_join(n_obs, by = c(groups = "group"), relationship = "many-to-many") |>
    group_by(site_id, contrast) |> summarise(eligible = all(n >= MVAL_MIN_OBS), .groups = "drop")
  info <- d$abund |> keep_je2() |>
    distinct(site_id, protein_key, locus_tag, product, residue, position, peptide = Peptide) |>
    left_join(loc, by = "site_id") |> left_join(best, by = "site_id") |>
    mutate(localised = !is.na(median_loc_prob) & median_loc_prob >= LOC_CUT,
           localised_best_only = !is.na(best_loc_prob) & best_loc_prob >= LOC_CUT)
  list(info = info, elig = elig,
       n_wide = n_obs |> mutate(group = paste0("n_", group)) |>
         tidyr::pivot_wider(names_from = group, values_from = n))
}
si <- purrr::map(purrr::set_names(SEARCHES), site_info)

residue_counts <- purrr::imap(si, \(x, s) {
  any_elig <- x$elig |> group_by(site_id) |> summarise(eligible_any_primary = any(eligible))
  x$info |> left_join(any_elig, by = "site_id") |>
    group_by(residue) |>
    summarise(all_sites = n(),
              detected_in_primary = sum(!is.na(n_primary_detected)),
              eligible_in_any_primary_contrast = sum(eligible_any_primary, na.rm = TRUE),
              localised_median_per_sample = sum(localised),
              # same denominator as the median rule: sites detected in >= 1 primary
              # sample (fixed 16 Nov: first version counted best_loc_prob over ALL sites)
              localised_best_loc_prob_only = sum(localised_best_only & !is.na(n_primary_detected)),
              eligible_and_localised = sum(eligible_any_primary & localised, na.rm = TRUE)) |>
    mutate(search = s, .before = 1)
}) |> bind_rows()
write_csv(residue_counts, file.path(out, "residue_counts_mval.csv"))
print(residue_counts)

loc_long <- purrr::imap(si, \(x, s) x$info |> mutate(search = s)) |> bind_rows() |>
  filter(!is.na(n_primary_detected)) |>
  select(search, site_id, residue, `median per-sample loc_prob (our rule)` = median_loc_prob,
         `best_loc_prob (best of 27 samples)` = best_loc_prob) |>
  tidyr::pivot_longer(c(-search, -site_id, -residue), names_to = "score", values_to = "loc_prob")
p <- ggplot(loc_long, aes(loc_prob, fill = score)) +
  geom_histogram(binwidth = 0.05, boundary = 0, position = "identity", alpha = 0.55) +
  geom_vline(xintercept = LOC_CUT, linetype = "dashed") +
  facet_grid(residue ~ search, scales = "free_y") +
  scale_fill_manual(values = c("#0072B2", "#D55E00")) +
  labs(x = "localisation probability", y = "sites",
       title = "Step 8: localisation, per-sample median vs best score",
       subtitle = "Sites detected in >= 1 of the 18 primary samples. Dashed = 0.75.",
       fill = NULL) +
  theme(legend.position = "bottom")
ggsave(file.path(out, "fig_localisation_mval.png"), p, width = 8, height = 8, dpi = 150)

# ---- 4. coverage funnel per primary contrast ------------------------------------
coverage <- purrr::imap(si, \(x, s) {
  x$elig |> left_join(select(x$info, site_id, residue, localised), by = "site_id") |>
    group_by(contrast) |>
    summarise(je2_sites = n(),
              eligible_2of3 = sum(eligible),
              eligible_localised_all_residues = sum(eligible & localised),
              supported_STY = sum(eligible & localised & residue %in% c("S", "T", "Y")),
              eligible_localised_H = sum(eligible & localised & residue == "H"),
              .groups = "drop") |>
    mutate(search = s, .before = 1)
}) |> bind_rows() |>
  mutate(contrast = factor(contrast, levels = PRIMARY$contrast)) |> arrange(search, contrast)
write_csv(coverage, file.path(out, "coverage_per_contrast_mval.csv"))
print(coverage, n = 20)

totals <- coverage |> group_by(search) |> summarise(supported = sum(supported_STY))
PRIMARY_SEARCH <- if (totals$supported[totals$search == "STY"] > totals$supported[totals$search == "STYH"]) "STY" else "STYH"
OTHER_SEARCH   <- setdiff(SEARCHES, PRIMARY_SEARCH)
note(sprintf("primary search by rule 6: %s (supported S/T/Y site x contrast pairs: STYH %d, STY %d)",
             PRIMARY_SEARCH, totals$supported[totals$search == "STYH"], totals$supported[totals$search == "STY"]))

p <- coverage |>
  select(search, contrast, `JE2 sites` = je2_sites, `eligible (2 of 3 per group)` = eligible_2of3,
         `supported (eligible + localised, S/T/Y)` = supported_STY) |>
  tidyr::pivot_longer(-c(search, contrast), names_to = "stage", values_to = "sites") |>
  mutate(stage = factor(stage, levels = unique(stage))) |>
  ggplot(aes(sites, contrast, fill = stage)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = sites), position = position_dodge(width = 0.8), hjust = -0.1, size = 2.4) +
  facet_wrap(~search) +
  scale_fill_manual(values = c("grey75", "#56B4E9", "#0072B2")) +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "phosphosites", y = NULL, fill = NULL,
       title = "Step 8: how many phosphosites can each primary contrast use?") +
  theme(legend.position = "bottom")
ggsave(file.path(out, "fig_coverage_funnel_mval.png"), p, width = 11, height = 6, dpi = 150)

# ---- 5. our refit of the site values (observed only, Step 5 model) --------------
# abundance_matrix() / fit_limma_treat() key rows by `protein_key`; here one row
# is one SITE, so we pass site_id in that column and rename it back afterwards.
fit_sites <- function(s, drop_samples = character()) {
  smp <- samples |> filter(!sample %in% drop_samples)
  a <- dat[[s]]$abund |> keep_je2() |> filter(sample %in% smp$sample) |>
    transmute(protein_key = site_id, sample, normalized_abundance)
  m <- abundance_matrix(a, smp)
  fit_limma_treat(m, smp) |> rename(site_id = protein_key) |>
    filter(contrast %in% PRIMARY$contrast)
}
add_supported_fdr <- function(res, s) {
  # NEW test family: BH within contrast over the SUPPORTED sites only (rule 5).
  sup <- si[[s]]$elig |> filter(eligible) |>
    inner_join(si[[s]]$info |> filter(localised, residue %in% c("S", "T", "Y")) |> select(site_id),
               by = "site_id") |> mutate(supported = TRUE)
  res |> left_join(sup |> select(site_id, contrast, supported), by = c("site_id", "contrast")) |>
    mutate(supported = coalesce(supported, FALSE),
           p.value = if_else(supported, p.value, NA_real_)) |>
    group_by(contrast) |> mutate(FDR = p.adjust(p.value, "BH")) |> ungroup() |>
    mutate(hit = supported & FDR < MVAL_FDR)
}
fits <- purrr::map(purrr::set_names(SEARCHES), \(s) {
  fit_sites(s) |> add_supported_fdr(s) |>
    left_join(si[[s]]$info, by = "site_id") |> left_join(meta, by = "locus_tag") |>
    left_join(select(PRIMARY, contrast, family), by = "contrast") |> mutate(search = s)
})
for (s in SEARCHES) {
  write_tsv(fits[[s]], sprintf("data/derived/limma_treat_phospho_%s_mval.tsv.gz", s))  # local only
  check(paste(s, ": tested == eligible (fit_limma_treat uses the same 2-of-3 rule)"),
        fits[[s]] |> inner_join(si[[s]]$elig, by = c("site_id", "contrast")) |>
          summarise(all(tested == eligible)) |> pull())
  check(paste(s, ": every hit has |diff| > 1 (treat lfc = 1)"),
        all(abs(fits[[s]]$diff[fits[[s]]$hit]) > MVAL_LFC))
}

# Guardrail: same sign as the supplied DEA where both are on observed values.
supplied <- purrr::map(purrr::set_names(SEARCHES), \(s)
  read_dea(paste0("phospho_", s), contrasts = PRIMARY$contrast) |> keep_je2() |>
    select(site_id, contrast, diff_supplied = diff, FDR_supplied = FDR, estimate_type, modelName))
for (s in SEARCHES) {
  cmp <- fits[[s]] |> filter(supported) |> inner_join(supplied[[s]], by = c("site_id", "contrast"))
  r <- cor(cmp$diff, cmp$diff_supplied)
  check(sprintf("%s: our site diff vs supplied diff on supported sites, r = %.3f (> 0.9, same sign convention)", s, r),
        r > 0.9)
}

# Supplied model audit (what the course tables say about themselves)
model_audit <- purrr::imap(supplied, \(x, s) x |> count(contrast, modelName, estimate_type) |>
                             mutate(search = s, .before = 1)) |> bind_rows()
write_csv(model_audit, file.path(out, "supplied_model_audit_mval.csv"))
note("STY supplied DEA: modelName ", paste(unique(supplied$STY$modelName), collapse = " / "),
     " but estimate_type is always ", paste(unique(supplied$STY$estimate_type), collapse = "/"),
     " -> in STY the label cannot separate imputed from observed estimates.")

# Hit counts: ours (primary search and the other) vs the supplied screen
hit_counts <- purrr::imap(fits, \(f, s) {
  f |> left_join(supplied[[s]], by = c("site_id", "contrast")) |>
    group_by(contrast) |>
    summarise(supported_sites = sum(supported),
              our_hits_up = sum(hit & diff > 0), our_hits_down = sum(hit & diff < 0),
              supplied_screen_hits_all_sites = sum(FDR_supplied < 0.1 & abs(diff_supplied) > 1, na.rm = TRUE),
              supplied_screen_hits_supported = sum(supported & FDR_supplied < 0.1 & abs(diff_supplied) > 1, na.rm = TRUE),
              .groups = "drop") |>
    mutate(search = s, .before = 1)
}) |> bind_rows() |>
  left_join(sheet |> select(contrast, meaning), by = "contrast") |>
  mutate(contrast = factor(contrast, levels = PRIMARY$contrast)) |> arrange(search, contrast)
write_csv(hit_counts, file.path(out, "phospho_hit_counts_mval.csv"))
print(hit_counts |> select(-meaning), n = 20)

# ---- 6. sensitivity checks -------------------------------------------------------
# (a) the two reinjected stk_pulseSN phospho samples only change the pooled
#     variance for the primary contrasts. Refit without them.
reinj <- samples$sample[samples$reinjected]
check("2 reinjected phospho samples (stk_pulseSN_1/2)", length(reinj) == 2 && all(str_detect(reinj, "stk_pulseSN")))
sens <- fit_sites(PRIMARY_SEARCH, drop_samples = reinj) |> add_supported_fdr(PRIMARY_SEARCH)
sens_cmp <- fits[[PRIMARY_SEARCH]] |> select(site_id, contrast, hit) |>
  inner_join(sens |> select(site_id, contrast, hit_noreinj = hit), by = c("site_id", "contrast"))
note(sprintf("reinjected samples removed (%s): hits %d -> %d, kept %d, lost %d, new %d",
             PRIMARY_SEARCH, sum(sens_cmp$hit), sum(sens_cmp$hit_noreinj),
             sum(sens_cmp$hit & sens_cmp$hit_noreinj), sum(sens_cmp$hit & !sens_cmp$hit_noreinj),
             sum(!sens_cmp$hit & sens_cmp$hit_noreinj)))

# (b) STY vs STYH on sites supported in both searches (same site_id string)
both <- fits$STYH |> filter(supported) |> select(site_id, contrast, locus_tag, residue, position, diff_STYH = diff, hit_STYH = hit) |>
  inner_join(fits$STY |> filter(supported) |> select(site_id, contrast, diff_STY = diff, hit_STY = hit),
             by = c("site_id", "contrast"))
r_both <- cor(both$diff_STYH, both$diff_STY)
note(sprintf("STY vs STYH: %d supported site x contrast pairs in both; r = %.3f; hits STYH %d, STY %d, both %d",
             nrow(both), r_both, sum(both$hit_STYH), sum(both$hit_STY), sum(both$hit_STYH & both$hit_STY)))
both <- both |> mutate(hit_class = case_when(hit_STYH & hit_STY ~ "hit in both",
                                              hit_STYH ~ "hit in STYH only", hit_STY ~ "hit in STY only",
                                              TRUE ~ "hit in neither"))
p <- ggplot(both, aes(diff_STYH, diff_STY, colour = hit_class)) +
  geom_abline(linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.6, size = 1.2) +
  facet_wrap(~ factor(contrast, levels = PRIMARY$contrast), ncol = 3) +
  scale_colour_manual(values = c(`hit in neither` = "grey65", `hit in STYH only` = "#0072B2",
                                 `hit in STY only` = "#E69F00", `hit in both` = "#D55E00")) +
  coord_equal() +
  labs(x = "site effect, STYH search (log2, A - B)", y = "site effect, STY search (log2, A - B)", colour = NULL,
       title = "Step 8: STY vs STYH, same site, same contrast (our refit)",
       subtitle = sprintf("Sites supported in both searches; r = %.2f; dashed = identity", r_both)) +
  theme(legend.position = "bottom", strip.text = element_text(size = 7))
ggsave(file.path(out, "fig_sty_vs_styh_mval.png"), p, width = 9, height = 10, dpi = 150)

# ---- 7. shortlist (primary search) and histidine list -----------------------------
f <- fits[[PRIMARY_SEARCH]]
shortlist <- f |> filter(supported) |>
  group_by(site_id, locus_tag, gene, product, residue, position, peptide, median_loc_prob, best_loc_prob) |>
  summarise(n_contrasts_supported = n(),
            hit_contrasts = paste(contrast[hit], collapse = ";"),
            .groups = "drop") |>
  left_join(si[[PRIMARY_SEARCH]]$n_wide |> select(site_id, any_of(paste0("n_", GROUPS[str_detect(GROUPS, "_(TSB|SN)$")]))),
            by = "site_id") |>
  mutate(search = PRIMARY_SEARCH, .before = 1) |>
  arrange(desc(hit_contrasts != ""), locus_tag, position)
write_csv(shortlist, file.path(out, "supported_sites_mval.csv"))
note(sprintf("shortlist (%s): %d supported S/T/Y sites on %d proteins; %d sites are a hit in >= 1 primary contrast",
             PRIMARY_SEARCH, nrow(shortlist), n_distinct(shortlist$locus_tag), sum(shortlist$hit_contrasts != "")))

his <- purrr::imap(si, \(x, s) x$info |> filter(residue == "H") |>
                     left_join(x$elig |> group_by(site_id) |> summarise(eligible_contrasts = sum(eligible)), by = "site_id") |>
                     left_join(meta, by = "locus_tag") |> mutate(search = s, .before = 1)) |> bind_rows()
write_csv(his, file.path(out, "histidine_sites_mval.csv"))

# ---- 8. raw values behind the hits (rule: plot before trusting a fold change) ----
top <- f |> filter(hit) |> arrange(FDR) |> distinct(site_id, .keep_all = TRUE) |> head(9)
if (nrow(top) > 0) {
  lab <- top |> transmute(site_id, label = paste0(coalesce(gene, locus_tag), " ", residue, position,
                                                  "\nbest in ", contrast, " (", sprintf("%+.1f", diff), ")"))
  raw <- dat[[PRIMARY_SEARCH]]$abund |> filter(site_id %in% top$site_id, medium %in% c("TSB", "SN")) |>
    inner_join(lab, by = "site_id") |> add_x_primary()
  ndet <- raw |> group_by(label, genotype, medium, x_pos) |>
    summarise(n = sum(!is.na(normalized_abundance)), .groups = "drop")
  p <- ggplot(raw, aes(x_pos, normalized_abundance, colour = genotype)) +
    geom_point(position = position_jitter(width = 0.04, height = 0, seed = 1), na.rm = TRUE, size = 1.8) +
    geom_text(data = ndet, aes(y = Inf, label = paste0(n, "/3")), vjust = 1.3, size = 2.5, show.legend = FALSE) +
    facet_wrap(~label, scales = "free_y") +
    primary_x_axis() +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    scale_colour_manual(values = GENOTYPE_COLOURS) +
    labs(x = "medium", y = "normalised site abundance (log2, value the model used)",
         title = sprintf("Step 8: per-sample values of the %d strongest supported site hits (%s)", nrow(top), PRIMARY_SEARCH),
         caption = "numbers = samples with a value (of 3); missing = not detected. Effect in the label: A - B of that contrast.") +
    theme(strip.text = element_text(size = 7))
  ggsave(file.path(out, "fig_supported_hits_raw_mval.png"), p, width = 10, height = 8, dpi = 150)
  p <- plot_sites(top$site_id, dat[[PRIMARY_SEARCH]]$sites) +
    labs(title = "Step 8: the same sites in the FragPipe site table (all 27 samples, not normalised)") +
    theme(strip.text = element_text(size = 7))
  ggsave(file.path(out, "fig_supported_hits_sites_mval.png"), p, width = 10, height = 8, dpi = 150)
}

writeLines(c(sprintf("Step 8 checks, %s", Sys.Date()), checks), file.path(out, "coverage_checks_mval.txt"))
message("done: ", out)
