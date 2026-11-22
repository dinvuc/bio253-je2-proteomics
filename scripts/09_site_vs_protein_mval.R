# Step 9 (PROJECT_WORKFLOW_mval.md): match supported phosphosites to their protein.
#
# Run from the repository root:   Rscript scripts/09_site_vs_protein_mval.R   (~30 s)
# Needs scripts/05_formal_test_mval.R and scripts/08_phospho_coverage_mval.R first.
#
# ---- Plan (coding rules, section 4) -----------------------------------------
# Aim:      A phosphosite signal goes up if the site is phosphorylated more OR if
#           there is simply more of the protein. For every supported site hit
#           from Step 8: does the site change follow the protein change, or is
#           something left over once the protein change is taken out?
# Input:    data/derived/limma_treat_phospho_<search>_mval.tsv.gz  (Step 8, our refit)
#           data/derived/limma_treat_total_mval.tsv.gz             (Step 5, our refit)
#           Both: limma ~0+group on all 27 samples, observed values only, so the
#           two estimates come from the same kind of model and the same contrast.
#           read_abundances("phospho_<search>" / "total") for the per-sample plots.
#           mval_medi/results/step08_phospho_coverage/supported_sites_mval.csv
#           (gives the primary search chosen in Step 8).
# Analysis: match site -> protein by the exact protein_key (+ locus_tag check) and EXACTLY the
#           same contrast name. Then, per supported site x primary contrast:
#             adjusted = diff_site - diff_protein   (log2; "relative protein-adjusted
#             phosphorylation signal", NOT occupancy/stoichiometry)
#           only when BOTH are observed estimates (site supported, protein tested).
#           Nothing imputed is ever subtracted (CLAUDE.md rule 5).
#           Uncertainty: se_adj = sqrt(se_site^2 + se_protein^2). This ASSUMES the
#           two estimates are independent (Cov = 0). The covariance term
#           -2 Cov(site, protein) cannot be estimated here: the phospho-enriched and
#           total runs are separate measurements and Step 2 did not establish that
#           sample "X_1" is the same culture aliquot in both layers. If they are
#           positively correlated, se_adj is too large (conservative).
#           NO p-value and NO FDR for the adjusted effect (rules section 5: do not
#           reuse the site FDR for a new quantity). Categories are DESCRIPTIVE.
# Categories, FIXED BEFORE looking at the matched values (proposed 13 Nov 2026,
# marked as checked by the group on 13 Nov 2026), for site hits (Step 8: treat lfc = 1, our BH-FDR < 0.05):
#   "consistent with protein amount": protein tested and |adjusted| < 1
#                                     (site change within 2-fold of the protein change)
#   "residual adjusted change":       protein tested, |adjusted| >= 1 AND
#                                     |adjusted| >= 2 * se_adj
#   "unresolved":                     protein not tested in this contrast (< 2 of 3
#                                     in a group), or |adjusted| >= 1 but < 2 * se_adj
#           These categories do NOT identify Stk1/Stp1 substrates (workflow Step 9).
# Output:   mval_medi/results/step09_site_vs_protein/
#           data/derived/site_protein_matched_mval.tsv.gz  (local, git-ignored: full table)
#                                             every supported site x primary contrast:
#                                             site + protein estimates, SE, adjusted, category
#             site_hit_categories_mval.csv    site hits only, with category
#             category_counts_mval.csv        per contrast (primary search + sensitivity)
#             unmatched_sites_mval.csv        sites whose protein is missing or not tested
#             site_vs_protein_tracking_mval.csv  per contrast: r(site, protein), median offsets
#             offset_sensitivity_mval.csv     residual hits after removing the global layer offset
#             fig_site_vs_protein_mval.png    Fig 4a: site effect vs protein effect
#             fig_adjusted_effects_mval.png   Fig 4b: adjusted effect +/- 2 se_adj, site hits
#             fig_paired_layers_mval.png      Fig 4c: per-sample site and protein values
#             figure_guide_mval.md, site_vs_protein_checks_mval.txt
# Check:    protein_key unique in the total table; contrast names identical; one example
#           recomputed by hand from group means; sample names identical in both layers;
#           category counts add up to the number of site hits.

suppressMessages(source("R/setup.R"))
out <- "mval_medi/results/step09_site_vs_protein"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
checks <- character()
check <- function(label, ok) {
  line <- sprintf("[%s] %s", if (isTRUE(ok)) "PASS" else "FAIL", label)
  message(line); checks <<- c(checks, line)
  invisible(ok)
}
note <- function(...) { line <- paste0("[NOTE] ", ...); message(line); checks <<- c(checks, line) }

ADJ_CUT <- 1   # |adjusted| below this = "consistent with protein amount" (2-fold)
SE_MULT <- 2   # residual must also be >= 2 approximate SE
lab_site <- \(d) paste0(coalesce(d$gene, d$locus_tag), " ", d$residue, d$position)

shortlist <- read_csv("mval_medi/results/step08_phospho_coverage/supported_sites_mval.csv", show_col_types = FALSE)
PRIMARY_SEARCH <- unique(shortlist$search)
stopifnot(length(PRIMARY_SEARCH) == 1)
OTHER_SEARCH <- setdiff(c("STYH", "STY"), PRIMARY_SEARCH)
note("primary search from Step 8: ", PRIMARY_SEARCH, "; sensitivity: ", OTHER_SEARCH)
sheet   <- read_csv("mval_medi/results/step05_formal_test/contrast_sheet_mval.csv", show_col_types = FALSE)
PRIMARY <- sheet |> filter(primary_18_samples) |> pull(contrast)

total <- read_tsv("data/derived/limma_treat_total_mval.tsv.gz", show_col_types = FALSE) |>
  filter(contrast %in% PRIMARY)
# Match by the EXACT protein_key (same search FASTA in both layers), not by locus
# tag: one locus tag (B7H15_13880) also has a 6-frame pseudogene/frameshift entry
# in the total proteome, so a locus-tag join could pick the wrong entry.
prot_keys <- distinct(total, protein_key, locus_tag)
check("total proteome: protein_key unique (exact match target)", !anyDuplicated(prot_keys$protein_key))
amb <- prot_keys |> group_by(locus_tag) |> filter(n() > 1) |> ungroup()
note("locus tags with > 1 total-proteome entry (not used for matching): ",
     paste(unique(amb$locus_tag), collapse = ", "), " -> ", paste(amb$protein_key, collapse = " / "))

# ---- 1. match site -> protein, same contrast --------------------------------------
match_search <- function(s) {
  sites <- read_tsv(sprintf("data/derived/limma_treat_phospho_%s_mval.tsv.gz", s), show_col_types = FALSE) |>
    filter(supported)
  check(paste(s, ": site table has only the 9 primary contrasts"), setequal(unique(sites$contrast), PRIMARY))
  m <- sites |>
    select(site_id, protein_key, locus_tag, gene, product, residue, position, peptide,
           contrast, family, diff_site = diff, se_site = std.error, FDR_site = FDR, hit_site = hit) |>
    left_join(total |> select(protein_key, locus_tag_prot = locus_tag, contrast,
                              diff_protein = diff, se_protein = std.error, FDR_protein = FDR,
                              tested_protein = tested, min_obs_protein = min_obs_in_groups),
              by = c("protein_key", "contrast"), relationship = "many-to-one") |>
    mutate(search = s, .before = 1,
           matched = !is.na(locus_tag_prot),
           tested_protein = coalesce(tested_protein, FALSE),
           adjusted = if_else(tested_protein, diff_site - diff_protein, NA_real_),
           se_adj = if_else(tested_protein, sqrt(se_site^2 + se_protein^2), NA_real_),
           category = case_when(
             !hit_site ~ "no site hit",
             !tested_protein ~ "unresolved (protein not tested)",
             abs(adjusted) < ADJ_CUT ~ "consistent with protein amount",
             abs(adjusted) >= SE_MULT * se_adj ~ "residual adjusted change",
             TRUE ~ "unresolved (adjusted effect uncertain)"))
  check(paste(s, ": every site matched to a total protein by exact protein_key"), all(m$matched))
  check(paste(s, ": locus_tag identical in both layers for every matched site"),
        all(m$locus_tag == m$locus_tag_prot, na.rm = TRUE))
  check(paste(s, ": one row per site x contrast after the join"), !anyDuplicated(m[c("site_id", "contrast")]))
  check(paste(s, ": no adjusted value uses an untested (sparse) protein estimate"),
        all(is.na(m$adjusted[!m$tested_protein])))
  m
}
matched <- purrr::map(purrr::set_names(c(PRIMARY_SEARCH, OTHER_SEARCH)), match_search)
mp <- matched[[PRIMARY_SEARCH]] |> mutate(contrast = factor(contrast, levels = PRIMARY))

# full table: local only (git-ignored), like Steps 5 and 7
write_tsv(mp |> arrange(contrast, category, desc(abs(diff_site))), "data/derived/site_protein_matched_mval.tsv.gz")
write_csv(mp |> filter(!matched | !tested_protein) |> select(search, site_id, locus_tag, gene, residue, position,
                                                            contrast, matched, tested_protein, min_obs_protein),
          file.path(out, "unmatched_sites_mval.csv"))
samples <- read_samples()
site_ab <- read_abundances(paste0("phospho_", PRIMARY_SEARCH)) |> keep_je2()
prot_ab <- read_abundances("total") |> keep_je2()
check("phospho and total layers use the same 27 sample names (join by name possible)",
      setequal(unique(site_ab$sample), samples$sample) && setequal(unique(prot_ab$sample), samples$sample))

# QC FLAG for each hit (added 13 Nov AFTER the first run showed B7H15_04940 T269
# at +10.8 resting on two stp_SN values near the detection floor). It does NOT
# change the categories or the hit rule (fixed before results); it shows how
# fragile a hit is: detections per group and the largest within-group range of
# the observed site values (log2) among the groups the contrast uses.
site_grp <- site_ab |> group_by(site_id, group) |>
  summarise(n = sum(!is.na(normalized_abundance)),
            range = if (n > 1) diff(range(normalized_abundance, na.rm = TRUE)) else NA_real_, .groups = "drop")
hit_qc <- mp |> filter(hit_site) |> distinct(site_id, contrast) |>
  mutate(groups = purrr::map(as.character(contrast), \(ct) contrast_groups(mval_contrasts()[[ct]]))) |>
  tidyr::unnest(groups) |>
  left_join(site_grp, by = c("site_id", groups = "group")) |>
  group_by(site_id, contrast) |>
  summarise(detected_per_group = paste0(groups, " ", n, "/3", collapse = "; "),
            max_within_group_range = max(range, na.rm = TRUE), .groups = "drop") |>
  mutate(fragile_spread_gt3 = max_within_group_range > 3)
hits <- mp |> filter(hit_site) |> left_join(hit_qc, by = c("site_id", "contrast")) |>
  arrange(contrast, category, desc(abs(adjusted)))
write_csv(hits, file.path(out, "site_hit_categories_mval.csv"))
note(sprintf("QC flag: %d of %d site hits have a within-group range > 3 log2 in a group of their contrast: %s",
             sum(hits$fragile_spread_gt3), nrow(hits),
             paste(unique(paste0(coalesce(hits$gene, hits$locus_tag), " ", hits$residue, hits$position)[hits$fragile_spread_gt3]),
                   collapse = ", ")))

counts <- bind_rows(matched) |> filter(hit_site) |>
  count(search, contrast, category) |>
  tidyr::pivot_wider(names_from = category, values_from = n, values_fill = 0) |>
  mutate(contrast = factor(contrast, levels = PRIMARY)) |> arrange(search, contrast)
write_csv(counts, file.path(out, "category_counts_mval.csv"))
print(counts, n = 30, width = 200)
check("category counts add up to the number of site hits (primary search)",
      sum(counts |> filter(search == PRIMARY_SEARCH) |> select(-search, -contrast)) == nrow(hits))

# How much does the protein change explain across ALL supported sites (not only hits)?
track <- mp |> filter(tested_protein) |> group_by(contrast) |>
  summarise(sites_with_both = n(), r_site_vs_protein = cor(diff_site, diff_protein),
            median_abs_site = median(abs(diff_site)), median_abs_adjusted = median(abs(adjusted)),
            # A median SIGNED adjusted effect far from 0 means a global offset between
            # the layers (e.g. different enrichment/normalisation in TSB and SN), not
            # site-specific regulation. Report it; do not correct for it.
            median_signed_adjusted = median(adjusted),
            .groups = "drop")
write_csv(track, file.path(out, "site_vs_protein_tracking_mval.csv"))
print(track)

# SENSITIVITY (added after seeing median_signed_adjusted of -0.2 to -0.65 in the
# TSB_vs_SN contrasts; in anc and stk the SN phospho samples detect ~2x more sites, Step 8). If that
# global offset is a layer artefact, how many "residual" hits remain once it is
# removed? Labelled sensitivity only: the main categories above are unchanged.
offset_sens <- hits |> filter(tested_protein) |>
  left_join(select(track, contrast, offset = median_signed_adjusted), by = "contrast") |>
  mutate(adj_minus_offset = adjusted - offset,
         is_res_main  = category == "residual adjusted change",
         is_res_after = abs(adj_minus_offset) >= ADJ_CUT & abs(adj_minus_offset) >= SE_MULT * se_adj,
         label = lab_site(pick(everything()))) |>
  group_by(contrast) |>
  # per-site flags first (lost/gained), counts last, with NEW names: reusing a
  # column name inside summarise() would replace the flag by its sum.
  summarise(lost   = paste(label[is_res_main & !is_res_after], collapse = "; "),
            gained = paste(label[!is_res_main & is_res_after], collapse = "; "),
            offset = first(offset),
            n_residual_main = sum(is_res_main), n_residual_after_offset = sum(is_res_after),
            .groups = "drop") |>
  relocate(lost, gained, .after = last_col())
write_csv(offset_sens, file.path(out, "offset_sensitivity_mval.csv"))
print(offset_sens)
note(sprintf("offset sensitivity: residual site hits %d (main) -> %d after removing the per-contrast median signed offset",
             sum(offset_sens$n_residual_main), sum(offset_sens$n_residual_after_offset)))

# Sensitivity: same site, same contrast, other search -> same category?
sens <- matched[[PRIMARY_SEARCH]] |> filter(hit_site) |> select(site_id, contrast, category) |>
  left_join(matched[[OTHER_SEARCH]] |> select(site_id, contrast, category_other = category),
            by = c("site_id", "contrast"))
note(sprintf("%s site hits: %d; same category in %s: %d; different: %d; site not supported in %s: %d",
             PRIMARY_SEARCH, nrow(sens), OTHER_SEARCH, sum(sens$category == sens$category_other, na.rm = TRUE),
             sum(sens$category != sens$category_other, na.rm = TRUE), OTHER_SEARCH, sum(is.na(sens$category_other))))

# ---- 2. one example by hand (coding rules section 12: check a number yourself) ----
ex <- hits |> filter(category != "unresolved (protein not tested)") |> slice_max(abs(adjusted), n = 1, with_ties = FALSE)
if (nrow(ex) == 1) {
  grp <- contrast_groups(mval_contrasts()[[as.character(ex$contrast)]])
  stopifnot(length(grp) == 2)   # the example is a simple A - B contrast
  gmean <- function(ab, id_col, id) ab |> filter(.data[[id_col]] == id, group %in% grp) |>
    group_by(group) |> summarise(m = mean(normalized_abundance, na.rm = TRUE)) |>
    tibble::deframe()
  ms <- gmean(site_ab, "site_id", ex$site_id); mpr <- gmean(prot_ab, "protein_key", ex$protein_key)
  by_hand <- (ms[[grp[1]]] - ms[[grp[2]]]) - (mpr[[grp[1]]] - mpr[[grp[2]]])
  note(sprintf("by hand: %s %s%d in %s: site %s - %s = %.3f, protein = %.3f, adjusted = %.3f (model: %.3f)",
               coalesce(ex$gene, ex$locus_tag), ex$residue, ex$position, ex$contrast, grp[1], grp[2],
               ms[[grp[1]]] - ms[[grp[2]]], mpr[[grp[1]]] - mpr[[grp[2]]], by_hand, ex$adjusted))
  check("hand-computed adjusted effect equals the model value", abs(by_hand - ex$adjusted) < 1e-6)
}

# ---- 3. figures ----------------------------------------------------------------------
cat_cols <- c("no site hit" = "grey75", "consistent with protein amount" = "#0072B2",
              "residual adjusted change" = "#D55E00", "unresolved (adjusted effect uncertain)" = "#CC79A7",
              "unresolved (protein not tested)" = "#E69F00")

# 4a: site vs protein effect, same contrast
pd <- mp |> filter(tested_protein) |> mutate(category = factor(category, levels = names(cat_cols)))
p <- ggplot(pd, aes(diff_protein, diff_site)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey40") +
  geom_abline(slope = 1, intercept = c(-ADJ_CUT, ADJ_CUT), colour = "grey40", linetype = "dotted") +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(data = filter(pd, category == "no site hit"), colour = "grey75", size = 1, alpha = 0.7) +
  geom_point(data = filter(pd, category != "no site hit"), aes(colour = category), size = 2) +
  ggrepel::geom_text_repel(data = filter(pd, category == "residual adjusted change"),
                           aes(label = lab_site(filter(pd, category == "residual adjusted change"))),
                           size = 2.4, max.overlaps = 20, min.segment.length = 0, seed = 1) +
  facet_wrap(~contrast, ncol = 3) +
  scale_colour_manual(values = cat_cols, drop = FALSE) +
  labs(x = "protein effect, total proteome (log2, A - B)", y = "phosphosite effect (log2, A - B)", colour = NULL,
       title = sprintf("Step 9: does the site change follow the protein change? (%s search)", PRIMARY_SEARCH),
       subtitle = "Solid line: site = protein (change fully explained by protein amount). Dotted: +/- 1 log2.\nGrey = supported site, not a site hit. Only sites whose protein was tested in the same contrast.") +
  theme(legend.position = "bottom", strip.text = element_text(size = 7)) +
  guides(colour = guide_legend(nrow = 2))
ggsave(file.path(out, "fig_site_vs_protein_mval.png"), p, width = 10, height = 11, dpi = 150)

# 4b: adjusted effect +/- 2 se_adj for every site hit (with a tested protein)
hd <- hits |> filter(tested_protein) |>
  mutate(label = paste0(lab_site(pick(everything())), "  [", contrast, "]"),
         label = forcats::fct_reorder(label, adjusted),
         category = factor(category, levels = names(cat_cols)))
if (nrow(hd) > 0) {
  p <- ggplot(hd, aes(adjusted, label, colour = category)) +
    geom_vline(xintercept = 0, colour = "grey40") +
    geom_vline(xintercept = c(-ADJ_CUT, ADJ_CUT), linetype = "dotted", colour = "grey40") +
    geom_errorbar(aes(xmin = adjusted - SE_MULT * se_adj, xmax = adjusted + SE_MULT * se_adj),
                  width = 0.3, orientation = "y") +
    geom_point(data = filter(hd, fragile_spread_gt3), aes(x = adjusted), shape = 4, colour = "black", size = 3.2) +
    geom_point(aes(x = diff_site), shape = 1, colour = "grey45", size = 1.8) +
    geom_point(size = 2.2) +
    scale_colour_manual(values = cat_cols, drop = TRUE) +
    labs(x = "log2 (A - B):  filled = adjusted (site - protein) +/- 2 approx. SE;  open = unadjusted site effect",
         y = NULL, colour = NULL,
         title = "Step 9: protein-adjusted phosphosite effects of all site hits",
         subtitle = "SE assumes independent site and protein estimates (covariance not estimable). Descriptive, no p-value.\nx = QC flag: a group of this contrast has a within-group range > 3 log2 (fragile).") +
    theme(legend.position = "bottom", axis.text.y = element_text(size = 6.5)) +
    guides(colour = guide_legend(nrow = 2))
  ggsave(file.path(out, "fig_adjusted_effects_mval.png"), p, width = 10, height = max(4, 1.5 + 0.17 * nrow(hd)), dpi = 150)
}

# 4c: per-sample values of BOTH layers for a small module: the 5 residual sites
# with the largest |adjusted| that are NOT QC-flagged (fragile ones are marked in
# 4b), the strongest "consistent with protein amount" example, and the Stk1 sites
# (link to H1: does the site loss follow the protein loss?).
stk1 <- KEY_PROTEINS$locus_tag[KEY_PROTEINS$short_name == "Stk1"]
stopifnot(length(stk1) == 1)
show <- bind_rows(
  hits |> filter(category == "residual adjusted change", !fragile_spread_gt3) |>
    distinct(site_id, .keep_all = TRUE) |> slice_max(abs(adjusted), n = 5, with_ties = FALSE),
  hits |> filter(category == "consistent with protein amount") |> slice_max(abs(diff_site), n = 1, with_ties = FALSE),
  hits |> filter(locus_tag == stk1) |> slice_max(abs(diff_site), n = 1, with_ties = FALSE)) |>
  distinct(site_id, .keep_all = TRUE)
if (nrow(show) > 0) {
  show <- show |> mutate(label = paste0(lab_site(pick(everything())), "\n", contrast, "\n", category))
  primary_samples <- \(x) filter(x, medium %in% c("TSB", "SN"))
  pl <- bind_rows(
    site_ab |> primary_samples() |> inner_join(select(show, site_id, label), by = "site_id") |>
      mutate(layer = "phosphosite"),
    # two shown sites can sit on the same protein: its values are then drawn under both
    prot_ab |> primary_samples() |>
      inner_join(select(show, protein_key, label), by = "protein_key", relationship = "many-to-many") |>
      mutate(layer = "total protein")) |>
    add_x_primary() |>
    mutate(label = factor(label, levels = show$label))
  ndet <- pl |> group_by(layer, label, genotype, medium, x_pos) |>
    summarise(n = sum(!is.na(normalized_abundance)), .groups = "drop")
  p <- ggplot(pl, aes(x_pos, normalized_abundance, colour = genotype)) +
    geom_point(position = position_jitter(width = 0.04, height = 0, seed = 1), na.rm = TRUE, size = 1.6) +
    geom_text(data = ndet, aes(y = Inf, label = paste0(n, "/3")), vjust = 1.3, size = 2.2, show.legend = FALSE) +
    # facet_wrap (not facet_grid) so EVERY panel gets its own y-range: row 1 =
    # phosphosite, row 2 = its protein. facet_grid would share y along a row.
    facet_wrap(~ layer + label, ncol = nrow(show), scales = "free_y",
               labeller = labeller(.multi_line = FALSE)) +
    primary_x_axis() +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    scale_colour_manual(values = GENOTYPE_COLOURS) +
    labs(x = "medium", y = "normalised abundance (log2)",
         title = "Step 9: phosphosite and its protein, per sample (18 primary samples)",
         caption = "Top row: phosphosite; bottom row: its total protein. Numbers = samples with a value (of 3). Every panel has its own y-scale; compare the PATTERN, not the level.") +
    theme(strip.text = element_text(size = 6.5))
  ggsave(file.path(out, "fig_paired_layers_mval.png"), p, width = max(8, 2.2 * nrow(show)), height = 6.5, dpi = 150)
}

writeLines(c(sprintf("Step 9 checks, %s", Sys.Date()), checks), file.path(out, "site_vs_protein_checks_mval.txt"))
message("done: ", out)
