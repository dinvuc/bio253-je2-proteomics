# Step 9c: refined phosphosite hits + N-to-C expression plots closer to Fig. 1 of
# Wolski, Dittmann, Panse, Kunz & Grossmann, Chapter 8 (ISBN 978-1-0716-5519-1).
#
# Run from the repository root:   Rscript scripts/09c_refined_hits_n_to_c_mval.R   (~1 min)
# Needs scripts/05, 08 and 09 first.
#
# ---- Plan (coding rules, section 4) -----------------------------------------
# Aim:      (1) Which of the 62 phosphosite hits (Step 8) are robust? (2) Re-draw the
#           N-to-C plots so they follow the chapter's figure more closely, for the
#           refined-hit proteins plus the key proteins of fig_paired_layers_mval.png.
# Input:    data/derived/limma_treat_phospho_STY_mval.tsv.gz / _STYH_ (Step 8 refits)
#           data/derived/limma_treat_total_mval.tsv.gz (Step 5)
#           read_abundances("phospho_STY") (for the leave-one-out refits), read_samples()
#           step09 site_hit_categories_mval.csv (category, QC flag), offset_sensitivity
# Refined hit (SENSITIVITY labels; the Step 8 hits themselves stay as they are). A Step 8
# site hit is "refined" when it passes ALL of:
#   a. leave-one-sample-out: still a hit (same model, treat lfc = 1, BH over supported
#      sites) in all 18 refits that each drop one primary sample (as Step 11 did for proteins)
#   b. search: also a hit in the STYH search (same site_id, same contrast)
#   c. QC: no group of the contrast with a within-group range > 3 log2 (Step 9 flag)
#   Reported next to it, not used as a filter: the Step 9 category (consistent with
#   protein / residual / unresolved) and whether a residual call survives the layer-offset
#   sensitivity.
# Figures (closer to the chapter): fixed 3 x 2 grid for every protein, so all proteins
#   are comparable (the chapter also shows the same contrasts for every protein):
#   left = genotype in TSB (stk_vs_anc, stp_vs_anc, stp_vs_stk), right = PASN response
#   within each strain (TSB_vs_SN at anc, stk, stp). Overall title "Protein: <id> (length: L)";
#   panel title "Prot : <id>; length: L; # sites:n; # not localized sites:m"; y "diff : <contrast>";
#   x "posInProtein"; legends "modAA", localisation status, "Rectangle: diff of protein";
#   plain lines (no end points); asterisks above the line end: ** FDR < 0.01, * FDR < 0.05.
#   Our data: solid = localised, dashed = not localised (the chapter: observed/imputed;
#   we do not impute). Refined hits get a black outline ring around the asterisks.
# Proteins: the proteins of fig_paired_layers_mval.png that also have >= 1 refined hit
#   (glyA, traP, B7H15_04690, Stk1; deoB and pgcA have no refined hit). The paired-layer
#   figure is copied into the folder. Plus, on group request (18 Nov), kdpD (B7H15_11550),
#   a Step 6 protein interaction hit (stk vs anc) without any site hit.
# Output:   mval_medi/results/step09_site_vs_protein/n_to_c_refined/
#             refined_hits_mval.csv          all 62 hits with checks a-c and the verdict
#             refined_summary_mval.csv       counts per contrast
#             fig_n_to_c_refined_<locus>_<gene>_mval.png, n_to_c_refined_all_mval.pdf
#             figure_guide_mval.md
# Check:    18 LOSO fits; full refit reproduces the Step 8 hits; every refined hit is drawn.

suppressMessages(source("R/setup.R"))
library(patchwork)
out <- "mval_medi/results/step09_site_vs_protein/n_to_c_refined"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
check <- function(label, ok) { message(sprintf("[%s] %s", if (isTRUE(ok)) "PASS" else "FAIL", label)); invisible(ok) }

sty  <- read_tsv("data/derived/limma_treat_phospho_STY_mval.tsv.gz", show_col_types = FALSE)
styh <- read_tsv("data/derived/limma_treat_phospho_STYH_mval.tsv.gz", show_col_types = FALSE)
prot <- read_tsv("data/derived/limma_treat_total_mval.tsv.gz", show_col_types = FALSE) |>
  select(protein_key, contrast, prot_diff = diff, prot_tested = tested, prot_FDR = FDR)
cat9 <- read_csv("mval_medi/results/step09_site_vs_protein/site_hit_categories_mval.csv", show_col_types = FALSE)
meta <- read_protein_metadata() |>
  transmute(locus_tag, aa_length, gene = if_else(ortholog_quality %in% "good", gene_symbol_best, NA_character_)) |>
  distinct(locus_tag, .keep_all = TRUE)
samples <- read_samples() |> arrange(genotype, medium, replicate)

# ---- 1. leave-one-sample-out on the site layer ---------------------------------------
supported <- sty |> filter(supported) |> distinct(site_id, contrast)
ab <- read_abundances("phospho_STY") |> keep_je2() |>
  transmute(protein_key = site_id, sample, normalized_abundance)   # one row = one site
m_all <- abundance_matrix(ab, samples)
loso_fit <- function(smp) {
  fit_limma_treat(m_all[, smp$sample], smp) |> rename(site_id = protein_key) |>
    inner_join(supported, by = c("site_id", "contrast")) |>      # same family as Step 8
    mutate(p = if_else(tested, p.value, NA_real_)) |>
    group_by(contrast) |> mutate(FDR2 = p.adjust(p, "BH")) |> ungroup() |>
    mutate(hit2 = !is.na(FDR2) & FDR2 < MVAL_FDR)
}
full <- loso_fit(samples)
check("full refit reproduces the 62 Step 8 hits",
      setequal(paste(full$site_id, full$contrast)[full$hit2], paste(cat9$site_id, cat9$contrast)))
prim <- samples$sample[samples$medium %in% c("TSB", "SN")]
loso <- purrr::map(prim, \(d) loso_fit(filter(samples, sample != d)) |>
                     semi_join(cat9, by = c("site_id", "contrast")) |>
                     transmute(site_id, contrast, dropped = d, hit2)) |> bind_rows()
check("18 leave-one-out refits", n_distinct(loso$dropped) == 18)
loso_s <- loso |> group_by(site_id, contrast) |>
  summarise(loso_still_hit = sum(hit2), loso_lost_when_dropping = paste(dropped[!hit2], collapse = "; "), .groups = "drop")

# ---- 2. assemble the checks ---------------------------------------------------------
off <- read_csv("mval_medi/results/step09_site_vs_protein/site_vs_protein_tracking_mval.csv", show_col_types = FALSE) |>
  select(contrast, offset = median_signed_adjusted)
refined <- cat9 |>
  select(site_id, protein_key, locus_tag, gene, product, residue, position, contrast, diff_site, FDR_site,
         diff_protein, adjusted, se_adj, category, fragile_spread_gt3, detected_per_group) |>
  left_join(loso_s, by = c("site_id", "contrast")) |>
  left_join(styh |> select(site_id, contrast, styh_supported = supported, styh_hit = hit), by = c("site_id", "contrast")) |>
  left_join(off, by = "contrast") |>
  mutate(styh_hit = coalesce(styh_hit, FALSE), styh_supported = coalesce(styh_supported, FALSE),
         check_loso = loso_still_hit == 18, check_styh = styh_hit, check_qc = !fragile_spread_gt3,
         refined = check_loso & check_styh & check_qc,
         failed = purrr::pmap_chr(list(check_loso, check_styh, check_qc, loso_still_hit, styh_supported), \(a, b, c, n, s)
           paste(c(if (!a) sprintf("leave-one-out %d/18", n), if (!b) (if (s) "not a hit in STYH" else "not supported in STYH"),
                   if (!c) "QC: within-group range > 3"), collapse = "; ")),
         residual_survives_offset = if_else(category == "residual adjusted change",
                                            abs(adjusted - offset) >= 1 & abs(adjusted - offset) >= 2 * se_adj, NA),
         direction_in_words = purrr::map2_chr(contrast, diff_site, \(ct, d) {
           g <- str_match(ct, "^([A-Za-z]+)_vs_([A-Za-z]+)_at_(.+)$")
           if (str_starts(ct, "interaction")) return(if (d > 0) "TSB-vs-SN difference larger in the first strain" else "larger in the second strain")
           paste0("site higher in ", if (d > 0) g[2] else g[3], " than ", if (d > 0) g[3] else g[2], " (", g[4], ")")
         })) |>
  arrange(desc(refined), contrast, desc(abs(diff_site)))
write_csv(refined, file.path(out, "refined_hits_mval.csv"))
summ <- refined |> group_by(contrast) |>
  # counts that use the per-site flag `refined` come FIRST; n_refined gets a new name
  # (reusing "refined" would replace the flag by its sum for the later lines)
  summarise(refined_residual = sum(refined & category == "residual adjusted change"),
            refined_consistent = sum(refined & category == "consistent with protein amount"),
            refined_unresolved = sum(refined & str_starts(category, "unresolved")),
            site_hits = n(), n_refined = sum(refined), lost_loso = sum(!check_loso), lost_styh = sum(!check_styh),
            lost_qc = sum(!check_qc), .groups = "drop") |>
  relocate(site_hits, n_refined, .after = contrast)
write_csv(summ, file.path(out, "refined_summary_mval.csv"))
print(summ, width = 200)
message(sprintf("refined: %d of %d site hits (%d sites, %d proteins)", sum(refined$refined), nrow(refined),
                n_distinct(refined$site_id[refined$refined]), n_distinct(refined$locus_tag[refined$refined])))
print(refined |> filter(refined) |> transmute(contrast, site = paste0(coalesce(gene, locus_tag), " ", residue, position),
                                               diff_site = round(diff_site, 2), category), n = 70)
print(refined |> filter(!refined) |> transmute(contrast, site = paste0(coalesce(gene, locus_tag), " ", residue, position), failed), n = 70)


# ---- 3. N-to-C plots following the chapter's Fig. 1 -----------------------------------
# Style taken from prophosqua (github.com/prolfqua/prophosqua, R/plotNtoC.R, MIT licence;
# the package itself could not be installed here: its dependencies need Rtools) and from
# the printed Fig. 1 of the chapter:
#   theme_minimal; black zero segment from 0 to the protein length; x limited to the protein;
#   yellow rectangle (alpha 0.3) from 0 to the protein log2FC = "diff of protein";
#   "N"/"C" written at the rectangle height (hjust = 0, vjust = 0);
#   one bar per site from 0 to the site log2FC, colour = modAA with prophosqua's colours
#   (S #0072B2, T #009E73, Y #D55E00, not-localised sites as modAA "NotLoc" #CC79A7);
#   linetype = imputation_status (we never impute, so every site is "observed", solid);
#   red asterisks at the bar tip: ** FDR < 0.01, * FDR < 0.05 (chapter caption);
#   panel title "Prot : <id>; length: L; # sites:n; # not localized sites:m";
#   y "diff : <contrast>", x "posInProtein"; a legend next to every panel; free y per panel;
#   no protein estimate -> red dashed box "No estimate for diff of protein" (as prophosqua).
# Our additions: refined hits keep red asterisks, hits that failed a robustness check get
# faded asterisks; sites too sparse to test are not drawn (prophosqua plots only estimated sites).
KEY_PAIRED <- c("B7H15_00775", "B7H15_11735", "B7H15_13850", "B7H15_10065", "B7H15_04690", "B7H15_06355")  # fig_paired_layers
GRID <- c("stk_vs_anc_at_TSB", "TSB_vs_SN_at_anc",      # row 1   (left: genotype in TSB,
          "stp_vs_anc_at_TSB", "TSB_vs_SN_at_stk",      # row 2    right: PASN response
          "stp_vs_stk_at_TSB", "TSB_vs_SN_at_stp")      # row 3    within each strain)
MOD_AA_COLORS <- c(S = "#0072B2", T = "#009E73", Y = "#D55E00", NotLoc = "#CC79A7")   # prophosqua
hit_key <- refined |> select(site_id, contrast, refined)

panel <- function(pk, ct, len, id) {
  d <- sty |> filter(protein_key == pk, contrast == ct, tested) |>
    left_join(hit_key, by = c("site_id", "contrast")) |>
    mutate(modAA = factor(if_else(supported, residue, "NotLoc"), levels = names(MOD_AA_COLORS)),
           imputation_status = "observed",
           significance = case_when(hit & FDR < 0.01 ~ "**", hit & FDR < 0.05 ~ "*", TRUE ~ ""),
           robust = if_else(coalesce(refined, FALSE), "refined", "failed a check"),
           vj = if_else(diff >= 0, -0.1, 1.1))
  p <- prot |> filter(protein_key == pk, contrast == ct)
  dprot <- if (nrow(p) == 1 && p$prot_tested) p$prot_diff else NA_real_
  title <- sprintf("Prot : %s; length: %d; # sites:%d; # not localized sites:%d",
                   id, len, nrow(d), sum(d$modAA == "NotLoc"))
  g <- ggplot(d) +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(size = 8.5), axis.title = element_text(size = 8),
          legend.title = element_text(size = 8), legend.text = element_text(size = 7),
          legend.key.size = unit(0.4, "cm"))
  # Layer order as in the printed figure: rectangle FIRST, bars on top (drawing the
  # semi-transparent rectangle over the bars tinted them, e.g. blue looked green).
  if (!is.na(dprot)) {
    g <- g +
      geom_rect(data = data.frame(xmin = 0, xmax = len, ymin = 0, ymax = dprot, fill = "diff of protein"),
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), alpha = 0.3, inherit.aes = FALSE) +
      scale_fill_manual(values = c("diff of protein" = "yellow"), name = "Rectangle") +
      annotate("text", x = 0, y = dprot, label = "N", vjust = 0, hjust = 0, size = 3) +
      annotate("text", x = len, y = dprot, label = "C", vjust = 0, hjust = 0, size = 3)
  } else {
    yext <- max(abs(d$diff), 1, na.rm = TRUE)
    g <- g +
      annotate("rect", xmin = 0, xmax = len, ymin = -yext / 2, ymax = yext / 2, fill = NA, colour = "red", linetype = "dashed") +
      annotate("text", x = len / 2, y = 0, label = "No estimate for diff of protein", colour = "red", angle = 20, size = 3.5)
  }
  g <- g +
    geom_segment(aes(x = position, xend = position, y = 0, yend = diff, colour = modAA, linetype = imputation_status),
                 linewidth = 1.1) +
    geom_text(data = \(x) filter(x, significance != ""),
              aes(position, diff, label = significance, alpha = robust, vjust = vj), colour = "red", size = 4.2,
              show.legend = FALSE) +
    annotate("segment", x = 0, xend = len, y = 0, yend = 0, colour = "black", linewidth = 0.8) +
    # legends list only what the panel contains (as in the printed figure: "S, T", "observed")
    scale_colour_manual(values = MOD_AA_COLORS, name = "modAA") +
    scale_linetype_manual(values = c(observed = "solid", imputed = "dashed"), name = "imputation_status") +
    scale_alpha_manual(values = c(refined = 1, `failed a check` = 0.4), guide = "none") +
    scale_x_continuous(limits = c(0, len)) +
    scale_y_continuous(expand = expansion(mult = 0.15)) +
    labs(x = "posInProtein", y = paste0("diff : ", ct), title = title) +
    guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 2),
           fill = guide_legend(order = 3))   # legend order as in the printed figure
  if (nrow(d) == 0)
    g <- g + annotate("text", x = len / 2, y = -0.25, label = "no site tested in this contrast", colour = "grey40", size = 3)
  g
}

# Plotted = IMPORTANT and SIGNIFICANT (group request, 17 Nov): proteins of
# fig_paired_layers_mval.png that also have >= 1 refined site hit. deoB and pgcA are in
# the paired-layer figure but have no refined hit, so they are not plotted any more.
# Added on request (group, 18 Nov): kdpD, a formal protein interaction hit of Step 6
# (fig_response_vs_anc_mval.png, interaction_stk_vs_anc_at_TSB_vs_SN). It has no site hit,
# so it is drawn for its protein rectangles; its phosphosites are mostly too sparse to test.
EXTRA <- tibble(locus_tag = "B7H15_11550",
                why = "Step 6 interaction hit (stk vs anc), group request")
proteins <- tibble(locus_tag = intersect(KEY_PAIRED, unique(refined$locus_tag[refined$refined])),
                   why = "paired-layer figure + refined hit") |>
  bind_rows(EXTRA) |>
  left_join(distinct(sty, locus_tag, protein_key, product), by = "locus_tag") |>   # product: also for proteins without a site hit
  left_join(meta, by = "locus_tag") |>
  arrange(locus_tag)
check("every plotted protein is in the STY site table", !anyNA(proteins$protein_key))
check("one protein_key per plotted locus tag", !anyDuplicated(proteins$locus_tag))
write_csv(proteins |> select(locus_tag, gene, product, aa_length, why), file.path(out, "plotted_proteins_mval.csv"))
# remove figures of proteins that are no longer selected (the folder holds only the current set)
unlink(list.files(out, pattern = "^fig_n_to_c_refined_.*\\.png$", full.names = TRUE))
# the paired-layer figure (Step 9) the selection is based on, next to the N-to-C plots
check("paired-layer figure copied", file.copy("mval_medi/results/step09_site_vs_protein/fig_paired_layers_mval.png",
                                              file.path(out, "fig_paired_layers_mval.png"), overwrite = TRUE))

pdf(file.path(out, "n_to_c_refined_all_mval.pdf"), width = 14, height = 11)
for (i in seq_len(nrow(proteins))) {
  pr <- proteins[i, ]
  id <- paste0(pr$locus_tag, if (!is.na(pr$gene)) paste0(" (", pr$gene, ")") else "")
  plots <- purrr::map(GRID, \(ct) panel(pr$protein_key, ct, pr$aa_length, id))
  n_ref <- sum(refined$refined & refined$locus_tag == pr$locus_tag)
  fig <- wrap_plots(plots, ncol = 2, byrow = TRUE) +
    plot_annotation(
      title = sprintf("Protein: %s (length: %d)", id, pr$aa_length),
      subtitle = sprintf("%s  |  %s  |  refined site hits: %d", str_trunc(coalesce(pr$product, ""), 60), pr$why, n_ref),
      caption = paste("Left: genotype difference in TSB (A - B as named). Right: PASN response TSB - SN within each strain (> 0 = lower in SN).",
                      "Yellow rectangle = protein log2FC (total proteome); bars = phosphosite log2FC (STY); NotLoc = site not localised (median loc_prob < 0.75).",
                      "** FDR < 0.01, * FDR < 0.05 (site test, treat |log2FC| > 1). Bright red = refined hit (leave-one-out 18/18, also a hit in STYH, no QC flag);",
                      "faded = hit that failed a check. Every estimate is observed (no imputation). Each panel has its own y-scale (as in the chapter).",
                      "Style after Wolski, Dittmann, Panse, Kunz & Grossmann, ch. 8, Fig. 1 and prophosqua (MIT).", sep = "\n"),
      theme = theme(plot.title = element_text(size = 12, face = "bold"), plot.subtitle = element_text(size = 9),
                    plot.caption = element_text(size = 7.5)))
  ggsave(file.path(out, sprintf("fig_n_to_c_refined_%s_%s_mval.png", pr$locus_tag, coalesce(pr$gene, "nogene"))),
         fig, width = 14, height = 11, dpi = 200, bg = "white")
  print(fig)
}
invisible(dev.off())
check(sprintf("%d protein figures written", nrow(proteins)),
      length(list.files(out, pattern = "^fig_n_to_c_refined_.*png$")) == nrow(proteins))
check("every refined hit lies in one of the 6 grid contrasts (so it is drawn)",
      all(refined$contrast[refined$refined] %in% GRID))
message("done: ", out)
