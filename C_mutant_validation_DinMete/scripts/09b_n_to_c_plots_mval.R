# Step 9b (PROJECT_WORKFLOW_mval.md, Step 9 figures): N-to-C expression plots of the
# proteins with a phosphosite hit, in the style of Fig. 1 of Wolski, Dittmann, Panse,
# Kunz & Grossmann, "Integrated Analysis of Posttranslational Modifications and Total
# Proteome: Methods for Distinguishing Abundance from Usage Changes" (Chapter 8,
# book ISBN 978-1-0716-5519-1; there made with prophosqua::n_to_c_expression_multicontrast).
#
# Run from the repository root:   Rscript C_mutant_validation_DinMete/scripts/09b_n_to_c_plots_mval.R   (~30 s)
# Needs C_mutant_validation_DinMete/scripts/05, 08 and 09 first.
#
# ---- Plan (coding rules, section 4) -----------------------------------------
# Aim:      For every protein in fig_adjusted_effects_mval.png (= every protein with a
#           phosphosite hit, Step 8/9), show along the protein sequence whether the site
#           changes go beyond the change of the protein itself, in every contrast in
#           which the protein has a site hit (all those contrasts together, one panel each).
# Input:    data/derived/limma_treat_phospho_STY_mval.tsv.gz  (Step 8 site refit, STY)
#           data/derived/limma_treat_total_mval.tsv.gz        (Step 5 protein refit)
#           C_mutant_validation_DinMete/results/step09_site_vs_protein/site_hit_categories_mval.csv (which proteins/contrasts)
#           read_protein_metadata(): aa_length (protein length), gene symbol (good orthologs only)
# Plot (as in the chapter's Fig. 1):
#           x = amino-acid position N -> C (0 .. protein length); y = log2 fold change (A - B).
#           Light-yellow rectangle over the whole protein, height = protein log2FC in the
#           same contrast ("diff of protein"); vertical lines = sites at their position,
#           height = site log2FC, colour = residue (S/T/Y); red asterisks: * FDR < 0.05,
#           ** FDR < 0.01 (our site test, Step 8: treat lfc = 1, BH over supported sites).
#           Adapted: the chapter uses solid/dashed for observed/imputed. Our refit has no
#           imputation, so solid = localised (supported, median per-sample loc_prob >= 0.75),
#           dashed = tested but not localised. Sites too sparse to test in that contrast
#           (< 2 of 3 in a group) = grey tick on the zero line. Same y-scale in all panels of
#           one protein. This is a DPA view (raw site change); the protein-adjusted (DPU)
#           view is fig_adjusted_effects_mval.png.
# Output:   C_mutant_validation_DinMete/results/step09_site_vs_protein/n_to_c/
#             fig_n_to_c_<locus_tag>_<gene>_mval.png  (one per protein)
#             n_to_c_all_proteins_mval.pdf             (all proteins, one page each)
# Check:    every protein of the hit table gets a figure; each panel's hit sites appear
#           with an asterisk; positions lie within 1 .. protein length.

suppressMessages(source("R/setup.R"))
library(patchwork)
out <- "C_mutant_validation_DinMete/results/step09_site_vs_protein/n_to_c"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
check <- function(label, ok) { message(sprintf("[%s] %s", if (isTRUE(ok)) "PASS" else "FAIL", label)); invisible(ok) }

hits  <- read_csv("C_mutant_validation_DinMete/results/step09_site_vs_protein/site_hit_categories_mval.csv", show_col_types = FALSE)
sites <- read_tsv("data/derived/limma_treat_phospho_STY_mval.tsv.gz", show_col_types = FALSE)
prot  <- read_tsv("data/derived/limma_treat_total_mval.tsv.gz", show_col_types = FALSE) |>
  select(protein_key, contrast, prot_diff = diff, prot_tested = tested, prot_FDR = FDR)
meta  <- read_protein_metadata() |>
  transmute(locus_tag, aa_length, gene = if_else(ortholog_quality %in% "good", gene_symbol_best, NA_character_)) |>
  distinct(locus_tag, .keep_all = TRUE)

panels <- hits |> distinct(protein_key, locus_tag, contrast) |>
  left_join(meta, by = "locus_tag")
check("every hit protein has a protein length", !anyNA(panels$aa_length))
PRIMARY_ORDER <- c("stk_vs_anc_at_TSB", "stp_vs_anc_at_TSB", "stp_vs_stk_at_TSB", "TSB_vs_SN_at_anc",
                   "TSB_vs_SN_at_stk", "TSB_vs_SN_at_stp", "interaction_stk_vs_anc_at_TSB_vs_SN",
                   "interaction_stp_vs_anc_at_TSB_vs_SN", "interaction_stp_vs_stk_at_TSB_vs_SN")

# direction in words for the panel label (CLAUDE.md rule 3)
in_words <- function(ct) {
  if (str_starts(ct, "interaction_")) {
    g <- str_match(ct, "^interaction_([a-z]+)_vs_([a-z]+)_at_")
    return(paste0("> 0: TSB-vs-SN difference larger in ", g[2], " than in ", g[3]))
  }
  g <- str_match(ct, "^([A-Za-z]+)_vs_([A-Za-z]+)_at_([A-Za-z]+)$")
  paste0("> 0: higher in ", g[2], " than in ", g[3], " (", g[4], ")")
}
RES_COL <- c(S = "#0072B2", T = "#E69F00", Y = "#009E73")

one_panel <- function(pk, ct, len, label, ylim) {
  s <- sites |> filter(protein_key == pk, contrast == ct)
  tested <- s |> filter(tested) |>
    mutate(status = if_else(supported, "localised (median loc_prob >= 0.75)", "not localised"),
           stars = case_when(supported & FDR < 0.01 ~ "**", supported & FDR < 0.05 ~ "*", TRUE ~ ""))
  sparse <- s |> filter(!tested)
  p <- prot |> filter(protein_key == pk, contrast == ct)
  rect_h <- if (nrow(p) == 1 && p$prot_tested) p$prot_diff else NA_real_
  n_sites <- n_distinct(s$site_id); n_notloc <- n_distinct(s$site_id[!s$localised])
  g <- ggplot() + geom_hline(yintercept = 0, colour = "grey60")
  if (!is.na(rect_h)) {
    g <- g + geom_rect(aes(xmin = 0, xmax = len, ymin = 0, ymax = rect_h, linetype = "diff of protein"),
                       fill = "#FFF2B3", colour = "#E0C200", alpha = 0.8) +
      annotate("text", x = 0, y = rect_h, label = "N", hjust = 1.4, vjust = 0.5, fontface = "bold") +
      annotate("text", x = len, y = rect_h, label = "C", hjust = -0.4, vjust = 0.5, fontface = "bold")
  }
  # Invisible dummy for every residue x line type, far outside the visible y-range
  # (coord_cartesian clips it): every panel then has the same legend keys, so
  # patchwork merges the legends into one instead of showing one per panel.
  far <- ylim[2] + 100 * diff(ylim)
  dummy <- tidyr::expand_grid(residue = names(RES_COL),
                              status = c("localised (median loc_prob >= 0.75)", "not localised")) |>
    mutate(position = len / 2, diff = far)
  g <- g +
    geom_segment(data = dummy, aes(x = position, xend = position, y = far, yend = diff + 1, colour = residue,
                                   linetype = status), linewidth = 0.9) +
    geom_point(data = dummy, aes(position, diff, colour = residue), size = 1.6) +
    geom_segment(data = tested, aes(x = position, xend = position, y = 0, yend = diff, colour = residue,
                                    linetype = status), linewidth = 0.9) +
    geom_point(data = tested, aes(position, diff, colour = residue), size = 1.6) +
    geom_text(data = filter(tested, stars != ""),
              aes(position, diff + sign(diff) * 0.04 * diff(ylim), label = stars),
              colour = "red", size = 6, fontface = "bold", vjust = 0.75) +
    geom_point(data = sparse, aes(position, 0), shape = 124, size = 3, colour = "grey45") +
    # fixed limits: every panel gets the SAME legend, so patchwork can merge them into one
    scale_colour_manual(values = RES_COL, limits = names(RES_COL), name = "residue (modAA)") +
    scale_linetype_manual(values = c("localised (median loc_prob >= 0.75)" = "solid", "not localised" = "dashed",
                                     "diff of protein" = "solid"),
                          limits = c("localised (median loc_prob >= 0.75)", "not localised", "diff of protein"),
                          breaks = c("localised (median loc_prob >= 0.75)", "not localised"), name = "site") +
    scale_x_continuous(limits = c(-0.06 * len, 1.06 * len), breaks = scales::pretty_breaks(5)) +
    coord_cartesian(ylim = ylim) +
    labs(x = "position in protein (amino acid)", y = paste0("log2 FC\n", ct),
         title = sprintf("Prot: %s; length: %d; # sites: %d; # not localized sites: %d", label, len, n_sites, n_notloc),
         subtitle = paste0(in_words(ct),
                           if (is.na(rect_h)) "   |   protein not tested in this contrast (no rectangle)"
                           else sprintf("   |   protein: %+.2f (FDR %s)", rect_h, signif(p$prot_FDR, 2)))) +
    theme_bw(base_size = 11) + theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 9))
  g
}

proteins <- panels |> distinct(protein_key, locus_tag, gene, aa_length) |> arrange(locus_tag)
pdf(file.path(out, "n_to_c_all_proteins_mval.pdf"), width = 12, height = 10)
n_done <- 0
for (i in seq_len(nrow(proteins))) {
  pr  <- proteins[i, ]
  cts <- panels |> filter(protein_key == pr$protein_key) |> pull(contrast) |> unique()
  cts <- cts[order(match(cts, PRIMARY_ORDER))]
  # common y-range over all panels of this protein: sites, protein, 0
  vals <- c(0, sites$diff[sites$protein_key == pr$protein_key & sites$contrast %in% cts & sites$tested],
            prot$prot_diff[prot$protein_key == pr$protein_key & prot$contrast %in% cts & prot$prot_tested])
  ylim <- range(vals, na.rm = TRUE) + c(-0.15, 0.15) * max(1, diff(range(vals, na.rm = TRUE)))
  label <- paste0(pr$locus_tag, if (!is.na(pr$gene)) paste0(" (", pr$gene, ")") else "")
  # positions inside the protein?
  pos <- sites$position[sites$protein_key == pr$protein_key]
  if (any(pos < 1 | pos > pr$aa_length)) message("[FAIL] site position outside protein: ", pr$locus_tag)
  # every hit of this protein must carry an asterisk
  h <- hits |> filter(protein_key == pr$protein_key)
  s_hit <- sites |> filter(protein_key == pr$protein_key, hit)
  if (!setequal(paste(h$site_id, h$contrast), paste(s_hit$site_id, s_hit$contrast)))
    message("[FAIL] hit sites do not match for ", pr$locus_tag)
  plots <- purrr::map(cts, \(ct) one_panel(pr$protein_key, ct, pr$aa_length, label, ylim))
  # one legend for the whole figure: keep it on the first panel only (it has every key
  # thanks to the dummy layer); patchwork would otherwise stack near-identical copies
  plots[-1] <- purrr::map(plots[-1], \(p) p + theme(legend.position = "none"))
  fig <- wrap_plots(plots, ncol = 1) + plot_layout(guides = "collect") +
    plot_annotation(title = paste0("N-to-C expression plot: ", label, " - ",
                                   str_trunc(unique(hits$product[hits$protein_key == pr$protein_key])[1], 60)),
                    caption = paste("Rectangle = protein log2FC (total proteome, same contrast); lines = phosphosite log2FC (STY search).",
                                    "* FDR < 0.05, ** FDR < 0.01 (site test, treat |log2FC| > 1).",
                                    "Solid = localised, dashed = not localised (the chapter uses solid/dashed for observed/imputed;",
                                    "we do not impute). Grey tick = site too sparse to test in that contrast. Style after Wolski et al., Fig. 1.",
                                    sep = "\n"))
  ggsave(file.path(out, sprintf("fig_n_to_c_%s_%s_mval.png", pr$locus_tag, coalesce(pr$gene, "nogene"))),
         fig, width = 11, height = 1.2 + 3.1 * length(cts), dpi = 200)
  print(fig)
  n_done <- n_done + 1
}
invisible(dev.off())
check(sprintf("one figure per hit protein (%d of %d)", n_done, nrow(proteins)), n_done == nrow(proteins))
message("done: ", out)
