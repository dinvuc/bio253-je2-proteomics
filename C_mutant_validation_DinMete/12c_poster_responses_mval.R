# Step 12c (poster): large-text, transparent versions of the Step 6 response scatter plots
# (fig_response_vs_anc_mval.png, fig_response_stk_vs_stp_mval.png) and one combined 3-panel figure.
# Run from the repository root:  Rscript C_mutant_validation_DinMete/12c_poster_responses_mval.R
# Needs data/derived/limma_treat_total_mval.tsv.gz (git-ignored; written by C_mutant_validation_DinMete/05_formal_test_mval.R).
# Aim:      PASN response in each mutant vs anc, and stp vs stk, at poster size.
# Input:    Step 5 limma table (TSB_vs_SN_at_* and the three interaction contrasts), read_protein_metadata()
# Analysis: none new; hit = treat lfc = 1, BH-FDR < 0.05; pink = exploratory 0.05 <= FDR < 0.1 (not a hit).
# Output:   C_mutant_validation_DinMete/results/step06_responses/fig_response_vs_anc_poster_mval.png,
#           fig_response_stk_vs_stp_poster_mval.png, fig_responses_combined_poster_mval.png (all transparent)
#           Combined figure: every panel exactly 3.6 x 3.6 in, square ranges (diagonal at 45 degrees).
# Check:    n / r / hits as in Step 6: 1,882 / 0.85 / 5 (stk), 1,931 / 0.86 / 1 (stp), 1,876 / 0.84 / 4 (stp vs stk).
suppressMessages(source("R/setup.R")); source("C_mutant_validation_DinMete/12c_poster_transparent_helper_mval.R")
out <- "C_mutant_validation_DinMete/results/step06_responses"
BS <- 17; TXT <- 5

res <- read_tsv("data/derived/limma_treat_total_mval.tsv.gz", show_col_types = FALSE)
meta <- read_protein_metadata() |> filter(ortholog_quality == "good") |>
  distinct(locus_tag, gene = gene_symbol_best)

wide <- res |>
  filter(contrast %in% c("TSB_vs_SN_at_anc", "TSB_vs_SN_at_stk", "TSB_vs_SN_at_stp",
                         "interaction_stk_vs_anc_at_TSB_vs_SN", "interaction_stp_vs_anc_at_TSB_vs_SN",
                         "interaction_stp_vs_stk_at_TSB_vs_SN")) |>
  select(protein_key, locus_tag, contrast, diff, FDR, tested) |>
  pivot_wider(names_from = contrast, values_from = c(diff, FDR, tested), names_glue = "{contrast}__{.value}")
resp <- wide |>
  transmute(protein_key, locus_tag,
            resp_anc = -`TSB_vs_SN_at_anc__diff`, resp_stk = -`TSB_vs_SN_at_stk__diff`,
            resp_stp = -`TSB_vs_SN_at_stp__diff`,
            int_stk_tested = `interaction_stk_vs_anc_at_TSB_vs_SN__tested`,
            int_stp_tested = `interaction_stp_vs_anc_at_TSB_vs_SN__tested`,
            int_stk_FDR = `interaction_stk_vs_anc_at_TSB_vs_SN__FDR`,
            int_stp_FDR = `interaction_stp_vs_anc_at_TSB_vs_SN__FDR`,
            int_direct_tested = `interaction_stp_vs_stk_at_TSB_vs_SN__tested`,
            int_direct_FDR = `interaction_stp_vs_stk_at_TSB_vs_SN__FDR`)

CL <- c("hit (FDR < 0.05)", "exploratory (FDR < 0.1)", "no evidence")
hit_class <- function(tested, fdr) factor(case_when(tested & fdr < MVAL_FDR ~ CL[1],
                                                    tested & fdr < 0.1 ~ CL[2], TRUE ~ CL[3]), CL)
CLASS_COLOURS <- setNames(c("#D55E00", "#CC79A7", "grey45"), CL)

count_text <- function(d) {
  d |> summarise(r = cor(x, y), n = n(), hits = sum(class == CL[1]), .groups = "drop") |>
    mutate(text = sprintf("r = %.2f\n%d %s of %s", r, hits, if_else(hits == 1, "hit", "hits"),
                          format(n, big.mark = ",")))
}

scatter <- function(d, x_lab, y_lab) {
  lab <- filter(d, class != "no evidence")
  ggplot(d, aes(x, y)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey35", linewidth = 0.7) +
    geom_abline(slope = 1, intercept = c(-MVAL_LFC, MVAL_LFC), colour = "grey60", linetype = 2, linewidth = 0.6) +
    geom_point(data = filter(d, class == "no evidence"), alpha = 0.3, size = 1.1, colour = "grey45") +
    geom_point(data = lab, aes(colour = class), size = 3.4) +
    geom_text(data = lab, aes(label = label, colour = class), size = TXT, vjust = -0.75,
              fontface = "bold", show.legend = FALSE) +
    geom_label(data = count_text(d) |> mutate(across(any_of("panel"), identity)),
               aes(x = -Inf, y = Inf, label = text), hjust = 0, vjust = 1, size = TXT,
               linewidth = 0, fill = NA, label.padding = unit(0.7, "lines"), lineheight = 0.95) +
    scale_colour_manual(values = CLASS_COLOURS, breaks = CL[1:2], drop = TRUE, name = NULL) +
    scale_x_continuous(breaks = seq(-8, 8, 4)) + scale_y_continuous(breaks = seq(-8, 8, 4)) +
    coord_equal(clip = "off") +
    labs(x = x_lab, y = y_lab) +
    guides(colour = guide_legend(override.aes = list(size = 4.5))) +
    theme_bw(base_size = BS) +
    theme(legend.position = "bottom", legend.text = element_text(size = BS),
          legend.margin = margin(-4, 0, 0, 0),
          axis.text = element_text(size = BS - 2, colour = "black"),
          strip.text = element_text(face = "bold", size = BS),
          panel.grid.minor = element_blank())
}

# Fig 2: each mutant against anc
sc <- bind_rows(
  resp |> filter(int_stk_tested) |> transmute(locus_tag, panel = "stk vs anc", x = resp_anc, y = resp_stk,
                                              class = hit_class(int_stk_tested, int_stk_FDR)),
  resp |> filter(int_stp_tested) |> transmute(locus_tag, panel = "stp vs anc", x = resp_anc, y = resp_stp,
                                              class = hit_class(int_stp_tested, int_stp_FDR))) |>
  left_join(meta, by = "locus_tag") |>
  mutate(label = coalesce(gene, str_remove(locus_tag, "B7H15_")))
cors <- sc |> group_by(panel) |> count_text()
message("vs anc: "); print(cors)
p_sc <- scatter(sc, "PASN response in anc (log2)", "PASN response in mutant (log2)") +
  facet_wrap(~panel)
p_sc$layers[[6]]$data <- cors          # per-panel r / hit counts
for (ext in "png") {
  ggsave(file.path(out, paste0("fig_response_vs_anc_poster_mval.", ext)), p_sc,
         width = 9.5, height = 5.9, dpi = 300, device = if (ext == "pdf") cairo_pdf else "png")
}

# stp against stk (direct comparison)
sd <- resp |> filter(int_direct_tested) |>
  transmute(locus_tag, x = resp_stk, y = resp_stp, class = hit_class(int_direct_tested, int_direct_FDR)) |>
  left_join(meta, by = "locus_tag") |>
  mutate(label = coalesce(gene, str_remove(locus_tag, "B7H15_")))
message("stp vs stk: "); print(count_text(sd))
p_d <- scatter(sd, "PASN response in stk (log2)", "PASN response in stp (log2)") +
  facet_wrap(~"stp vs stk")
for (ext in "png") {
  ggsave(file.path(out, paste0("fig_response_stk_vs_stp_poster_mval.", ext)), p_d,
         width = 6.3, height = 6.9, dpi = 300, device = if (ext == "pdf") cairo_pdf else "png")
}
message("Done")
# Square panels: same span on x and y within each plot, and aspect ratio 1, so the
# diagonal stays at 45 degrees and all three panels get the same size.
sq <- function(d) range(c(d$x, d$y), na.rm = TRUE)
PANEL_IN <- 3.6                                   # every panel exactly 3.6 x 3.6 inches
SPACING_PT <- 5.5                                 # gap between facet panels (theme default)
p_sc_sq <- p_sc + coord_cartesian(xlim = sq(sc), ylim = sq(sc), clip = "off") +
  theme(panel.spacing = unit(SPACING_PT, "pt"), plot.margin = margin(5.5, SPACING_PT / 2, 5.5, 5.5))
# third panel: y axis (ticks + title) on the right, so the three panels sit side by side
p_d_sq <- p_d + coord_cartesian(xlim = sq(sd), ylim = sq(sd), clip = "off") +
  scale_y_continuous(breaks = seq(-8, 8, 4), position = "right") +
  theme(plot.margin = margin(5.5, 5.5, 5.5, SPACING_PT / 2))
# fixed panel-area sizes: left plot = 2 panels + facet gap, right plot = 1 panel; same height
# the gap between plot 2 and 3 equals the facet gap, so all three boxes are evenly spaced
fig <- (p_sc_sq + p_d_sq) +
  plot_layout(widths = unit(c(2 * PANEL_IN + SPACING_PT / 72.27, PANEL_IN), "in"),
              heights = unit(PANEL_IN, "in"), guides = "collect") &
  theme(legend.position = "bottom")
ggsave(file.path(out, "fig_responses_combined_poster_mval.png"), fig,
       width = 13, height = 5.6, dpi = 300)
message("Done combined")
