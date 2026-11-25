# Step 12 (PROJECT_WORKFLOW_mval.md): the four poster figures, 300 dpi, from saved results.
#
# Run from the repository root:   Rscript scripts/12_poster_figures_mval.R   (~20 s)
# Needs Steps 3, 4b, 5, 6, 7, 10b and 11 first (it reads their outputs; nothing is refitted).
#
# ---- Plan (coding rules, section 4) -----------------------------------------
# Aim:      One figure per poster question, readable on an A0 poster, every number
#           traceable to a step script.
#             Fig 1 (H1)  Stp1/Stk1 per sample, peptides per sample, ErmB control
#             Fig 2 (H2)  mutant PASN response vs the anc response; interaction hits
#                         marked with their Step 11 stability
#             Fig 3 (H3)  stp effects vs stk effects (baseline, response), noise-corrected r
#             Fig 4       pathway level (Step 10) + genome-position scan (Step 11)
# Input:    step03 enzyme_per_sample, step04b ermb_per_sample, data/derived/
#           limma_treat_total_mval.tsv.gz, data/derived/effect_pairs_step07_mval.tsv.gz,
#           step07 similarity_stats, step10 summary/robust_terms, step11
#           loso_hit_stability + genome_top_windows (thresholds) + genome_windows (blocks)
# Analysis: no new statistics. The genome-scan curve (mean of 10 neighbouring proteins)
#           is recomputed exactly as in Step 11; its thresholds are read from Step 11.
#           18 primary samples (TSB, SN) in every per-sample panel.
# Output:   mval_medi/results/step12_poster/fig1_h1_enzymes_mval.png ... fig4_pathways_genome_mval.png
#           + fig3_h3_shared_distinct_compact_mval.png (Fig 3, narrow)
#           + fig3_h3_shared_distinct_compact_transparent_mval.png (same, no background)
# Check:    18 samples, 3 per group in the per-sample panels; scan curve reproduces the
#           Step 11 top window of stk_vs_anc_at_TSB (+1.510).

suppressMessages(source("R/setup.R"))
library(patchwork)
out <- "mval_medi/results/step12_poster"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 13))
check <- function(label, ok) { message(sprintf("[%s] %s", if (isTRUE(ok)) "PASS" else "FAIL", label)); invisible(ok) }
prim <- \(x) filter(x, medium %in% c("TSB", "SN")) |> mutate(medium = factor(medium, levels = c("TSB", "SN")))
GCOL <- GENOTYPE_COLOURS

# ---- Fig 1: H1 ----------------------------------------------------------------------
enz <- read_csv("mval_medi/results/step03_validation/enzyme_per_sample_mval.csv", show_col_types = FALSE) |>
  prim() |> mutate(genotype = factor(genotype, levels = GENOTYPES))
check("Fig 1: 18 samples x 2 enzymes, 3 per group",
      nrow(enz) == 36 && all(table(enz$short_name, enz$group)[, unique(enz$group)] == 3))
floor_y <- min(enz$normalized_abundance, na.rm = TRUE) - 1.5
p1a <- ggplot(enz |> add_x_primary(), aes(x_pos, normalized_abundance, colour = genotype)) +
  geom_point(size = 2.6, position = position_jitter(width = 0.05, height = 0, seed = 1), na.rm = TRUE) +
  geom_point(data = \(d) filter(d, is.na(normalized_abundance)), aes(y = floor_y), shape = 4, size = 2.6,
             position = position_jitter(width = 0.05, height = 0, seed = 1), show.legend = FALSE) +
  facet_wrap(~short_name) + primary_x_axis() + scale_colour_manual(values = GCOL) +
  labs(x = NULL, y = "protein abundance (log2)", title = "a  Kinase and phosphatase per sample",
       subtitle = "x at the bottom = not detected (no value; not a zero)")
p1b <- ggplot(enz |> add_x_primary(), aes(x_pos, coalesce(n_peptides, 0), colour = genotype)) +
  geom_point(size = 2.6, position = position_jitter(width = 0.05, height = 0.15, seed = 1)) +
  facet_wrap(~short_name) + primary_x_axis() + scale_colour_manual(values = GCOL) +
  labs(x = NULL, y = "peptides quantified", title = "b  Peptides behind each value",
       subtitle = "0 = protein not quantified in that sample")
erm <- read_csv("mval_medi/results/step04b_transposon/ermb_per_sample_mval.csv", show_col_types = FALSE) |>
  prim() |> mutate(genotype = factor(genotype, levels = GENOTYPES))
p1c <- ggplot(erm |> add_x_primary(), aes(x_pos, spectral_count, colour = genotype)) +
  geom_point(size = 2.6, position = position_jitter(width = 0.05, height = 0, seed = 1)) +
  primary_x_axis() + scale_colour_manual(values = GCOL) + expand_limits(y = 0) +
  labs(x = NULL, y = "ErmB spectra", title = "c  Transposon marker ErmB",
       subtitle = "positive control for the insertion")
p1 <- (p1a / (p1b | p1c + theme(legend.position = "none")) + plot_layout(guides = "collect", heights = c(1, 1))) +
  plot_annotation(title = "H1: what do the mutants actually lack?",
                  caption = "18 samples (TSB, SN), n = 3 per group. Steps 3 and 4b.")
ggsave(file.path(out, "fig1_h1_enzymes_mval.png"), p1, width = 11, height = 9, dpi = 300)

# ---- Fig 2: H2 ----------------------------------------------------------------------
res  <- read_tsv("data/derived/limma_treat_total_mval.tsv.gz", show_col_types = FALSE) |> keep_je2()
meta <- read_protein_metadata() |> filter(ortholog_quality == "good") |> distinct(locus_tag, gene = gene_symbol_best)
stab <- read_csv("mval_medi/results/step11_robustness/loso_hit_stability_mval.csv", show_col_types = FALSE)
resp <- function(mut) {
  ct_int <- paste0("interaction_", mut, "_vs_anc_at_TSB_vs_SN")
  w <- res |> filter(contrast %in% c("TSB_vs_SN_at_anc", paste0("TSB_vs_SN_at_", mut), ct_int)) |>
    select(protein_key, locus_tag, contrast, diff, tested, FDR) |>
    mutate(contrast = case_when(contrast == "TSB_vs_SN_at_anc" ~ "anc", str_starts(contrast, "interaction") ~ "int",
                                TRUE ~ "mut")) |>
    tidyr::pivot_wider(names_from = contrast, values_from = c(diff, tested, FDR))
  w |> filter(tested_int) |> mutate(mutant = mut, hit = FDR_int < MVAL_FDR) |>
    left_join(meta, by = "locus_tag") |>
    left_join(stab |> filter(contrast == ct_int) |> select(locus_tag, status, in_own_block), by = "locus_tag") |>
    mutate(hit_class = case_when(!hit ~ "no altered response",
                                 !is.na(in_own_block) ~ "altered; lies in a shifted genomic block",
                                 status == "stable (18/18)" ~ "altered; stable (18/18 leave-one-out)",
                                 TRUE ~ "altered; fragile (lost in some leave-one-out refits)"))
}
r2 <- bind_rows(resp("stk"), resp("stp")) |>
  mutate(panel = paste0(mutant, " (", sum(hit), " altered of ", n(), ")"), .by = mutant)
rr <- r2 |> summarise(r = cor(diff_anc, diff_mut), .by = mutant)
cls <- c("no altered response" = "grey75", "altered; stable (18/18 leave-one-out)" = "#009E73",
         "altered; fragile (lost in some leave-one-out refits)" = "#E69F00",
         "altered; lies in a shifted genomic block" = "#CC79A7")
p2 <- ggplot(r2, aes(diff_anc, diff_mut)) +
  geom_abline(colour = "grey40") + geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(data = \(d) filter(d, !hit), colour = "grey75", size = 1, alpha = 0.6) +
  geom_point(data = \(d) filter(d, hit), aes(fill = hit_class), shape = 21, size = 3.4) +
  ggrepel::geom_text_repel(data = \(d) filter(d, hit), aes(label = coalesce(gene, locus_tag)), size = 4, seed = 1,
                           min.segment.length = 0, box.padding = 0.5) +
  geom_text(data = rr |> left_join(distinct(r2, mutant, panel), by = "mutant"),
            aes(x = -Inf, y = Inf, label = sprintf("r = %.2f", r)), hjust = -0.2, vjust = 1.5, size = 4.5) +
  facet_wrap(~panel) + coord_equal() +
  scale_fill_manual(values = cls, drop = TRUE, name = NULL) +
  labs(x = "PASN response of anc: log2 (TSB - SN)\n> 0 = lower in SN",
       y = "PASN response of the mutant: log2 (TSB - SN)",
       title = "H2: do the mutants respond differently to P. aeruginosa supernatant?",
       subtitle = "Each point = one protein measured in both media in both strains. Diagonal = same response as anc.\nAltered = interaction (mutant vs anc) with limma treat |log2| > 1, BH-FDR < 0.05.",
       caption = "Steps 5, 6, 11. n = 3 per group.") +
  theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 2))
ggsave(file.path(out, "fig2_h2_responses_mval.png"), p2, width = 12, height = 7.5, dpi = 300)

# ---- Fig 3: H3 ----------------------------------------------------------------------
pairs <- read_tsv("data/derived/effect_pairs_step07_mval.tsv.gz", show_col_types = FALSE) |> filter(eligible_both)
sim <- read_csv("mval_medi/results/step07_shared_distinct/similarity_stats_mval.csv", show_col_types = FALSE)
lab3 <- sim |> transmute(layer, txt = sprintf("n = %d proteins\nnoise-corrected r = %.2f (95%% CI %.2f-%.2f)\n(shared anc noise alone gives r = %.2f)",
                                             n, r_true_corrected, r_true_low, r_true_high, r_noise_only))
pairs <- pairs |> mutate(pt = case_when(direct_hit ~ "stp and stk differ directly",
                                        class == "hit in both, same direction" ~ "altered in both, same direction",
                                        TRUE ~ "other"),
                         layer = factor(layer, levels = sim$layer))
lay_lab <- c(setNames(c("a  Baseline in TSB: mutant - anc", "b  PASN response: interaction mutant vs anc"), sim$layer))
p3 <- ggplot(pairs, aes(stk_diff, stp_diff)) +
  geom_abline(colour = "grey40") + geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(data = \(d) filter(d, pt == "other"), colour = "grey70", size = 1, alpha = 0.6) +
  geom_point(data = \(d) filter(d, pt != "other"), aes(fill = pt), shape = 21, size = 3.2) +
  ggrepel::geom_text_repel(data = \(d) filter(d, pt != "other"), aes(label = coalesce(gene, locus_tag)), size = 3.4,
                           seed = 1, max.overlaps = 15, min.segment.length = 0) +
  geom_text(data = lab3, aes(x = -Inf, y = Inf, label = txt), hjust = -0.05, vjust = 1.2, size = 3.8, lineheight = 0.9) +
  facet_wrap(~layer, labeller = labeller(layer = lay_lab)) + coord_equal() +
  scale_fill_manual(values = c("stp and stk differ directly" = "#D55E00", "altered in both, same direction" = "#0072B2"), name = NULL) +
  labs(x = "stk vs anc (log2)", y = "stp vs anc (log2)",
       title = "H3: does stp (which also lacks most Stk1) look like stk?",
       subtitle = "Diagonal = identical effect in both mutants. Both axes subtract the same anc values, so noise alone correlates them;\nthe noise-corrected r removes that part.",
       caption = "Steps 5, 7. Direct stp-vs-stk differences: limma treat |log2| > 1, BH-FDR < 0.05.") +
  theme(legend.position = "bottom")
ggsave(file.path(out, "fig3_h3_shared_distinct_mval.png"), p3, width = 12, height = 7.5, dpi = 300)

# ---- Fig 3, compact ---------------------------------------------------------------------
# stk effects spread wider than stp effects, so stk goes on the y axis: the panels become narrow and
# the scales stay equal (diagonal at 45 degrees). Named: proteins with a gene symbol, all 'both'
# proteins and all hits of the response; the other locus tags made the centre unreadable.
named <- pairs |> filter(pt != "other") |>
  left_join(read_protein_metadata() |> distinct(locus_tag, annot_product = product), by = "locus_tag") |>
  mutate(label = case_when(!is.na(gene) ~ gene,
                           pt == "stp and stk differ directly" & layer == sim$layer[1] ~ "",
                           str_detect(annot_product, "hypothetical") ~ locus_tag,
                           TRUE ~ str_replace(str_remove(annot_product, " protein$"), "-binding", "-\nbinding")))
XL <- range(pairs$stp_diff) + c(-0.4, 0.4)
YL <- range(pairs$stk_diff) + c(-0.4, 0.4)
compact3 <- function(i, title, meaning) {
  s <- sim[i, ]
  ggplot(filter(pairs, layer == s$layer), aes(stp_diff, stk_diff)) +
    geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
    geom_abline(colour = "grey40", linetype = 2, linewidth = 0.4) +
    annotate("text", x = 2.75, y = 3.3, label = "same in both", angle = 45, size = 2.6, colour = "grey35") +
    geom_point(data = \(d) filter(d, pt == "other"), colour = "grey70", size = 0.6, alpha = 0.5) +
    geom_point(data = \(d) filter(d, pt != "other"), aes(fill = pt), shape = 21, size = 2.2) +
    ggrepel::geom_text_repel(data = filter(named, layer == s$layer), aes(label = label), size = 2.9, seed = 1,
                             lineheight = 0.85, box.padding = 0.5, point.padding = 0.15, force = 4, max.overlaps = Inf,
                             min.segment.length = 0, segment.colour = "grey40", segment.size = 0.25,
                             bg.color = "white", bg.r = 0.15) +
    coord_equal(xlim = XL, ylim = YL, expand = FALSE) +
    scale_x_continuous(breaks = c(-2, 0, 2, 4)) + scale_y_continuous(breaks = c(-4, -2, 0, 2, 4)) +
    scale_fill_manual(values = c("stp and stk differ directly" = "#D55E00", "altered in both, same direction" = "#0072B2"),
                      limits = c("altered in both, same direction", "stp and stk differ directly"), name = NULL) +
    labs(title = title,
         subtitle = paste(sprintf("r = %.2f (95%% CI %.2f-%.2f), n = %d", s$r_true_corrected, s$r_true_low, s$r_true_high, s$n),
                          meaning, sep = "\n"),
         x = "stp vs anc (log2)", y = if (i == 1) "stk vs anc (log2)") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 8.2),
          panel.grid = element_blank())
}
p3c <- compact3(1, "a  Baseline in TSB", "> 0 = higher in the mutant") +
  (compact3(2, "b  PASN response", "> 0 = response lower in the mutant") + guides(fill = "none")) +
  plot_layout(guides = "collect") +
  plot_annotation(caption = paste("r is corrected for the shared anc noise (0 = unrelated, 1 = identical). n = 3 per group.",
                                  "Coloured = limma treat |log2| > 1, BH-FDR < 0.05. Only some coloured proteins are named.",
                                  sep = "\n"),
                  theme = theme(plot.caption = element_text(hjust = 0, size = 8))) &
  theme(legend.position = "bottom")
check("Fig 3 compact: same 20 + 4 direct and 3 shared proteins",
      identical(as.integer(table(named$layer[named$direct_hit])), c(20L, 4L)) &&
        sum(named$pt == "altered in both, same direction") == 3)
ggsave(file.path(out, "fig3_h3_shared_distinct_compact_mval.png"), p3c, width = 7, height = 6.15, dpi = 400)
save_transparent(p3c, file.path(out, "fig3_h3_shared_distinct_compact_transparent_mval.png"), width = 7, height = 6.15)

# ---- Fig 4: pathway level + genome scan ------------------------------------------------
rob <- read_csv("mval_medi/results/step10_enrichment/summary/robust_terms_mval.csv", show_col_types = FALSE)
cnt <- rob |> mutate(side = if_else(auc > 0.5, "up", "down"),
                     layer = factor(layer, levels = c("baseline TSB", "baseline SN", "SN response (interaction)"))) |>
  count(mutant, layer, side) |> mutate(n_signed = if_else(side == "up", n, -n))
p4a <- ggplot(cnt, aes(layer, n_signed, fill = mutant)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) + geom_hline(yintercept = 0) +
  scale_fill_manual(values = GCOL[c("stk", "stp")]) +
  labs(x = NULL, y = "shifted functional terms\n(up: higher in mutant / response lower in mutant)",
       title = "a  How many functional terms shift?",
       subtitle = "Robust = rank test and cameraPR both FDR < 0.05. Terms overlap.") +
  theme(axis.text.x = element_text(angle = 12, hjust = 1))
top <- read_csv("mval_medi/results/step11_robustness/genome_top_windows_mval.csv", show_col_types = FALSE)
pos <- read_protein_metadata(cds_only = TRUE) |> filter(!is.na(start)) |> distinct(locus_tag, start)
WIN <- 10
win_mean <- function(v) as.numeric(stats::filter(c(v, v[seq_len(WIN - 1)]), rep(1 / WIN, WIN), sides = 1))[WIN:(length(v) + WIN - 1)]
SC <- c("stk_vs_anc_at_TSB", "stp_vs_anc_at_TSB")
scan <- purrr::map(SC, \(ct) {
  d <- res |> filter(contrast == ct, tested) |> inner_join(pos, by = "locus_tag") |> arrange(start)
  tibble(contrast = ct, start = d$start, w = win_mean(d$diff))
}) |> bind_rows() |> left_join(select(top, contrast, threshold), by = "contrast")
check("scan curve reproduces the Step 11 top window of stk_vs_anc_at_TSB (+1.510)",
      abs(max(scan$w[scan$contrast == "stk_vs_anc_at_TSB"]) - top$window_mean[top$contrast == "stk_vs_anc_at_TSB"]) < 0.001)
# labels sit LEFT of the peak/dip, inside the panel (they were clipped at the edges)
ann <- tibble(contrast = SC, x = c(2.14, 2.06), y = c(1.35, -1.75),
              txt = c("SigB / mazEF / kdp block\n(B7H15_11470-11590)", "prophage block\n(B7H15_10910-11150)"))
stk1_pos <- pos$start[pos$locus_tag == "B7H15_06355"] / 1e6
p4b <- ggplot(scan, aes(start / 1e6, w)) +
  geom_hline(aes(yintercept = threshold), linetype = 2, colour = "grey50") +
  geom_hline(aes(yintercept = -threshold), linetype = 2, colour = "grey50") +
  geom_vline(xintercept = stk1_pos, colour = "#CC79A7") +
  geom_line(aes(colour = contrast), linewidth = 0.5, show.legend = FALSE) +
  expand_limits(y = c(-2.5, 2)) +
  geom_text(data = ann, aes(x, y, label = txt), size = 3.6, hjust = 1.05, lineheight = 0.9) +
  facet_wrap(~contrast, ncol = 1, labeller = labeller(contrast = c(stk_vs_anc_at_TSB = "stk vs anc, TSB",
                                                                   stp_vs_anc_at_TSB = "stp vs anc, TSB"))) +
  scale_colour_manual(values = c(stk_vs_anc_at_TSB = unname(GCOL["stk"]), stp_vs_anc_at_TSB = unname(GCOL["stp"]))) +
  labs(x = "genome position (Mb)", y = "mean log2 effect,\n10 neighbouring proteins",
       title = "b  Are effects clustered along the genome?",
       subtitle = "Dashed = genome-wide 5% threshold (permutations). Pink line = stp1/stk1 locus.")
p4 <- (p4a | p4b) + plot_layout(widths = c(1, 1.3)) +
  plot_annotation(title = "Pathway level: broad shifts in both mutants, partly from two genomic regions",
                  caption = "Steps 10 and 11. A block of neighbouring proteins changing together can be co-regulation or a strain difference (e.g. copy number, prophage); DNA data are needed.")
ggsave(file.path(out, "fig4_pathways_genome_mval.png"), p4, width = 15, height = 7, dpi = 300)
message("done: ", out)
