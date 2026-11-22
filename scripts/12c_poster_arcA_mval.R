# Step 12c (poster): large-text, transparent version of Step 6 fig_neighbourhood_arcA_14675_mval.png.
# Run from the repository root:  Rscript scripts/12c_poster_arcA_mval.R
# Needs data/derived/limma_treat_total_mval.tsv.gz (git-ignored; written by scripts/05_formal_test_mval.R).
# Aim:      arcA and its 10 protein-coding neighbours: (a) gene map, (b) every sample, (c) Step 5 interactions.
# Input:    JE2_genome.gff, read_protein_metadata(), read_abundances("total"), read_samples(), Step 5 limma table
# Analysis: none new (Step 5 estimates +/- 2 SE; filled = formal hit, treat lfc = 1, BH-FDR < 0.05).
# Output:   mval_medi/results/step06_responses/hit_neighbourhoods/fig_neighbourhood_arcA_poster_mval.png (transparent)
# Check:    18 samples, 11 genes in the window; same estimates as fig_neighbourhood_arcA_14675_mval.png.
suppressMessages(source("R/setup.R")); source("scripts/12c_poster_transparent_helper_mval.R")
library(patchwork)
out <- "mval_medi/results/step06_responses/hit_neighbourhoods"
h <- "B7H15_14675"            # arcA
N_SIDE <- 5; RUN_GAP <- 150
INTER <- c("stk vs anc" = "interaction_stk_vs_anc_at_TSB_vs_SN",
           "stp vs anc" = "interaction_stp_vs_anc_at_TSB_vs_SN",
           "stp vs stk" = "interaction_stp_vs_stk_at_TSB_vs_SN")
COMP_COLOURS <- c("stk vs anc" = unname(GENOTYPE_COLOURS["stk"]),
                  "stp vs anc" = unname(GENOTYPE_COLOURS["stp"]),
                  "stp vs stk" = "#CC79A7")
BS <- 17                       # base font size, as in the other poster figures
TXT <- 5                       # geom_text size (mm) ~ 14 pt

gff <- read_tsv("data/raw/resources/JE2_genome.gff", comment = "#", col_names = FALSE,
                col_types = "cccddcccc", progress = FALSE) |>
  filter(X3 == "gene") |>
  transmute(locus_tag = str_match(X9, "locus_tag=([^;]+)")[, 2],
            biotype = str_match(X9, "gene_biotype=([^;]+)")[, 2],
            start = X4, end = X5, strand = X7) |>
  arrange(start)
coding <- gff |> filter(biotype == "protein_coding") |> mutate(idx = row_number())
meta <- read_protein_metadata() |>
  transmute(locus_tag, product,
            gene = if_else(ortholog_quality == "good", gene_symbol_best, NA_character_))
coding <- coding |> left_join(distinct(meta, locus_tag, .keep_all = TRUE), by = "locus_tag") |>
  mutate(tag = str_remove(locus_tag, "B7H15_"), short = coalesce(gene, tag))

samples <- read_samples()
abund <- read_abundances("total") |> keep_je2() |>
  select(locus_tag, sample, normalized_abundance) |>
  inner_join(select(samples, sample, genotype, medium), by = "sample") |>
  filter(medium %in% c("TSB", "SN")) |>
  mutate(medium = factor(medium, c("TSB", "SN")))
stopifnot(n_distinct(abund$sample) == 18)

stats <- read_tsv("data/derived/limma_treat_total_mval.tsv.gz", show_col_types = FALSE) |>
  filter(contrast %in% INTER) |>
  transmute(locus_tag, comparison = names(INTER)[match(contrast, INTER)],
            diff, std.error, FDR, tested, hit = tested & FDR < MVAL_FDR)
all_hits <- unique(stats$locus_tag[stats$hit])

arrow_poly <- function(start, end, strand, y = 0, h = 0.22, head_bp) {
  hd <- min(head_bp, 0.45 * (end - start))
  x <- if (strand == "+") c(start, end - hd, end - hd, end, end - hd, end - hd, start)
       else c(end, start + hd, start + hd, start, start + hd, start + hd, end)
  tibble(x = x, y = y + c(-h, -h, -1.7 * h, 0, 1.7 * h, h, h))
}
same_strand_run <- function(i) {
  gi <- which(gff$locus_tag == coding$locus_tag[coding$idx == i]); run <- gi
  for (dir in c(-1, 1)) {
    j <- gi
    repeat {
      k <- j + dir
      if (k < 1 || k > nrow(gff) || gff$strand[k] != gff$strand[gi]) break
      gap <- if (dir > 0) gff$start[k] - gff$end[j] - 1 else gff$start[j] - gff$end[k] - 1
      if (gap >= RUN_GAP) break
      run <- c(run, k); j <- k
    }
  }
  gff$locus_tag[sort(run)]
}

i0 <- coding$idx[coding$locus_tag == h]
win <- coding |> filter(idx %in% (i0 - N_SIDE):(i0 + N_SIDE)) |> mutate(offset = idx - i0)
stopifnot(nrow(win) == 11)
run <- same_strand_run(i0)
det <- abund |> filter(locus_tag %in% win$locus_tag) |>
  group_by(locus_tag) |> summarise(n_detected = sum(!is.na(normalized_abundance)), .groups = "drop")
win <- win |> left_join(det, by = "locus_tag") |>
  mutate(n_detected = coalesce(n_detected, 0L),
         role = case_when(offset == 0 ~ "the hit (arcA)",
                          locus_tag %in% all_hits ~ "also a hit",
                          n_detected == 0 ~ "not detected",
                          TRUE ~ "neighbour"),
         x_lab = factor(short, levels = short))

# (a) gene map
span <- c(min(win$start), max(win$end)); pad <- 0.02 * diff(span); head_bp <- 0.014 * diff(span)
map_genes <- gff |> filter(start >= span[1], end <= span[2]) |>
  left_join(select(win, locus_tag, role, short, offset), by = "locus_tag") |>
  mutate(role = coalesce(role, "other gene"))
polys <- map_genes |> mutate(id = row_number()) |> group_by(id, locus_tag, role) |>
  reframe(arrow_poly(start, end, strand, head_bp = head_bp))
labs_map <- map_genes |> filter(!is.na(short)) |> arrange(start) |>
  mutate(mid = (start + end) / 2, ylab = if_else(row_number() %% 2 == 1, -0.62, -1.02))
run_box <- gff |> filter(locus_tag %in% run)
run_x <- c(max(min(run_box$start), span[1] - pad), min(max(run_box$end), span[2] + pad))
fill_vals <- c("the hit (arcA)" = "#D55E00", "also a hit" = "#E69F00",
               "neighbour" = "grey55", "not detected" = "grey88", "other gene" = "white")
fill_vals <- fill_vals[names(fill_vals) %in% polys$role]
p_map <- ggplot() +
  annotate("rect", xmin = run_x[1], xmax = run_x[2], ymin = -0.45, ymax = 0.7, fill = "#DCEBF7") +
  annotate("text", x = run_x[1] + pad / 3, y = 0.68, hjust = 0, vjust = 1, size = TXT * 0.85,
           colour = "#2B6A99", label = "same-strand run: co-transcription possible") +
  geom_polygon(data = polys, aes(x, y, group = id, fill = role), colour = "grey25", linewidth = 0.3) +
  geom_segment(data = labs_map, aes(x = mid, xend = mid, y = -0.42, yend = ylab + 0.13),
               colour = "grey70", linewidth = 0.3) +
  geom_text(data = labs_map, aes(mid, ylab, label = short,
                                 fontface = if_else(offset == 0, "bold", "plain")), size = TXT) +
  scale_fill_manual(values = fill_vals, breaks = names(fill_vals), name = NULL) +
  scale_x_continuous(limits = span + c(-pad, pad), breaks = NULL, name = NULL) +
  scale_y_continuous(limits = c(-1.15, 0.72), breaks = NULL, name = NULL) +
  theme_bw(base_size = BS) +
  theme(legend.position = "top", legend.text = element_text(size = BS - 2),
        legend.key.size = unit(0.5, "cm"), legend.margin = margin(0, 0, -4, 0),
        panel.grid = element_blank())

# (b) every sample, TSB and PASN
ab <- abund |> filter(locus_tag %in% win$locus_tag) |>
  left_join(select(win, locus_tag, x_lab), by = "locus_tag")
ab_det <- ab |> group_by(x_lab, medium, genotype) |>
  summarise(n_det = sum(!is.na(normalized_abundance)), .groups = "drop")
y_low <- floor(min(ab$normalized_abundance, na.rm = TRUE)) - 1.5
nd <- win |> filter(n_detected == 0) |> select(x_lab) |>
  tidyr::crossing(medium = factor(c("TSB", "SN"), c("TSB", "SN")))
dodge <- position_dodge(width = 0.75)
p_ab <- ggplot(ab |> filter(!is.na(normalized_abundance)),
               aes(x_lab, normalized_abundance, colour = genotype, group = genotype)) +
  annotate("rect", xmin = N_SIDE + 0.5, xmax = N_SIDE + 1.5, ymin = -Inf, ymax = Inf, fill = "grey90") +
  geom_point(size = 2.4, position = dodge) +
  geom_text(data = ab_det, aes(x_lab, y_low + 0.6, group = genotype,
                               label = if_else(n_det < 3 & !x_lab %in% nd$x_lab, as.character(n_det), "")),
            position = dodge, size = TXT, fontface = "bold", colour = "firebrick", inherit.aes = FALSE) +
  geom_text(data = nd, aes(x_lab, y_low + 3, label = "not\ndetected"), size = TXT * 0.85,
            colour = "grey45", lineheight = 0.85, inherit.aes = FALSE) +
  facet_grid(medium ~ ., labeller = labeller(medium = c(TSB = "TSB", SN = "PASN"))) +
  scale_x_discrete(drop = FALSE) +
  scale_colour_manual(values = GENOTYPE_COLOURS, name = NULL) +
  scale_y_continuous(limits = c(y_low, NA), breaks = c(9, 12, 15, 18)) +
  labs(x = NULL, y = "log2 abundance") +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  theme_bw(base_size = BS) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = BS - 2, colour = "black"),
        strip.text = element_text(face = "bold"),
        legend.position = "top", legend.text = element_text(size = BS),
        legend.margin = margin(0, 0, -4, 0), panel.grid.minor = element_blank())

# (c) interactions, estimate +/- 2 SE
st <- stats |> filter(locus_tag %in% win$locus_tag) |>
  left_join(select(win, locus_tag, x_lab), by = "locus_tag") |>
  mutate(comparison = factor(comparison, names(INTER)), shown = tested & !is.na(diff))
nt_gene <- st |> filter(!shown) |> count(x_lab) |> filter(n == length(INTER))
pd <- position_dodge(width = 0.7)
p_int <- ggplot(st |> filter(shown), aes(x_lab, diff, colour = comparison, group = comparison)) +
  annotate("rect", xmin = N_SIDE + 0.5, xmax = N_SIDE + 1.5, ymin = -Inf, ymax = Inf, fill = "grey90") +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_hline(yintercept = c(-MVAL_LFC, MVAL_LFC), linetype = 2, colour = "grey60") +
  geom_errorbar(aes(ymin = diff - 2 * std.error, ymax = diff + 2 * std.error),
                width = 0, position = pd, linewidth = 0.7) +
  geom_point(aes(shape = hit), size = 3, stroke = 1.1, position = pd, fill = "white") +
  geom_text(data = nt_gene, aes(x_lab, 0.45, label = "n.t."), size = TXT * 0.85, colour = "grey45",
            inherit.aes = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_colour_manual(values = COMP_COLOURS, name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                     labels = c(`TRUE` = "hit", `FALSE` = "not a hit"), name = NULL) +
  labs(x = NULL, y = "interaction (log2)") +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 3.5)),
         shape = guide_legend(order = 2, override.aes = list(size = 3.5))) +
  theme_bw(base_size = BS) +
  theme(axis.text.x = element_text(size = BS - 2.5, colour = "black"),
        axis.text.y = element_text(size = BS - 2, colour = "black"),
        legend.position = "bottom", legend.text = element_text(size = BS - 1),
        legend.margin = margin(-6, 0, 0, 0), panel.grid.minor = element_blank())

fig <- (p_map / p_ab / p_int) +
  plot_layout(heights = c(0.9, 2.2, 1.3)) +
  plot_annotation(tag_levels = "a",
                  caption = paste0("red number = samples detected (of 3), shown when < 3;  ",
                                   "n.t. = not testable\n",
                                   "interaction > 0 = PASN response lower in the first strain;  error bars = \u00b1 2 SE"),
                  theme = theme(plot.caption = element_text(size = BS - 3, colour = "grey25", hjust = 0))) &
  theme(plot.tag = element_text(size = BS + 2, face = "bold"))

for (ext in "png") {
  ggsave(file.path(out, paste0("fig_neighbourhood_arcA_poster_mval.", ext)), fig,
         width = 10, height = 11, dpi = 300, device = if (ext == "pdf") cairo_pdf else "png")
}
message("Done")
