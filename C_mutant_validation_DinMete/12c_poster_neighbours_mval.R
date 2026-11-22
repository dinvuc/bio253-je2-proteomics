# Step 12c (poster): large-text, transparent version of Step 4b fig_neighbours_by_genotype_mval.png.
# Run from the repository root:  Rscript C_mutant_validation_DinMete/12c_poster_neighbours_mval.R   (needs Step 4 output)
# Aim:      the stp1/stk1 neighbourhood, one panel per strain, at poster size.
# Input:    step04_polar_effect/gene_context_mval.csv; read_abundances("total") (log2, NA = not detected)
# Analysis: none (raw per-sample values; red digit = samples detected of 3, shown when < 3).
# Output:   C_mutant_validation_DinMete/results/step04_polar_effect/fig_neighbours_by_genotype_poster_mval.png (transparent)
# Check:    27 samples, 11 proteins; same values as fig_neighbours_by_genotype_mval.png.
suppressMessages(source("R/setup.R")); source("C_mutant_validation_DinMete/12c_poster_transparent_helper_mval.R")

out_dir <- "C_mutant_validation_DinMete/results/step04_polar_effect"
context_file <- "C_mutant_validation_DinMete/results/step04_polar_effect/gene_context_mval.csv"
STP1 <- "B7H15_06350"; STK1 <- "B7H15_06355"

context <- read_csv(context_file, show_col_types = FALSE) |>
  filter(n_samples_detected > 0) |>
  arrange(start) |>
  mutate(tag = str_remove(locus_tag, "B7H15_"),
         gene = label,   # gene name, or the locus number when the gene has no name
         gene = factor(gene, levels = gene))

abund <- read_abundances("total") |> keep_je2() |>
  filter(locus_tag %in% context$locus_tag) |>
  select(locus_tag, sample, genotype, medium, normalized_abundance) |>
  left_join(context |> select(locus_tag, gene), by = "locus_tag")

stopifnot(!anyDuplicated(abund[c("locus_tag", "sample")]),
          n_distinct(abund$sample) == 27, n_distinct(abund$locus_tag) == 11)

detect <- abund |>
  group_by(gene, locus_tag, genotype, medium) |>
  summarise(n_det = sum(!is.na(normalized_abundance)), .groups = "drop")

enzyme_pos <- which(context$locus_tag %in% c(STP1, STK1))
band <- tibble(xmin = min(enzyme_pos) - 0.5, xmax = max(enzyme_pos) + 0.5)
dodge <- position_dodge(width = 0.85)
y_low <- floor(min(abund$normalized_abundance, na.rm = TRUE)) - 1.5

p <- ggplot(abund |> filter(!is.na(normalized_abundance)),
            aes(gene, normalized_abundance, colour = genotype, shape = medium, group = medium)) +
  geom_rect(data = band, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey88") +
  geom_point(size = 2.6, position = dodge) +
  geom_text(data = detect, aes(gene, y_low + 0.6, group = medium,
                               label = if_else(n_det < 3, as.character(n_det), "")),
            position = dodge, size = 5.5, colour = "firebrick", fontface = "bold",
            inherit.aes = FALSE) +
  facet_wrap(~ genotype, ncol = 1,
             labeller = as_labeller(c(anc = "anc (wild type)",
                                      stk = "stk (stk1::Tn)",
                                      stp = "stp (stp1::Tn)"))) +
  scale_colour_manual(values = GENOTYPE_COLOURS, guide = "none") +
  scale_shape_manual(values = c(TSB = 16, SNpulse = 17, SN = 15), name = NULL,
                     labels = c(TSB = "TSB", SNpulse = "PASN pulse", SN = "PASN")) +
  scale_y_continuous(limits = c(y_low, NA), breaks = c(10, 14, 18)) +
  labs(x = "genes in genome order (grey = stp1, stk1)",
       y = "log2 abundance",
       caption = "red number = samples detected (of 3), shown when < 3") +
  theme_bw(base_size = 17) +
  theme(legend.position = "top",
        legend.text = element_text(size = 17),
        legend.margin = margin(0, 0, -6, 0),
        strip.text = element_text(face = "bold", size = 17),
        axis.text.x = element_text(size = 15, colour = "black"),
        plot.caption = element_text(size = 15, colour = "firebrick", hjust = 0),
        plot.caption.position = "plot",
        axis.text.y = element_text(size = 15, colour = "black"),
        panel.grid.minor = element_blank())

for (ext in "png") {
  ggsave(file.path(out_dir, paste0("fig_neighbours_by_genotype_poster_mval.", ext)), p,
         width = 8.5, height = 7, dpi = 300)
}
message("Done")
