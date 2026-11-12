# Step 4b (add-on to scripts/04_polar_effect_mval.R): the neighbourhood of
# stp1/stk1 seen one strain at a time.
#
# Run from the repository root:   Rscript scripts/04b_neighbours_by_genotype_mval.R   (~10 s)
# Needs the Step 4 output gene_context_mval.csv: run scripts/04_polar_effect_mval.R first.
#
# ---- Plan (coding rules, section 4) -----------------------------------------
# Aim:      fig_neighbours_mval.png shows one panel per protein. Here we show one
#           panel per STRAIN (anc / stk / stp), with the proteins along the x-axis in
#           genome order and stp1 + stk1 included where they sit in the genome.
#           That way you can read along the chromosome and see where abundance drops.
# Input:    gene_context_mval.csv (Step 4: locus tags, start, label, n_samples_detected)
#           read_abundances("total")  normalised log2 abundance, NA = not detected
# Analysis: none: raw per-sample values only (CLAUDE.md rule 4). NA stays NA, never 0
#           (rule 6); instead we print "n/3 detected" where a group has < 3 values.
# Output:   mval_medi/results/step04_polar_effect/fig_neighbours_by_genotype_mval.png
# Check:    27 samples, 3 per group; 11 proteins (9 neighbours + stp1 + stk1);
#           detection counts for stp1/stk1 match the Step 3/4 tables.

suppressMessages(source("R/setup.R"))

out_dir <- "mval_medi/results/step04_polar_effect"
context_file <- file.path(out_dir, "gene_context_mval.csv")
if (!file.exists(context_file)) stop("Run scripts/04_polar_effect_mval.R first.")

STP1 <- "B7H15_06350"; STK1 <- "B7H15_06355"

# Genes in the window that were measured at all (unmeasured ones have no dots to show),
# in genome order. stp1/stk1 are kept even though they are missing in some samples.
context <- read_csv(context_file, show_col_types = FALSE) |>
  filter(n_samples_detected > 0) |>
  arrange(start) |>
  mutate(tag = str_remove(locus_tag, "B7H15_"),
         gene = if_else(label == tag, tag, paste0(label, "\n(", tag, ")")),
         gene = factor(gene, levels = gene))

abund <- read_abundances("total") |> keep_je2() |>
  filter(locus_tag %in% context$locus_tag) |>
  select(locus_tag, sample, genotype, medium, normalized_abundance) |>
  left_join(context |> select(locus_tag, gene), by = "locus_tag")   # join by locus tag, not position

# Checks (coding rules section 6): every sample once per protein, 3 samples per group.
stopifnot(!anyDuplicated(abund[c("locus_tag", "sample")]),
          n_distinct(abund$sample) == 27,
          all(count(distinct(abund, sample, genotype, medium), genotype, medium)$n == 3))
message("CHECK proteins plotted: ", n_distinct(abund$locus_tag),
        " (", paste(levels(context$gene) |> str_replace("\n", " "), collapse = ", "), ")")

# How many of the 3 replicates were detected per protein x strain x medium.
detect <- abund |>
  group_by(gene, locus_tag, genotype, medium) |>
  summarise(n_det = sum(!is.na(normalized_abundance)), .groups = "drop")
gaps <- detect |> filter(n_det < 3)
message("CHECK groups with < 3/3 detected:")
print(gaps |> select(locus_tag, genotype, medium, n_det))

# Grey band behind stp1 and stk1 so the two enzymes stand out in every panel.
enzyme_pos <- which(context$locus_tag %in% c(STP1, STK1))
band <- tibble(xmin = min(enzyme_pos) - 0.5, xmax = max(enzyme_pos) + 0.5)

dodge <- position_dodge(width = 0.6)
y_low <- floor(min(abund$normalized_abundance, na.rm = TRUE)) - 1   # room for the n/3 labels

p <- ggplot(abund |> filter(!is.na(normalized_abundance)),
            aes(gene, normalized_abundance, colour = genotype, shape = medium, group = medium)) +
  geom_rect(data = band, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey90") +
  geom_point(size = 1.9, position = dodge) +
  # Use ALL groups (blank label at 3/3): dodging needs every medium present at each
  # gene, otherwise a lone label is centred and lands under the wrong medium.
  geom_text(data = detect, aes(gene, y_low + 0.4, group = medium,
                               label = if_else(n_det < 3, paste0(n_det, "/3"), "")),
            position = dodge, size = 2.4, colour = "firebrick", inherit.aes = FALSE) +
  facet_wrap(~ genotype, ncol = 1,
             labeller = as_labeller(c(anc = "anc (wild type)",
                                      stk = "stk (stk1::Tn, NE217)",
                                      stp = "stp (stp1::Tn, NE1919; also low Stk1 = polar effect)"))) +
  scale_colour_manual(values = GENOTYPE_COLOURS, guide = "none") +
  scale_shape_manual(values = c(TSB = 16, SNpulse = 17, SN = 15), name = "medium") +
  scale_y_continuous(limits = c(y_low, NA)) +
  labs(x = "measured genes in genome order, left to right (JE2, CP020619.1)",
       y = "normalised log2 abundance",
       title = "Proteins around stp1/stk1, one panel per strain (each dot = one sample)",
       subtitle = paste("Grey band = stp1 and stk1. Red n/3 = replicates detected when fewer than 3;",
                        "missing samples are not plotted (NA, not 0).")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"))

ggsave(file.path(out_dir, "fig_neighbours_by_genotype_mval.png"), p,
       width = 10, height = 9, dpi = 300)
message("Done: ", file.path(out_dir, "fig_neighbours_by_genotype_mval.png"))

# ---- Figure 2: rows = medium (TSB, SN), columns = strain -> 6 panels ----------
# Same data, but every panel holds ONE strain in ONE medium, so the medium effect
# (read down a column) and the strain effect (read along a row) are separated.
# SNpulse is left out on purpose (group request: TSB vs SN only).
two_media <- c("TSB", "SN")
abund2 <- abund |> filter(medium %in% two_media) |> mutate(medium = factor(medium, levels = two_media))
detect2 <- detect |> filter(medium %in% two_media, n_det < 3) |>
  mutate(medium = factor(medium, levels = two_media))
stopifnot(n_distinct(abund2$sample) == 18)          # 3 strains x 2 media x 3 replicates
message("CHECK figure 2: ", n_distinct(abund2$sample), " samples, ",
        nrow(detect2), " panel x protein cells with < 3/3 detected")

p2 <- ggplot(abund2 |> filter(!is.na(normalized_abundance)),
             aes(gene, normalized_abundance, colour = genotype)) +
  geom_rect(data = band, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey90") +
  geom_point(size = 1.8, position = position_jitter(width = 0.12, height = 0, seed = 1)) +
  geom_text(data = detect2, aes(gene, y_low + 0.4, label = paste0(n_det, "/3")),
            size = 2.6, colour = "firebrick", inherit.aes = FALSE) +
  facet_grid(medium ~ genotype,
             labeller = labeller(genotype = c(anc = "anc (wild type)", stk = "stk (stk1::Tn)",
                                               stp = "stp (stp1::Tn; also low Stk1)"))) +
  scale_colour_manual(values = GENOTYPE_COLOURS, guide = "none") +
  scale_x_discrete(labels = \(x) str_replace(x, "\n", " ")) +
  scale_y_continuous(limits = c(y_low, NA)) +
  labs(x = "measured genes in genome order, left to right (JE2, CP020619.1)",
       y = "normalised log2 abundance",
       title = "Proteins around stp1/stk1 per strain (columns) and medium (rows); each dot = one sample",
       subtitle = paste("Grey band = stp1 and stk1. Red n/3 = replicates detected when fewer than 3;",
                        "missing samples are not plotted (NA, not 0). SNpulse not shown.")) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1), strip.text = element_text(face = "bold"))

ggsave(file.path(out_dir, "fig_neighbours_by_genotype_medium_mval.png"), p2,
       width = 14, height = 7, dpi = 300)
message("Done: ", file.path(out_dir, "fig_neighbours_by_genotype_medium_mval.png"))
