# Step 12c (poster): large-text, transparent version of Step 3 Fig 1b (fig_enzyme_peptides_mval.png).
# Run from the repository root:  Rscript C_mutant_validation_DinMete/12c_poster_enzyme_peptides_mval.R
# Aim:      the Stk1/Stp1 peptide-evidence figure at poster size (strain on x, medium as symbol).
# Input:    C_mutant_validation_DinMete/results/step03_validation/enzyme_per_sample_mval.csv (peptides per sample; NA = not detected)
# Analysis: none (counts only; NA drawn as an open symbol at 0, never used in a log2 calculation).
# Output:   C_mutant_validation_DinMete/results/step03_validation/fig_enzyme_peptides_poster_mval.png (transparent)
# Check:    54 rows (2 proteins x 27 samples); same counts as fig_enzyme_peptides_mval.png.
suppressMessages(source("R/setup.R")); source("C_mutant_validation_DinMete/12c_poster_transparent_helper_mval.R")

out_dir <- "C_mutant_validation_DinMete/results/step03_validation"
per_sample <- read_csv("C_mutant_validation_DinMete/results/step03_validation/enzyme_per_sample_mval.csv",
                       show_col_types = FALSE)
stopifnot(nrow(per_sample) == 54, n_distinct(per_sample$sample) == 27)

d <- per_sample |>
  mutate(detected = !is.na(n_peptides),
         n_plot = replace_na(n_peptides, 0),          # count scale only, never logged
         medium = factor(medium, levels = c("TSB", "SNpulse", "SN")),
         genotype = factor(genotype, levels = c("anc", "stk", "stp")),
         short_name = factor(short_name, levels = c("Stk1", "Stp1")),
         # filled shape when detected, open shape at 0 when not detected
         shp = paste(medium, detected),
         # fixed, evenly spaced x positions: strain -> medium block -> replicate
         x = as.integer(genotype) + (as.integer(medium) - 2) * 0.29 + (replicate - 2) * 0.085)

shapes <- c("TSB TRUE" = 16, "SNpulse TRUE" = 17, "SN TRUE" = 15,
            "TSB FALSE" = 1, "SNpulse FALSE" = 2, "SN FALSE" = 0)

p <- ggplot(d, aes(x, n_plot, colour = genotype, shape = shp)) +
  geom_point(size = 2.8, stroke = 1.2) +
  scale_x_continuous(breaks = 1:3, labels = c("anc", "stk", "stp"), limits = c(0.5, 3.5),
                     expand = c(0, 0)) +
  facet_wrap(~ short_name, scales = "free_y") +
  scale_colour_manual(values = GENOTYPE_COLOURS, guide = "none") +
  scale_shape_manual(values = shapes, name = NULL,
                     breaks = c("TSB TRUE", "SNpulse TRUE", "SN TRUE"),
                     labels = c("TSB", "PASN pulse", "PASN")) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.06, 0.08))) +
  labs(x = NULL, y = "peptides per sample",
       caption = "open symbol at 0 = protein not detected in that sample") +
  theme_bw(base_size = 17) +
  theme(legend.position = "top",
        legend.text = element_text(size = 17),
        legend.margin = margin(0, 0, -6, 0),
        strip.text = element_text(face = "bold", size = 17),
        axis.text.x = element_text(size = 17, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 15, colour = "black"),
        plot.caption = element_text(size = 15, colour = "grey25", hjust = 0),
        plot.caption.position = "plot",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

for (ext in "png") {
  ggsave(file.path(out_dir, paste0("fig_enzyme_peptides_poster_mval.", ext)), p,
         width = 8.5, height = 5, dpi = 300)
}
message("Done")
