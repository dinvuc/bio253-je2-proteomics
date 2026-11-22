# Helper for the C_mutant_validation_DinMete/scripts/12c_poster_*_mval.R scripts: poster figures with transparent backgrounds.
# Plot, panel and legend backgrounds are transparent (the strip bars and grey bands stay),
# so a poster background can show through. Sourced after R/setup.R; overrides ggsave() locally.
library(patchwork)
TRANSP <- theme(plot.background = element_rect(fill = "transparent", colour = NA),
                panel.background = element_rect(fill = "transparent", colour = NA),
                legend.background = element_rect(fill = "transparent", colour = NA),
                legend.box.background = element_rect(fill = "transparent", colour = NA),
                legend.key = element_rect(fill = "transparent", colour = NA))
ggsave <- function(filename, plot, ...) {
  ggplot2::ggsave(filename, plot & TRANSP, ..., bg = "transparent")
}
