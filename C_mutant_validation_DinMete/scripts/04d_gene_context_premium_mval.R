# Step 4d of our workflow (C_mutant_validation_DinMete/project_outlines/PROJECT_WORKFLOW_mval.md):
# a cleaner, presentation-quality version of the Step 4 gene map (fig_gene_context_mval.png).
#
# Run from the repository root:   Rscript C_mutant_validation_DinMete/scripts/04d_gene_context_premium_mval.R   (~5 s)
# Needs the Step 4 output gene_context_mval.csv (run C_mutant_validation_DinMete/scripts/04_polar_effect_mval.R first).
# Does NOT read anything from data/, so it also runs on a fresh clone.
#
# ---- Plan (coding rules, section 4) -----------------------------------------
# Aim:      Same content as fig_gene_context_mval.png, drawn so the polar-effect argument is
#           readable at a glance: (A) the genomic window, (B) a zoom into stp1-stk1 with both
#           insertions, their codon, and the 4 bp overlap.
# Input:    C_mutant_validation_DinMete/results/step04_polar_effect/gene_context_mval.csv (coordinates, strand,
#           labels, products, in_aureowiki_unit, n_samples_detected; written by Step 4)
#           TRANSPOSON_INSERTIONS, TRANSPOSON_LENGTH_BP, GENOTYPE_COLOURS (R/constants.R)
# Analysis: none (drawing only). Codons are recomputed from the coordinates, not typed in.
# Output:   C_mutant_validation_DinMete/results/step04_polar_effect/fig_gene_context_premium_mval.{png,pdf}
# Check:    codons 189 (stp1) and 46 (stk1); protein lengths 247 / 664 aa; 4 bp overlap;
#           12 genes in the window (the same as Step 4).
#
# Colours: the insertions use the genotype colours of every other figure (stk orange,
# stp green), so "the green triangle" always means the stp strain.
# The bent arrow at the start of the unit marks where the unit begins; its exact promoter
# position is schematic (AureoWiki / Mäder et al. 2016 give the unit, not the base).

suppressMessages(source("R/setup.R"))

out_dir <- "C_mutant_validation_DinMete/results/step04_polar_effect"
ctx_file <- file.path(out_dir, "gene_context_mval.csv")
if (!file.exists(ctx_file)) stop("Run C_mutant_validation_DinMete/scripts/04_polar_effect_mval.R first.")
check <- function(label, ok) message(sprintf("[%s] %s", if (isTRUE(ok)) "PASS" else "FLAG", label))

STP1 <- "B7H15_06350"; STK1 <- "B7H15_06355"
FONT <- "Arial"
INK <- "#1F2933"; MUTED <- "#6B7785"
FILLS <- c("stp1" = unname(GENOTYPE_COLOURS["stp"]),
           "stk1" = unname(GENOTYPE_COLOURS["stk"]),
           "same transcription unit" = "#9AAABB",
           "flanking gene" = "#D5DBE1",
           "protein not detected in our data" = "white")

ctx <- read_csv(ctx_file, show_col_types = FALSE) |>
  mutate(category = case_when(locus_tag == STP1 ~ "stp1",
                              locus_tag == STK1 ~ "stk1",
                              n_samples_detected == 0 ~ "protein not detected in our data",
                              in_aureowiki_unit ~ "same transcription unit",
                              TRUE ~ "flanking gene"),
         category = factor(category, levels = names(FILLS)),
         # gene symbols in italics, bare locus numbers upright
         face = if_else(str_detect(label, "^[0-9]"), "plain", "italic"))
check("12 genes in the window (as Step 4)", nrow(ctx) == 12)
check("stp1 and stk1 overlap by 4 bp", ctx$gap_to_next_bp[ctx$locus_tag == STP1] == -4)

ins <- TRANSPOSON_INSERTIONS |>
  mutate(codon = (insertion_pos - cds_start) %/% 3 + 1,
         n_aa  = (cds_end - cds_start + 1) / 3 - 1,             # minus the stop codon
         # how far into the gene, by nucleotide (CLAUDE.md: stp1 76 %)
         pct_in = round(100 * (insertion_pos - cds_start + 1) / (cds_end - cds_start + 1)),
         strain = if_else(short_name == "Stp1", "stp", "stk"))
check("insertion codons 189 (stp1) and 46 (stk1)", identical(as.numeric(ins$codon), c(189, 46)))
check("protein lengths 247 / 664 aa", identical(as.numeric(ins$n_aa), c(247, 664)))

# ---- Geometry -------------------------------------------------------------------------
# One ggplot with two tracks sharing one x range. Track A (y = 3) is in genome bp. Track B
# (y = 0) is a linear zoom of stp1-stk1 stretched over the full width; zx() maps bp into it.
XL   <- c(1226700, 1238150)
ZOOM <- c(min(ins$cds_start) - 120, max(ins$cds_end) + 120)
zx   <- function(x) XL[1] + (x - ZOOM[1]) * diff(XL) / diff(ZOOM)
YA <- 3; YB <- 0

# gene arrow as a polygon: body half-height h, head half-height H, head length `head` (x units)
gene_poly <- function(start, end, strand, y, h, H, head) {
  head <- min(head, 0.5 * (end - start))
  if (strand == "+") {
    x <- c(start, end - head, end - head, end, end - head, end - head, start)
  } else {
    x <- c(end, start + head, start + head, start, start + head, start + head, end)
  }
  tibble(x = x, y = y + c(-h, -h, -H, 0, H, h, h))
}

polys_a <- ctx |>
  mutate(poly = purrr::pmap(list(start, end, strand), \(s, e, st) gene_poly(s, e, st, YA, 0.17, 0.30, 230))) |>
  select(locus_tag, category, poly) |> unnest(poly)

zoom_genes <- ctx |> filter(locus_tag %in% c(STP1, STK1)) |>
  left_join(select(ins, locus_tag, n_aa), by = "locus_tag") |>
  mutate(zs = zx(start), ze = zx(end))
polys_b <- zoom_genes |>
  mutate(poly = purrr::pmap(list(zs, ze, strand), \(s, e, st) gene_poly(s, e, st, YB, 0.30, 0.50, 520))) |>
  select(locus_tag, category, poly) |> unnest(poly)

unit <- ctx |> filter(in_aureowiki_unit) |> summarise(s = min(start), e = max(end))
pair <- ctx |> filter(locus_tag %in% c(STP1, STK1)) |> summarise(s = min(start), e = max(end))
overlap_x <- zx(mean(c(ctx$start[ctx$locus_tag == STK1], ctx$end[ctx$locus_tag == STP1])))
ins <- ins |> mutate(zi = zx(insertion_pos))

# light "lens" connecting the stp1-stk1 stretch in track A to its zoom in track B
# (starts below the gene labels so it does not run through them)
lens <- tibble(x = c(pair$s, pair$e, XL[2], XL[1]), y = c(YA - 0.72, YA - 0.72, YB + 1.75, YB + 1.75))
lens_stem <- tibble(xmin = pair$s, xmax = pair$e, ymin = YA - 0.72, ymax = YA - 0.36)

ticks <- seq(1228000, 1238000, 2000)
txt <- function(pt) pt / .pt

p <- ggplot() +
  # ---- track A: genomic window -------------------------------------------------------
  annotate("segment", x = XL[1], xend = XL[2], y = YA + 1.15, yend = YA + 1.15, colour = "grey80", linewidth = 0.3) +
  annotate("segment", x = ticks, xend = ticks, y = YA + 1.15, yend = YA + 1.07, colour = "grey60", linewidth = 0.3) +
  annotate("text", x = ticks, y = YA + 1.24, label = paste0(scales::comma(ticks / 1000), " kb"),
           family = FONT, size = txt(7.5), colour = MUTED, vjust = 0) +
  annotate("text", x = XL[1], y = YA + 1.24, label = "JE2 (CP020619.1)", hjust = 0, vjust = 0,
           family = FONT, size = txt(7.5), colour = MUTED, fontface = "bold") +
  # transcription unit: bent arrow at its start, line over the unit, arrow head at its end
  annotate("path", x = c(unit$s, unit$s, unit$e), y = c(YA + 0.34, YA + 0.62, YA + 0.62),
           colour = INK, linewidth = 0.55,
           arrow = arrow(length = unit(1.8, "mm"), type = "closed")) +
  annotate("text", x = unit$s + 120, y = YA + 0.72, hjust = 0, vjust = 0, family = FONT, size = txt(8),
           colour = INK, label = "one transcription unit (annotated in NCTC 8325: AureoWiki; Mäder et al. 2016)") +
  geom_polygon(data = lens, aes(x, y), fill = "#EEF1F4") +
  geom_rect(data = lens_stem, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "#EEF1F4") +
  geom_polygon(data = polys_a, aes(x, y, group = locus_tag, fill = category),
               colour = "white", linewidth = 0.4) +
  geom_polygon(data = filter(polys_a, category == "protein not detected in our data"),
               aes(x, y, group = locus_tag), fill = NA, colour = "grey55", linewidth = 0.35, linetype = "22") +
  geom_text(data = ctx, aes((start + end) / 2, YA - 0.52, label = label, fontface = face),
            family = FONT, size = txt(8.5), colour = INK) +
  # ---- track B: zoom into stp1-stk1 -----------------------------------------------------
  geom_polygon(data = polys_b, aes(x, y, group = locus_tag, fill = category), colour = "white", linewidth = 0.6) +
  geom_text(data = zoom_genes, aes(zs + 180, YB, label = label), hjust = 0, family = FONT,
            fontface = "bold.italic", size = txt(13), colour = "white") +
  geom_text(data = zoom_genes, aes(zs + 820, YB,
                                   label = sprintf("%s  ·  %d aa", product, n_aa)),
            hjust = 0, family = FONT, size = txt(8.5), colour = "white") +
  # insertions: stem, triangle, two-line label (stp label to the left, stk label to the right)
  geom_segment(data = ins, aes(x = zi, xend = zi, y = YB + 0.52, yend = YB + 0.88), colour = INK, linewidth = 0.5) +
  geom_point(data = ins, aes(zi, YB + 0.95, fill = strain), shape = 25, size = 4.2, colour = INK, stroke = 0.6) +
  geom_text(data = ins, aes(zi + if_else(strain == "stp", -160, 160), YB + 1.25,
                            hjust = if_else(strain == "stp", 1, 0), colour = strain,
                            label = sprintf("%s mutant  (%s)", strain, ntml_id)),
            family = FONT, fontface = "bold", size = txt(9.5)) +
  geom_text(data = ins, aes(zi + if_else(strain == "stp", -160, 160), YB + 0.98,
                            hjust = if_else(strain == "stp", 1, 0),
                            label = sprintf("Tn in codon %d of %d (%d %% into the gene)", codon, n_aa, pct_in)),
            family = FONT, size = txt(8.5), colour = INK) +
  # the 4 bp overlap
  annotate("segment", x = overlap_x, xend = overlap_x, y = YB - 0.52, yend = YB - 0.86, colour = INK, linewidth = 0.4) +
  annotate("text", x = overlap_x, y = YB - 0.95, vjust = 1, family = FONT, size = txt(8.5), colour = INK,
           label = "stp1 and stk1 overlap by 4 bp: stk1 starts before stp1 ends") +
  annotate("text", x = XL[1], y = YB + 1.62, hjust = 0, family = FONT, size = txt(7.5), colour = MUTED,
           fontface = "bold", label = "zoom: stp1-stk1") +
  scale_fill_manual(values = c(FILLS, GENOTYPE_COLOURS[c("stk", "stp")]), breaks = names(FILLS), name = NULL,
                    labels = c("stp1 (Ser/Thr phosphatase)", "stk1 (Ser/Thr kinase)",
                               "same transcription unit", "flanking gene", "protein not detected in our data")) +
  scale_colour_manual(values = GENOTYPE_COLOURS, guide = "none") +
  guides(fill = guide_legend(nrow = 1, override.aes = list(colour = c(NA, NA, NA, NA, "grey55"),
                                                           linetype = c(1, 1, 1, 1, 2)))) +
  coord_cartesian(xlim = XL, ylim = c(-1.2, YA + 1.45), expand = FALSE, clip = "off") +
  labs(title = "stp1 and stk1 lie in one transcription unit",
       subtitle = paste0("The stp insertion sits upstream of the whole stk1 gene, so it can also reduce stk1 ",
                         "expression (polar effect, a hypothesis to test).\nBoth mutants carry the same ",
                         scales::comma(TRANSPOSON_LENGTH_BP), " bp bursa aurealis transposon (ermB, gfp)."),
       caption = paste0("Gene coordinates: JE2 genome CP020619.1 (course GFF). Insertions: NTML Tn insertion map. ",
                        "Transcription unit: AureoWiki SAOUHSC_01182-01187 (NCTC 8325; Mäder et al. 2016, PMID 27035918),\n",
                        "placed on JE2 by gene order and UniProt ortholog (an assumption for the first four genes). ",
                        "Promoter position schematic. 'Not detected' = no value in any of the 27 total-proteome samples.")) +
  theme_void(base_family = FONT) +
  theme(plot.title = element_text(face = "bold", size = 15, colour = INK, margin = margin(b = 4)),
        plot.subtitle = element_text(size = 10, colour = MUTED, lineheight = 1.15, margin = margin(b = 14)),
        plot.caption = element_text(size = 7, colour = MUTED, hjust = 0, lineheight = 1.15, margin = margin(t = 12)),
        plot.title.position = "plot", plot.caption.position = "plot",
        legend.position = "bottom", legend.text = element_text(size = 8.5, colour = INK),
        legend.key.size = unit(3.2, "mm"), legend.margin = margin(t = 10),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(14, 18, 10, 18))

png_file <- file.path(out_dir, "fig_gene_context_premium_mval.png")
ggsave(png_file, p, width = 11, height = 5.6, dpi = 300, device = ragg::agg_png, bg = "white")
ggsave(sub("png$", "pdf", png_file), p, width = 11, height = 5.6, device = cairo_pdf)
message("saved ", png_file, " (+ .pdf)")

# ---- Wide version: one row, true genome scale, no zoom track ----------------------------
# For a wide slot (slide or poster banner). Gene names above the arrows, stp1/stk1 named
# inside their arrows, insertions below the genes so they do not cross the unit arrow.
XW <- c(1226750, 1238100)
polys_w <- ctx |>
  mutate(poly = purrr::pmap(list(start, end, strand), \(s, e, st) gene_poly(s, e, st, 0, 0.24, 0.40, 210))) |>
  select(locus_tag, category, poly) |> unnest(poly)
inside <- ctx |> filter(locus_tag %in% c(STP1, STK1)) |>
  left_join(select(ins, locus_tag, n_aa), by = "locus_tag") |>
  # stp1 is too short for more than its name; its length is in the insertion label below
  mutate(text = if_else(locus_tag == STK1, sprintf("%s  ·  %d aa", label, n_aa), label))
ins_w <- ins |> mutate(side = if_else(strain == "stp", -1, 1), hj = if_else(strain == "stp", 1, 0))

p_wide <- ggplot() +
  # genome axis
  annotate("segment", x = XW[1], xend = XW[2], y = 1.45, yend = 1.45, colour = "grey80", linewidth = 0.3) +
  annotate("segment", x = ticks, xend = ticks, y = 1.45, yend = 1.38, colour = "grey60", linewidth = 0.3) +
  annotate("text", x = ticks, y = 1.52, label = paste0(scales::comma(ticks / 1000), " kb"),
           family = FONT, size = txt(7.5), colour = MUTED, vjust = 0) +
  annotate("text", x = XW[1], y = 1.52, label = "JE2 (CP020619.1)", hjust = 0, vjust = 0,
           family = FONT, size = txt(7.5), colour = MUTED, fontface = "bold") +
  # transcription unit
  annotate("path", x = c(unit$s, unit$s, unit$e), y = c(0.72, 0.92, 0.92), colour = INK, linewidth = 0.55,
           arrow = arrow(length = unit(1.8, "mm"), type = "closed")) +
  annotate("text", x = unit$s + 110, y = 1.0, hjust = 0, vjust = 0, family = FONT, size = txt(8), colour = INK,
           label = "one transcription unit (NCTC 8325: AureoWiki; Mäder et al. 2016)") +
  # genes
  geom_polygon(data = polys_w, aes(x, y, group = locus_tag, fill = category), colour = "white", linewidth = 0.5) +
  geom_polygon(data = filter(polys_w, category == "protein not detected in our data"),
               aes(x, y, group = locus_tag), fill = NA, colour = "grey55", linewidth = 0.35, linetype = "22") +
  geom_text(data = filter(ctx, !locus_tag %in% c(STP1, STK1)),
            aes((start + end) / 2, 0.58, label = label, fontface = face),
            family = FONT, size = txt(8.5), colour = INK) +
  geom_text(data = inside, aes(start + 90, 0, label = text), hjust = 0, family = FONT,
            fontface = "bold.italic", size = txt(8.5), colour = "white") +
  # insertions below the genes
  geom_segment(data = ins_w, aes(x = insertion_pos, xend = insertion_pos, y = -0.27, yend = -0.56),
               colour = INK, linewidth = 0.5) +
  geom_point(data = ins_w, aes(insertion_pos, -0.62, fill = strain), shape = 24, size = 3.8,
             colour = INK, stroke = 0.6) +
  geom_text(data = ins_w, aes(insertion_pos + side * 110, -0.62, hjust = hj, colour = strain,
                              label = sprintf("%s mutant (%s)", strain, ntml_id)),
            family = FONT, fontface = "bold", size = txt(9)) +
  geom_text(data = ins_w, aes(insertion_pos + side * 110, -0.82, hjust = hj,
                              label = sprintf("Tn in codon %d of %d (%d %% into the gene)", codon, n_aa, pct_in)),
            family = FONT, size = txt(8), colour = INK) +
  scale_fill_manual(values = c(FILLS, GENOTYPE_COLOURS[c("stk", "stp")]), breaks = names(FILLS), name = NULL,
                    labels = c("stp1 (Ser/Thr phosphatase)", "stk1 (Ser/Thr kinase)",
                               "same transcription unit", "flanking gene", "protein not detected\nin our data")) +
  scale_colour_manual(values = GENOTYPE_COLOURS, guide = "none") +
  guides(fill = guide_legend(ncol = 1, override.aes = list(colour = c(NA, NA, NA, NA, "grey55"),
                                                           linetype = c(1, 1, 1, 1, 2)))) +
  coord_cartesian(xlim = XW, ylim = c(-0.95, 1.75), expand = FALSE, clip = "off") +
  labs(title = "stp1 and stk1 lie in one transcription unit",
       subtitle = paste0("stk1 starts 4 bp before stp1 ends; the stp insertion sits upstream of the whole stk1 gene, ",
                         "so it can also reduce stk1 expression (polar effect, a hypothesis to test). Both mutants carry the ",
                         scales::comma(TRANSPOSON_LENGTH_BP), " bp bursa aurealis transposon (ermB, gfp)."),
       caption = paste0("Coordinates: JE2 genome CP020619.1 (course GFF). Insertions: NTML Tn insertion map. ",
                        "Transcription unit: AureoWiki SAOUHSC_01182-01187 (NCTC 8325; Mäder et al. 2016, PMID 27035918), ",
                        "placed on JE2 by gene order and UniProt ortholog (an assumption for the first four genes).\n",
                        "Promoter position schematic. 'Not detected' = no value in any of the 27 total-proteome samples.")) +
  theme_void(base_family = FONT) +
  theme(plot.title = element_text(face = "bold", size = 14, colour = INK, margin = margin(b = 3)),
        plot.subtitle = element_text(size = 9.5, colour = MUTED, margin = margin(b = 10)),
        plot.caption = element_text(size = 7, colour = MUTED, hjust = 0, lineheight = 1.15, margin = margin(t = 8)),
        plot.title.position = "plot", plot.caption.position = "plot",
        legend.position = "right", legend.text = element_text(size = 8.5, colour = INK),
        legend.key.size = unit(3.2, "mm"), legend.key.spacing.y = unit(1.5, "mm"),
        legend.margin = margin(l = 14),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(12, 16, 8, 16))

wide_file <- file.path(out_dir, "fig_gene_context_premium_wide_mval.png")
ggsave(wide_file, p_wide, width = 14, height = 3.4, dpi = 300, device = ragg::agg_png, bg = "white")
ggsave(sub("png$", "pdf", wide_file), p_wide, width = 14, height = 3.4, device = cairo_pdf)
message("saved ", wide_file, " (+ .pdf)")
