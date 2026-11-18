# Step 12b (PROJECT_WORKFLOW_mval.md): poster versions of the Step 12 figures for the A0 poster.
#
# Run from the repository root:   Rscript scripts/12b_poster_a0_figures_mval.R   (~25 s)
# Then build the poster:          quarto typst compile mval_medi/results/step12_poster/poster_a0_mval.typ
# Needs the same inputs as Step 12 (Steps 3, 4b, 5, 6, 7, 10b, 11). Nothing is refitted.
#
# ---- Plan (coding rules, section 4) -----------------------------------------
# Aim:      The Step 12 figures are report figures (11-15 in wide, ~10 pt text, legends,
#           gridlines). On an A0 poster they would be shrunk or blown up and become unreadable.
#           Here every figure is drawn at its PRINTED size (one poster column = 243 mm), so the
#           sizes below are the real sizes on paper. Changes follow the lecture
#           "Academic posters" (BIO253, Nov 2026):
#             - text in figures >= 22-24 pt, Arial (same font as the poster text)
#             - explanations written on the figure (direct labels), no legends where possible
#             - no generic titles; the message sits in the poster heading above the figure
#             - no unnecessary gridlines; thick lines; filled shapes instead of thin lines
#               (genome scan drawn as a filled area)
#             - vector PDF for printing + a 300 dpi PNG preview
#             - fewer labels: only what supports the message (omit the rest)
# Input:    as Step 12, plus TRANSPOSON_INSERTIONS and KEY_PROTEINS (R/constants.R) for the operon
#           sketch.
# Analysis: no new statistics. Numbers on the figures (fold change, r, counts) are read from
#           the saved step results, never typed in.
# Output:   mval_medi/results/step12_poster/poster/p0_operon_mval.{pdf,png} ... p4_pathways_genome_mval.{pdf,png},
#           qr_repo_mval.svg
# Check:    18 samples / 3 per group in the per-sample panels; the Stk1 fold change printed on
#           Fig 1 equals the Step 5 estimate (log2 -5.42); same hit counts as Step 12 (5 / 1).

suppressMessages(source("R/setup.R"))
library(patchwork)
out <- "mval_medi/results/step12_poster/poster"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
check <- function(label, ok) { message(sprintf("[%s] %s", if (isTRUE(ok)) "PASS" else "FAIL", label)); invisible(ok) }

# ---- Poster geometry and style --------------------------------------------------------
# One poster column is 243 mm wide (A0 = 841 mm; margins 2 x 35 mm, gutters 2 x 20 mm;
# see poster_a0_mval.typ). Change both files together.
COL_W_MM <- 243
FONT     <- "Arial"
PT       <- 24                       # smallest text on the poster (lecture: >= 24 pt body text)
txt      <- function(pt = PT) pt / .pt   # geom_text() sizes are in mm; this converts from points
GCOL     <- GENOTYPE_COLOURS          # anc blue, stk orange, stp green: same colours in every figure
theme_set(
  theme_classic(base_size = 26, base_family = FONT) +
    theme(axis.text = element_text(size = PT, colour = "black"),
          axis.title = element_text(size = 26),
          axis.line = element_line(linewidth = 0.9), axis.ticks = element_line(linewidth = 0.9),
          axis.ticks.length = unit(2.5, "mm"),
          strip.background = element_blank(),
          strip.text = element_text(size = 26, face = "bold", hjust = 0),
          plot.title = element_text(size = 28, face = "bold", margin = margin(b = 3, unit = "mm")),
          plot.title.position = "plot",
          legend.position = "none",
          panel.spacing = unit(8, "mm"),
          plot.margin = margin(3, 6, 3, 3, "mm"))
)
save_fig <- function(p, name, h_mm) {
  w <- COL_W_MM / 25.4; h <- h_mm / 25.4
  ggsave(file.path(out, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(out, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  message("saved ", name, sprintf(" (%d x %d mm)", COL_W_MM, h_mm))
}
prim <- \(x) filter(x, medium %in% c("TSB", "SN")) |>
  mutate(medium = factor(medium, levels = c("TSB", "SN")), genotype = factor(genotype, levels = GENOTYPES))

res <- read_tsv("data/derived/limma_treat_total_mval.tsv.gz", show_col_types = FALSE) |> keep_je2()

# ---- p0: the stp1-stk1 operon and the two insertions ----------------------------------------
# Drawn from TRANSPOSON_INSERTIONS (NTML map), not from memory. Both genes are on the + strand.
ti <- TRANSPOSON_INSERTIONS |>
  left_join(select(KEY_PROTEINS, short_name, role), by = "short_name") |>
  mutate(off = min(cds_start),
         x0 = (cds_start - off) / 1000, x1 = (cds_end - off) / 1000, xi = (insertion_pos - off) / 1000,
         codon = (insertion_pos - cds_start) %/% 3 + 1,
         n_aa = (cds_end - cds_start + 1) / 3 - 1,             # minus the stop codon
         gene = tolower(short_name),
         strain = c(Stp1 = "stp", Stk1 = "stk")[short_name],
         enzyme = if_else(str_detect(role, "phosphatase"), "phosphatase", "kinase"))
check("p0: insertion codons are 189 (stp1) and 46 (stk1), as in CLAUDE.md",
      identical(as.numeric(ti$codon), c(189, 46)))
arrow_poly <- function(x0, x1, head = 0.13, h = 0.32, H = 0.5)
  tibble(x = c(x0, x1 - head, x1 - head, x1, x1 - head, x1 - head, x0), y = c(-h, -h, -H, 0, H, h, h))
genes <- ti |> select(gene, strain, x0, x1) |> purrr::pmap(\(gene, strain, x0, x1)
  arrow_poly(x0, x1) |> mutate(gene = gene, strain = strain)) |> bind_rows()
p0 <- ggplot() +
  geom_polygon(data = genes, aes(x, y, group = gene, fill = strain), colour = "white", linewidth = 0.6) +
  geom_text(data = ti, aes((x0 + x1) / 2 - 0.05, 0, label = gene), family = FONT, fontface = "bold.italic",
            colour = "white", size = txt(26)) +
  geom_text(data = ti, aes((x0 + x1) / 2, -0.85, label = sprintf("%s, %d aa", enzyme, n_aa)),
            family = FONT, size = txt(22)) +
  geom_point(data = ti, aes(xi, 0.72, fill = strain), shape = 25, size = 9, colour = "black", stroke = 0.8) +
  # two label rows so the two insertion labels (only 0.3 kb apart) do not overlap
  geom_segment(data = ti, aes(x = xi, xend = xi, y = 0.9, yend = if_else(strain == "stp", 1.42, 1.0), colour = strain),
               linewidth = 0.9) +
  geom_text(data = ti, aes(if_else(strain == "stp", 0, xi + 0.06), if_else(strain == "stp", 1.62, 1.05),
                           label = sprintf("%s mutant: Tn in codon %d", strain, codon), colour = strain),
            hjust = 0, family = FONT, fontface = "bold", size = txt(22)) +
  scale_fill_manual(values = GCOL) + scale_colour_manual(values = GCOL) +
  coord_cartesian(xlim = c(-0.05, max(ti$x1) + 0.05), ylim = c(-1.05, 1.8), clip = "off") +
  theme_void(base_family = FONT) + theme(legend.position = "none", plot.margin = margin(2, 4, 2, 4, "mm"))
save_fig(p0, "p0_operon_mval", 80)

# ---- p1: H1, what the mutants lack --------------------------------------------------------
enz <- read_csv("mval_medi/results/step03_validation/enzyme_per_sample_mval.csv", show_col_types = FALSE) |>
  prim() |> mutate(short_name = factor(short_name, levels = c("Stk1", "Stp1")))
check("Fig 1: 18 samples x 2 enzymes, 3 per group",
      nrow(enz) == 36 && all(count(enz, short_name, group)$n == 3))
stk1_fc <- res |> filter(contrast == "stp_vs_anc_at_TSB", locus_tag == "B7H15_06355") |> pull(diff)
check(sprintf("Fig 1: Stk1 stp_vs_anc_at_TSB = %.2f (Step 5: -5.42)", stk1_fc), abs(stk1_fc + 5.42) < 0.01)
ND <- 8   # y position of 'not detected' (below every measured value, ~10.4 is the lowest)
check("Fig 1: every measured value lies above the n.d. row", min(enz$normalized_abundance, na.rm = TRUE) > ND + 1.5)
jit <- position_jitter(width = 0.14, height = 0, seed = 1)
# protein in columns, medium in rows: 4 side-by-side facets were too narrow for 'anc stk stp' at 24 pt
fac <- facet_grid(medium ~ short_name)
p1a <- ggplot(enz, aes(genotype, coalesce(normalized_abundance, ND), colour = genotype)) +
  geom_point(aes(shape = is.na(normalized_abundance)), size = 5, stroke = 1.6, position = jit) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4)) +
  geom_text(data = tibble(short_name = factor("Stk1", levels(enz$short_name)), medium = factor("TSB", c("TSB", "SN")),
                          x = 2.25, y = 14.3, label = sprintf("%.0f× lower", 2^-stk1_fc)),
            aes(x = x, y = y, label = label), inherit.aes = FALSE, colour = GCOL[["stp"]], family = FONT, fontface = "bold",
            size = txt(24)) +
  scale_y_continuous(breaks = c(ND, seq(10, 18, 2)), labels = c("n.d.", seq(10, 18, 2)), limits = c(ND - 0.5, 18.8)) +
  fac + scale_colour_manual(values = GCOL) +
  labs(x = NULL, y = "abundance (log2)", title = "a  Protein per sample   (n.d. = not detected)")
p1b <- ggplot(enz, aes(genotype, coalesce(n_peptides, 0), colour = genotype)) +
  geom_point(size = 5, position = position_jitter(width = 0.14, height = 0.2, seed = 1)) +
  fac + scale_colour_manual(values = GCOL) + scale_y_continuous(breaks = c(0, 10, 20, 30)) +
  labs(x = NULL, y = "peptides", title = "b  Peptides behind each value")
erm <- read_csv("mval_medi/results/step04b_transposon/ermb_per_sample_mval.csv", show_col_types = FALSE) |> prim()
check("Fig 1c: ErmB in 18 samples, 0 spectra in every anc sample, > 0 in every mutant sample",
      nrow(erm) == 18 && all(erm$spectral_count[erm$genotype == "anc"] == 0) &&
        all(erm$spectral_count[erm$genotype != "anc"] > 0))
p1c <- ggplot(erm, aes(genotype, spectral_count, colour = genotype)) +
  geom_point(size = 5, position = jit) + facet_grid(~ medium) +
  scale_colour_manual(values = GCOL) + expand_limits(y = 0) +
  labs(x = NULL, y = "ErmB spectra", title = "c  Transposon marker ErmB: only in mutants")
p1 <- (p1a / p1b / p1c + plot_layout(heights = c(1.5, 1.1, 0.6))) &
  theme(strip.text.y = element_text(angle = 0, hjust = 0))
save_fig(p1, "p1_h1_enzymes_mval", 510)

# ---- p2: H2, PASN response of each mutant vs anc --------------------------------------------
meta <- read_protein_metadata() |> filter(ortholog_quality == "good") |> distinct(locus_tag, gene = gene_symbol_best)
stab <- read_csv("mval_medi/results/step11_robustness/loso_hit_stability_mval.csv", show_col_types = FALSE)
resp <- function(mut) {
  ct_int <- paste0("interaction_", mut, "_vs_anc_at_TSB_vs_SN")
  res |> filter(contrast %in% c("TSB_vs_SN_at_anc", paste0("TSB_vs_SN_at_", mut), ct_int)) |>
    select(protein_key, locus_tag, contrast, diff, tested, FDR) |>
    mutate(contrast = case_when(contrast == "TSB_vs_SN_at_anc" ~ "anc", str_starts(contrast, "interaction") ~ "int",
                                TRUE ~ "mut")) |>
    tidyr::pivot_wider(names_from = contrast, values_from = c(diff, tested, FDR)) |>
    filter(tested_int) |> mutate(mutant = mut, hit = FDR_int < MVAL_FDR) |>
    left_join(meta, by = "locus_tag") |>
    left_join(stab |> filter(contrast == ct_int) |> select(locus_tag, status, in_own_block), by = "locus_tag") |>
    # robust = kept in every leave-one-out refit AND not part of a block that shifts as a unit
    mutate(robust = hit & status %in% "stable (18/18)" & is.na(in_own_block))
}
r2 <- bind_rows(resp("stk"), resp("stp")) |>
  mutate(panel = sprintf("%s: %d of %s proteins altered", mutant, sum(hit), format(n(), big.mark = ",")),
         .by = mutant)
n_hit <- r2 |> summarise(n = sum(hit), .by = mutant)
check("Fig 2: 5 altered responses in stk, 1 in stp (as Step 12)", identical(n_hit$n, c(5L, 1L)))
rr <- r2 |> summarise(r = cor(diff_anc, diff_mut), .by = c(mutant, panel))
lim <- c(-6.5, 6)
key <- tibble(panel = unique(r2$panel)[2], x = lim[2], y = lim[1] + 0.3,
              label = "● robust (all leave-one-out refits)\n○ fragile or in a shifted block")
p2 <- ggplot(r2, aes(diff_anc, diff_mut)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.6) + geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.6) +
  geom_abline(colour = "grey30", linewidth = 0.9) +
  geom_point(data = \(d) filter(d, !hit), colour = "grey65", size = 1.6, alpha = 0.55) +
  geom_point(data = \(d) filter(d, hit), aes(colour = mutant, fill = if_else(robust, mutant, "open")),
             shape = 21, size = 7, stroke = 2) +
  ggrepel::geom_text_repel(data = \(d) filter(d, hit), aes(label = coalesce(gene, locus_tag)), family = FONT,
                           size = txt(24), seed = 3, min.segment.length = 0, box.padding = 1, point.padding = 0.6,
                           force = 4) +
  geom_text(data = rr, aes(x = lim[1], y = lim[2], label = sprintf("r = %.2f", r)), hjust = 0, vjust = 1,
            family = FONT, fontface = "bold", size = txt(28)) +
  annotate("text", x = 3.2, y = 5.6, label = "same response\nas anc", hjust = 1, vjust = 1, family = FONT,
           size = txt(22), colour = "grey30", lineheight = 0.85) +
  geom_text(data = key, aes(x, y, label = label), hjust = 1, vjust = 0, family = FONT, size = txt(22), lineheight = 0.9) +
  facet_wrap(~panel, ncol = 1) + coord_equal(xlim = lim, ylim = lim) +
  scale_colour_manual(values = GCOL) + scale_fill_manual(values = c(GCOL, open = "white")) +
  labs(x = "response of anc: log2 (TSB / PASN)\n> 0 = lower in PASN", y = "response of the mutant: log2 (TSB / PASN)")
save_fig(p2, "p2_h2_responses_mval", 445)

# ---- p3: H3, stp effects vs stk effects ------------------------------------------------------
pairs <- read_tsv("data/derived/effect_pairs_step07_mval.tsv.gz", show_col_types = FALSE) |> filter(eligible_both)
sim <- read_csv("mval_medi/results/step07_shared_distinct/similarity_stats_mval.csv", show_col_types = FALSE)
lay_lab <- setNames(c("Baseline in TSB", "Response to PASN"), sim$layer)
pairs <- pairs |>
  mutate(pt = case_when(direct_hit ~ "direct", class == "hit in both, same direction" ~ "both", TRUE ~ "other"),
         layer = factor(layer, levels = sim$layer))
n3 <- pairs |> summarise(n_direct = sum(pt == "direct"), n_both = sum(pt == "both"), .by = layer)
check("Fig 3: 20 + 4 direct stp-vs-stk hits (Step 7)", identical(n3$n_direct, c(20L, 4L)))
# only the noise-corrected r goes on the panel; 'noise alone gives 0.50' is said in the poster text
lab3 <- sim |> transmute(layer = factor(layer, levels = sim$layer), txt = sprintf("r = %.2f", r_true_corrected)) |>
  left_join(n3, by = "layer") |>
  mutate(key = sprintf("%d differ directly", n_direct))
# label the 'same direction in both' proteins and the few direct hits of the response layer; the 20 baseline
# direct hits stay unlabelled (count only), so the panel stays readable from 1.5 m
lab_pts <- pairs |> filter(pt == "both" | (pt == "direct" & layer == sim$layer[2]))
p3 <- ggplot(pairs, aes(stk_diff, stp_diff)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.6) + geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.6) +
  geom_abline(colour = "grey30", linewidth = 0.9) +
  geom_point(data = \(d) filter(d, pt == "other"), colour = "grey65", size = 1.4, alpha = 0.55) +
  geom_point(data = \(d) filter(d, pt != "other"), aes(fill = pt), shape = 21, size = 5, stroke = 0.8) +
  ggrepel::geom_text_repel(data = lab_pts, aes(label = coalesce(gene, locus_tag)), family = FONT, size = txt(22),
                           seed = 1, min.segment.length = 0, box.padding = 0.7, force = 3) +
  geom_text(data = lab3, aes(x = -5, y = 5, label = txt), hjust = 0, vjust = 1, family = FONT,
            size = txt(26), fontface = "bold") +
  geom_text(data = lab3, aes(x = 5, y = -5, label = key), hjust = 1, vjust = 0, family = FONT,
            size = txt(22), colour = "#D55E00", fontface = "bold") +
  facet_wrap(~layer, ncol = 1, labeller = labeller(layer = lay_lab)) + coord_equal(xlim = c(-5, 5), ylim = c(-5, 5)) +
  scale_x_continuous(breaks = c(-4, 0, 4)) + scale_y_continuous(breaks = c(-4, 0, 4)) +
  scale_fill_manual(values = c(direct = "#D55E00", both = "black")) +
  labs(x = "stk vs anc (log2)", y = "stp vs anc (log2)")
save_fig(p3, "p3_h3_shared_distinct_mval", 295)

# ---- p4: pathway counts + genome scan --------------------------------------------------------
rob <- read_csv("mval_medi/results/step10_enrichment/summary/robust_terms_mval.csv", show_col_types = FALSE)
cnt <- rob |> count(mutant, layer) |>
  mutate(layer = factor(layer, levels = c("baseline TSB", "baseline SN", "SN response (interaction)"),
                        labels = c("in TSB", "in PASN", "response")))
check("Fig 4a: robust terms stk 35/143/44, stp 192/168/8 (Step 10 guide)",
      identical(arrange(cnt, mutant, layer)$n, c(35L, 143L, 44L, 192L, 168L, 8L)))
# horizontal bars: the layer names are readable as y labels (as x labels they collided)
p4a <- ggplot(cnt, aes(n, forcats::fct_rev(layer), fill = mutant)) +
  geom_col(width = 0.75) + geom_text(aes(label = n), hjust = -0.2, family = FONT, size = txt(22)) +
  facet_wrap(~mutant, labeller = labeller(mutant = c(stk = "stk vs anc", stp = "stp vs anc"))) +
  scale_fill_manual(values = GCOL) + scale_x_continuous(expand = expansion(mult = c(0, 0.3)), breaks = c(0, 100, 200)) +
  labs(x = "shifted functional terms", y = NULL, title = "a  Pathway level")
top <- read_csv("mval_medi/results/step11_robustness/genome_top_windows_mval.csv", show_col_types = FALSE)
pos <- read_protein_metadata(cds_only = TRUE) |> filter(!is.na(start)) |> distinct(locus_tag, start)
WIN <- 10   # exactly as Step 11 and Step 12: mean of 10 neighbouring proteins, circular genome
win_mean <- function(v) as.numeric(stats::filter(c(v, v[seq_len(WIN - 1)]), rep(1 / WIN, WIN), sides = 1))[WIN:(length(v) + WIN - 1)]
SC <- c("stk_vs_anc_at_TSB", "stp_vs_anc_at_TSB")
scan <- purrr::map(SC, \(ct) {
  d <- res |> filter(contrast == ct, tested) |> inner_join(pos, by = "locus_tag") |> arrange(start)
  tibble(contrast = ct, start = d$start / 1e6, w = win_mean(d$diff))
}) |> bind_rows() |> left_join(select(top, contrast, threshold), by = "contrast")
check("Fig 4b: scan reproduces the Step 11 top window of stk_vs_anc_at_TSB (+1.510)",
      abs(max(scan$w[scan$contrast == SC[1]]) - top$window_mean[top$contrast == SC[1]]) < 0.001)
peaks <- bind_rows(scan |> filter(contrast == SC[1]) |> slice_max(w, n = 1),
                   scan |> filter(contrast == SC[2]) |> slice_min(w, n = 1)) |>
  # both block labels sit LEFT of their peak/dip (no room to the right); the stp1/stk1 label goes to
  # the empty bottom-left, the threshold label to the empty bottom-right of the stk track
  mutate(txt = c("SigB /\nmazEF /\nkdp block", "prophage\nblock"),
         x = start - 0.05, y = c(2.5, -1.75), vjust = c(1, 0.5))
stk1_pos <- pos$start[pos$locus_tag == "B7H15_06355"] / 1e6
thr_lab <- top |> filter(contrast == SC[1]) |> transmute(contrast, x = 2.95, y = -threshold - 0.45)
p4b <- ggplot(scan, aes(start, w)) +
  geom_ribbon(aes(ymin = 0, ymax = w, fill = contrast), outline.type = "full", colour = NA) +
  geom_hline(aes(yintercept = threshold), linetype = 2, linewidth = 0.8, colour = "grey35") +
  geom_hline(aes(yintercept = -threshold), linetype = 2, linewidth = 0.8, colour = "grey35") +
  geom_vline(xintercept = stk1_pos, linewidth = 0.8) +
  annotate("text", x = stk1_pos - 0.03, y = -2.35, label = "stp1/stk1", hjust = 1, family = FONT,
           fontface = "italic", size = txt(22)) +
  geom_text(data = thr_lab, aes(x, y, label = "5 % threshold"), hjust = 1, family = FONT, size = txt(22), colour = "grey35") +
  geom_text(data = peaks, aes(x, y, label = txt, vjust = vjust), hjust = 1, family = FONT, fontface = "bold",
            size = txt(22), lineheight = 0.85) +
  facet_wrap(~contrast, ncol = 1, labeller = labeller(contrast = c(stk_vs_anc_at_TSB = "stk vs anc (TSB)",
                                                                   stp_vs_anc_at_TSB = "stp vs anc (TSB)"))) +
  scale_fill_manual(values = c(stk_vs_anc_at_TSB = unname(GCOL["stk"]), stp_vs_anc_at_TSB = unname(GCOL["stp"]))) +
  scale_y_continuous(limits = c(-2.6, 2.6), breaks = c(-2, -1, 0, 1, 2)) +
  labs(x = "genome position (Mb)", y = "mean log2 effect\n(10 neighbours)",
       title = "b  Effects along the genome")
p4 <- p4a / p4b + plot_layout(heights = c(0.55, 1.45))
save_fig(p4, "p4_pathways_genome_mval", 470)

# ---- QR code to the repository (code and reports only, no raw data) ------------------------------
qrcode::generate_svg(qrcode::qr_code("https://github.com/dinvuc/bio253-je2-proteomics"),
                     filename = file.path(out, "qr_repo_mval.svg"), show = FALSE)
message("done: ", out)
