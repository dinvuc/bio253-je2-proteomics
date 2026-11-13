# Step 6 add-on: one figure per hit with its genomic neighbourhood.
# Do the genes around a hit change too, or only the hit itself?
#
# Run from the repository root:   Rscript scripts/06c_hit_neighbourhoods_mval.R   (~30 s)
# Needs scripts/05_formal_test_mval.R, 06_altered_responses_mval.R and
# 06b_hit_genome_context_mval.R first.
#
# ---- Plan (coding rules, section 4) -----------------------------------------
# Aim:      For each of the 8 Step 6 hits, show the hit in the middle with the 5
#           protein-coding genes before and after it (genome order), and ask
#           whether the neighbours show the same strain effect.
# Input:    data/raw/resources/JE2_genome.gff          genes, strand, biotype
#           read_protein_metadata()                      gene symbols (good orthologs), products
#           read_abundances("total"), read_samples()    per-sample log2 values, joined by name
#           data/derived/limma_treat_total_mval.tsv.gz  interaction estimates, SE, treat FDR
#           step06_responses/hits_genome_context_mval.csv  hit list, distances to the insertions
# Analysis: none new (descriptive). Per gene: raw values in TSB and SN, and the three
#           TSB-vs-SN interactions from Step 5 (estimate +/- 2 SE; filled = formal hit,
#           treat lfc = 1, BH-FDR < 0.05; n.t. = not testable, < 2/3 in a group).
#           Same-strand run = neighbours on the hit's strand with gaps < 150 bp:
#           a hint of co-transcription only, not a proven operon.
# Output:   mval_medi/results/step06_responses/hit_neighbourhoods/
#             fig_neighbourhood_<gene>_mval.png   one per hit (map / abundance / interactions)
#             neighbourhood_genes_mval.csv        every gene in every window with its estimates
#             neighbourhood_summary_mval.csv      per hit: do the neighbours change too?
# Check:    18 samples, 3 per group, joined by name; 11 protein-coding genes per window
#           with the hit at position 0; hit list = Step 6 hits.

suppressMessages(source("R/setup.R"))
library(patchwork)
out <- "mval_medi/results/step06_responses/hit_neighbourhoods"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
N_SIDE <- 5
RUN_GAP <- 150
INTER <- c("stk vs anc" = "interaction_stk_vs_anc_at_TSB_vs_SN",
           "stp vs anc" = "interaction_stp_vs_anc_at_TSB_vs_SN",
           "stp vs stk" = "interaction_stp_vs_stk_at_TSB_vs_SN")
COMP_COLOURS <- c("stk vs anc" = unname(GENOTYPE_COLOURS["stk"]),
                  "stp vs anc" = unname(GENOTYPE_COLOURS["stp"]),
                  "stp vs stk" = "#CC79A7")

# ---- 1. hits and genome -----------------------------------------------------------
hits_ctx <- read_csv("mval_medi/results/step06_responses/hits_genome_context_mval.csv", show_col_types = FALSE)
hits <- hits_ctx$locus_tag
stopifnot(length(hits) == 8)

gff <- read_tsv("data/raw/resources/JE2_genome.gff", comment = "#", col_names = FALSE,
                col_types = "cccddcccc", progress = FALSE) |>
  filter(X3 == "gene") |>
  transmute(locus_tag = str_match(X9, "locus_tag=([^;]+)")[, 2],
            biotype = str_match(X9, "gene_biotype=([^;]+)")[, 2],
            start = X4, end = X5, strand = X7) |>
  arrange(start)
stopifnot(!anyNA(gff$locus_tag), !anyNA(gff$biotype), all(hits %in% gff$locus_tag))
coding <- gff |> filter(biotype == "protein_coding") |> mutate(idx = row_number())

meta <- read_protein_metadata() |>
  transmute(locus_tag, product,
            gene = if_else(ortholog_quality == "good", gene_symbol_best, NA_character_))
coding <- coding |> left_join(distinct(meta, locus_tag, .keep_all = TRUE), by = "locus_tag") |>
  mutate(tag = str_remove(locus_tag, "B7H15_"),
         short = coalesce(gene, tag),
         x_lab = if_else(is.na(gene), tag, paste0(gene, "\n", tag)))

# ---- 2. measurements and statistics -------------------------------------------------
samples <- read_samples()
abund <- read_abundances("total") |> keep_je2() |>
  select(locus_tag, sample, normalized_abundance) |>
  inner_join(select(samples, sample, genotype, medium), by = "sample") |>   # by NAME (pitfall 2)
  filter(medium %in% c("TSB", "SN")) |>
  mutate(medium = factor(medium, c("TSB", "SN")))
stopifnot(n_distinct(abund$sample) == 18,
          all(count(distinct(abund, sample, genotype, medium), genotype, medium)$n == 3))

stats <- read_tsv("data/derived/limma_treat_total_mval.tsv.gz", show_col_types = FALSE) |>
  filter(contrast %in% INTER) |>
  transmute(locus_tag, comparison = names(INTER)[match(contrast, INTER)],
            diff, std.error, FDR, tested, hit = tested & FDR < MVAL_FDR)
hit_comps <- stats |> filter(hit, locus_tag %in% hits) |>
  group_by(locus_tag) |> summarise(comps = list(comparison), .groups = "drop")
stopifnot(setequal(hit_comps$locus_tag, hits))
all_hits <- unique(stats$locus_tag[stats$hit])   # any interaction hit, genome-wide

# ---- 3. drawing helpers ---------------------------------------------------------------
# gene arrow as a polygon (body + head), so arrows keep their real length
arrow_poly <- function(start, end, strand, y = 0, h = 0.22, head_bp) {
  hd <- min(head_bp, 0.45 * (end - start))
  if (strand == "+") {
    x <- c(start, end - hd, end - hd, end, end - hd, end - hd, start)
  } else {
    x <- c(end, start + hd, start + hd, start, start + hd, start + hd, end)
  }
  tibble(x = x, y = y + c(-h, -h, -1.7 * h, 0, 1.7 * h, h, h))
}
same_strand_run <- function(i) {          # indices (in gff) of the hit's same-strand run
  gi <- which(gff$locus_tag == coding$locus_tag[coding$idx == i])
  run <- gi
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
fmt_kb <- \(d) paste0(if_else(d >= 0, "+", "−"), formatC(abs(d) / 1000, format = "f", digits = 1, big.mark = ","), " kb")

# ---- 4. one figure per hit ----------------------------------------------------------------
gene_rows <- list(); summary_rows <- list()
for (h in hits) {
  i0 <- coding$idx[coding$locus_tag == h]
  win <- coding |> filter(idx %in% (i0 - N_SIDE):(i0 + N_SIDE)) |>
    mutate(offset = idx - i0)
  stopifnot(nrow(win) == 2 * N_SIDE + 1, sum(win$offset == 0) == 1)
  hc <- hits_ctx |> filter(locus_tag == h)
  me <- win |> filter(offset == 0)
  comps <- hit_comps$comps[[which(hit_comps$locus_tag == h)]]
  run <- same_strand_run(i0)

  # per-gene status for colours and the tables
  det <- abund |> filter(locus_tag %in% win$locus_tag) |>
    group_by(locus_tag) |> summarise(n_detected = sum(!is.na(normalized_abundance)), .groups = "drop")
  win <- win |> left_join(det, by = "locus_tag") |>
    mutate(n_detected = coalesce(n_detected, 0L),
           role = case_when(offset == 0 ~ "the hit",
                            locus_tag %in% all_hits ~ "also an interaction hit",
                            n_detected == 0 ~ "protein not detected",
                            TRUE ~ "measured neighbour"),
           in_run = locus_tag %in% run,
           x_lab = factor(x_lab, levels = x_lab))

  # (a) gene map, to scale, including RNA genes in the window
  span <- c(min(win$start), max(win$end))
  pad <- 0.02 * diff(span)
  head_bp <- 0.012 * diff(span)
  # genes fully inside the window (a gene just outside can overlap the edge gene by a few bp)
  map_genes <- gff |> filter(start >= span[1], end <= span[2]) |>
    left_join(select(win, locus_tag, role, short, offset), by = "locus_tag") |>
    mutate(role = coalesce(role, paste("other gene:", biotype)))
  polys <- map_genes |> mutate(id = row_number()) |>
    group_by(id, locus_tag, role) |>
    reframe(arrow_poly(start, end, strand, head_bp = head_bp))
  labs_map <- map_genes |> filter(!is.na(short)) |> arrange(start) |>
    mutate(mid = (start + end) / 2, ylab = if_else(row_number() %% 2 == 1, -0.62, -0.95))
  run_box <- gff |> filter(locus_tag %in% run)
  # a long run can reach beyond the window: clip the box and say so
  run_cut <- min(run_box$start) < span[1] || max(run_box$end) > span[2]
  run_x <- c(max(min(run_box$start), span[1] - pad), min(max(run_box$end), span[2] + pad))
  fill_vals <- c("the hit" = "#D55E00", "also an interaction hit" = "#E69F00",
                 "measured neighbour" = "grey55", "protein not detected" = "grey88")
  other_roles <- setdiff(unique(polys$role), names(fill_vals))
  fill_vals <- c(fill_vals, setNames(rep("white", length(other_roles)), other_roles))
  p_map <- ggplot() +
    { if (length(run) > 1)
        list(annotate("rect", xmin = run_x[1], xmax = run_x[2], ymin = -0.5, ymax = 0.62,
                      fill = "#DCEBF7"),
             annotate("text", x = run_x[1], y = 0.6, hjust = 0, vjust = 1, size = 2.6,
                      colour = "#2B6A99",
                      label = sprintf("same-strand run with the hit (%d genes, gaps < %d bp%s): co-transcription possible, not shown",
                                      length(run), RUN_GAP,
                                      if (run_cut) "; continues beyond this window" else ""))) } +
    geom_polygon(data = polys, aes(x, y, group = id, fill = role), colour = "grey25", linewidth = 0.25) +
    geom_text(data = labs_map, aes(mid, ylab, label = short,
                                   fontface = if_else(offset == 0, "bold", "plain")), size = 2.8) +
    geom_segment(data = labs_map, aes(x = mid, xend = mid, y = -0.42, yend = ylab + 0.12),
                 colour = "grey70", linewidth = 0.25) +
    scale_fill_manual(values = fill_vals, name = NULL) +
    scale_x_continuous(labels = scales::label_comma(), limits = span + c(-pad, pad),
                       name = "JE2 genome position (CP020619.1); arrows point in the direction of transcription (strand)") +
    scale_y_continuous(limits = c(-1.1, 0.65), breaks = NULL, name = NULL) +
    theme_bw(base_size = 9) +
    theme(legend.position = "right", panel.grid.minor = element_blank())

  # (b) per-sample abundance, TSB and SN, strains side by side
  ab <- abund |> filter(locus_tag %in% win$locus_tag) |>
    left_join(select(win, locus_tag, x_lab), by = "locus_tag")
  ab_det <- ab |> group_by(x_lab, medium, genotype) |>
    summarise(n_det = sum(!is.na(normalized_abundance)), .groups = "drop")
  y_low <- floor(min(ab$normalized_abundance, na.rm = TRUE)) - 1
  nd <- win |> filter(n_detected == 0) |> select(x_lab) |> tidyr::crossing(medium = factor(c("TSB", "SN"), c("TSB", "SN")))
  dodge <- position_dodge(width = 0.65)
  p_ab <- ggplot(ab |> filter(!is.na(normalized_abundance)),
                 aes(x_lab, normalized_abundance, colour = genotype, group = genotype)) +
    annotate("rect", xmin = N_SIDE + 0.5, xmax = N_SIDE + 1.5, ymin = -Inf, ymax = Inf, fill = "grey90") +
    geom_point(size = 1.5, position = dodge) +
    geom_text(data = ab_det, aes(x_lab, y_low + 0.35, group = genotype,
                                 # 0/3 is shown too, unless the protein is never detected
                                 # (then the grey "not detected" text covers it)
                                 label = if_else(n_det < 3 & !x_lab %in% nd$x_lab, paste0(n_det, "/3"), "")),
              position = dodge, size = 2.2, colour = "firebrick", inherit.aes = FALSE) +
    geom_text(data = nd, aes(x_lab, y_low + 2, label = "not\ndetected"), size = 2.4, colour = "grey45",
              inherit.aes = FALSE) +
    facet_grid(medium ~ ., labeller = labeller(medium = c(TSB = "TSB (control)", SN = "SN (PASN)"))) +
    scale_x_discrete(drop = FALSE) +
    scale_colour_manual(values = GENOTYPE_COLOURS, name = "strain") +
    scale_y_continuous(limits = c(y_low, NA)) +
    labs(x = NULL, y = "normalised log2 abundance") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "right")

  # (c) the three interactions per gene (estimate +/- 2 SE)
  st <- stats |> filter(locus_tag %in% win$locus_tag) |>
    left_join(select(win, locus_tag, x_lab), by = "locus_tag") |>
    mutate(comparison = factor(comparison, names(INTER)),
           shown = tested & !is.na(diff))
  # not testable: one label per gene when no comparison is testable, otherwise per comparison
  nt <- st |> filter(!shown) |> group_by(x_lab) |> mutate(all_nt = n() == length(INTER)) |> ungroup()
  nt_gene <- nt |> filter(all_nt) |> distinct(x_lab)
  nt <- nt |> filter(!all_nt)
  p_int <- ggplot(st |> filter(shown), aes(x_lab, diff, colour = comparison, group = comparison)) +
    annotate("rect", xmin = N_SIDE + 0.5, xmax = N_SIDE + 1.5, ymin = -Inf, ymax = Inf, fill = "grey90") +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_hline(yintercept = c(-MVAL_LFC, MVAL_LFC), linetype = 2, colour = "grey60") +
    geom_errorbar(aes(ymin = diff - 2 * std.error, ymax = diff + 2 * std.error),
                  width = 0, position = position_dodge(width = 0.6), linewidth = 0.4) +
    geom_point(aes(shape = hit), size = 2, position = position_dodge(width = 0.6), fill = "white") +
    geom_text(data = nt, aes(x_lab, 0, label = "n.t.", group = comparison),
              position = position_dodge(width = 0.6), size = 2, colour = "grey45", inherit.aes = FALSE) +
    geom_text(data = nt_gene, aes(x_lab, 0, label = "not\ntestable"), size = 2.3, colour = "grey45",
              lineheight = 0.85, inherit.aes = FALSE) +
    scale_x_discrete(drop = FALSE) +
    scale_colour_manual(values = COMP_COLOURS, name = "interaction") +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), labels = c(`TRUE` = "formal hit", `FALSE` = "not a hit"),
                       name = NULL) +
    labs(x = sprintf("protein-coding genes in genome order: the hit (grey band) with %d before and %d after", N_SIDE, N_SIDE),
         y = "interaction (log2)\n> 0 = SN response lower in\nthe first strain") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(size = 7.5, lineheight = 0.9), legend.position = "right")

  title <- sprintf("%s (%s): %s", coalesce(me$gene, me$tag), h, me$product)
  subtitle <- paste0(
    "Hit in: ", paste(comps, collapse = ", "), ".  Strand: ",
    if_else(me$strand == "+", "forward (+)", "reverse (−)"),
    ".  Start codon to stk1 insertion ", fmt_kb(hc$dist_to_stk1_Tn_bp), ", to stp1 insertion ",
    fmt_kb(hc$dist_to_stp1_Tn_bp), if_else(hc$path_to_Tn == "direct", "", " (shortest way via position 1)"), ".\n",
    "(a) gene map  (b) every sample, TSB and SN (red n/3 = replicates detected if < 3; missing = not plotted)  ",
    "(c) Step 5 interactions, estimate ± 2 SE; filled = formal hit (treat lfc = 1, BH-FDR < 0.05); n.t. = not testable.")
  fig <- (p_map / p_ab / p_int) +
    plot_layout(heights = c(1.15, 2.2, 1.35)) +
    plot_annotation(title = title, subtitle = subtitle, tag_levels = "a",
                    theme = theme(plot.title = element_text(face = "bold", size = 12),
                                  plot.subtitle = element_text(size = 8.5)))
  fname <- sprintf("fig_neighbourhood_%s_%s_mval.png", coalesce(me$gene, "x"), me$tag) |> str_replace("_x_", "_")
  ggsave(file.path(out, fname), fig, width = 12, height = 10.5, dpi = 250)
  message("written ", fname)

  # tables
  st_w <- st |> select(locus_tag, comparison, diff, FDR, tested) |>
    mutate(comparison = str_replace(comparison, " vs ", "_vs_")) |>
    pivot_wider(names_from = comparison, values_from = c(diff, FDR, tested), names_glue = "{comparison}_{.value}")
  gene_rows[[h]] <- win |>
    transmute(hit = h, hit_gene = coalesce(me$gene, me$tag), offset, locus_tag, gene, product, strand,
              start, end, role, in_same_strand_run_with_hit = in_run, n_detected_of_18 = n_detected) |>
    left_join(st_w, by = "locus_tag")
  nb <- st |> filter(locus_tag != h, as.character(comparison) %in% comps)
  summary_rows[[h]] <- tibble(
    hit = h, hit_gene = coalesce(me$gene, me$tag), hit_in = paste(comps, collapse = "; "),
    neighbours_detected = sum(win$offset != 0 & win$n_detected > 0),
    neighbours_testable = n_distinct(nb$locus_tag[nb$shown]),
    neighbour_formal_hits = paste(unique(win$short[win$locus_tag %in% nb$locus_tag[nb$hit]]), collapse = ", "),
    neighbours_estimate_above_1 = paste(unique(win$short[win$locus_tag %in% nb$locus_tag[nb$shown & abs(nb$diff) > MVAL_LFC]]),
                                        collapse = ", "),
    same_strand_run = paste(coding$short[match(run, coding$locus_tag)], collapse = ", "))
}

genes_tbl <- bind_rows(gene_rows)
summary_tbl <- bind_rows(summary_rows)
stopifnot(all(count(genes_tbl, hit)$n == 2 * N_SIDE + 1))
write_csv(genes_tbl, file.path(out, "neighbourhood_genes_mval.csv"))
write_csv(summary_tbl, file.path(out, "neighbourhood_summary_mval.csv"))
print(summary_tbl, width = 250)
