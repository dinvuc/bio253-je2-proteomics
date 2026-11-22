# Step 11c (PROJECT_WORKFLOW_mval.md, after Step 11b): genome context of regions R7 and R6
# at baseline TSB. Do the changed proteins sit in one operon, many operons, or a whole
# stretch of DNA? Can a polar effect of the transposon explain them?
#
# Run from the repository root:   Rscript C_mutant_validation_DinMete/scripts/11c_R6_R7_context_mval.R   (~20 s)
# Needs C_mutant_validation_DinMete/scripts/05 (limma fit) and 11b (region table) first.
#
# ---- Plan (coding rules, section 4) -----------------------------------------
# Aim:      R7 (sigB, mazEF, kdp; flagged in stk_vs_anc_at_TSB) and R6 (prophage at
#           2.09 Mb; flagged in stp_vs_anc_at_TSB) are the two baseline-TSB blocks of the
#           Step 11 genome scan. Ask: (1) which putative transcription units (TUs) do the
#           changed proteins belong to? (2) Genome-wide: do stk / stp effects follow TUs,
#           i.e. do neighbours in the same TU change together more than neighbours with their
#           own promoters? (3) Could a polar effect of the stk1/stp1 transposon reach R6/R7?
#           (4) Is R7 "more SigB activity" or "more of the whole DNA stretch"?
# Input:    data/derived/limma_treat_total_mval.tsv.gz (Step 5 fit: diff, tested)
#           read_protein_metadata(cds_only = TRUE): start, end, strand, gene symbol
#           read_abundances("total"): per-sample normalised log2 abundance (NA = not seen)
#           C_mutant_validation_DinMete/results/step11_robustness/genome_regions_mval.csv (Step 11b regions)
# Analysis: Putative TU = run of neighbouring CDS on the same strand with gaps < 150 bp
#           (same rule as Steps 6b/6c). This is a HINT of co-transcription from the genome
#           sequence only, not a measured operon. (a) gene maps + per-gene diffs for R7 and R6;
#           (b) per-sample block score = mean over the block's proteins of
#           log2(sample) - mean log2(anc TSB), so replicates and outliers are visible;
#           (c) Pearson r of the diffs of adjacent gene pairs, split by pair type
#           (same TU / same strand but gap >= 150 bp / convergent / divergent); divergent
#           pairs (two separate promoters) are the control; repeated without R6 and R7;
#           (d) SigB-operon proteins (inside R7) vs SigB-regulon reporters (outside R7).
# Output:   C_mutant_validation_DinMete/results/step11_robustness/
#             fig_R6_R7_genome_context_mval.png   the figure
#             fig_R6_R7_genome_context_compact_mval.png   panels A + B only, narrow (poster)
#             R6_R7_genes_mval.csv                gene table: TU, strand, diffs, estimate type
#             R6_R7_tu_effects_mval.csv           TU means per contrast
#             operon_coherence_mval.csv           r per pair type and contrast
#             R6_R7_context_checks_mval.txt
# Check:    R6 and R7 tag ranges equal Step 11b; R7 mean in stk TSB reproduces Step 11b (1.12);
#           hlb is split by the R6 prophage (hlb on both sides); every block protein used in
#           the per-sample score is measured in >= 1 anc TSB replicate.

suppressMessages(source("R/setup.R"))
library(patchwork)
out <- "C_mutant_validation_DinMete/results/step11_robustness"
checks <- character()
check <- function(label, ok) {
  line <- sprintf("[%s] %s", if (isTRUE(ok)) "PASS" else "FAIL", label)
  message(line); checks <<- c(checks, line); invisible(ok)
}
note <- function(...) { line <- paste0("[NOTE] ", ...); message(line); checks <<- c(checks, line) }

RUN_GAP <- 150
COL <- c(stk = "#E69F00", stp = "#009E73")
CTS <- c("stk_vs_anc_at_TSB", "stp_vs_anc_at_TSB", "stk_vs_anc_at_SN", "stp_vs_anc_at_SN")
tag_in <- \(tag, a, b) tag >= a & tag <= b    # B7H15_nnnnn: text order = genome order

# ---- 1. inputs + putative transcription units ------------------------------------------
regions <- read_csv(file.path(out, "genome_regions_mval.csv"), show_col_types = FALSE)
R <- regions |> filter(region %in% c("R6", "R7")) |> select(region, first_tag, last_tag)
check("R6 / R7 ranges as in Step 11b",
      identical(R$first_tag, c("B7H15_10910", "B7H15_11470")) && identical(R$last_tag, c("B7H15_11150", "B7H15_11590")))

fit <- read_tsv("data/derived/limma_treat_total_mval.tsv.gz", show_col_types = FALSE) |>
  filter(str_starts(protein_key, "lcl_"), contrast %in% CTS)   # drop the frameshift re-entries
check("one fit row per locus tag and contrast", !anyDuplicated(fit[c("locus_tag", "contrast")]))

gff <- read_protein_metadata(cds_only = TRUE) |> filter(!is.na(start)) |> arrange(start) |>
  transmute(locus_tag, start, end, strand,
            gene = if_else(ortholog_quality %in% "good", gene_symbol_best, NA_character_),
            sym = gene_symbol_best, product) |>
  mutate(gap = start - lag(end) - 1,
         tu = cumsum(is.na(gap) | strand != lag(strand) | gap >= RUN_GAP))
note(sprintf("%d CDS in %d putative TUs (%d with > 1 gene); rule: same strand, gap < %d bp",
             nrow(gff), max(gff$tu), sum(table(gff$tu) > 1), RUN_GAP))

wide <- fit |> select(locus_tag, contrast, diff, tested, estimate_type) |>
  mutate(diff = if_else(tested, diff, NA_real_)) |>
  pivot_wider(names_from = contrast, values_from = c(diff, tested, estimate_type))
g <- gff |> left_join(wide, by = "locus_tag")

# ---- 2. the two regions: gene table + TU means --------------------------------------------
VIEW <- tribble(          # view = region plus flanks, so the edges of the effect are visible
  ~region, ~from,         ~to,           ~focus,
  "R7",    "B7H15_11380", "B7H15_11700", "stk",
  "R6",    "B7H15_10815", "B7H15_11185", "stp")
genes <- purrr::pmap(VIEW, \(region, from, to, focus) g |> filter(tag_in(locus_tag, from, to)) |>
  mutate(view = region, in_region = tag_in(locus_tag, R$first_tag[R$region == region], R$last_tag[R$region == region]))) |>
  bind_rows()
write_csv(genes |> select(view, in_region, locus_tag, start, end, strand, gap, tu, gene, product,
                          starts_with("diff_"), starts_with("estimate_type_")) |>
            mutate(across(starts_with("diff_"), \(x) round(x, 2))),
          file.path(out, "R6_R7_genes_mval.csv"))

r7_mean <- mean(genes$diff_stk_vs_anc_at_TSB[genes$view == "R7" & genes$in_region], na.rm = TRUE)
check(sprintf("R7 mean in stk_vs_anc_at_TSB reproduces Step 11b (%.2f)", r7_mean), abs(r7_mean - 1.12) < 0.01)
r6_mean <- mean(genes$diff_stp_vs_anc_at_TSB[genes$view == "R6" & genes$in_region], na.rm = TRUE)
check(sprintf("R6 mean in stp_vs_anc_at_TSB reproduces Step 11b (%.2f)", r6_mean), abs(r6_mean - (-1.09)) < 0.01)
check("R6 prophage sits inside hlb (hlb fragment on both sides)",
      sum(genes$sym[genes$view == "R6"] %in% "hlb") == 2)
check("immune-evasion genes scn, chp present at the R6 prophage edge",
      all(c("scn", "chp") %in% genes$sym[genes$view == "R6"]))
note("R6 = beta-haemolysin-converting prophage (phiSa3-like): integrase + scn/chp/staphylokinase, inserted in hlb")
lod <- genes |> filter(in_region) |> select(starts_with("estimate_type_")) |> unlist()
note(sprintf("estimate_type in R6/R7 (TSB + SN, tested proteins): %s",
             paste(names(table(lod)), table(lod), sep = " = ", collapse = ", ")))

tu_eff <- genes |> group_by(view, tu) |>
  summarise(first_tag = first(locus_tag), last_tag = last(locus_tag), strand = first(strand),
            n_genes = n(), in_region = any(in_region),
            genes = paste(na.omit(gene), collapse = ","),
            across(starts_with("diff_"), list(n = \(x) sum(!is.na(x)), mean = \(x) mean(x, na.rm = TRUE))),
            .groups = "drop") |>
  mutate(across(ends_with("_mean"), \(x) round(if_else(is.nan(x), NA_real_, x), 2)))
write_csv(tu_eff, file.path(out, "R6_R7_tu_effects_mval.csv"))
r7_tus <- tu_eff |> filter(view == "R7", in_region, diff_stk_vs_anc_at_TSB_n > 0)
note(sprintf("R7: %d putative TUs with tested proteins, %d on '+' and %d on '-' strand; TU means in stk TSB %.2f to %.2f",
             nrow(r7_tus), sum(r7_tus$strand == "+"), sum(r7_tus$strand == "-"),
             min(r7_tus$diff_stk_vs_anc_at_TSB_mean), max(r7_tus$diff_stk_vs_anc_at_TSB_mean)))
r7_bp <- genes |> filter(view == "R7", in_region) |> summarise(bp = max(end) - min(start) + 1) |> pull(bp)
note(sprintf("R7 spans %.1f kb", r7_bp / 1000))

# polar effect: distance of each region to the transposon insertions
tn <- TRANSPOSON_INSERTIONS
dist_kb <- genes |> filter(in_region) |> group_by(view) |>
  summarise(start = min(start), end = max(end)) |>
  mutate(kb_to_stk1_Tn = round((start - tn$insertion_pos[tn$short_name == "Stk1"]) / 1000),
         kb_to_stp1_Tn = round((start - tn$insertion_pos[tn$short_name == "Stp1"]) / 1000))
for (i in seq_len(nrow(dist_kb))) note(sprintf("%s starts %d kb from the stk1 transposon, %d kb from the stp1 transposon",
             dist_kb$view[i], dist_kb$kb_to_stk1_Tn[i], dist_kb$kb_to_stp1_Tn[i]))
gap_tn <- gff |> filter(start > tn$cds_end[2], start < min(dist_kb$start)) |>
  summarise(n_cds = n(), n_tu = n_distinct(tu))
note(sprintf("between stk1 and R6 lie %d CDS in %d putative TUs: a polar effect only reaches the downstream genes of the SAME transcript",
             gap_tn$n_cds, gap_tn$n_tu))

# ---- 3. per-sample block score (plot raw values before trusting a fold change) ---------------
ab <- read_abundances("total") |> filter(str_starts(protein_key, "lcl_"))
block_score <- purrr::map(R$region, \(rg) {
  tags <- genes |> filter(view == rg, in_region) |> pull(locus_tag)
  x <- ab |> filter(locus_tag %in% tags)
  ref <- x |> filter(group == "anc_TSB") |> group_by(locus_tag) |>
    summarise(ref = mean(normalized_abundance, na.rm = TRUE), n_ref = sum(!is.na(normalized_abundance)))
  x |> inner_join(ref |> filter(n_ref >= 1), by = "locus_tag") |>
    group_by(sample, genotype, medium, replicate, reinjected) |>
    summarise(score = mean(normalized_abundance - ref, na.rm = TRUE),
              n_prot = sum(!is.na(normalized_abundance)), .groups = "drop") |>
    mutate(region = rg)
}) |> bind_rows()
check("per-sample score built from >= 10 proteins in every sample", min(block_score$n_prot) >= 10)
print(block_score |> filter(medium == "TSB") |> select(region, sample, score, n_prot) |> arrange(region, sample))

# ---- 4. genome-wide: do effects follow putative TUs? --------------------------------------------
pairs <- tibble(a = head(g$locus_tag, -1), b = g$locus_tag[-1],
                type = case_when(g$tu[-1] == head(g$tu, -1) ~ "same TU",
                                 g$strand[-1] == head(g$strand, -1) ~ "same strand, gap >= 150 bp",
                                 head(g$strand, -1) == "+" ~ "convergent  → ←",
                                 TRUE ~ "divergent  ← →  (2 promoters)"))
in_r67 <- \(t) tag_in(t, R$first_tag[1], R$last_tag[1]) | tag_in(t, R$first_tag[2], R$last_tag[2])
coh <- tidyr::expand_grid(contrast = CTS, set = c("whole genome", "without R6 and R7")) |>
  mutate(res = purrr::map2(contrast, set, \(ct, s) {
    d <- setNames(g[[paste0("diff_", ct)]], g$locus_tag)
    p <- pairs |> mutate(da = d[a], db = d[b]) |> filter(!is.na(da), !is.na(db))
    if (s != "whole genome") p <- p |> filter(!in_r67(a), !in_r67(b))
    p |> group_by(type) |> summarise(n_pairs = n(), r = cor(da, db),
                                      lo = cor.test(da, db)$conf.int[1], hi = cor.test(da, db)$conf.int[2], .groups = "drop")
  })) |> tidyr::unnest(res) |>
  mutate(mutant = str_sub(contrast, 1, 3), medium = str_extract(contrast, "(TSB|SN)$"))
write_csv(coh |> mutate(across(c(r, lo, hi), \(x) round(x, 3))), file.path(out, "operon_coherence_mval.csv"))
print(coh |> filter(set == "without R6 and R7"), n = 40)
st <- \(m, md, tp, s = "without R6 and R7") coh |> filter(mutant == m, medium == md, type == tp, set == s)
check("stk TSB: same-TU pairs more alike than divergent pairs (95% CIs do not overlap), without R6/R7",
      st("stk", "TSB", "same TU")$lo > st("stk", "TSB", "divergent  ← →  (2 promoters)")$hi)

# ---- 5. R7 = more SigB activity, or more of the whole stretch? ---------------------------------
# SigB-dependent reporters outside R7: asp23 (alkaline shock protein 23) and the crt carotenoid
# operon (staphyloxanthin) are classic SigB targets; the sigB operon (rsbU-V-W-sigB) is inside R7.
sigb <- g |> filter(locus_tag %in% c("B7H15_11470", "B7H15_11475", "B7H15_11480", "B7H15_11485") |
                    gene %in% c("asp23", "crtM", "crtN")) |>
  mutate(name = case_when(locus_tag == "B7H15_11470" ~ "sigB", locus_tag == "B7H15_11485" ~ "rsbU",
                          TRUE ~ gene),
         set = if_else(tag_in(locus_tag, R$first_tag[2], R$last_tag[2]),
                       "sigB operon (inside R7)", "SigB targets (outside R7)")) |>
  select(locus_tag, name, set, diff_stk_vs_anc_at_TSB, diff_stp_vs_anc_at_TSB)
check("sigB-operon and 3 SigB reporters found", nrow(sigb) == 7)
print(sigb)
note(sprintf("stk TSB: sigB operon mean %+.2f, SigB reporters (asp23, crtM, crtN) mean %+.2f",
             mean(sigb$diff_stk_vs_anc_at_TSB[str_starts(sigb$set, "sigB")]),
             mean(sigb$diff_stk_vs_anc_at_TSB[str_starts(sigb$set, "SigB")])))

# ---- 6. figure ---------------------------------------------------------------------------------
arrow_poly <- function(start, end, strand, y = 0, h = 0.22, head_bp) {   # as in 06c
  hd <- min(head_bp, 0.45 * (end - start))
  x <- if (strand == "+") c(start, end - hd, end - hd, end, end - hd, end - hd, start)
       else               c(end, start + hd, start + hd, start, start + hd, start + hd, end)
  tibble(x = x, y = y + c(-h, -h, -1.7 * h, 0, 1.7 * h, h, h))
}
LIM <- 3

region_panel <- function(rg, focus, title, label_tags, block_label = paste0(" ", rg, " (Step 11 block)")) {
  d <- genes |> filter(view == rg) |> mutate(mid = (start + end) / 2 / 1e6,
    fdiff = .data[[paste0("diff_", focus, "_vs_anc_at_TSB")]])
  span <- range(c(d$start, d$end)) / 1e6
  hb <- 0.012 * diff(span)
  polys <- d |> group_by(locus_tag, fdiff) |>
    reframe(arrow_poly(start / 1e6, end / 1e6, strand, head_bp = hb))
  reg <- d |> filter(in_region) |> summarise(x0 = min(start) / 1e6, x1 = max(end) / 1e6)
  # putative TUs with > 1 gene: bracket under the arrows
  tus <- d |> group_by(tu) |> filter(n() > 1) |>
    summarise(x0 = min(start) / 1e6, x1 = max(end) / 1e6, .groups = "drop")
  labs_g <- d |> filter(locus_tag %in% label_tags | (sym %in% c("hlb", "scn", "chp"))) |>
    mutate(lab = case_when(locus_tag == "B7H15_11470" ~ "sigB", locus_tag == "B7H15_11485" ~ "rsbU",
                           locus_tag == "B7H15_10870" ~ "sak", locus_tag == "B7H15_11170" ~ "int",
                           locus_tag == "B7H15_10965" ~ "capsid", locus_tag == "B7H15_10975" ~ "portal",
                           locus_tag == "B7H15_10980" ~ "terminase", locus_tag == "B7H15_10880" ~ "holin",
                           TRUE ~ coalesce(gene, sym)),
           ylab = c(0.5, 0.72, 0.94)[(row_number() - 1) %% 3 + 1])
  p_map <- ggplot() +
    annotate("rect", xmin = reg$x0, xmax = reg$x1, ymin = -0.62, ymax = 1.05, fill = "#FBE3D6", alpha = 0.7) +
    annotate("text", x = reg$x0, y = -0.6, label = block_label, hjust = 0, vjust = 0,
             size = 3, colour = "#A0400B", fontface = "bold") +
    geom_segment(data = tus, aes(x = x0, xend = x1, y = -0.48, yend = -0.48), colour = "#2B6A99", linewidth = 0.8) +
    geom_polygon(data = polys, aes(x, y, group = locus_tag, fill = pmax(pmin(fdiff, LIM), -LIM)),
                 colour = "grey25", linewidth = 0.2) +
    geom_text(data = labs_g, aes(mid, ylab, label = lab), size = 2.7, fontface = "italic") +
    geom_segment(data = labs_g, aes(x = mid, xend = mid, y = 0.4, yend = ylab - 0.1), colour = "grey70", linewidth = 0.2) +
    scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", limits = c(-LIM, LIM), na.value = "grey85",
                         name = sprintf("%s vs anc, TSB\nlog2 diff (grey =\nnot tested)", focus)) +
    scale_x_continuous(limits = span, expand = expansion(mult = 0.01)) +
    scale_y_continuous(limits = c(-0.62, 1.05), breaks = NULL) +
    labs(x = NULL, y = NULL, title = title) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(),
          plot.title = element_text(size = 10.5, face = "bold"))
  pts <- d |> select(locus_tag, mid, stk = diff_stk_vs_anc_at_TSB, stp = diff_stp_vs_anc_at_TSB) |>
    pivot_longer(c(stk, stp), names_to = "mutant", values_to = "diff") |> filter(!is.na(diff))
  tu_m <- pts |> left_join(d |> select(locus_tag, tu, start, end), by = "locus_tag") |>
    group_by(tu, mutant) |> summarise(x0 = min(start) / 1e6, x1 = max(end) / 1e6, m = mean(diff), .groups = "drop")
  p_dot <- ggplot(pts, aes(mid, diff, colour = mutant)) +
    annotate("rect", xmin = reg$x0, xmax = reg$x1, ymin = -Inf, ymax = Inf, fill = "#FBE3D6", alpha = 0.7) +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_segment(data = tu_m, aes(x = x0, xend = x1, y = m, yend = m, colour = mutant), linewidth = 1.6, alpha = 0.45,
                 inherit.aes = FALSE) +
    geom_point(size = 1.6) +
    scale_colour_manual(values = COL, labels = c(stk = "stk vs anc", stp = "stp vs anc"), name = "TSB") +
    scale_x_continuous(limits = span, expand = expansion(mult = 0.01)) +
    labs(x = "JE2 genome position (Mb)", y = "log2 diff\n(> 0: higher in mutant)") +
    theme(legend.position = "right")
  p_map / p_dot + plot_layout(heights = c(1, 1.1))
}

pA <- region_panel("R7", "stk",
  sprintf("A  R7 (2.19 Mb): in stk %d of %d putative TUs (both strands) are higher by > 0.5 log2 (sigB, mazEF, kdp, murF-ddl); stp is flat",
    sum(r7_tus$diff_stk_vs_anc_at_TSB_mean > 0.5), nrow(r7_tus)),
  c("B7H15_11470", "B7H15_11475", "B7H15_11480", "B7H15_11485", "B7H15_11495", "B7H15_11500", "B7H15_11505",
    "B7H15_11530", "B7H15_11540", "B7H15_11550", "B7H15_11555", "B7H15_11565", "B7H15_11570", "B7H15_11390", "B7H15_11620"))
pB <- region_panel("R6", "stp",
  "B  R6 (2.09 Mb): a prophage inside hlb (phiSa3-like). Most phage proteins (several '-' strand TUs) are lower in stp, not in stk; a few are higher",
  c("B7H15_10870", "B7H15_10880", "B7H15_10965", "B7H15_10975", "B7H15_10980", "B7H15_11170"))

bs <- block_score |> mutate(medium = factor(medium, levels = c("TSB", "SN", "SNpulse")),
                            genotype = factor(genotype, levels = c("anc", "stk", "stp")),
                            region = factor(region, levels = c("R7", "R6"),
                                            labels = c("R7", "R6")))
pC <- ggplot(bs, aes(genotype, score, colour = genotype)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_point(aes(shape = reinjected), size = 2, position = position_jitter(width = 0.12, seed = 1)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.3, colour = "grey20") +
  facet_grid(region ~ medium) +
  scale_colour_manual(values = c(anc = "grey40", COL), guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
  labs(x = NULL, y = "block score per sample\n(mean log2 vs anc TSB)",
       title = "C  Raw values: every replicate, all media",
       subtitle = "R7 stays up in stk in SN and SNpulse too") +
  theme(plot.title = element_text(size = 10.5, face = "bold"), plot.subtitle = element_text(size = 8.5),
        strip.text = element_text(size = 8))

cd <- coh |> filter(medium == "TSB") |>
  mutate(type = factor(type, levels = rev(c("same TU", "same strand, gap >= 150 bp", "convergent  → ←",
                                             "divergent  ← →  (2 promoters)"))))
pD <- ggplot(cd, aes(r, type, colour = mutant, shape = set)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6), size = 0.3) +
  scale_colour_manual(values = COL, labels = c(stk = "stk vs anc", stp = "stp vs anc"), name = "TSB") +
  scale_shape_manual(values = c(`whole genome` = 16, `without R6 and R7` = 1), name = NULL) +
  labs(x = "correlation r of the two genes' log2 diffs (95% CI)", y = "adjacent gene pair",
       title = "D  Genome-wide: do effects follow transcription units?",
       subtitle = "Same-TU neighbours change together; neighbours with their own promoters do not") +
  theme(plot.title = element_text(size = 10.5, face = "bold"), plot.subtitle = element_text(size = 8.5),
        legend.position = "bottom", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0))

sg <- sigb |> pivot_longer(starts_with("diff_"), names_to = "mutant", values_to = "diff") |>
  mutate(mutant = str_sub(str_remove(mutant, "diff_"), 1, 3),
         name = factor(name, levels = c("rsbU", "rsbV", "rsbW", "sigB", "asp23", "crtM", "crtN")))
pE <- ggplot(sg, aes(name, diff, fill = mutant)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  facet_grid(~set, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = COL, labels = c(stk = "stk vs anc", stp = "stp vs anc"), name = "TSB") +
  labs(x = NULL, y = "log2 diff (> 0: higher in mutant)",
       title = "E  More SigB protein, but not more SigB activity",
       subtitle = "In stk the operon is ~2x, its targets are unchanged") +
  theme(plot.title = element_text(size = 10.5, face = "bold"), plot.subtitle = element_text(size = 8.5),
        axis.text.x = element_text(face = "italic"), strip.text = element_text(size = 8), legend.position = "bottom")

cap <- paste0(
  "Putative TU (blue bar) = neighbouring genes on the same strand with gaps < 150 bp: a hint of co-transcription from the sequence, not a measured operon. ",
  "Faint bars in the dot tracks = TU means. All proteins in R6/R7 are observed in both groups (no LOD imputation). n = 3 per group; Step 11 scan is permutation-tested genome-wide; ",
  "the TU-level numbers here are descriptive. Polar effect: the transposons sit in stp1/stk1 at 1.23 Mb, ",
  sprintf("%d kb (%d putative TUs) before R6 and %d kb before R7. ", dist_kb$kb_to_stk1_Tn[dist_kb$view == "R6"], gap_tn$n_tu,
          dist_kb$kb_to_stk1_Tn[dist_kb$view == "R7"]),
  "A polar effect only reaches genes downstream in the SAME transcript, so it cannot reach R6 or R7. ",
  "Note: stp also lacks most Stk1 (polar effect inside the stp1-stk1 operon), yet R7 is flat in stp, so R7 does not simply follow Stk1 loss. ",
  "Source for reduced SigB activity in a USA300 pknB/stk1 mutant: Tamber et al. 2010, doi:10.1128/IAI.00296-10.")

bottom <- (pC | pD | pE) + plot_layout(widths = c(1.15, 1, 0.9))
fig <- wrap_elements(full = pA) / wrap_elements(full = pB) / wrap_elements(full = bottom) +
  plot_layout(heights = c(1.2, 1.2, 1.25)) +
  plot_annotation(
    title = "Step 11c: R7 and R6 at baseline TSB: one operon, many operons, or a stretch of DNA?",
    subtitle = paste0("R7: in stk the whole ~", round(r7_bp / 1000), " kb block is up across ", nrow(r7_tus), " putative TUs on both strands, including the sigB operon, but SigB targets are not up (E): ",
                      "this looks like a REGIONAL effect (e.g. a copy-number change or a local regulator) rather than more SigB activity.\n",
                      "R6: in stp, most proteins of the hlb prophage are lower (capsid, portal, replication module), not in stk: a prophage-level change, e.g. less spontaneous induction (hypothesis). Neither block can be a polar effect of the transposon (> 850 kb away)."),
    caption = str_wrap(cap, 230),
    theme = theme(plot.title = element_text(size = 14, face = "bold"), plot.subtitle = element_text(size = 9.5),
                  plot.caption = element_text(hjust = 0, size = 8)))
ggsave(file.path(out, "fig_R6_R7_genome_context_mval.png"), fig, width = 16, height = 15, dpi = 200)

# ---- 7. compact A + B for the poster -----------------------------------------------------------
# Same data and genome windows, drawn narrow: fewer gene labels (they would collide), shorter
# titles and one shared legend below, so no width is lost on the right.
compact <- function(rg, focus, title, label_tags) {
  p <- region_panel(rg, focus, title, label_tags, block_label = "")
  list(p[[1]] + labs(y = sprintf("colour:\n%s vs anc", focus)) +
         guides(fill = guide_colourbar(title = "log2 diff (grey = not tested)", order = 1,
                                       theme = theme(legend.key.width = unit(2.6, "cm"),
                                                     legend.key.height = unit(0.3, "cm")))),
       p[[2]] + labs(x = if (rg == "R6") "JE2 genome position (Mb)",
                     y = "log2 diff\n(> 0: higher\nin mutant)"))
}
cA <- compact("R7", "stk",
  sprintf("A  R7 (2.19 Mb): higher in stk, flat in stp\n     %d of %d putative operons up by > 0.5 log2, on both strands",
          sum(r7_tus$diff_stk_vs_anc_at_TSB_mean > 0.5), nrow(r7_tus)),
  c("B7H15_11470", "B7H15_11500", "B7H15_11540", "B7H15_11565"))
cB <- compact("R6", "stp",
  "B  R6 (2.09 Mb): prophage in hlb, mostly lower in stp, not in stk\n     a few phage proteins are higher",
  c("B7H15_10870", "B7H15_10965", "B7H15_10975", "B7H15_11170"))
fig_compact <- wrap_plots(c(cA, cB), ncol = 1, heights = c(1, 1, 1, 1), guides = "collect") +
  plot_annotation(caption = "TSB, n = 3. Shaded = block. Blue bar = putative operon. Faint bar = operon mean.",
                  theme = theme(plot.caption = element_text(hjust = 0, size = 8.5))) &
  theme(legend.position = "bottom", legend.box = "horizontal", legend.margin = margin(0, 4, 0, 4),
        legend.title = element_text(size = 9), axis.title.y = element_text(size = 8.5),
        plot.margin = margin(2, 6, 2, 2))
ggsave(file.path(out, "fig_R6_R7_genome_context_compact_mval.png"), fig_compact, width = 7, height = 6.4, dpi = 400)

writeLines(c(sprintf("Step 11c checks, %s", Sys.Date()), checks), file.path(out, "R6_R7_context_checks_mval.txt"))
message("done: ", out)
