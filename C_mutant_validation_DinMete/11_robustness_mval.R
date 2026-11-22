# Step 11 (PROJECT_WORKFLOW_mval.md): robustness checks of the headline conclusions.
#
# Run from the repository root:   Rscript C_mutant_validation_DinMete/11_robustness_mval.R   (~3-5 min)
# Needs C_mutant_validation_DinMete/05 (formal test), 06, 07 and 10b first.
#
# ---- Plan (coding rules, section 4) -----------------------------------------
# Aim:      Do the headline conclusions of Steps 3-10 survive reasonable checks, or do
#           they rest on one replicate, a marginal threshold, imputed values or a
#           single genomic region?
# Input:    read_abundances("total") |> keep_je2(), read_samples() (joined by name)
#           data/derived/limma_treat_total_mval.tsv.gz (Step 5, main fit)
#           read_protein_metadata() (genome start, strand)
# Analysis (all labelled SENSITIVITY; the main results and hit rule stay as they are):
#   1. Leave-one-sample-out (LOSO): drop each of the 18 primary samples in turn,
#      refit the Step 5 model (limma ~0+group, treat lfc = 1, BH within contrast,
#      >= 2 of 3 per group). For every main formal hit in the H1-H3 contrasts:
#      in how many of the 18 refits is it still a hit? Which drop removes it?
#      "stable" = hit in all 18; "fragile" = lost in >= 1 refit.
#      (With n = 3, dropping a sample leaves 2: the test loses power, so some loss
#      is expected even for real effects. Read together with the estimate.)
#   2. Thresholds fixed in the workflow as alternatives, not as replacements:
#      (a) treat lfc = log2(1.5) with BH-FDR < 0.05, (b) treat lfc = 1 with FDR < 0.1.
#      Question: does the SHAPE of the H2/H3 answer change (stk vs stp hit counts)?
#   3. Genome-position scan (prompted by Step 10: the SigB/mazEF block
#      B7H15_11465-11570 is higher in stk): tested proteins ordered by genome start;
#      mean diff in windows of 10 consecutive proteins; null = 1000 random
#      PERMUTATIONS of the diffs over the positions (keeps the value distribution,
#      breaks neighbourhood). (Fixed 17 Nov: a first version used circular shifts,
#      which keep every neighbourhood intact, so the null max always equalled the
#      observed max and nothing could be flagged.) Operons and regulons are real
#      neighbourhood structure, so a flagged block means "neighbours more alike than
#      random", NOT "artefact". Window is flagged when |mean| exceeds the 95th
#      percentile of the genome-wide MAXIMUM |window mean| under the null
#      (family-wise 5 % per contrast). Windows of 10 are a fixed choice, stated here.
#   4. Headline table: conclusion, evidence, checks, status (stable / tentative /
#      unresolved), written in figure_guide_mval.md from the numbers produced here.
# Output:   C_mutant_validation_DinMete/results/step11_robustness/
#             loso_hit_stability_mval.csv     per main hit: n refits still a hit, dropping sample that kills it
#             threshold_sensitivity_mval.csv  hits per contrast under main and 2 alternative rules
#             genome_windows_mval.csv         flagged windows (proteins, positions, mean diff)
#             genome_top_windows_mval.csv     most extreme window per contrast + family-wise p
#             genome_block_members_mval.csv   proteins inside the flagged blocks, with their effects
#             fig_loso_stability_mval.png, fig_genome_scan_mval.png
#             figure_guide_mval.md, robustness_checks_mval.txt
# Check:    18 LOSO fits with 26 samples each, every group >= 2 samples; the full refit
#           reproduces the Step 5 hits exactly; windows contain 10 proteins.

suppressMessages(source("R/setup.R"))
out <- "C_mutant_validation_DinMete/results/step11_robustness"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(20261117)
checks <- character()
check <- function(label, ok) {
  line <- sprintf("[%s] %s", if (isTRUE(ok)) "PASS" else "FAIL", label)
  message(line); checks <<- c(checks, line)
  invisible(ok)
}
note <- function(...) { line <- paste0("[NOTE] ", ...); message(line); checks <<- c(checks, line) }

H_CONTRASTS <- c("stk_vs_anc_at_TSB", "stp_vs_anc_at_TSB", "stp_vs_stk_at_TSB",
                 "interaction_stk_vs_anc_at_TSB_vs_SN", "interaction_stp_vs_anc_at_TSB_vs_SN",
                 "interaction_stp_vs_stk_at_TSB_vs_SN")
PRIMARY <- read_csv("C_mutant_validation_DinMete/results/step05_formal_test/contrast_sheet_mval.csv", show_col_types = FALSE) |>
  filter(primary_18_samples) |> pull(contrast)

samples <- read_samples() |> arrange(genotype, medium, replicate)
abund   <- read_abundances("total") |> keep_je2()
check("27 samples, 3 per group", nrow(samples) == 27 && all(table(samples$group) == 3))
m_all <- abundance_matrix(abund, samples)
meta  <- read_protein_metadata() |> filter(ortholog_quality == "good") |> distinct(locus_tag, gene = gene_symbol_best)
key2tag <- distinct(abund, protein_key, locus_tag)

main <- read_tsv("data/derived/limma_treat_total_mval.tsv.gz", show_col_types = FALSE) |>
  filter(contrast %in% PRIMARY)
hit_of <- \(x, fdr = MVAL_FDR) x$tested & !is.na(x$FDR) & x$FDR < fdr

# Guardrail: refitting the full data reproduces the Step 5 hits.
full <- fit_limma_treat(m_all, samples) |> filter(contrast %in% PRIMARY)
same <- full |> select(protein_key, contrast, h1 = FDR) |>
  inner_join(main |> select(protein_key, contrast, h0 = FDR), by = c("protein_key", "contrast"))
check("full refit reproduces Step 5 FDR (max abs diff < 1e-10)", max(abs(same$h1 - same$h0), na.rm = TRUE) < 1e-10)

# ---- 1. leave-one-sample-out --------------------------------------------------------
main_hits <- main |> filter(contrast %in% H_CONTRASTS, hit_of(main)) |>
  select(protein_key, contrast, diff, FDR)
prim_samples <- samples$sample[samples$medium %in% c("TSB", "SN")]
loso <- purrr::map(prim_samples, \(drop) {
  smp <- samples |> filter(sample != drop)
  stopifnot(all(table(smp$group) >= 2))
  fit_limma_treat(m_all[, smp$sample], smp) |>
    filter(contrast %in% H_CONTRASTS) |>
    semi_join(main_hits, by = c("protein_key", "contrast")) |>
    transmute(protein_key, contrast, dropped = drop, still_hit = tested & !is.na(FDR) & FDR < MVAL_FDR,
              diff_loso = diff)
}) |> bind_rows()
check("18 LOSO refits done", n_distinct(loso$dropped) == 18)
stab <- loso |> group_by(protein_key, contrast) |>
  summarise(refits_still_hit = sum(still_hit), n_refits = n(),
            lost_when_dropping = paste(dropped[!still_hit], collapse = "; "),
            diff_range_loso = paste(sprintf("%.2f", range(diff_loso, na.rm = TRUE)), collapse = " to "),
            .groups = "drop") |>
  inner_join(main_hits, by = c("protein_key", "contrast")) |>
  left_join(key2tag, by = "protein_key") |> left_join(meta, by = "locus_tag") |>
  mutate(status = if_else(refits_still_hit == n_refits, "stable (18/18)", "fragile"),
         contrast = factor(contrast, levels = H_CONTRASTS)) |>
  arrange(contrast, refits_still_hit)
write_csv(stab, file.path(out, "loso_hit_stability_mval.csv"))
print(stab |> count(contrast, status) |> tidyr::pivot_wider(names_from = status, values_from = n, values_fill = 0))
note(sprintf("LOSO: %d of %d main hits in the H1-H3 contrasts stay a hit in all 18 refits",
             sum(stab$status == "stable (18/18)"), nrow(stab)))
# which samples are behind most losses?
print(loso |> filter(!still_hit) |> count(dropped, sort = TRUE) |> head(8))

pl <- stab |> mutate(label = paste0(coalesce(gene, locus_tag), "  ", sprintf("%+.1f", diff)),
                     label = forcats::fct_reorder(label, refits_still_hit))
p <- ggplot(pl, aes(refits_still_hit, label, fill = status)) +
  geom_col() + geom_vline(xintercept = 18, linetype = 2) +
  facet_wrap(~contrast, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c(`stable (18/18)` = "#009E73", fragile = "#E69F00")) +
  labs(x = "refits (of 18, one primary sample left out each) in which the protein is still a formal hit",
       y = NULL, fill = NULL,
       title = "Step 11: does each formal hit survive leaving out one sample?",
       subtitle = "Main rule: limma treat lfc = 1, BH-FDR < 0.05. Number after the name = main log2 effect (A - B).") +
  theme(axis.text.y = element_text(size = 6.5), strip.text = element_text(size = 7), legend.position = "bottom")
ggsave(file.path(out, "fig_loso_stability_mval.png"), p, width = 10, height = 10, dpi = 150)

# ---- 2. threshold sensitivity --------------------------------------------------------
alt15 <- fit_limma_treat(m_all, samples, lfc = log2(1.5)) |> filter(contrast %in% PRIMARY)
thr <- bind_rows(
  main  |> mutate(rule = "main: treat lfc 1, FDR < 0.05", hit = hit_of(main)),
  alt15 |> mutate(rule = "alt: treat lfc log2(1.5), FDR < 0.05", hit = hit_of(alt15)),
  main  |> mutate(rule = "alt: treat lfc 1, FDR < 0.10", hit = hit_of(main, 0.10))) |>
  group_by(contrast, rule) |> summarise(hits = sum(hit), up = sum(hit & diff > 0), down = sum(hit & diff < 0), .groups = "drop") |>
  mutate(contrast = factor(contrast, levels = PRIMARY)) |> arrange(contrast, rule)
write_csv(thr, file.path(out, "threshold_sensitivity_mval.csv"))
print(thr |> select(contrast, rule, hits) |> tidyr::pivot_wider(names_from = rule, values_from = hits), width = 200)

# ---- 3. genome-position scan --------------------------------------------------------
pos <- read_protein_metadata(cds_only = TRUE) |> filter(!is.na(start)) |> distinct(locus_tag, start, strand)
WIN <- 10; NPERM <- 1000
SCAN <- c("stk_vs_anc_at_TSB", "stk_vs_anc_at_SN", "stp_vs_anc_at_TSB", "stp_vs_anc_at_SN",
          "interaction_stk_vs_anc_at_TSB_vs_SN", "interaction_stp_vs_anc_at_TSB_vs_SN")
win_mean <- function(v) as.numeric(stats::filter(c(v, v[seq_len(WIN - 1)]), rep(1 / WIN, WIN), sides = 1))[WIN:(length(v) + WIN - 1)]
all5 <- read_tsv("data/derived/limma_treat_total_mval.tsv.gz", show_col_types = FALSE) |>
  filter(contrast %in% SCAN)   # the SN baselines are not among the 9 "primary" contrasts
check("all 6 scanned contrasts present in the Step 5 table", setequal(unique(all5$contrast), SCAN))
scan <- purrr::map(SCAN, \(ct) {
  d <- all5 |> filter(contrast == ct, tested) |>       # carries locus_tag already
    inner_join(pos, by = "locus_tag") |> arrange(start)
  w <- win_mean(d$diff)                   # window i = proteins i .. i+WIN-1 (circular genome)
  null_max <- replicate(NPERM, max(abs(win_mean(sample(d$diff)))))   # random permutation
  cut <- quantile(null_max, 0.95)
  tibble(contrast = ct, i = seq_along(w), first_tag = d$locus_tag,
         last_tag = d$locus_tag[(seq_along(w) + WIN - 2) %% length(w) + 1],
         start = d$start, window_mean = w, threshold = unname(cut), flagged = abs(w) > cut,
         # family-wise p of this window: share of null genomes whose BEST window is at least as extreme
         p_familywise = purrr::map_dbl(abs(w), \(a) mean(null_max >= a)))
}) |> bind_rows()
# the most extreme window of each contrast, whether flagged or not
top_win <- scan |> group_by(contrast) |> slice_max(abs(window_mean), n = 1, with_ties = FALSE) |> ungroup() |>
  transmute(contrast, first_tag, last_tag, start, window_mean = round(window_mean, 3),
            threshold = round(threshold, 3), p_familywise)
write_csv(top_win, file.path(out, "genome_top_windows_mval.csv"))
print(top_win, width = 200)
flag <- scan |> filter(flagged)
# merge overlapping flagged windows into blocks
blocks <- flag |> arrange(contrast, i) |> group_by(contrast) |>
  mutate(block = cumsum(c(1, diff(i) > 1))) |> group_by(contrast, block) |>
  summarise(first_tag = first(first_tag), last_tag = last(last_tag), start = first(start),
            n_windows = n(), max_abs_mean = max(abs(window_mean)) * sign(window_mean[which.max(abs(window_mean))]),
            threshold = first(threshold), .groups = "drop")
write_csv(blocks, file.path(out, "genome_windows_mval.csv"))
print(blocks, width = 200)
note(sprintf("genome scan: %d flagged blocks over %d contrasts (window %d proteins, %d permutations)",
             nrow(blocks), length(SCAN), WIN, NPERM))
for (ct in c("stk_vs_anc_at_TSB", "stk_vs_anc_at_SN", "stp_vs_anc_at_TSB")) {
  s <- scan |> filter(contrast == ct, first_tag >= "B7H15_11465", first_tag <= "B7H15_11570")
  b <- s |> slice_max(abs(window_mean), n = 1, with_ties = FALSE)
  note(sprintf("%s, SigB/mazEF region (windows starting B7H15_11465-11570): best window %s-%s mean %.3f, threshold %.3f, family-wise p = %.3f, rank %d of %d windows",
               ct, b$first_tag, b$last_tag, b$window_mean, b$threshold, b$p_familywise,
               sum(abs(scan$window_mean[scan$contrast == ct]) >= abs(b$window_mean)), sum(scan$contrast == ct)))
}

tn <- tibble(tag = c("B7H15_06355", "B7H15_06350"), what = c("stk1 (NE217)", "stp1 (NE1919)")) |>
  left_join(pos, by = c(tag = "locus_tag"))
p <- ggplot(scan, aes(start / 1e6, window_mean)) +
  geom_hline(aes(yintercept = threshold), linetype = 2, colour = "grey50") +
  geom_hline(aes(yintercept = -threshold), linetype = 2, colour = "grey50") +
  geom_vline(xintercept = tn$start[1] / 1e6, colour = "#CC79A7") +
  geom_line(colour = "grey40", linewidth = 0.3) +
  geom_point(data = flag, colour = "#D55E00", size = 0.8) +
  facet_wrap(~factor(contrast, levels = SCAN), ncol = 2) +
  labs(x = "genome position (Mb, JE2 CP020619.1; window start)", y = "mean log2 effect of 10 neighbouring proteins",
       title = "Step 11: are effects clustered along the genome?",
       subtitle = "Dashed = genome-wide 5% threshold (1000 random permutations). Orange = flagged window. Pink = stp1/stk1 locus.") +
  theme(strip.text = element_text(size = 7))
ggsave(file.path(out, "fig_genome_scan_mval.png"), p, width = 11, height = 8, dpi = 150)

# ---- 3b. what is inside the flagged blocks, and which formal hits lie in one? ----------
# Locus tags share one format (B7H15_nnnnn), so text order = genome order within a block.
tag_in <- \(tag, a, b) tag >= a & tag <= b
annot <- read_protein_metadata(cds_only = TRUE) |>
  transmute(locus_tag, start, strand, gene = if_else(ortholog_quality %in% "good", gene_symbol_best, NA_character_), product)
diffs4 <- all5 |> filter(contrast %in% c("stk_vs_anc_at_TSB", "stk_vs_anc_at_SN", "stp_vs_anc_at_TSB", "stp_vs_anc_at_SN")) |>
  transmute(protein_key, locus_tag, contrast, diff = if_else(tested, round(diff, 2), NA_real_)) |>
  tidyr::pivot_wider(names_from = contrast, values_from = diff)
members <- purrr::pmap(blocks |> select(contrast, block, first_tag, last_tag), \(contrast, block, first_tag, last_tag)
  annot |> filter(tag_in(locus_tag, first_tag, last_tag)) |> mutate(block_contrast = contrast, block = block)) |>
  bind_rows() |> left_join(diffs4, by = "locus_tag") |> arrange(block_contrast, block, start)
write_csv(members, file.path(out, "genome_block_members_mval.csv"))

mutants_of <- \(ct) if (str_detect(ct, "stp_vs_stk")) c("stk", "stp") else str_extract(ct, "stk|stp")
stab <- stab |> mutate(in_own_block = purrr::map2_chr(as.character(contrast), locus_tag, \(ct, tag) {
  b <- blocks |> filter(str_extract(contrast, "stk|stp") %in% mutants_of(ct), tag_in(tag, first_tag, last_tag))
  if (nrow(b) == 0) NA_character_ else paste(unique(paste0(b$first_tag, "-", b$last_tag)), collapse = "; ")
}))
write_csv(stab, file.path(out, "loso_hit_stability_mval.csv"))   # rewritten with the block flag
note(sprintf("formal hits (H1-H3 contrasts) lying in a flagged block of their own mutant: %d of %d (%s)",
             sum(!is.na(stab$in_own_block)), nrow(stab),
             paste(unique(paste0(coalesce(stab$gene, stab$locus_tag), " [", stab$contrast, "]")[!is.na(stab$in_own_block)]), collapse = ", ")))

writeLines(c(sprintf("Step 11 checks, %s", Sys.Date()), checks), file.path(out, "robustness_checks_mval.txt"))
message("done: ", out)
