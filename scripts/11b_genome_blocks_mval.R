# Step 11b (PROJECT_WORKFLOW_mval.md, after Step 11): what changes inside the flagged
# genome blocks, shown in the same layout as the Step 10 pathway-shift figure (Fig 4b).
#
# Run from the repository root:   Rscript scripts/11b_genome_blocks_mval.R   (~20 s)
# Needs scripts/05, 10, 10b and 11 first.
#
# ---- Plan (coding rules, section 4) -----------------------------------------
# Aim:      The Step 11 genome scan says WHERE neighbouring proteins change together, but
#           not WHAT they are. Name the content of each flagged block, show its effect in
#           the same 6 contrasts as Fig 4b (stk / stp x baseline TSB, baseline SN, SN
#           response), and count how many members of each Fig 4b pathway term lie inside
#           a flagged block. Question: how much of a "pathway shift" is one genomic region?
# Input:    mval_medi/results/step11_robustness/genome_windows_mval.csv (flagged blocks,
#             per-contrast thresholds)
#           data/derived/limma_treat_total_mval.tsv.gz (Step 5 fit: diff, tested, statistic)
#           read_protein_metadata(cds_only = TRUE) (start, gene symbol, product)
#           read_protein_terms(); step10 rank_enrichment_mval.csv + summary/robust_terms_mval.csv
# Analysis: No new test. (1) Merge flagged blocks that overlap across contrasts into
#           regions R1..Rn (genome order). (2) Per region and contrast: mean log2 diff of the
#           TESTED proteins in the region (not the window mean; different number), and
#           whether a flagged window of that contrast lies in it. (3) Rebuild the 16 Fig 4b
#           terms with the 10b rule and count their tested members inside each mutant's
#           flagged regions. Region labels are written by hand FROM the annotation
#           (gene symbols / products in the member table), and the script checks that every
#           gene named in a label is really in that region.
# Output:   mval_medi/results/step11_robustness/
#             genome_regions_mval.csv          one row per region: range, content, where flagged
#             genome_region_effects_mval.csv   region x contrast: mean diff, n tested, flagged, words
#             block_term_overlap_mval.csv      Fig 4b term x mutant: members inside flagged regions
#             fig_genome_blocks_mval.png       Fig 4b-style dot plot of the regions
#             fig_genome_scan_labelled_mval.png the Step 11 scan with region labels
#             fig_pathway_terms_in_blocks_mval.png  share of each Fig 4b term inside the blocks
#             genome_blocks_checks_mval.txt
# Check:    recomputed windows reproduce the Step 11 flags exactly; Fig 4b rebuilt with 16
#           terms; every gene named in a label lies in its region.

suppressMessages(source("R/setup.R"))
out <- "mval_medi/results/step11_robustness"
checks <- character()
check <- function(label, ok) {
  line <- sprintf("[%s] %s", if (isTRUE(ok)) "PASS" else "FAIL", label)
  message(line); checks <<- c(checks, line)
  invisible(ok)
}
note <- function(...) { line <- paste0("[NOTE] ", ...); message(line); checks <<- c(checks, line) }

# The 6 contrasts of the Step 11 scan = the 6 columns of Fig 4b.
CT <- tibble(
  contrast = c("stk_vs_anc_at_TSB", "stk_vs_anc_at_SN", "interaction_stk_vs_anc_at_TSB_vs_SN",
               "stp_vs_anc_at_TSB", "stp_vs_anc_at_SN", "interaction_stp_vs_anc_at_TSB_vs_SN"),
  mutant = rep(c("stk", "stp"), each = 3),
  layer = factor(rep(c("baseline TSB", "baseline SN", "SN response (interaction)"), 2),
                 levels = c("baseline TSB", "baseline SN", "SN response (interaction)")))
tag_in <- \(tag, a, b) tag >= a & tag <= b    # B7H15_nnnnn: text order = genome order

# ---- 1. flagged blocks -> regions ---------------------------------------------------
blocks <- read_csv(file.path(out, "genome_windows_mval.csv"), show_col_types = FALSE)
check("Step 11 blocks cover the 6 contrasts", setequal(unique(blocks$contrast), CT$contrast))
# merge tag ranges that overlap in ANY contrast (same stretch of DNA flagged twice = one region)
b <- blocks |> arrange(first_tag)
reg <- list(); cur <- b[1, c("first_tag", "last_tag")]
for (i in seq_len(nrow(b))[-1]) {
  if (b$first_tag[i] <= cur$last_tag) cur$last_tag <- max(cur$last_tag, b$last_tag[i])
  else { reg <- c(reg, list(cur)); cur <- b[i, c("first_tag", "last_tag")] }
}
regions <- bind_rows(c(reg, list(cur))) |> mutate(region = paste0("R", row_number()))

# Hand-written labels, read off the member table (gene symbols / products). `genes` must
# all lie in the region; the check below fails if Step 11 ever moves a block.
LABELS <- tribble(
  ~first_tag,    ~label,                                              ~genes,
  "B7H15_01055", "ptsG, murQ, azoR (sugar uptake)",                   "ptsG|murQ|azoR",
  "B7H15_01445", "rbs operon, lytM, esxA/esaA (type VII secretion)",  "rbsK|rbsD|lytM|esxA|esaA",
  "B7H15_04460", "cold-shock protein, gcvH",                          "gcvH",
  "B7H15_07835", "prophage at 1.56 Mb (with leukocidin genes)",       "lukDv|hlgC",
  "B7H15_09880", "splA-F proteases, lantibiotic genes, lukD",         "splB|splF",
  "B7H15_10910", "prophage at 2.09 Mb",                               "",
  "B7H15_11470", "sigB, mazEF, alr1, kdp, murF, ddl",                 "rsbW|mazE|mazF|alr1|kdpA|murF|ddl",
  "B7H15_11640", "sceD, fabZ, murA, ATP synthase (atpA/C/D/G/H)",     "sceD|fabZ|murA|atpA|atpC|atpD|atpG|atpH")
check("one hand-written label per region (regions unchanged since labelling)",
      nrow(regions) == nrow(LABELS) && all(regions$first_tag == LABELS$first_tag))
regions <- regions |> left_join(LABELS, by = "first_tag")

annot <- read_protein_metadata(cds_only = TRUE) |>
  transmute(locus_tag, start, gene = if_else(ortholog_quality %in% "good", gene_symbol_best, NA_character_), product)
members <- purrr::pmap(regions |> select(region, first_tag, last_tag), \(region, first_tag, last_tag)
  annot |> filter(tag_in(locus_tag, first_tag, last_tag)) |> mutate(region = region)) |> bind_rows()
label_ok <- purrr::map_lgl(seq_len(nrow(regions)), \(i) {
  g <- str_split(regions$genes[i], "\\|")[[1]]; g <- g[g != ""]
  all(g %in% members$gene[members$region == regions$region[i]])
})
check("every gene named in a region label lies in that region", all(label_ok))
phage_words <- "phage|terminase|capsid|portal|holin"
regions <- regions |> left_join(
  members |> group_by(region) |> summarise(
    start_mb = round(min(start) / 1e6, 2), n_cds = n(),
    n_phage_products = sum(str_detect(product, phage_words)),
    gene_symbols = paste(na.omit(unique(gene)), collapse = ", "), .groups = "drop"), by = "region")
check("both prophage regions carry phage structural genes",
      all(regions$n_phage_products[str_detect(regions$label, "prophage")] >= 3))

# ---- 2. region effect per contrast -------------------------------------------------
fit <- read_tsv("data/derived/limma_treat_total_mval.tsv.gz", show_col_types = FALSE) |>
  filter(contrast %in% CT$contrast)
eff <- tidyr::expand_grid(regions |> select(region, first_tag, last_tag), CT) |>
  mutate(
    d = purrr::pmap(list(first_tag, last_tag, contrast), \(a, z, ct)
      fit |> filter(contrast == ct, tested, tag_in(locus_tag, a, z)) |> pull(diff)),
    n_tested = lengths(d),
    mean_diff = purrr::map_dbl(d, \(x) if (length(x)) mean(x) else NA_real_),
    flagged = purrr::pmap_lgl(list(first_tag, last_tag, contrast), \(a, z, ct)
      any(blocks$contrast == ct & blocks$first_tag <= z & blocks$last_tag >= a)),
    in_words = case_when(
      is.na(mean_diff) ~ "no tested protein",
      str_starts(contrast, "interaction") & mean_diff > 0 ~ paste0("SN response lower in ", mutant, " than anc"),
      str_starts(contrast, "interaction") ~ paste0("SN response higher in ", mutant, " than anc"),
      mean_diff > 0 ~ paste0("higher in ", mutant, " than anc"),
      TRUE ~ paste0("lower in ", mutant, " than anc"))) |>
  select(-d)
check("every Step 11 block lies inside exactly one region",
      all(purrr::map2_int(blocks$first_tag, blocks$last_tag, \(a, z)
        sum(regions$first_tag <= a & regions$last_tag >= z)) == 1))
regions <- regions |> left_join(
  eff |> filter(flagged) |> group_by(region) |> summarise(flagged_in = paste(contrast, collapse = "; ")), by = "region")
write_csv(regions |> select(region, first_tag, last_tag, start_mb, label, n_cds, n_phage_products,
                            gene_symbols, flagged_in), file.path(out, "genome_regions_mval.csv"))
write_csv(eff |> select(region, contrast, mutant, layer, n_tested, mean_diff, flagged, in_words) |>
            mutate(mean_diff = round(mean_diff, 2)), file.path(out, "genome_region_effects_mval.csv"))
print(regions |> select(region, first_tag, last_tag, start_mb, n_cds, label, flagged_in), width = 200)

# ---- 3. Fig 4b-style dot plot of the regions -----------------------------------------
pd <- eff |> left_join(regions |> select(region, label, start_mb), by = "region") |>
  mutate(row = factor(sprintf("%s  %s  [%s-%s, %.2f Mb]", region, label,
                              str_remove(first_tag, "B7H15_"), str_remove(last_tag, "B7H15_"), start_mb)),
         shape = if_else(flagged, "flagged block in this contrast", "not flagged"))
pd$row <- factor(pd$row, levels = rev(sort(unique(as.character(pd$row)))))   # R1 on top
lim <- max(abs(pd$mean_diff), na.rm = TRUE)
p <- ggplot(pd, aes(layer, row)) +
  geom_point(aes(fill = mean_diff, shape = shape), size = 6, colour = "grey30") +
  geom_text(aes(label = sprintf("%+.1f", mean_diff)), nudge_x = 0.32, size = 3, colour = "grey20") +
  facet_wrap(~mutant, labeller = labeller(mutant = c(stk = "stk vs anc", stp = "stp vs anc"))) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, limits = c(-lim, lim),
                       name = "mean log2 diff of the\ntested proteins in the region\n(> 0: higher in mutant /\nSN response lower in mutant)") +
  scale_shape_manual(values = c(`flagged block in this contrast` = 21, `not flagged` = 24), name = NULL) +
  labs(x = NULL, y = NULL,
       title = "Step 11b: what changes inside the flagged genome blocks?",
       subtitle = "Rows = regions flagged by the Step 11 genome scan in at least one contrast (content from the annotation).
Same columns as the pathway figure (Fig 4b). Number = mean log2 diff of the tested proteins in the region.",
       caption = "Circle = a window inside the region passes the genome-wide 5% threshold of that contrast. No new test.") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1), axis.text.y = element_text(size = 8),
        panel.grid.major.y = element_line(colour = "grey92"))
ggsave(file.path(out, "fig_genome_blocks_mval.png"), p, width = 12, height = 5.2, dpi = 300)

# ---- 4. the Step 11 scan again, now with region labels --------------------------------
# Window means are deterministic, so they are recomputed; the thresholds come from Step 11
# (they need the permutations). Check: the recomputed flags equal Step 11's.
WIN <- 10
pos <- read_protein_metadata(cds_only = TRUE) |> filter(!is.na(start)) |> distinct(locus_tag, start)
win_mean <- function(v) as.numeric(stats::filter(c(v, v[seq_len(WIN - 1)]), rep(1 / WIN, WIN), sides = 1))[WIN:(length(v) + WIN - 1)]
thr <- blocks |> distinct(contrast, threshold)
scan <- purrr::map(CT$contrast, \(ct) {
  d <- fit |> filter(contrast == ct, tested) |> inner_join(pos, by = "locus_tag") |> arrange(start)
  tibble(contrast = ct, start = d$start, first_tag = d$locus_tag, window_mean = win_mean(d$diff))
}) |> bind_rows() |> left_join(thr, by = "contrast") |> mutate(flagged = abs(window_mean) > threshold)
check("recomputed windows reproduce the Step 11 flagged-window count",
      sum(scan$flagged) == sum(blocks$n_windows))
lab <- scan |> filter(flagged) |>
  mutate(region = purrr::map_chr(first_tag, \(t) regions$region[tag_in(t, regions$first_tag, regions$last_tag)][1])) |>
  group_by(contrast, region) |> slice_max(abs(window_mean), n = 1, with_ties = FALSE) |> ungroup()
check("every flagged window falls in a region", !anyNA(lab$region))
tn <- annot |> filter(locus_tag == "B7H15_06355")   # stk1; stp1 is 4 bp upstream
facet_lab <- setNames(paste0(CT$mutant, " vs anc: ", CT$layer), CT$contrast)
# one dotted line per region in EVERY panel, so a block flagged in one contrast can be
# compared with the other five. Neighbouring regions (R6/R7/R8) get staggered labels.
reg_lines <- members |> group_by(region) |> summarise(mid = (min(start) + max(start)) / 2e6) |>
  arrange(mid) |> mutate(y_lab = c(2.55, 2.25)[(row_number() - 1) %% 2 + 1])
p2 <- ggplot(scan, aes(start / 1e6, window_mean)) +
  geom_vline(data = reg_lines, aes(xintercept = mid), linetype = 3, colour = "#56B4E9", linewidth = 0.4) +
  geom_text(data = reg_lines, aes(x = mid, y = y_lab, label = region), size = 2.3, colour = "#1F6E9C") +
  geom_hline(aes(yintercept = threshold), linetype = 2, colour = "grey50") +
  geom_hline(aes(yintercept = -threshold), linetype = 2, colour = "grey50") +
  geom_vline(xintercept = tn$start / 1e6, colour = "#CC79A7") +
  geom_line(colour = "grey40", linewidth = 0.3) +
  geom_point(data = filter(scan, flagged), colour = "#D55E00", size = 0.8) +
  geom_label(data = lab, aes(label = region, vjust = if_else(window_mean > 0, -0.3, 1.3)),
             size = 3, label.size = 0.2, label.padding = unit(0.12, "lines"), fill = "white") +
  facet_wrap(~factor(contrast, levels = CT$contrast), ncol = 3, labeller = labeller(.default = \(x) facet_lab[x])) +
  scale_y_continuous(expand = expansion(mult = 0.18)) +
  labs(x = "genome position (Mb, JE2 CP020619.1; window start)", y = "mean log2 effect of 10 neighbouring proteins",
       title = "Step 11: are effects clustered along the genome? (regions named in Step 11b)",
       subtitle = "Dashed = genome-wide 5% threshold (1000 random permutations). Orange = flagged window. Pink = stp1/stk1 locus. Blue dotted = regions R1-R8, drawn in every panel.",
       caption = paste(sprintf("%s %s", regions$region, regions$label), collapse = " | ") |> str_wrap(170)) +
  theme(strip.text = element_text(size = 8), plot.caption = element_text(hjust = 0, size = 8))
ggsave(file.path(out, "fig_genome_scan_labelled_mval.png"), p2, width = 13, height = 7, dpi = 200)

# ---- 5. how much of each Fig 4b pathway term lies in a flagged region? -----------------
# Rebuild the 16 terms exactly as scripts/10b does (robust stk terms, size <= 200, by
# cameraPR FDR, max 6 per contrast, overlap coefficient < 0.8).
terms <- read_protein_terms() |> distinct(protein_key, category, term, term_description)
term_members <- terms |> group_by(category, term) |> summarise(m = list(unique(protein_key)), .groups = "drop")
cand <- read_csv("mval_medi/results/step10_enrichment/summary/robust_terms_mval.csv", show_col_types = FALSE) |>
  filter(mutant == "stk", term_size <= 200) |> arrange(camera_FDR) |>
  left_join(term_members, by = c("category", "term"))
keep <- list(); keep_rows <- integer()
for (i in seq_len(nrow(cand))) {
  if (sum(cand$contrast[keep_rows] == cand$contrast[i]) >= 6) next
  ov <- purrr::map_dbl(keep, \(k) length(intersect(k, cand$m[[i]])) / min(length(k), length(cand$m[[i]])))
  if (length(ov) && any(ov >= 0.8)) next
  keep <- c(keep, list(cand$m[[i]])); keep_rows <- c(keep_rows, i)
}
shown <- cand[keep_rows, ] |> distinct(category, term, term_description)
check("Fig 4b rebuilt with the same 16 terms", nrow(shown) == 16)

tested_keys <- fit |> filter(tested) |> distinct(protein_key, locus_tag)
own_regions <- \(mut) regions |> filter(str_detect(coalesce(flagged_in, ""), paste0(mut, "_vs_anc")))
overlap <- tidyr::expand_grid(shown, mutant = c("stk", "stp")) |>
  mutate(res = purrr::pmap(list(category, term, mutant), \(cat, tm, mut) {
    keys <- terms |> filter(category == cat, term == tm) |> pull(protein_key)
    tg <- tested_keys |> filter(protein_key %in% keys) |> distinct(locus_tag) |> pull(locus_tag)
    rg <- own_regions(mut)
    inside <- tg[purrr::map_lgl(tg, \(t) any(tag_in(t, rg$first_tag, rg$last_tag)))]
    which_reg <- unique(unlist(purrr::map(inside, \(t) rg$region[tag_in(t, rg$first_tag, rg$last_tag)])))
    nm <- annot |> filter(locus_tag %in% inside) |> mutate(n = coalesce(gene, locus_tag)) |> pull(n)
    tibble(n_tested_members = length(tg), n_in_flagged_region = length(inside),
           regions = paste(sort(which_reg), collapse = ", "), members_in_region = paste(nm, collapse = ", "))
  })) |> tidyr::unnest(res) |>
  mutate(share_in_region = round(n_in_flagged_region / n_tested_members, 2))
write_csv(overlap, file.path(out, "block_term_overlap_mval.csv"))
print(overlap |> select(term, mutant, n_tested_members, n_in_flagged_region, regions, members_in_region), n = 40, width = 200)

# same row order as Fig 4b: mean stk AUC over the 3 stk contrasts
auc_stk <- purrr::map(CT$contrast[CT$mutant == "stk"], \(ct)
  read_csv(file.path("mval_medi/results/step10_enrichment", ct, "rank_enrichment_mval.csv"), show_col_types = FALSE)) |>
  bind_rows() |> semi_join(shown, by = c("category", "term")) |> group_by(category, term) |>
  summarise(o = mean(auc), .groups = "drop")
po <- overlap |> left_join(auc_stk, by = c("category", "term")) |>
  mutate(label = paste0(str_trunc(coalesce(term_description, term), 45), " [", term, "]"),
         label = reorder(label, o),
         txt = if_else(n_in_flagged_region > 0,
                       sprintf("%d/%d  %s", n_in_flagged_region, n_tested_members, str_trunc(members_in_region, 28)),
                       sprintf("0/%d", n_tested_members)))
p3 <- ggplot(po, aes(share_in_region, label)) +
  geom_col(aes(fill = mutant), width = 0.65) +
  geom_text(aes(label = txt), hjust = -0.05, size = 2.7, colour = "grey20") +
  facet_wrap(~mutant, labeller = labeller(mutant = c(stk = "inside a flagged stk region", stp = "inside a flagged stp region"))) +
  scale_fill_manual(values = c(stk = "#E69F00", stp = "#009E73"), guide = "none") +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1.9), breaks = c(0, 0.5, 1)) +
  labs(x = "share of the term's tested proteins", y = NULL,
       title = "Step 11b: how much of each pathway-figure term lies in a flagged genome block?",
       subtitle = "Rows = the 16 terms of the pathway figure (Fig 4b), same order. Text = members inside / tested members, then their names.") +
  theme(axis.text.y = element_text(size = 7))
ggsave(file.path(out, "fig_pathway_terms_in_blocks_mval.png"), p3, width = 12, height = 5.5, dpi = 300)

writeLines(c(sprintf("Step 11b checks, %s", Sys.Date()), checks), file.path(out, "genome_blocks_checks_mval.txt"))
message("done: ", out)
