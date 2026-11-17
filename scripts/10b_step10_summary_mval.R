# Step 10 (PROJECT_WORKFLOW_mval.md): pathway-level summary, candidate annotation,
# phospho enrichment.
#
# Run from the repository root:   Rscript scripts/10b_step10_summary_mval.R   (~1 min)
# Needs: scripts/05, 07, 08, 09 and
#   Rscript scripts/10_enrichment_mval.R stk_vs_anc_at_TSB stk_vs_anc_at_SN \
#     interaction_stk_vs_anc_at_TSB_vs_SN stp_vs_anc_at_TSB stp_vs_anc_at_SN \
#     interaction_stp_vs_anc_at_TSB_vs_SN
#
# ---- Plan (coding rules, section 4) -----------------------------------------
# Aim:      (1) "Stk1 broadly involved" (presentation feedback): HOW MANY functional
#           categories shift in stk vs anc, at pathway level rather than single
#           proteins? (2) Do the same pathways shift in stp (H3 at pathway level)?
#           (3) A source-linked annotation table of our candidates. (4) Phospho
#           enrichment against the phospho background.
# Input:    mval_medi/results/step10_enrichment/<contrast>/  (ORA + rank results)
#           data/derived/limma_treat_total_mval.tsv.gz (Step 5: moderated t of the
#           TESTED proteins = the Step 5 eligible background)
#           read_protein_terms(), read_protein_metadata()
#           data/derived/limma_treat_phospho_STY_mval.tsv.gz (Step 8)
# Analysis: Rules fixed before this summary was run (17 Nov 2026, group to confirm):
#   * The rank test of 10_enrichment (Wilcoxon) treats proteins as independent.
#     Proteins of one pathway are correlated, which makes it anti-conservative
#     (Wu & Smyth 2012, Nucleic Acids Res 40:e133, doi:10.1093/nar/gks461).
#     Sensitivity: limma::cameraPR on the same moderated t, same terms (size
#     5-500), default inter-gene correlation 0.01, BH within category.
#     ROBUST pathway shift = FDR < 0.05 in BOTH tests, same direction.
#   * ORA terms with < 3 hit proteins are labelled "driven by 1-2 proteins": with
#     1-20 hits per contrast they restate the hit proteins, not a pathway.
#   * Display only: the dot plot shows the top robust terms (size <= 200, so that
#     very general terms do not crowd it) and drops a term whose members overlap
#     >= 80% (overlap coefficient, i.e. nested terms) with a term already shown.
#   * Phospho: protein-level ORA, foreground = proteins with >= 1 site hit (Step 8),
#     background = the proteins with >= 1 supported site; each protein counted once.
# Output:   mval_medi/results/step10_enrichment/summary/
#             term_counts_mval.csv           significant terms per contrast x method
#             robust_terms_mval.csv          terms significant in Wilcoxon AND cameraPR
#             candidate_annotation_mval.csv  hits of the H1-H3 contrasts + phospho, UniProt-linked
#             phospho_ora_mval.csv           phospho protein-level ORA
#             string_upload_candidates_mval.txt, string_background_tested_mval.txt
#             fig_pathway_counts_mval.png    Fig 4a: how many terms shift (stk vs stp)
#             fig_pathway_shifts_mval.png    Fig 4b: direction of top terms in all 6 contrasts
#             figure_guide_mval.md, summary_checks_mval.txt
# Check:    cameraPR and Wilcoxon agree in direction for robust terms; contrast
#           directions in words; candidate locus tags all found in the metadata.

suppressMessages(source("R/setup.R"))
out <- "mval_medi/results/step10_enrichment/summary"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
checks <- character()
check <- function(label, ok) {
  line <- sprintf("[%s] %s", if (isTRUE(ok)) "PASS" else "FAIL", label)
  message(line); checks <<- c(checks, line)
  invisible(ok)
}
note <- function(...) { line <- paste0("[NOTE] ", ...); message(line); checks <<- c(checks, line) }

CT <- tibble(
  contrast = c("stk_vs_anc_at_TSB", "stk_vs_anc_at_SN", "interaction_stk_vs_anc_at_TSB_vs_SN",
               "stp_vs_anc_at_TSB", "stp_vs_anc_at_SN", "interaction_stp_vs_anc_at_TSB_vs_SN"),
  mutant = rep(c("stk", "stp"), each = 3),
  layer = rep(c("baseline TSB", "baseline SN", "SN response (interaction)"), 2))
# "top" in the rank test = higher statistic = positive diff:
#   baseline: higher in the mutant; interaction: SN response LOWER in the mutant than in anc.
CT <- CT |> mutate(top_means = if_else(str_starts(contrast, "interaction"),
                                       paste0("SN response lower in ", mutant, " than anc"),
                                       paste0("higher in ", mutant, " than anc")),
                   bottom_means = if_else(str_starts(contrast, "interaction"),
                                          paste0("SN response higher in ", mutant, " than anc"),
                                          paste0("higher in anc than ", mutant)))
res_dir <- "mval_medi/results/step10_enrichment"
check("enrichment output exists for all 6 contrasts",
      all(file.exists(file.path(res_dir, CT$contrast, "rank_enrichment_mval.csv"))))

total <- read_tsv("data/derived/limma_treat_total_mval.tsv.gz", show_col_types = FALSE) |>
  keep_je2() |> filter(tested)
terms <- read_protein_terms() |> distinct(protein_key, category, term, term_description)

# ---- 1. cameraPR sensitivity next to the Wilcoxon rank test ------------------------
camera_one <- function(ct) {
  st <- total |> filter(contrast == ct)
  stat <- setNames(st$statistic, st$protein_key)
  tt <- terms |> filter(protein_key %in% names(stat))
  sets <- tt |> group_by(category, term) |> summarise(members = list(protein_key), n = n(), .groups = "drop") |>
    filter(n >= 5, n <= 500)
  idx <- setNames(sets$members, paste(sets$category, sets$term, sep = "||"))
  cam <- limma::cameraPR(stat, idx, inter.gene.cor = 0.01, sort = FALSE)
  tibble(key = rownames(cam), cam_direction = cam$Direction, cam_p = cam$PValue) |>
    tidyr::separate_wider_delim(key, "||", names = c("category", "term")) |>
    group_by(category) |> mutate(cam_FDR = p.adjust(cam_p, "BH")) |> ungroup() |>
    mutate(contrast = ct)
}
rank_all <- purrr::map(CT$contrast, \(ct)
  read_csv(file.path(res_dir, ct, "rank_enrichment_mval.csv"), show_col_types = FALSE) |>
    mutate(contrast = ct)) |> bind_rows()
cam_all <- purrr::map(CT$contrast, camera_one) |> bind_rows()
both <- rank_all |> inner_join(cam_all, by = c("contrast", "category", "term")) |>
  mutate(wil_dir = if_else(auc > 0.5, "Up", "Down"),
         robust = FDR < 0.05 & cam_FDR < 0.05 & wil_dir == cam_direction) |>
  left_join(CT, by = "contrast")
check("every rank-tested term also got a cameraPR result", nrow(both) == nrow(rank_all))
check("robust terms: Wilcoxon and cameraPR agree in direction (by construction, recounted)",
      all(both$wil_dir[both$robust] == both$cam_direction[both$robust]))

# ---- 2. counts: how many categories shift? ----------------------------------------
read_ora <- function(ct, d) {   # a direction with 0 hits gives a header-only file
  x <- read_csv(file.path(res_dir, ct, paste0("ora_", d, "_mval.csv")), show_col_types = FALSE)
  if (nrow(x) == 0) NULL else mutate(x, contrast = ct)
}
ora_all <- purrr::map(CT$contrast, \(ct) purrr::map(c("up", "down"), \(d) read_ora(ct, d)) |> bind_rows()) |>
  bind_rows()
counts <- bind_rows(
  ora_all |> filter(FDR < 0.05) |> count(contrast, method = if_else(hits_in_term >= 3,
    "ORA, >= 3 hit proteins", "ORA, driven by 1-2 proteins")),
  both |> filter(FDR < 0.05) |> count(contrast, method = "rank, Wilcoxon"),
  both |> filter(cam_FDR < 0.05) |> count(contrast, method = "rank, cameraPR"),
  both |> filter(robust) |> count(contrast, method = "rank, robust (both)")) |>
  tidyr::complete(contrast = CT$contrast, method, fill = list(n = 0)) |>
  left_join(CT |> select(contrast, mutant, layer), by = "contrast") |>
  arrange(mutant, layer, method)
write_csv(counts, file.path(out, "term_counts_mval.csv"))
print(counts |> select(-contrast) |> tidyr::pivot_wider(names_from = method, values_from = n), width = 200)
hits_n <- read_csv("mval_medi/results/step05_formal_test/hit_counts_mval.csv", show_col_types = FALSE) |>
  filter(contrast %in% CT$contrast) |> select(contrast, n_tested, formal_hits)
print(hits_n)

robust <- both |> filter(robust) |>
  transmute(contrast, mutant, layer, category, term, term_description, term_size, auc = round(auc, 3),
            wilcoxon_FDR = signif(FDR, 3), camera_FDR = signif(cam_FDR, 3),
            direction_in_words = if_else(auc > 0.5, top_means, bottom_means)) |>
  arrange(contrast, camera_FDR)
write_csv(robust, file.path(out, "robust_terms_mval.csv"))

# Fig 4a: robust term counts per contrast and category, stk next to stp
pc <- both |> filter(robust) |>
  mutate(category = str_remove(category, " \\(.*"),
         side = if_else(auc > 0.5, "top", "bottom")) |>
  count(mutant, layer, category, side) |>
  mutate(n_signed = if_else(side == "top", n, -n),
         layer = factor(layer, levels = c("baseline TSB", "baseline SN", "SN response (interaction)")))
p <- ggplot(pc, aes(layer, n_signed, fill = category)) +
  geom_col() + geom_hline(yintercept = 0) +
  facet_wrap(~mutant, labeller = labeller(mutant = c(stk = "stk vs anc", stp = "stp vs anc"))) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = NULL, y = "robust shifted terms\n(above 0: higher in mutant / SN response lower in mutant;\nbelow 0: higher in anc / SN response higher in mutant)",
       fill = "annotation source",
       title = "Step 10: how many functional categories shift? (rank-based, total proteome)",
       subtitle = "Robust = FDR < 0.05 in the Wilcoxon rank test AND in cameraPR (correlation-aware), same direction.\nBackground = proteins tested in that contrast. Many terms overlap (GO hierarchy): counts are not independent pathways.") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1), legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 3))
ggsave(file.path(out, "fig_pathway_counts_mval.png"), p, width = 10, height = 6.5, dpi = 300)

# Fig 4b: top robust terms of the stk contrasts, their direction in all 6 contrasts
members <- terms |> group_by(category, term) |> summarise(m = list(unique(protein_key)), .groups = "drop")
cand <- both |> filter(robust, mutant == "stk", term_size <= 200) |>
  arrange(cam_FDR) |> left_join(members, by = c("category", "term"))
keep <- list(); keep_rows <- integer()
for (i in seq_len(nrow(cand))) {
  if (sum(cand$contrast[keep_rows] == cand$contrast[i]) >= 6) next
  # overlap coefficient = shared members / size of the SMALLER term, so a term nested
  # in one already shown is dropped (Jaccard missed nested terms; changed 17 Nov)
  ov <- purrr::map_dbl(keep, \(k) length(intersect(k, cand$m[[i]])) / min(length(k), length(cand$m[[i]])))
  if (length(ov) && any(ov >= 0.8)) next
  keep <- c(keep, list(cand$m[[i]])); keep_rows <- c(keep_rows, i)
}
shown <- cand[keep_rows, ] |> distinct(category, term, term_description)
note(sprintf("Fig 4b shows %d terms (top robust stk terms, size <= 200, overlap coefficient < 0.8)", nrow(shown)))
pd <- both |> semi_join(shown, by = c("category", "term")) |>
  # Label by term ID, NOT by term_size: the size depends on which proteins were
  # tested in each contrast, so a size in the label split one term into several
  # rows (fixed 17 Nov).
  mutate(label = paste0(str_trunc(coalesce(term_description, term), 45), " [", term, "]"),
         layer = factor(layer, levels = c("baseline TSB", "baseline SN", "SN response (interaction)")),
         sig = case_when(robust ~ "robust (both tests)", FDR < 0.05 | cam_FDR < 0.05 ~ "one test only",
                         TRUE ~ "not significant"))
ord <- pd |> filter(mutant == "stk") |> group_by(label) |> summarise(o = mean(auc)) |> arrange(o)
pd$label <- factor(pd$label, levels = ord$label)
p <- ggplot(pd, aes(layer, label)) +
  geom_point(aes(fill = auc - 0.5, size = -log10(cam_FDR), shape = sig), colour = "grey30") +
  facet_wrap(~mutant, labeller = labeller(mutant = c(stk = "stk vs anc", stp = "stp vs anc"))) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0,
                       name = "AUC - 0.5\n(> 0: higher in mutant /\nSN response lower in mutant)") +
  scale_shape_manual(values = c(`robust (both tests)` = 21, `one test only` = 22, `not significant` = 24), name = NULL) +
  scale_size_continuous(range = c(1.5, 6), name = "-log10 cameraPR FDR") +
  labs(x = NULL, y = NULL,
       title = "Step 10: do the stk pathway shifts also occur in stp?",
       subtitle = "Rows = top robust terms in the stk contrasts (term ID in brackets; nested/overlapping terms removed).
The same terms are shown for stp.") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1), axis.text.y = element_text(size = 7))
ggsave(file.path(out, "fig_pathway_shifts_mval.png"), p, width = 11, height = max(5, 1.5 + 0.28 * nrow(shown)), dpi = 300)

# Is the stk pattern shared by stp at pathway level? Correlation of AUCs over all terms.
auc_w <- both |> select(layer, mutant, category, term, auc) |>
  tidyr::pivot_wider(names_from = mutant, values_from = auc)
auc_r <- auc_w |> group_by(layer) |> summarise(n_terms = n(), r_auc_stk_stp = cor(stk, stp, use = "complete.obs"))
print(auc_r)
note(paste(sprintf("AUC correlation stk vs stp, %s: r = %.2f over %d terms", auc_r$layer, auc_r$r_auc_stk_stp, auc_r$n_terms),
           collapse = "; "))
write_csv(auc_r, file.path(out, "auc_correlation_stk_stp_mval.csv"))

# ---- 3. candidate annotation table (source-linked) -----------------------------------
meta <- read_protein_metadata()
H_CONTRASTS <- c("stk_vs_anc_at_TSB", "stp_vs_anc_at_TSB", "stp_vs_stk_at_TSB",
                 "interaction_stk_vs_anc_at_TSB_vs_SN", "interaction_stp_vs_anc_at_TSB_vs_SN",
                 "interaction_stp_vs_stk_at_TSB_vs_SN")
prot_c <- read_tsv("data/derived/limma_treat_total_mval.tsv.gz", show_col_types = FALSE) |>
  keep_je2() |> filter(contrast %in% H_CONTRASTS, tested, FDR < MVAL_FDR) |>
  transmute(layer = "total protein", contrast, locus_tag, feature = "protein",
            diff = round(diff, 2), se = round(std.error, 2), FDR = signif(FDR, 2))
phos_c <- read_csv("mval_medi/results/step09_site_vs_protein/site_hit_categories_mval.csv", show_col_types = FALSE) |>
  filter(!str_starts(contrast, "TSB_vs_SN")) |>
  transmute(layer = "phosphosite", contrast, locus_tag, feature = paste0(residue, position, " (", category, ")"),
            diff = round(diff_site, 2), se = round(se_site, 2), FDR = signif(FDR_site, 2))
cands <- bind_rows(prot_c, phos_c) |>
  left_join(meta |> select(locus_tag, any_of(c("gene_symbol_best", "up_protein_name", "product", "uniprot_acc",
                                               "ortholog_quality", "up_ec", "kegg", "uniprot_keywords"))),
            by = "locus_tag", relationship = "many-to-one") |>
  mutate(direction_in_words = purrr::map2_chr(contrast, diff, \(ct, d) {
           if (str_starts(ct, "interaction_")) {
             g <- str_match(ct, "^interaction_([a-z]+)_vs_([a-z]+)_at_")
             paste0("SN response lower in ", if (d > 0) g[2] else g[3], " than in ", if (d > 0) g[3] else g[2])
           } else {
             g <- str_match(ct, "^([a-z]+)_vs_([a-z]+)_at_([A-Za-z]+)")
             paste0("higher in ", if (d > 0) g[2] else g[3], " than ", if (d > 0) g[3] else g[2], " (", g[4], ")")
           }}),
         gene_symbol_best = if_else(ortholog_quality == "good", gene_symbol_best, NA_character_),
         uniprot_link = if_else(!is.na(uniprot_acc), paste0("https://www.uniprot.org/uniprotkb/", uniprot_acc), NA_character_),
         across(any_of(c("kegg", "uniprot_keywords")), \(x) str_trunc(x, 120))) |>
  arrange(layer, contrast, desc(abs(diff)))
check("all candidate locus tags found in the metadata", !anyNA(cands$product))
write_csv(cands, file.path(out, "candidate_annotation_mval.csv"))
note(sprintf("candidates: %d rows, %d proteins (%d with a good UniProt ortholog)",
             nrow(cands), n_distinct(cands$locus_tag),
             n_distinct(cands$locus_tag[cands$ortholog_quality %in% "good"])))
write_string_upload(cands$locus_tag, file.path(out, "string_upload_candidates_mval.txt"))
write_string_upload(unique(total$locus_tag), file.path(out, "string_background_tested_mval.txt"))

# ---- 4. phospho enrichment, protein level, phospho background ------------------------
ph <- read_tsv("data/derived/limma_treat_phospho_STY_mval.tsv.gz", show_col_types = FALSE) |> filter(supported)
ph_bg  <- unique(ph$protein_key)
ph_fg  <- unique(ph$protein_key[ph$hit])
ph_ora <- run_ora(ph_fg, ph_bg, read_protein_terms())
write_csv(select(ph_ora, -hit_keys), file.path(out, "phospho_ora_mval.csv"))
note(sprintf("phospho ORA: %d hit proteins of %d background proteins; terms with FDR < 0.05: %d (with >= 3 hits: %d)",
             length(ph_fg), length(ph_bg), sum(ph_ora$FDR < 0.05), sum(ph_ora$FDR < 0.05 & ph_ora$hits_in_term >= 3)))
print(ph_ora |> filter(FDR < 0.05) |> select(category, term_description, hits_in_term, term_size, FDR) |> head(10))

writeLines(c(sprintf("Step 10 summary checks, %s", Sys.Date()), checks), file.path(out, "summary_checks_mval.txt"))
message("done: ", out)
