# Step 10 tool (mval_medi/project_outlines/PROJECT_WORKFLOW_mval.md):
# functional enrichment of one or more total-proteome contrasts, done locally in R.
#
# Why locally and not only on the STRING website?
#   - we control the BACKGROUND: only proteins we actually measured (and that
#     were tested in that contrast), not the whole genome. Otherwise categories
#     that are easy to detect by MS look "enriched".
#   - it is reproducible: same input -> same table, re-run any time.
#   - STRING (organism STRG0B09NBK, "JE2_Bio253_HS26") stays useful for network
#     pictures: this script writes ready-to-paste locus-tag lists.
#
# Two methods (lecture 06, slide 25):
#   ORA  = formal hits from step 5 (limma treat, H0 |log2FC| <= 1, BH-FDR < 0.05),
#          split up/down, Fisher test. No post-hoc |log2FC| filter: the threshold
#          is part of the test. The formal test uses observed values only, so
#          presence/absence (lod_imputed) proteins are not hits here; they are
#          counted separately in step 5 (hit_counts_mval.csv).
#   rank = GSEA-like: all TESTED proteins ranked by the ordinary moderated t (eBayes),
#          Wilcoxon test per term. Only OBSERVED estimates are ranked, because an
#          imputed fold change is not a precise number (CLAUDE.md rule 5).
#
# Run from the repository root:
#   Rscript scripts/10_enrichment_mval.R                         (default: stk_vs_anc)
#   Rscript scripts/10_enrichment_mval.R stk_vs_anc_at_TSB stp_vs_anc_at_TSB
#
# Output per contrast: mval_medi/results/step10_enrichment/<contrast>/
#   hits_up_mval.csv / hits_down_mval.csv        selected proteins (+ imputed flag)
#   ora_up_mval.csv / ora_down_mval.csv          ORA results, all categories
#   rank_enrichment_mval.csv                     rank-based results
#   string_upload_up_mval.txt / _down_mval.txt   locus tags for STRING
#   fig_top_terms_mval.png                       top terms of both methods

suppressMessages(source("R/setup.R"))

contrasts <- commandArgs(trailingOnly = TRUE)
if (length(contrasts) == 0) contrasts <- "stk_vs_anc"
FDR_CUT <- MVAL_FDR  # 0.05, group decision 2026-11-12 (rules section 5)
LFC_CUT <- 0         # only splits up/down: |log2FC| > 1 is already tested by treat()

# Needs step 5 first: Rscript scripts/05_formal_test_mval.R
formal_file <- "data/derived/limma_treat_total_mval.tsv.gz"
if (!file.exists(formal_file)) stop("Run scripts/05_formal_test_mval.R first.")
# background = proteins actually TESTED in the contrast (enough observed values)
dea   <- read_tsv(formal_file, show_col_types = FALSE) |> filter(tested) |> keep_je2()
terms <- read_protein_terms()
stopifnot(all(contrasts %in% dea$contrast))

# Terms reported for stk_vs_anc on the lecture slides (06_web_resources, slides
# 36-38). We look them up in our results as a plausibility check of the method.
LECTURE_TERMS <- c("GO:0051715", "GO:0031640", "GO:0090729", "KW-0843", "GO:0005576", "GO:0006825")

for (ct in contrasts) {
  message("\n==== ", ct, " ====")
  print(explain_contrast(ct))
  out <- file.path("mval_medi/results/step10_enrichment", ct)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  tested <- filter(dea, contrast == ct)
  background <- unique(tested$protein_key)

  # ---- ORA, separately for up and down -----------------------------------------
  ora <- list()
  for (dir in c("up", "down")) {
    hits <- hit_list(dea, ct, fdr = FDR_CUT, lfc = LFC_CUT, direction = dir)
    message(sprintf("%-4s hits: %d (of which lod_imputed: %d) out of %d tested proteins",
                    dir, nrow(hits), sum(hits$imputed), length(background)))
    write_csv(hits, file.path(out, paste0("hits_", dir, "_mval.csv")))
    write_string_upload(hits$locus_tag, file.path(out, paste0("string_upload_", dir, "_mval.txt")))

    res <- run_ora(hits$protein_key, background, terms) |>
      mutate(hits_imputed = purrr::map_int(strsplit(hit_keys, ";"),
                                           \(k) sum(k %in% hits$protein_key[hits$imputed])),
             direction = dir)
    write_csv(select(res, -hit_keys), file.path(out, paste0("ora_", dir, "_mval.csv")))
    ora[[dir]] <- res
  }

  # ---- rank-based (GSEA-like), observed estimates only ---------------------------
  obs <- filter(tested, estimate_type == "observed")
  stats <- setNames(obs$statistic, obs$protein_key)
  message("rank test on ", length(stats), " observed proteins (",
          nrow(tested) - nrow(obs), " imputed left out)")
  rank_res <- run_rank_enrichment(stats, terms)
  write_csv(rank_res, file.path(out, "rank_enrichment_mval.csv"))

  # ---- short console summary -------------------------------------------------------
  first_group  <- str_split_i(str_remove(ct, "_at_.*"), "_vs_", 1)
  second_group <- str_split_i(str_remove(ct, "_at_.*"), "_vs_", 2)
  show <- \(x) x |> filter(FDR < 0.05) |> slice_head(n = 8) |>
    transmute(category = str_trunc(category, 28), term, description = str_trunc(term_description, 45),
              n = paste0(hits_in_term, "/", term_size), imputed = hits_imputed, FDR = signif(FDR, 2))
  message("\nORA up (higher in ", first_group, "), FDR < 0.05:");   print(show(ora$up), n = 8)
  message("\nORA down (higher in ", second_group, "), FDR < 0.05:"); print(show(ora$down), n = 8)
  message("\nRank-based, FDR < 0.05 (top = higher in ", first_group, "):")
  print(rank_res |> filter(FDR < 0.05) |> slice_head(n = 10) |>
          transmute(term, description = str_trunc(term_description, 45), size = term_size,
                    auc = round(auc, 2), direction, FDR = signif(FDR, 2)), n = 10)

  message("\nPlausibility check against the lecture slides (stk_vs_anc):")
  print(bind_rows(
    bind_rows(ora) |> filter(term %in% LECTURE_TERMS) |>
      transmute(method = paste("ORA", direction), term, term_description, n = paste0(hits_in_term, "/", term_size), FDR = signif(FDR, 2)),
    rank_res |> filter(term %in% LECTURE_TERMS) |>
      transmute(method = paste("rank", direction), term, term_description, n = as.character(term_size), FDR = signif(FDR, 2))) |>
      arrange(term), n = 30)

  # ---- figure: top terms per method --------------------------------------------------
  top <- bind_rows(
    bind_rows(ora) |> group_by(direction) |> slice_min(FDR, n = 8, with_ties = FALSE) |> ungroup() |>
      transmute(panel = paste0("ORA ", direction), label = str_trunc(term_description, 50), FDR),
    rank_res |> slice_min(FDR, n = 10, with_ties = FALSE) |>
      transmute(panel = paste0("rank (", if_else(str_starts(direction, "top"), "top", "bottom"), ")"),
                label = str_trunc(term_description, 50), FDR))
  p <- ggplot(top, aes(-log10(FDR), reorder(label, -log10(FDR)))) +
    geom_col() + geom_vline(xintercept = -log10(0.05), linetype = 2) +
    facet_wrap(~panel, scales = "free_y", ncol = 1) +
    labs(y = NULL, x = "-log10(FDR)  (dashed line: FDR = 0.05)",
         title = paste("Enrichment:", ct),
         subtitle = paste0("up = higher in ", first_group, ", down = higher in ", second_group,
                           "; background = ", length(background), " measured JE2 proteins"))
  ggsave(file.path(out, "fig_top_terms_mval.png"), p, width = 8, height = 9, dpi = 150)
  message("written to ", out)
}
