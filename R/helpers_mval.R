# helpers_mval.R -- group mval helpers (loaded automatically by R/setup.R).
#
# 1. keep_je2()             one agreed filter for "real JE2 proteins"
# 2. hit_list()             up/down hits of one contrast, imputed ones labelled
# 3. run_ora()              over-representation analysis (Fisher test) against
#                           the MEASURED background, not the whole genome
# 4. run_rank_enrichment()  GSEA-like test on the full ranked list (Wilcoxon)
# 5. write_string_upload()  locus-tag list to paste into STRING (STRG0B09NBK)
#
# Why our own filter? The default `drop_contaminants = TRUE` misses 68
# contaminants (and in abundance tables 24 reversed decoys `rev_...`), see
# notebook 2026-11-11. The lecture's `grep 'ARG'` trick drops 15 real JE2
# proteins whose key has no "ARG". A locus tag is the cleanest test.

#' Keep only real JE2 proteins: has a locus tag and is not a reversed decoy.
keep_je2 <- function(x) {
  filter(x, !is.na(locus_tag), !str_detect(protein_key, "^rev_"))
}

#' Hits of one contrast. direction "up" = diff > lfc (higher in the FIRST group
#' of the contrast name), "down" = diff < -lfc (higher in the SECOND group).
#' lod_imputed rows are kept (presence/absence evidence) but labelled.
hit_list <- function(dea, contrast, fdr = 0.1, lfc = 1, direction = c("up", "down")) {
  direction <- match.arg(direction)
  sel <- dea |> keep_je2() |> filter(.data$contrast == !!contrast, FDR < fdr)
  sel <- if (direction == "up") filter(sel, diff > lfc) else filter(sel, diff < -lfc)
  sel |>
    mutate(imputed = estimate_type == "lod_imputed") |>
    select(protein_key, locus_tag, product, contrast, diff, FDR, statistic, estimate_type, imputed) |>
    arrange(FDR)
}

#' Over-representation analysis with a one-sided Fisher exact test.
#'
#' hits:       protein_keys of the selected proteins
#' background: protein_keys of ALL proteins tested in that contrast
#' terms:      read_protein_terms() (protein_key, category, term, term_description)
#' Per category the universe is the background proteins that have at least
#' one term in that category; FDR (Benjamini-Hochberg) is computed per
#' category, the way STRING reports it.
run_ora <- function(hits, background, terms, min_size = 5, max_size = 500) {
  hits <- intersect(unique(hits), background)
  terms <- terms |> filter(protein_key %in% background) |>
    distinct(protein_key, category, term, term_description)

  terms |>
    group_by(category) |>
    mutate(universe = n_distinct(protein_key),
           hits_in_category = n_distinct(protein_key[protein_key %in% hits])) |>
    group_by(category, term, term_description, universe, hits_in_category) |>
    summarise(term_size = n(), hits_in_term = sum(protein_key %in% hits),
              hit_keys = paste(protein_key[protein_key %in% hits], collapse = ";"),
              .groups = "drop") |>
    filter(term_size >= min_size, term_size <= max_size, hits_in_term > 0) |>
    mutate(
      expected = term_size * hits_in_category / universe,
      fold_enrichment = hits_in_term / expected,
      # 2x2 table: in term & hit | in term & not hit | not in term & hit | neither
      p_value = purrr::pmap_dbl(list(hits_in_term, term_size, hits_in_category, universe),
                         \(a, s, h, u) fisher.test(matrix(c(a, s - a, h - a, u - s - h + a), 2),
                                                   alternative = "greater")$p.value)) |>
    group_by(category) |> mutate(FDR = p.adjust(p_value, "BH")) |> ungroup() |>
    arrange(FDR, p_value)
}

#' GSEA-like rank test: are the members of a term shifted towards the top or the
#' bottom of the ranked list? Wilcoxon rank-sum test, members vs non-members.
#'
#' stats: named numeric vector, names = protein_key, values = ranking statistic
#'        (we use the moderated t `statistic`: it combines size and certainty).
#' auc:   probability that a random member ranks above a random non-member
#'        (0.5 = no shift, > 0.5 = towards high values, < 0.5 = towards low).
run_rank_enrichment <- function(stats, terms, min_size = 5, max_size = 500) {
  stats <- stats[!is.na(stats)]
  terms <- terms |> filter(protein_key %in% names(stats)) |>
    distinct(protein_key, category, term, term_description)

  terms |>
    group_by(category, term, term_description) |>
    summarise(members = list(protein_key), term_size = n(), .groups = "drop") |>
    filter(term_size >= min_size, term_size <= max_size) |>
    mutate(test = purrr::map(members, \(m) {
      inside <- stats[names(stats) %in% m]; outside <- stats[!names(stats) %in% m]
      w <- wilcox.test(inside, outside, exact = FALSE)
      tibble(auc = unname(w$statistic) / (length(inside) * length(outside)),
             median_stat_in_term = median(inside), p_value = w$p.value)
    })) |>
    select(-members) |> tidyr::unnest(test) |>
    mutate(direction = if_else(auc > 0.5, "top (high statistic)", "bottom (low statistic)")) |>
    group_by(category) |> mutate(FDR = p.adjust(p_value, "BH")) |> ungroup() |>
    arrange(FDR, p_value)
}

#' Write locus tags, one per line, for STRING "Multiple proteins" with the
#' course organism STRG0B09NBK (JE2_Bio253_HS26). Returns the path invisibly.
write_string_upload <- function(locus_tags, path) {
  writeLines(sort(unique(na.omit(locus_tags))), path)
  invisible(path)
}
