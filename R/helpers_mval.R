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

# ==============================================================================
# 6. Formal differential-abundance test (step 5, decided 2026-11-12)
# ==============================================================================
# Group decision (rules section 5): a protein is a HIT when limma::treat() with
# lfc = 1 gives BH-FDR < 0.05. treat() tests H0: |log2FC| <= 1, so the
# fold-change threshold is PART OF THE TEST. We do not filter |log2FC| afterwards.
# The supplied prolfquapp screen (FDR < 0.1, then |diff| > 1) stays as a
# labelled comparison only.
MVAL_FDR <- 0.05   # BH alpha
MVAL_LFC <- 1      # treat() threshold on |log2 fold change|
MVAL_MIN_OBS <- 2  # observed values needed in EACH group of a contrast to test it

#' Wide matrix (proteins x samples) of normalised log2 abundances.
#' Columns are joined BY SAMPLE NAME and ordered like `samples` (pitfall 2).
#' NA stays NA (CLAUDE.md rule 6).
abundance_matrix <- function(abund, samples) {
  stopifnot(setequal(unique(abund$sample), samples$sample),
            !anyDuplicated(abund[c("protein_key", "sample")]))
  wide <- abund |>
    select(protein_key, sample, normalized_abundance) |>
    pivot_wider(names_from = sample, values_from = normalized_abundance)
  m <- as.matrix(wide[, samples$sample])
  rownames(m) <- wide$protein_key
  stopifnot(identical(colnames(m), samples$sample))
  m
}

#' All course contrasts as limma expressions over cell means "<genotype>_<medium>".
#' Naming rule: A_vs_B = A - B; a contrast without "_at_" averages over the
#' other factor; interaction_X_vs_Y_at_A_vs_B = (X_A - X_B) - (Y_A - Y_B).
mval_contrasts <- function(genotypes = c("anc", "stk", "stp"),
                           media = c("TSB", "SN", "SNpulse")) {
  g <- \(geno, med) paste0(geno, "_", med)
  avg <- \(x) paste0("(", paste(x, collapse = " + "), ")/", length(x))
  # pair orders exactly as in the course tables (first - second)
  geno_pairs <- list(c("stk", "anc"), c("stp", "anc"), c("stp", "stk"))
  med_pairs  <- list(c("TSB", "SN"), c("TSB", "SNpulse"), c("SNpulse", "SN"))
  out <- c()
  for (p in med_pairs) {
    out[paste0(p[1], "_vs_", p[2])] <- paste(avg(g(genotypes, p[1])), "-", avg(g(genotypes, p[2])))
    for (ge in genotypes)
      out[paste0(p[1], "_vs_", p[2], "_at_", ge)] <- paste(g(ge, p[1]), "-", g(ge, p[2]))
  }
  for (p in geno_pairs) {
    out[paste0(p[1], "_vs_", p[2])] <- paste(avg(g(p[1], media)), "-", avg(g(p[2], media)))
    for (me in media)
      out[paste0(p[1], "_vs_", p[2], "_at_", me)] <- paste(g(p[1], me), "-", g(p[2], me))
    for (q in med_pairs)
      out[paste0("interaction_", p[1], "_vs_", p[2], "_at_", q[1], "_vs_", q[2])] <-
        sprintf("(%s - %s) - (%s - %s)", g(p[1], q[1]), g(p[1], q[2]), g(p[2], q[1]), g(p[2], q[2]))
  }
  out
}

#' Groups (cell means) that a contrast expression uses.
contrast_groups <- function(expr) unique(str_extract_all(expr, "[a-z]+_[A-Za-z]+")[[1]])

#' limma fit with model ~ 0 + group (one mean per genotype x medium cell, no
#' intercept, so no reference level can flip a sign: every contrast is written
#' out as A - B by mval_contrasts()). Returns one row per protein x contrast with
#'   diff        log2 fold change A - B
#'   statistic   ordinary moderated t (eBayes, H0: diff = 0) -> for ranking
#'   treat_t, p.value, FDR   treat() test of H0: |diff| <= lfc, BH within contrast
#'   tested      FALSE if a group of the contrast has < min_obs observed values
#'               (then p.value / FDR are NA and the protein is NOT in the BH family)
fit_limma_treat <- function(m, samples, lfc = MVAL_LFC, min_obs = MVAL_MIN_OBS) {
  stopifnot(identical(colnames(m), samples$sample))
  group  <- factor(samples$group)
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  contr  <- mval_contrasts()
  cm     <- limma::makeContrasts(contrasts = contr, levels = design)
  colnames(cm) <- names(contr)

  fit  <- limma::lmFit(m, design)
  fit2 <- limma::contrasts.fit(fit, cm)
  eb   <- limma::eBayes(fit2)
  tr   <- limma::treat(fit2, lfc = lfc)

  # observed values per protein x group
  n_obs <- sapply(levels(group), \(gr) rowSums(!is.na(m[, group == gr, drop = FALSE])))

  purrr::imap(contr, \(expr, ct) {
    ok <- apply(n_obs[, contrast_groups(expr), drop = FALSE] >= min_obs, 1, all)
    p  <- ifelse(ok, tr$p.value[, ct], NA_real_)
    tibble(protein_key = rownames(m), contrast = ct,
           diff = eb$coefficients[, ct], std.error = eb$stdev.unscaled[, ct] * sqrt(eb$s2.post),
           statistic = eb$t[, ct], treat_t = tr$t[, ct],
           p.value = p, FDR = p.adjust(p, "BH"),   # BH over tested proteins only
           min_obs_in_groups = apply(n_obs[, contrast_groups(expr), drop = FALSE], 1, min),
           tested = ok)
  }) |>
    bind_rows() |>
    mutate(modelName = paste0("limma ~0+group, treat(lfc=", lfc, ")"),
           estimate_type = "observed")   # no imputation in this model
}
