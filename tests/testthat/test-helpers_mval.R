# Tests for R/helpers_mval.R (group mval). Toy data only -- no course data needed.

test_that("keep_je2 drops contaminants (no locus tag) and reversed decoys", {
  x <- tibble(protein_key = c("lcl_CP020619.1_prot_1159", "zz_CON_K1C9",
                              "rev_tr|lcl_CP020619.1_prot_ARG44775.1_262", "lcl_CP020619.1_prot_ARG45702.1_1"),
              locus_tag = c("B7H15_05930", NA, "B7H15_00001", "B7H15_06350"))
  expect_equal(keep_je2(x)$locus_tag, c("B7H15_05930", "B7H15_06350"))
})

test_that("hit_list splits by direction and labels imputed estimates", {
  dea <- tibble(protein_key = paste0("p", 1:4), locus_tag = paste0("L", 1:4), product = "x",
                contrast = "stk_vs_anc", diff = c(2, -3, 0.5, 4), FDR = c(0.01, 0.01, 0.01, 0.5),
                statistic = c(5, -6, 1, 2), estimate_type = c("observed", "lod_imputed", "observed", "observed"))
  up <- hit_list(dea, "stk_vs_anc", direction = "up")
  down <- hit_list(dea, "stk_vs_anc", direction = "down")
  expect_equal(up$protein_key, "p1")          # p4 fails FDR, p3 fails |diff|
  expect_equal(down$protein_key, "p2")
  expect_true(down$imputed)
})

test_that("run_ora gives the Fisher p-value of the 2x2 table and uses the background", {
  bg <- paste0("p", 1:100)
  terms <- bind_rows(
    tibble(protein_key = paste0("p", 1:10), term = "T1"),     # term of size 10
    tibble(protein_key = paste0("p", 11:100), term = "T2")) |>
    mutate(category = "GO", term_description = term)
  hits <- c(paste0("p", 1:8), paste0("p", 11:12), "not_in_background")
  res <- run_ora(hits, bg, terms, min_size = 1)
  t1 <- filter(res, term == "T1")
  expect_equal(t1$hits_in_term, 8)
  expect_equal(t1$hits_in_category, 10)        # "not_in_background" is ignored
  expect_equal(t1$p_value,
               fisher.test(matrix(c(8, 2, 2, 88), 2), alternative = "greater")$p.value)
  expect_lt(t1$FDR, 0.001)
})

test_that("run_rank_enrichment detects a term shifted to the top", {
  stats <- setNames(c(10:1, rep(0, 40)), paste0("p", 1:50))
  terms <- tibble(protein_key = paste0("p", 1:5), category = "GO", term = "TOP",
                  term_description = "top term")
  res <- run_rank_enrichment(stats, terms, min_size = 1)
  expect_gt(res$auc, 0.9)
  expect_equal(res$direction, "top (high statistic)")
  expect_lt(res$p_value, 0.01)
})

test_that("mval_contrasts writes A_vs_B as A - B", {
  ct <- mval_contrasts()
  expect_length(ct, 33)
  expect_equal(unname(ct["TSB_vs_SN_at_anc"]), "anc_TSB - anc_SN")
  expect_equal(unname(ct["interaction_stp_vs_anc_at_TSB_vs_SN"]),
               "(stp_TSB - stp_SN) - (anc_TSB - anc_SN)")
})

test_that("fit_limma_treat: sign, treat threshold and min_obs on toy data", {
  set.seed(1)
  samples <- tidyr::expand_grid(genotype = c("anc", "stk", "stp"),
                                medium = c("TSB", "SN", "SNpulse"), replicate = 1:3) |>
    mutate(group = paste0(genotype, "_", medium), sample = paste0(group, "_", replicate))
  m <- matrix(rnorm(200 * 27, 20, 0.1), 200, 27, dimnames = list(paste0("p", 1:200), samples$sample))
  m["p1", samples$group == "anc_TSB"] <- m["p1", samples$group == "anc_TSB"] + 3    # big effect
  m["p2", samples$group == "anc_TSB"] <- m["p2", samples$group == "anc_TSB"] + 0.5  # below lfc 1
  m["p3", samples$group == "anc_SN"][1:2] <- NA                                     # only 1 obs
  res <- fit_limma_treat(m, samples) |> filter(contrast == "TSB_vs_SN_at_anc")
  p1 <- filter(res, protein_key == "p1"); p2 <- filter(res, protein_key == "p2")
  expect_gt(p1$diff, 2.5)              # TSB higher -> positive
  expect_lt(p1$FDR, 0.05)
  expect_gt(p2$FDR, 0.05)              # real but < 1 log2: not a hit under treat
  expect_false(filter(res, protein_key == "p3")$tested)
  expect_true(is.na(filter(res, protein_key == "p3")$FDR))
})
