# Step 5, part 2 (PROJECT_WORKFLOW_mval.md): contrast definitions, direction
# check, eligibility and one-sided signals. Part 1 (the formal test itself) is
# scripts/05_formal_test_mval.R and must be run first.
#
# Run from the repository root:   Rscript scripts/05b_contrast_sheet_mval.R   (~20 s)
#
# ---- Plan (coding rules, section 4) -----------------------------------------
# Aim:      Fix WHAT each contrast means before Steps 6-7 read results:
#           (1) sort the 33 contrasts into families (baseline, exposure,
#               interaction ...); (2) prove the direction of every estimate by
#               recomputing it from group means; (3) count eligible proteins per
#               contrast and keep the excluded ones (one-sided on/off, sparse)
#               in a separate table; (4) document the supplied FDR scope.
# Input:    data/derived/limma_treat_total_mval.tsv.gz (step 5 part 1)
#           read_abundances("total") |> keep_je2()   normalised log2 (not changed)
#           read_samples()                           join by sample name
#           read_dea("total", drop_contaminants = FALSE)   supplied lm_impute DEA
# Analysis: no new test. Eligibility rule (fixed in the workflow before results):
#           >= 2 of 3 observed values in EVERY group the contrast uses. This is
#           exactly `tested` of the formal test. Excluded proteins are classed
#             on/off  = >= 2/3 in one group, 0/3 in another (presence/absence)
#             sparse  = everything else (some group 1/3, or never >= 2/3)
# Output:   mval_medi/results/step05_formal_test/
#             contrast_sheet_mval.csv        33 contrasts: family, formula, meaning,
#                                            eligible/excluded counts, hits
#             direction_check_mval.csv       limma diff vs hand-computed group means
#             one_sided_signals_mval.csv     on/off proteins of the 9 primary contrasts
#             fig_eligibility_mval.png       eligible / on-off / sparse per contrast
#             contrast_sheet_checks_mval.txt
#           methods sheet (text): mval_medi/results/step05_formal_test/methods_sheet_mval.md
# Check:    limma diff == group-mean formula for every tested row (|dev| < 1e-8);
#           no eligible row is lod_imputed in the supplied DEA; supplied FDR =
#           BH within contrast over all supplied rows.

suppressMessages(source("R/setup.R"))
out <- "mval_medi/results/step05_formal_test"
formal_file <- "data/derived/limma_treat_total_mval.tsv.gz"
if (!file.exists(formal_file)) stop("Run scripts/05_formal_test_mval.R first.")
checks <- character()
check <- function(label, ok) {
  line <- sprintf("[%s] %s", if (isTRUE(ok)) "PASS" else "FAIL", label)
  message(line); checks <<- c(checks, line)
  invisible(ok)
}

PRIMARY <- c("stk_vs_anc_at_TSB", "stp_vs_anc_at_TSB", "stp_vs_stk_at_TSB",
             "TSB_vs_SN_at_anc", "TSB_vs_SN_at_stk", "TSB_vs_SN_at_stp",
             "interaction_stk_vs_anc_at_TSB_vs_SN", "interaction_stp_vs_anc_at_TSB_vs_SN",
             "interaction_stp_vs_stk_at_TSB_vs_SN")

res     <- read_tsv(formal_file, show_col_types = FALSE)
samples <- read_samples()
abund   <- read_abundances("total") |> keep_je2() |>
  select(protein_key, sample, normalized_abundance) |>
  inner_join(select(samples, sample, group), by = "sample")   # labels BY NAME (pitfall 2)
check("27 samples, 3 per group after the join",
      n_distinct(abund$sample) == 27 && all(table(distinct(abund, sample, group)$group) == 3))
contr <- mval_contrasts()

# ---- 1. contrast families ------------------------------------------------------
family_of <- function(ct) case_when(
  str_starts(ct, "interaction_")                     ~ "interaction: does the medium response differ between genotypes?",
  str_detect(ct, "^(stk|stp)_vs_(anc|stk)_at_TSB$")  ~ "baseline: genotype difference in TSB",
  str_detect(ct, "^(stk|stp)_vs_(anc|stk)_at_")      ~ "genotype difference in a PASN medium",
  str_detect(ct, "^(stk|stp)_vs_(anc|stk)$")         ~ "genotype difference averaged over media",
  str_detect(ct, "_at_(anc|stk|stp)$")               ~ "exposure: medium response within one genotype",
  TRUE                                               ~ "exposure averaged over genotypes")

# ---- 2. per protein x group: observed values and mean of observed values ---------
per_group <- abund |>
  group_by(protein_key, group) |>
  summarise(n_obs = sum(!is.na(normalized_abundance)),
            mean_obs = if (n_obs > 0) mean(normalized_abundance, na.rm = TRUE) else NA_real_,
            .groups = "drop")
stopifnot(nrow(per_group) == n_distinct(abund$protein_key) * 9)
n_obs_wide  <- per_group |> select(-mean_obs) |> pivot_wider(names_from = group, values_from = n_obs)
mean_wide   <- per_group |> select(-n_obs)    |> pivot_wider(names_from = group, values_from = mean_obs)

# ---- 3. direction check: recompute EVERY tested diff from group means -------------
# With the cell-means model ~ 0 + group, limma's coefficient for a group is the
# mean of that group's observed values, so a contrast estimate must equal the same
# formula evaluated on the plain group means. This checks the sign of all 33.
recompute <- purrr::imap(contr, \(expr, ct) {
  tibble(protein_key = mean_wide$protein_key, contrast = ct,
         diff_by_hand = eval(parse(text = expr), envir = mean_wide))
}) |> bind_rows()
dir_check <- res |> filter(tested) |>
  select(protein_key, locus_tag, contrast, diff) |>
  inner_join(recompute, by = c("protein_key", "contrast")) |>
  mutate(abs_dev = abs(diff - diff_by_hand))
dir_summary <- dir_check |> group_by(contrast) |>
  summarise(n_tested = n(), max_abs_dev = max(abs_dev), .groups = "drop")
write_csv(dir_summary, file.path(out, "direction_check_mval.csv"))
check(sprintf("direction: limma diff equals the group-mean formula in all %d tested rows (max |dev| %.1e)",
              nrow(dir_check), max(dir_check$abs_dev)), max(dir_check$abs_dev) < 1e-8)

# one worked example in words (Stk1, stp_vs_anc_at_TSB)
ex <- mean_wide |> filter(protein_key == res$protein_key[res$locus_tag == "B7H15_06355"][1])
ex_diff <- res |> filter(locus_tag == "B7H15_06355", contrast == "stp_vs_anc_at_TSB") |> pull(diff)
example <- c(explain_contrast("stp_vs_anc_at_TSB"),
             sprintf("Worked example, Stk1 (B7H15_06355): mean stp_TSB = %.2f, mean anc_TSB = %.2f, difference = %.2f; limma diff = %.2f -> Stk1 is HIGHER in anc = LOWER in stp.",
                     ex$stp_TSB, ex$anc_TSB, ex$stp_TSB - ex$anc_TSB, ex_diff),
             explain_contrast("interaction_stp_vs_anc_at_TSB_vs_SN"))
cat(example, sep = "\n")
check("worked example: Stk1 stp - anc in TSB equals limma diff", abs(ex$stp_TSB - ex$anc_TSB - ex_diff) < 1e-8)

# ---- 4. eligibility classes per protein x contrast ---------------------------------
classes <- purrr::imap(contr, \(expr, ct) {
  g <- contrast_groups(expr)
  n <- as.matrix(n_obs_wide[, g])
  tibble(protein_key = n_obs_wide$protein_key, contrast = ct,
         min_obs = apply(n, 1, min), max_obs = apply(n, 1, max),
         obs_pattern = apply(n, 1, \(r) paste(paste0(g, "=", r, "/3"), collapse = " ")))
}) |> bind_rows() |>
  mutate(eligibility = case_when(min_obs >= MVAL_MIN_OBS ~ "eligible",
                                 min_obs == 0 & max_obs >= MVAL_MIN_OBS ~ "on/off (one-sided)",
                                 TRUE ~ "sparse"))
classes <- classes |> left_join(select(res, protein_key, contrast, tested, diff, std.error, FDR),
                                by = c("protein_key", "contrast"))
check("eligibility rule == `tested` of the formal test (same rule, no drift)",
      all((classes$eligibility == "eligible") == classes$tested))

supplied_all <- read_dea("total", drop_contaminants = FALSE)
# supplied FDR scope: is it BH within each contrast over ALL supplied rows?
scope <- supplied_all |> group_by(contrast) |>
  summarise(n_rows = n(), max_dev = max(abs(p.adjust(p.value, "BH") - FDR)), .groups = "drop")
check(sprintf("supplied FDR = BH within each contrast over all %d supplied rows (incl. contaminants/decoys)",
              unique(scope$n_rows)), all(scope$max_dev < 1e-10) && n_distinct(scope$n_rows) == 1)

supplied <- supplied_all |> keep_je2() |>
  select(protein_key, contrast, diff_supplied = diff, FDR_supplied = FDR, estimate_type)
classes <- classes |> left_join(supplied, by = c("protein_key", "contrast"))
# Finding (13 Nov): the supplied `estimate_type` is a PER-PROTEIN label, not per
# contrast. If any of the 9 groups is 0/3, prolfquapp labels all 33 contrasts of
# that protein `lod_imputed`, even contrasts whose own groups are fully observed.
# So "eligible here" and "lod_imputed there" can both be true; our formal test
# uses only observed values for those rows. Check that this IS the rule:
label_rule <- supplied |>
  group_by(protein_key) |>
  summarise(n_labels = n_distinct(estimate_type), lod = any(estimate_type == "lod_imputed"), .groups = "drop") |>
  inner_join(n_obs_wide |> transmute(protein_key, any_group_empty = apply(across(-protein_key), 1, min) == 0),
             by = "protein_key")
check("supplied estimate_type is constant over the 33 contrasts of each protein", all(label_rule$n_labels == 1))
check(sprintf("supplied lod_imputed <=> at least one of the 9 groups is 0/3 (%d proteins)", sum(label_rule$lod)),
      all(label_rule$lod == label_rule$any_group_empty))

# ---- 5. contrast sheet ---------------------------------------------------------------
sheet <- classes |>
  group_by(contrast) |>
  summarise(n_proteins = n(),
            eligible = sum(eligibility == "eligible"),
            excluded_on_off = sum(eligibility == "on/off (one-sided)"),
            excluded_sparse = sum(eligibility == "sparse"),
            # eligible here, but the supplied row is labelled lod_imputed because
            # some OTHER group of the protein is empty (see finding above)
            eligible_with_supplied_lod_label = sum(eligibility == "eligible" & estimate_type == "lod_imputed", na.rm = TRUE),
            formal_hits_up = sum(tested & FDR < MVAL_FDR & diff > 0, na.rm = TRUE),
            formal_hits_down = sum(tested & FDR < MVAL_FDR & diff < 0, na.rm = TRUE),
            supplied_screen_hits = sum(FDR_supplied < FDR_CUTOFF & abs(diff_supplied) > LFC_CUTOFF, na.rm = TRUE),
            supplied_screen_lod_imputed = sum(FDR_supplied < FDR_CUTOFF & abs(diff_supplied) > LFC_CUTOFF &
                                                estimate_type == "lod_imputed", na.rm = TRUE),
            .groups = "drop") |>
  mutate(family = family_of(contrast),
         primary_18_samples = contrast %in% PRIMARY,
         formula = unname(contr[contrast]),
         meaning = purrr::map_chr(contrast, explain_contrast)) |>
  relocate(contrast, family, primary_18_samples, formula, meaning) |>
  arrange(desc(primary_18_samples), family, contrast)
stopifnot(nrow(sheet) == 33, all(sheet$eligible + sheet$excluded_on_off + sheet$excluded_sparse == sheet$n_proteins))
write_csv(sheet, file.path(out, "contrast_sheet_mval.csv"))
print(select(sheet, contrast, eligible, excluded_on_off, excluded_sparse, formal_hits_up, formal_hits_down,
             supplied_screen_hits) |> filter(contrast %in% PRIMARY))

# ---- 6. one-sided (on/off) signals of the primary contrasts, kept separately ----------
annot <- distinct(res, protein_key, locus_tag, product)
one_sided <- classes |>
  filter(contrast %in% PRIMARY, eligibility == "on/off (one-sided)") |>
  left_join(annot, by = "protein_key") |>
  mutate(supplied_screen = FDR_supplied < FDR_CUTOFF & abs(diff_supplied) > LFC_CUTOFF,
         # "higher in" only exists for simple A_vs_B contrasts; an interaction is a
         # difference of differences, read it with explain_contrast()
         higher_in = case_when(str_starts(contrast, "interaction_") ~ NA_character_,
                               diff_supplied > 0 ~ str_match(contrast, "^([A-Za-z]+)_vs_")[, 2],
                               TRUE ~ str_match(contrast, "^[A-Za-z]+_vs_([A-Za-z]+)")[, 2])) |>
  select(contrast, locus_tag, product, obs_pattern, diff_supplied, FDR_supplied, estimate_type,
         supplied_screen, higher_in) |>
  arrange(contrast, desc(supplied_screen), FDR_supplied)
write_csv(one_sided, file.path(out, "one_sided_signals_mval.csv"))
message("one-sided rows in primary contrasts: ", nrow(one_sided),
        "; of these in the supplied screen: ", sum(one_sided$supplied_screen, na.rm = TRUE))

# ---- 7. figure ------------------------------------------------------------------------------
elig_long <- sheet |> filter(primary_18_samples) |>
  select(contrast, eligible, `on/off (one-sided)` = excluded_on_off, sparse = excluded_sparse) |>
  pivot_longer(-contrast, names_to = "class", values_to = "n") |>
  mutate(contrast = factor(contrast, rev(PRIMARY)),
         class = factor(class, c("sparse", "on/off (one-sided)", "eligible")))
p <- ggplot(elig_long, aes(n, contrast, fill = class)) +
  geom_col() +
  geom_text(data = filter(elig_long, class == "eligible"), aes(label = n), hjust = 1.1, colour = "white", size = 3) +
  scale_fill_manual(values = c(eligible = "#0072B2", `on/off (one-sided)` = "#D55E00", sparse = "grey70")) +
  labs(x = "JE2 proteins (of 2232)", y = NULL, fill = NULL,
       title = "Which proteins can each primary contrast test?",
       subtitle = "eligible = >= 2 of 3 observed in every group the contrast uses (formal test);\non/off = >= 2/3 in one group, 0/3 in another (presence/absence, kept separately)") +
  theme(legend.position = "bottom")
ggsave(file.path(out, "fig_eligibility_mval.png"), p, width = 8, height = 4.5, dpi = 300)

writeLines(c(paste("Step 5 part 2: contrast sheet, run", Sys.Date()), checks, "", example),
           file.path(out, "contrast_sheet_checks_mval.txt"))
if (any(grepl("^\\[FAIL", checks))) stop("At least one check FAILED, see ", out)
message("done, written to ", out)
