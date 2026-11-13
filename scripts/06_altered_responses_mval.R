# Step 6 (PROJECT_WORKFLOW_mval.md): altered total-protein responses to PASN (H2).
#
# Run from the repository root:   Rscript scripts/06_altered_responses_mval.R   (~30 s)
# Needs scripts/05_formal_test_mval.R and scripts/05b_contrast_sheet_mval.R first.
#
# ---- Plan (coding rules, section 4) -----------------------------------------
# Aim:      Does a mutation change the RESPONSE to PASN (SN vs TSB), beyond any
#           difference already present in TSB? Every response is framed relative
#           to the ancestor (presentation feedback, 12 Nov).
# Input:    data/derived/limma_treat_total_mval.tsv.gz  formal test (Step 5)
#           mval_medi/results/step05_formal_test/one_sided_signals_mval.csv
#           read_abundances("total") (raw per-sample values for candidate plots)
#           read_protein_metadata()  gene symbols (ortholog_quality == "good" only)
# Analysis: no new test. Main results = the two interaction contrasts
#             interaction_stk_vs_anc_at_TSB_vs_SN, interaction_stp_vs_anc_at_TSB_vs_SN
#           formal hit = treat(lfc = 1) BH-FDR < 0.05 (Step 5, eligible proteins only).
#           For reading, a "response" is written as SN minus TSB:
#             response_X = log2(X_SN) - log2(X_TSB) = -diff(TSB_vs_SN_at_X)
#           positive = higher in SN. Then
#             interaction_X_vs_anc = response_anc - response_X
#           so a POSITIVE interaction = the mutant's SN response is lower (less up /
#           more down) than the ancestor's. Checked numerically below.
# Output:   mval_medi/results/step06_responses/
#             response_summary_mval.csv       per primary contrast: eligible, hits, effect spread
#             interaction_candidates_mval.csv formal interaction hits with responses, baseline, SE, FDR
#             one_sided_responses_mval.csv    on/off patterns in the interactions (exploratory)
#             fig_effect_distributions_mval.png
#             fig_response_vs_anc_mval.png    Fig 2: mutant response vs anc response
#             fig_candidates_raw_mval.png     per-sample values of every candidate
#             figure_guide_mval.md, responses_checks_mval.txt
# Check:    interaction == response_anc - response_mutant for every eligible protein;
#           candidate count equals Step 5 hit counts; raw plots inspected.

suppressMessages(source("R/setup.R"))
out <- "mval_medi/results/step06_responses"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
checks <- character()
check <- function(label, ok) {
  line <- sprintf("[%s] %s", if (isTRUE(ok)) "PASS" else "FAIL", label)
  message(line); checks <<- c(checks, line)
  invisible(ok)
}

res <- read_tsv("data/derived/limma_treat_total_mval.tsv.gz", show_col_types = FALSE)
one_sided <- read_csv("mval_medi/results/step05_formal_test/one_sided_signals_mval.csv", show_col_types = FALSE)
meta <- read_protein_metadata() |>
  filter(ortholog_quality == "good") |>
  distinct(locus_tag, gene = gene_symbol_best)

BASELINE  <- c("stk_vs_anc_at_TSB", "stp_vs_anc_at_TSB", "stp_vs_stk_at_TSB")
AT_SN     <- c("stk_vs_anc_at_SN", "stp_vs_anc_at_SN", "stp_vs_stk_at_SN")   # only for reading candidates
EXPOSURE  <- c("TSB_vs_SN_at_anc", "TSB_vs_SN_at_stk", "TSB_vs_SN_at_stp")
INTER     <- c("interaction_stk_vs_anc_at_TSB_vs_SN", "interaction_stp_vs_anc_at_TSB_vs_SN")
INTER_ALL <- c(INTER, "interaction_stp_vs_stk_at_TSB_vs_SN")
is_hit <- \(d) d$tested & d$FDR < MVAL_FDR

# ---- 1. summary: baseline, exposure and interaction contrasts separately ---------
summary_tbl <- res |>
  filter(contrast %in% c(BASELINE, EXPOSURE, INTER_ALL), tested) |>
  group_by(contrast) |>
  summarise(eligible = n(),
            hits_up = sum(FDR < MVAL_FDR & diff > 0), hits_down = sum(FDR < MVAL_FDR & diff < 0),
            median_abs_diff = median(abs(diff)), q95_abs_diff = quantile(abs(diff), 0.95),
            n_abs_diff_gt_1_estimate = sum(abs(diff) > 1),   # descriptive, NOT a hit rule
            median_se = median(std.error), .groups = "drop") |>
  mutate(family = case_when(contrast %in% BASELINE ~ "baseline (TSB)",
                            contrast %in% EXPOSURE ~ "exposure within genotype",
                            TRUE ~ "interaction"),
         meaning = purrr::map_chr(contrast, explain_contrast)) |>
  arrange(factor(contrast, c(BASELINE, EXPOSURE, INTER_ALL)))
write_csv(summary_tbl, file.path(out, "response_summary_mval.csv"))
print(select(summary_tbl, contrast, eligible, hits_up, hits_down, median_abs_diff, q95_abs_diff, n_abs_diff_gt_1_estimate))

# ---- 2. responses per genotype, written as SN minus TSB ---------------------------
wide <- res |>
  filter(contrast %in% c(BASELINE, AT_SN, EXPOSURE, INTER_ALL)) |>
  select(protein_key, locus_tag, product, contrast, diff, std.error, FDR, tested) |>
  pivot_wider(names_from = contrast, values_from = c(diff, std.error, FDR, tested), names_glue = "{contrast}__{.value}")
resp <- wide |>
  transmute(protein_key, locus_tag, product,
            resp_anc = -`TSB_vs_SN_at_anc__diff`, resp_stk = -`TSB_vs_SN_at_stk__diff`,
            resp_stp = -`TSB_vs_SN_at_stp__diff`,
            int_stk = `interaction_stk_vs_anc_at_TSB_vs_SN__diff`, int_stp = `interaction_stp_vs_anc_at_TSB_vs_SN__diff`,
            int_stk_tested = `interaction_stk_vs_anc_at_TSB_vs_SN__tested`,
            int_stp_tested = `interaction_stp_vs_anc_at_TSB_vs_SN__tested`,
            int_stk_FDR = `interaction_stk_vs_anc_at_TSB_vs_SN__FDR`, int_stp_FDR = `interaction_stp_vs_anc_at_TSB_vs_SN__FDR`)
dev <- with(resp, c(abs(int_stk - (resp_anc - resp_stk))[int_stk_tested],
                    abs(int_stp - (resp_anc - resp_stp))[int_stp_tested]))
check(sprintf("interaction == response_anc - response_mutant for all %d eligible rows (max |dev| %.1e)",
              length(dev), max(dev)), max(dev) < 1e-8)

# ---- 3. candidates: every formal hit of the three interactions ---------------------
cand <- res |>
  filter(contrast %in% INTER_ALL, tested, FDR < MVAL_FDR) |>
  select(contrast, protein_key, locus_tag, product, staa8_uniprot, interaction_diff = diff,
         interaction_se = std.error, interaction_FDR = FDR, min_obs_in_groups, modelName, estimate_type) |>
  left_join(meta, by = "locus_tag") |>
  left_join(select(resp, protein_key, resp_anc, resp_stk, resp_stp), by = "protein_key") |>
  left_join(wide |> select(protein_key,
                           stk_vs_anc_TSB = `stk_vs_anc_at_TSB__diff`, stk_vs_anc_SN = `stk_vs_anc_at_SN__diff`,
                           stp_vs_anc_TSB = `stp_vs_anc_at_TSB__diff`, stp_vs_anc_SN = `stp_vs_anc_at_SN__diff`,
                           stp_vs_stk_TSB = `stp_vs_stk_at_TSB__diff`, stp_vs_stk_SN = `stp_vs_stk_at_SN__diff`),
            by = "protein_key") |>
  # Where does the strain difference sit? The interaction = (difference in TSB) -
  # (difference in SN). Descriptive label on the ESTIMATES (|log2| > 1), not a test.
  mutate(pair = str_match(contrast, "^interaction_(st[kp]_vs_[a-z]+)_at")[, 2],
         diff_in_TSB = case_when(pair == "stk_vs_anc" ~ stk_vs_anc_TSB, pair == "stp_vs_anc" ~ stp_vs_anc_TSB,
                                 TRUE ~ stp_vs_stk_TSB),
         diff_in_SN  = case_when(pair == "stk_vs_anc" ~ stk_vs_anc_SN, pair == "stp_vs_anc" ~ stp_vs_anc_SN,
                                 TRUE ~ stp_vs_stk_SN),
         where_difference = case_when(
           is.na(diff_in_TSB) | is.na(diff_in_SN) ~ "not computable (group empty)",
           abs(diff_in_TSB) > 1 & abs(diff_in_SN) <= 1 ~ "mainly in TSB (baseline), gone in SN",
           abs(diff_in_TSB) <= 1 & abs(diff_in_SN) > 1 ~ "mainly in SN (response)",
           abs(diff_in_TSB) > 1 & abs(diff_in_SN) > 1 & sign(diff_in_TSB) != sign(diff_in_SN) ~ "both, opposite sign",
           abs(diff_in_TSB) > 1 & abs(diff_in_SN) > 1 ~ "both, same sign",
           TRUE ~ "small in both (< 1 each)")) |>
  mutate(in_words = case_when(
    str_detect(contrast, "stk_vs_anc") & interaction_diff > 0 ~ "SN response LOWER in stk than in anc",
    str_detect(contrast, "stk_vs_anc") ~ "SN response HIGHER in stk than in anc",
    str_detect(contrast, "stp_vs_anc") & interaction_diff > 0 ~ "SN response LOWER in stp than in anc",
    str_detect(contrast, "stp_vs_anc") ~ "SN response HIGHER in stp than in anc",
    interaction_diff > 0 ~ "SN response LOWER in stp than in stk",
    TRUE ~ "SN response HIGHER in stp than in stk")) |>
  arrange(contrast, interaction_FDR)
step5_counts <- read_csv("mval_medi/results/step05_formal_test/hit_counts_mval.csv", show_col_types = FALSE) |>
  filter(contrast %in% INTER_ALL)
check("candidate count equals the Step 5 formal hit counts", nrow(cand) == sum(step5_counts$formal_hits))
write_csv(cand, file.path(out, "interaction_candidates_mval.csv"))
print(select(cand, contrast, locus_tag, gene, interaction_diff, interaction_se, interaction_FDR,
             resp_anc, resp_stk, resp_stp, min_obs_in_groups), width = 200)
print(select(cand, pair, locus_tag, gene, diff_in_TSB, diff_in_SN, where_difference), width = 200)

# ---- 4. one-sided (on/off) signals in the interactions: exploratory --------------
one_sided_int <- one_sided |>
  filter(contrast %in% INTER_ALL) |>
  mutate(empty_groups = str_extract_all(obs_pattern, "[a-z]+_[A-Za-z]+(?==0/3)") |> purrr::map_chr(paste, collapse = "+")) |>
  left_join(meta, by = "locus_tag")
# Which on/off patterns could mean a DIFFERENT response? Only an on/off that
# happens in one strain and one medium. Both strains empty in the same medium =
# shared response; one strain empty in both media = absent from that strain.
classify_empty <- function(empty) {
  parts <- str_split(empty, "\\+")[[1]]
  geno <- unique(str_extract(parts, "^[a-z]+")); med <- unique(str_extract(parts, "[A-Za-z]+$"))
  if (length(geno) == 1 && length(med) == 1) paste0("on/off in ", geno, " only (", med, " empty): candidate altered response")
  else if (length(geno) == 2 && length(med) == 1 && length(parts) == 2) paste0("shared: both strains empty in ", med)
  else if (length(geno) == 1 && length(med) == 2) paste0("absent from ", geno, " in both media (not a response)")
  else "mixed (3+ empty groups)"
}
one_sided_int <- one_sided_int |> mutate(pattern_class = purrr::map_chr(empty_groups, classify_empty))
one_sided_summary <- one_sided_int |> count(contrast, pattern_class, name = "n_proteins") |> arrange(contrast, desc(n_proteins))
write_csv(one_sided_int, file.path(out, "one_sided_responses_mval.csv"))
print(one_sided_summary, n = 40)

# ---- 5. figures --------------------------------------------------------------------
# 5a effect distributions: reference response next to the two interactions
dist_df <- res |> filter(contrast %in% c("TSB_vs_SN_at_anc", INTER), tested) |>
  mutate(contrast = factor(contrast, c("TSB_vs_SN_at_anc", INTER)))
p_dist <- ggplot(dist_df, aes(diff)) +
  geom_histogram(binwidth = 0.1, fill = "grey55") +
  geom_vline(xintercept = c(-MVAL_LFC, MVAL_LFC), linetype = 2, colour = "#D55E00") +
  facet_wrap(~contrast, ncol = 1, scales = "free_y") +
  coord_cartesian(xlim = c(-5, 5)) +
  labs(x = "diff (log2)", y = "eligible proteins",
       title = "Effect sizes: reference SN response in anc vs the two interactions",
       subtitle = "TSB_vs_SN_at_anc > 0 = higher in TSB (lower in SN). Interaction > 0 = SN response lower in the mutant than in anc.\nDashed: |log2| = 1, the treat threshold. Eligible proteins only.") +
  theme(plot.subtitle = element_text(size = 8))
ggsave(file.path(out, "fig_effect_distributions_mval.png"), p_dist, width = 8, height = 7, dpi = 300)

# 5b Fig 2: response in the mutant against the response in anc
sc <- bind_rows(
  resp |> filter(int_stk_tested) |> transmute(locus_tag, mutant = "stk", resp_anc, resp_mut = resp_stk,
                                              hit = int_stk_FDR < MVAL_FDR),
  resp |> filter(int_stp_tested) |> transmute(locus_tag, mutant = "stp", resp_anc, resp_mut = resp_stp,
                                              hit = int_stp_FDR < MVAL_FDR)) |>
  left_join(meta, by = "locus_tag") |>
  mutate(label = coalesce(gene, locus_tag),
         mutant = factor(mutant, c("stk", "stp"),
                         c("stk vs anc  (interaction_stk_vs_anc_at_TSB_vs_SN)",
                           "stp vs anc  (interaction_stp_vs_anc_at_TSB_vs_SN)\nstp also lacks most Stk1")))
cors <- sc |> group_by(mutant) |>
  summarise(r = cor(resp_anc, resp_mut), n = n(), hits = sum(hit), .groups = "drop") |>
  mutate(text = sprintf("n = %d eligible\nPearson r = %.2f\nformal hits = %d", n, r, hits))
p_sc <- ggplot(sc, aes(resp_anc, resp_mut)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey40") +
  geom_abline(slope = 1, intercept = c(-MVAL_LFC, MVAL_LFC), colour = "grey60", linetype = 2) +
  geom_point(data = filter(sc, !hit), alpha = 0.25, size = 0.7, colour = "grey45") +
  geom_point(data = filter(sc, hit), colour = "#D55E00", size = 2) +
  geom_text(data = filter(sc, hit), aes(label = label), colour = "#D55E00", size = 3, vjust = -0.8) +
  geom_text(data = cors, aes(x = -Inf, y = Inf, label = text), hjust = -0.1, vjust = 1.2, size = 3) +
  facet_wrap(~mutant) +
  coord_equal() +
  labs(x = "SN response in anc: log2(SN) - log2(TSB)", y = "SN response in the mutant: log2(SN) - log2(TSB)",
       title = "Is the response to PASN different in the mutants? (relative to the ancestor)",
       subtitle = paste0("On the diagonal = same response as anc. Below = response lower in the mutant (interaction > 0).\n",
                         "Orange = formal interaction hit (treat lfc = 1, BH-FDR < 0.05). Dashed = diagonal +/- 1.")) +
  theme(plot.subtitle = element_text(size = 8), strip.text = element_text(size = 8))
ggsave(file.path(out, "fig_response_vs_anc_mval.png"), p_sc, width = 10, height = 5.8, dpi = 300)

# 5c raw per-sample values of every candidate (CLAUDE.md rule 4)
cand_tags <- unique(cand$locus_tag)
p_raw <- plot_protein(cand_tags, abundances = read_abundances("total"), max_proteins = 12) +
  labs(title = "Every formal interaction hit: per-sample values, all 27 samples",
       subtitle = "Interactions use TSB and SN only; SNpulse shown for context. Read: does the TSB->SN change differ between strains?")
ggsave(file.path(out, "fig_candidates_raw_mval.png"), p_raw, width = 11, height = 8, dpi = 200)

writeLines(c(paste("Step 6 altered responses, run", Sys.Date()), checks),
           file.path(out, "responses_checks_mval.txt"))
if (any(grepl("^\\[FAIL", checks))) stop("At least one check FAILED, see ", out)
message("done, written to ", out)
