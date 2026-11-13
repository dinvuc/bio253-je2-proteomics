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
#                                             (pink = exploratory 0.05 <= FDR < 0.1, not a hit)
#             fig_response_stk_vs_stp_mval.png  stp response vs stk response (direct interaction)
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
            int_stk_FDR = `interaction_stk_vs_anc_at_TSB_vs_SN__FDR`, int_stp_FDR = `interaction_stp_vs_anc_at_TSB_vs_SN__FDR`,
            int_direct = `interaction_stp_vs_stk_at_TSB_vs_SN__diff`,
            int_direct_tested = `interaction_stp_vs_stk_at_TSB_vs_SN__tested`,
            int_direct_FDR = `interaction_stp_vs_stk_at_TSB_vs_SN__FDR`)
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
# A hit in one mutant says nothing about the other mutant unless that one was
# testable too. Add, per candidate, the detections in all 6 TSB/SN groups and
# the OTHER mutant's interaction: our estimate if eligible, else the supplied
# lm_impute estimate (presence/absence evidence only, labelled).
samples <- read_samples()
obs <- read_abundances("total") |> keep_je2() |>
  filter(protein_key %in% cand$protein_key) |>
  select(protein_key, sample, normalized_abundance) |>
  inner_join(select(samples, sample, genotype, medium), by = "sample") |>   # by NAME
  filter(medium %in% c("TSB", "SN")) |>
  group_by(protein_key, genotype, medium) |>
  summarise(n = sum(!is.na(normalized_abundance)), .groups = "drop") |>
  arrange(protein_key, genotype, desc(medium == "TSB")) |>
  group_by(protein_key) |>
  summarise(detected_TSB_SN = paste(sprintf("%s %s %d/3", genotype, medium, n), collapse = ", "), .groups = "drop")
supplied_int <- read_dea("total") |> keep_je2() |>
  filter(contrast %in% INTER) |>
  select(protein_key, contrast, other_supplied_diff = diff, other_supplied_FDR = FDR, other_supplied_type = estimate_type)
other_of <- c(interaction_stk_vs_anc_at_TSB_vs_SN = "interaction_stp_vs_anc_at_TSB_vs_SN",
              interaction_stp_vs_anc_at_TSB_vs_SN = "interaction_stk_vs_anc_at_TSB_vs_SN")
cand <- cand |>
  left_join(obs, by = "protein_key") |>
  mutate(other_contrast = unname(other_of[contrast])) |>
  left_join(res |> select(protein_key, other_contrast = contrast, other_tested = tested,
                          other_diff = diff, other_FDR = FDR),
            by = c("protein_key", "other_contrast")) |>
  left_join(rename(supplied_int, other_contrast = contrast), by = c("protein_key", "other_contrast")) |>
  mutate(other_mutant = case_when(
    is.na(other_contrast) ~ "not applicable (direct stp_vs_stk contrast)",
    other_tested & other_FDR < MVAL_FDR ~ "also a formal hit in the other mutant",
    other_tested ~ sprintf("testable in the other mutant, not a hit (diff %.2f, FDR %.2f)", other_diff, other_FDR),
    TRUE ~ sprintf("NOT testable in the other mutant (< 2/3 in a group); supplied %s estimate %.2f, FDR %.3f",
                   other_supplied_type, other_supplied_diff, other_supplied_FDR)))
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

# Hit classes for the figures. The hit rule stays BH-FDR < MVAL_FDR (0.05, fixed
# 12 Nov before results). On request (13 Nov) proteins with 0.05 <= FDR < 0.1 are
# SHOWN as a separate, labelled exploratory layer. They are not hits, are not
# counted as hits anywhere, and the threshold was chosen after seeing results.
FDR_EXPLORATORY <- 0.1
hit_class <- function(tested, fdr) {
  factor(case_when(tested & fdr < MVAL_FDR ~ "formal hit (FDR < 0.05)",
                   tested & fdr < FDR_EXPLORATORY ~ "exploratory (0.05 <= FDR < 0.1)",
                   TRUE ~ "no evidence"),
         c("no evidence", "exploratory (0.05 <= FDR < 0.1)", "formal hit (FDR < 0.05)"))
}
CLASS_COLOURS <- c("formal hit (FDR < 0.05)" = "#D55E00", "exploratory (0.05 <= FDR < 0.1)" = "#CC79A7",
                   "no evidence" = "grey45")
draw_response_scatter <- function(d, cors, x_lab, y_lab, title, subtitle) {
  ggplot(d, aes(x, y)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey40") +
    geom_abline(slope = 1, intercept = c(-MVAL_LFC, MVAL_LFC), colour = "grey60", linetype = 2) +
    geom_point(data = filter(d, class == "no evidence"), alpha = 0.25, size = 0.7, colour = "grey45") +
    geom_point(data = filter(d, class != "no evidence"), aes(colour = class), size = 2) +
    geom_text(data = filter(d, class != "no evidence"), aes(label = label, colour = class),
              size = 3, vjust = -0.8, show.legend = FALSE) +
    # geom_label with hjust = 0 keeps all lines left-aligned; padding keeps it off the border
    geom_label(data = cors, aes(x = -Inf, y = Inf, label = text), hjust = 0, vjust = 1, size = 3,
               label.size = 0, fill = NA, label.padding = unit(0.6, "lines")) +
    scale_colour_manual(values = CLASS_COLOURS, drop = TRUE) +
    coord_equal() +
    labs(x = x_lab, y = y_lab, colour = NULL, title = title, subtitle = subtitle) +
    theme(plot.subtitle = element_text(size = 8), strip.text = element_text(size = 8), legend.position = "bottom")
}
count_text <- function(d, x, y) {
  d |> summarise(r = cor({{ x }}, {{ y }}), n = n(),
                 hits = sum(class == "formal hit (FDR < 0.05)"),
                 expl = sum(class == "exploratory (0.05 <= FDR < 0.1)"), .groups = "drop") |>
    mutate(text = sprintf("n = %d eligible\nPearson r = %.2f\nformal hits (FDR < 0.05) = %d\nexploratory (FDR < 0.1) = +%d",
                          n, r, hits, expl))
}

# 5b Fig 2: response in the mutant against the response in anc
sc <- bind_rows(
  resp |> filter(int_stk_tested) |> transmute(locus_tag, mutant = "stk", x = resp_anc, y = resp_stk,
                                              class = hit_class(int_stk_tested, int_stk_FDR)),
  resp |> filter(int_stp_tested) |> transmute(locus_tag, mutant = "stp", x = resp_anc, y = resp_stp,
                                              class = hit_class(int_stp_tested, int_stp_FDR))) |>
  left_join(meta, by = "locus_tag") |>
  mutate(label = coalesce(gene, locus_tag),
         mutant = factor(mutant, c("stk", "stp"),
                         c("stk vs anc  (interaction_stk_vs_anc_at_TSB_vs_SN)",
                           "stp vs anc  (interaction_stp_vs_anc_at_TSB_vs_SN)\nstp also lacks most Stk1")))
p_sc <- draw_response_scatter(
  sc, sc |> group_by(mutant) |> count_text(x, y),
  "SN response in anc: log2(SN) - log2(TSB)", "SN response in the mutant: log2(SN) - log2(TSB)",
  "Is the response to PASN different in the mutants? (relative to the ancestor)",
  paste0("On the diagonal = same response as anc. Below = response lower in the mutant (interaction > 0). Dashed = diagonal +/- 1.\n",
         "Hit rule: treat lfc = 1, BH-FDR < 0.05 (fixed before results). Pink = exploratory only (0.05 <= FDR < 0.1, chosen after results; not a hit).")) +
  facet_wrap(~mutant)
ggsave(file.path(out, "fig_response_vs_anc_mval.png"), p_sc, width = 10, height = 6.2, dpi = 300)

# 5b' separate figure: stp response against stk response (direct comparison).
# interaction_stp_vs_stk_at_TSB_vs_SN = (stp_TSB - stp_SN) - (stk_TSB - stk_SN)
#                                     = response_stk - response_stp
# -> positive = below the diagonal = SN response lower in stp than in stk.
# No anc term here, so the correlation has no shared-anc-reference problem.
dev_direct <- with(filter(resp, int_direct_tested), abs(int_direct - (resp_stk - resp_stp)))
check(sprintf("interaction_stp_vs_stk == response_stk - response_stp for all %d eligible rows (max |dev| %.1e)",
              length(dev_direct), max(dev_direct)), max(dev_direct) < 1e-8)
sc_direct <- resp |> filter(int_direct_tested) |>
  transmute(locus_tag, x = resp_stk, y = resp_stp, class = hit_class(int_direct_tested, int_direct_FDR)) |>
  left_join(meta, by = "locus_tag") |>
  mutate(label = coalesce(gene, locus_tag))
p_direct <- draw_response_scatter(
  sc_direct, count_text(sc_direct, x, y),
  "SN response in stk: log2(SN) - log2(TSB)", "SN response in stp: log2(SN) - log2(TSB)",
  "Do the two mutants respond to PASN in the same way?\nstp vs stk (interaction_stp_vs_stk_at_TSB_vs_SN)",
  paste0("On the diagonal = same response in stp and stk. Below = response lower in stp than in stk\n",
         "(interaction > 0). Dashed = diagonal +/- 1. Hit rule: treat lfc = 1, BH-FDR < 0.05.\n",
         "Pink = exploratory only (0.05 <= FDR < 0.1, chosen after results; not a hit). stp also lacks most Stk1."))
ggsave(file.path(out, "fig_response_stk_vs_stp_mval.png"), p_direct, width = 7.5, height = 8, dpi = 300)
print(sc_direct |> filter(class != "no evidence") |> select(locus_tag, gene, x, y, class))

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
