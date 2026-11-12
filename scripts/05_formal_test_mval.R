# Step 5 (PROJECT_WORKFLOW_mval.md): formal differential-abundance test, total proteome.
#
# Aim      Replace the post-hoc screen "FDR < 0.1, then |log2FC| > 1" by a test
#          whose null hypothesis already contains the fold-change threshold.
#          Thresholds fixed by the group on 2026-11-12, BEFORE seeing results.
# Input    read_abundances("total"): normalized_abundance (already log2, NOT
#          re-normalised), read_samples(); read_dea("total") only for comparison.
# Analysis limma, model  ~ 0 + group  (group = genotype x medium, 9 cell means,
#          3 replicates each, all 27 samples so the variance estimate uses all
#          data, as in the supplied model). Contrasts written out as A - B
#          (mval_contrasts()), so the sign cannot depend on a reference level.
#          treat(lfc = 1): H0 |log2FC| <= 1. Hit = BH-FDR < 0.05 (alpha 0.05),
#          BH within each contrast over the TESTED proteins only.
#          A contrast is tested for a protein only if every group in it has
#          >= 2 observed values. No imputation: NA stays NA.
#          No batch term: all 27 total-proteome runs were acquired on one day
#          (20260612) and `reinjected` is FALSE for all of them. Run order is in
#          blocks by group, so drift could not be separated from group anyway
#          (limitation, not fixable by a model).
# Output   data/derived/limma_treat_total_mval.tsv.gz      full table (git-ignored)
#          mval_medi/results/step05_formal_test/
#            hit_counts_mval.csv     formal hits vs supplied screen, per contrast
#            hits_formal_mval.csv    all formal hits (all contrasts)
#            sign_check_mval.csv     our diff vs supplied diff, per contrast
#            fig_volcano_primary_mval.png, fig_raw_top_hits_mval.png
#            formal_test_checks_mval.txt
# Check    27 samples / 3 per group; our diff correlates > 0.95 with supplied
#          diff on observed rows (same sign convention); every formal hit has
#          |diff| > 1; Stp1 drops in stp and Stk1 in stk (known answer).
#
# Run from the repository root:  Rscript scripts/05_formal_test_mval.R

suppressMessages(source("R/setup.R"))
out <- "mval_medi/results/step05_formal_test"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
checks <- character()
check <- function(label, ok) {
  line <- sprintf("[%s] %s", if (isTRUE(ok)) "PASS" else "FAIL", label)
  message(line); checks <<- c(checks, line)
  invisible(ok)
}

# ---- 1. data, joined by sample name --------------------------------------------
samples <- read_samples() |> arrange(genotype, medium, replicate)
abund   <- read_abundances("total") |> keep_je2()
check("27 samples, 3 per group", nrow(samples) == 27 && all(table(samples$group) == 3))
check("total proteome: no reinjected samples (no batch covariate needed)", !any(abund$reinjected))

m <- abundance_matrix(abund, samples)
message(nrow(m), " JE2 proteins x ", ncol(m), " samples")

# ---- 2. fit ----------------------------------------------------------------------
res <- fit_limma_treat(m, samples) |>
  left_join(distinct(abund, protein_key, locus_tag, product, staa8_uniprot), by = "protein_key")
stopifnot(nrow(res) == nrow(m) * 33, !anyDuplicated(res[c("protein_key", "contrast")]))
write_tsv(res, "data/derived/limma_treat_total_mval.tsv.gz")

# ---- 3. guardrails ------------------------------------------------------------------
supplied <- read_dea("total") |> keep_je2() |>
  select(protein_key, contrast, diff_supplied = diff, FDR_supplied = FDR, estimate_type)
check("same 33 contrast names as the supplied DEA", setequal(unique(res$contrast), unique(supplied$contrast)))

sign_check <- res |> filter(tested) |>
  inner_join(filter(supplied, estimate_type == "observed"), by = c("protein_key", "contrast")) |>
  group_by(contrast) |>
  summarise(n = n(), cor_diff = cor(diff, diff_supplied),
            median_abs_dev = median(abs(diff - diff_supplied)), .groups = "drop")
write_csv(sign_check, file.path(out, "sign_check_mval.csv"))
check(sprintf("sign convention: cor(our diff, supplied diff) > 0.95 in all contrasts (min %.3f)",
              min(sign_check$cor_diff)), all(sign_check$cor_diff > 0.95))

hits <- res |> filter(tested, FDR < MVAL_FDR)
check(sprintf("every formal hit has |diff| > %g (treat property)", MVAL_LFC), all(abs(hits$diff) > MVAL_LFC))

# known answer: the mutated protein should be lower in its own mutant
known <- res |>
  filter(locus_tag %in% KEY_PROTEINS$locus_tag[1:2],
         contrast %in% c("stp_vs_anc_at_TSB", "stk_vs_anc_at_TSB")) |>
  left_join(select(KEY_PROTEINS, locus_tag, short_name), by = "locus_tag") |>
  select(short_name, contrast, diff, FDR, tested, min_obs_in_groups)
print(known)
stp1 <- filter(known, short_name == "Stp1", contrast == "stp_vs_anc_at_TSB")
stk1 <- filter(known, short_name == "Stk1", contrast == "stk_vs_anc_at_TSB")
check("known answer: Stp1 lower in stp (TSB) or not measurable there",
      !stp1$tested || stp1$diff < 0)
check("known answer: Stk1 lower in stk (TSB) or not measurable there",
      !stk1$tested || stk1$diff < 0)

# ---- 4. counts: formal test vs supplied screen --------------------------------------
screen <- supplied |>
  group_by(contrast) |>
  summarise(screen_hits_supplied = sum(FDR_supplied < FDR_CUTOFF & abs(diff_supplied) > LFC_CUTOFF, na.rm = TRUE),
            screen_hits_lod_imputed = sum(FDR_supplied < FDR_CUTOFF & abs(diff_supplied) > LFC_CUTOFF &
                                            estimate_type == "lod_imputed", na.rm = TRUE),
            .groups = "drop")
counts <- res |>
  group_by(contrast) |>
  summarise(n_proteins = n(), n_tested = sum(tested),
            formal_up = sum(tested & FDR < MVAL_FDR & diff > 0),
            formal_down = sum(tested & FDR < MVAL_FDR & diff < 0), .groups = "drop") |>
  mutate(formal_hits = formal_up + formal_down) |>
  left_join(screen, by = "contrast")
write_csv(counts, file.path(out, "hit_counts_mval.csv"))
print(counts, n = 33)

hits |>
  select(contrast, locus_tag, product, staa8_uniprot, diff, std.error, treat_t, FDR, min_obs_in_groups) |>
  arrange(contrast, FDR) |>
  write_csv(file.path(out, "hits_formal_mval.csv"))

# ---- 5. figures ------------------------------------------------------------------------
# primary 18-sample contrasts (TSB and SN only), lines at the SAME thresholds as the test
primary <- c("TSB_vs_SN_at_anc", "TSB_vs_SN_at_stk", "TSB_vs_SN_at_stp",
             "stk_vs_anc_at_TSB", "stp_vs_anc_at_TSB", "stp_vs_stk_at_TSB",
             "interaction_stk_vs_anc_at_TSB_vs_SN", "interaction_stp_vs_anc_at_TSB_vs_SN",
             "interaction_stp_vs_stk_at_TSB_vs_SN")
vol <- res |> filter(contrast %in% primary, tested) |>
  mutate(contrast = factor(contrast, primary),
         class = case_when(FDR < MVAL_FDR & diff > 0 ~ "hit, higher in A",
                           FDR < MVAL_FDR & diff < 0 ~ "hit, higher in B",
                           TRUE ~ "no evidence |log2FC| > 1"))
p <- ggplot(vol, aes(diff, -log10(FDR), colour = class)) +
  geom_point(alpha = 0.6, size = 0.9) +
  geom_vline(xintercept = c(-MVAL_LFC, MVAL_LFC), linetype = 2, colour = "grey50") +
  geom_hline(yintercept = -log10(MVAL_FDR), linetype = 2, colour = "grey50") +
  scale_colour_manual(values = c("hit, higher in A" = "#D55E00", "hit, higher in B" = "#0072B2",
                                 "no evidence |log2FC| > 1" = "grey70")) +
  facet_wrap(~contrast, ncol = 3, labeller = label_wrap_gen(28)) +
  labs(x = "diff = log2(A) - log2(B)   (contrast A_vs_B)", y = "-log10(treat BH-FDR)", colour = NULL,
       title = "Formal test: limma treat, H0 |log2FC| <= 1, BH-FDR < 0.05",
       subtitle = "Model ~ 0 + group, 27 samples, n = 3 per group; observed values only (no imputation)") +
  theme(legend.position = "bottom", strip.text = element_text(size = 7))
ggsave(file.path(out, "fig_volcano_primary_mval.png"), p, width = 9, height = 9, dpi = 150)

# CLAUDE.md rule 4: look at the raw per-sample values behind the strongest hits
top <- hits |> filter(contrast == "stk_vs_anc_at_TSB") |> slice_min(FDR, n = 6, with_ties = FALSE)
if (nrow(top) > 0) {
  ggsave(file.path(out, "fig_raw_top_hits_mval.png"),
         plot_protein(top$locus_tag, abundances = read_abundances("total")) +
           labs(title = "Top formal hits stk_vs_anc_at_TSB: raw per-sample values"),
         width = 10, height = 7, dpi = 150)
}

writeLines(c(paste("Step 5 formal test, run", Sys.Date()),
             sprintf("Model ~ 0 + group; limma treat(lfc = %g); BH-FDR < %g; min %d obs per group",
                     MVAL_LFC, MVAL_FDR, MVAL_MIN_OBS),
             checks), file.path(out, "formal_test_checks_mval.txt"))
if (any(grepl("^\\[FAIL", checks))) stop("At least one check FAILED, see ", out)
message("done, written to ", out)
