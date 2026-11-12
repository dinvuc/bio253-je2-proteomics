# Step 3 of our workflow (mval_medi/project_outlines/PROJECT_WORKFLOW_mval.md):
# validate Stp1 and Stk1 first (H1).
#
# Run from the repository root:   Rscript scripts/03_mutant_validation_mval.R   (~30 s)
#
# ---- Plan (coding rules, section 4) -----------------------------------------
# Aim:      Do the mutants lack the enzyme they should lack? Does `stp` also lose
#           Stk1? How much of each "detected" value rests on real peptide evidence?
# Input:    read_abundances("total")       normalised log2 protein abundance, NA = not detected
#           DE_total.xlsx, sheet normalized_abundances:
#             nr_children_protein_Id = peptides quantified per protein PER SAMPLE
#             nrPeptides             = peptides of that protein in the whole data set
#           read_dea("total")              supplied prolfquapp contrasts (lm_impute)
#           limma_treat_total_mval.tsv.gz  our formal test (Step 5)
#           read_sites("STYH"/"STY")       Stk1 phosphosites, peptide, position, loc_prob
#           JE2_proteins.fasta             to check that peptides are unique to Stk1
# Analysis: counting detections and peptides per sample (no new test); the
#           supplied and our formal contrasts are only looked up, not refitted.
# Output:   mval_medi/results/step03_validation/ (small summaries, no raw data)
#             fig_enzyme_abundance_mval.png    Fig 1a: per-sample log2 abundance, label = n peptides
#             fig_enzyme_peptides_mval.png     Fig 1b: peptides quantified per sample
#             fig_stk1_phosphosite_map_mval.png Fig 1c: Stk1 sites, position x sample
#             enzyme_per_sample_mval.csv       one row per protein x sample
#             enzyme_per_group_mval.csv        detections, peptides, mean per group
#             enzyme_contrasts_mval.csv        supplied DEA + our treat test, mutant vs anc
#             stk1_sites_mval.csv              site detections per group, both searches
#             hiq_peptide_check_mval.csv       the HIQDSVPNVTtDVR peptide, per sample
#             evidence_table_mval.csv          rule-based call per enzyme x strain
#             validation_checks_mval.txt       printed CHECK lines
# Check:    27 samples, 3 per group after the join; Stk1 absent in stk (prediction);
#           peptide counts match the abundance NA pattern (NA abundance <=> NA peptides).
#
# ---- Predictions, written BEFORE looking at this script's output -------------
#   anc: Stp1 present, Stk1 present (both in every sample, many peptides).
#   stk: Stk1 strongly reduced / absent (insertion at codon 46 of 664);
#        Stp1 unchanged (it is upstream of the insertion).
#   stp: Stp1 strongly reduced (insertion at codon 189 of 247);
#        Stk1 ALSO reduced if the insertion is polar (stk1 is directly downstream,
#        4 bp overlap) -- this is what Step 4 examines.
#
# Missing values stay NA. In the peptide-count figure, NA is drawn as
# "not detected" at 0 peptides: a count of 0 is a real count, but it is never
# put into a log2 calculation.

suppressMessages(source("R/setup.R"))

out_dir <- "mval_medi/results/step03_validation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
save_fig <- function(p, name, w = 9, h = 5) {
  ggsave(file.path(out_dir, name), p, width = w, height = h, dpi = 300)
}
checks <- character()
check <- function(label, ok, detail = "") {
  line <- sprintf("[%s] %s %s", if (isTRUE(ok)) "PASS" else "FLAG", label, detail)
  checks <<- c(checks, line)
  message(line)
}

enzymes <- KEY_PROTEINS |> filter(short_name %in% c("Stp1", "Stk1")) |>
  select(short_name, locus_tag, protein_id)
samples <- read_samples()

# Insertion codons, recomputed from the genome coordinates (as in notebook 07)
# instead of typed in: codon = ceiling(nt into the CDS / 3).
insertions <- TRANSPOSON_INSERTIONS |>
  mutate(cds_length  = cds_end - cds_start + 1,
         protein_aa  = cds_length %/% 3 - 1,          # last codon is the stop codon
         nt_into_cds = insertion_pos - cds_start + 1,
         codon       = (nt_into_cds + 2) %/% 3)
check("insertion codons are Stk1 46/664 and Stp1 189/247",
      identical(insertions$codon, c(189, 46)) && identical(insertions$protein_aa, c(247, 664)))

# ---- 1. Protein abundance per sample -----------------------------------------
abund <- read_abundances("total") |> keep_je2() |>
  inner_join(enzymes, by = "locus_tag") |>
  select(short_name, locus_tag, sample, genotype, medium, group, replicate, normalized_abundance)

# ---- 2. Peptides quantified per sample (substitute for the missing peptide table)
# Why the Excel: nr_children_protein_Id is not in our derived tables. We identify
# the protein by its protein_id in the FASTA key (never by its peptide count) and
# the sample by parse_sample_name(), which also repairs truncated names.
xl <- suppressWarnings(readxl::read_excel("data/raw/total/DE_total.xlsx",
                                          sheet = "normalized_abundances"))
pep <- xl |>
  mutate(protein_id = str_extract(protein_Id, "ARG\\d+\\.\\d")) |>
  filter(protein_id %in% enzymes$protein_id) |>
  mutate(sample = parse_sample_name(raw.file)$sample) |>
  transmute(protein_id, sample, n_peptides = nr_children_protein_Id,
            n_peptides_dataset = nrPeptides,
            abundance_excel = normalized_abundance)
rm(xl)

per_sample <- abund |>
  left_join(enzymes |> select(locus_tag, protein_id), by = "locus_tag") |>
  left_join(pep, by = c("protein_id", "sample"))   # join BY NAME (pitfall 2)

check("per-sample table: 2 proteins x 27 samples, no duplicates",
      nrow(per_sample) == 54 && !anyDuplicated(per_sample[c("short_name", "sample")]))
check("every sample name in the Excel matches the sample sheet",
      setequal(unique(pep$sample), samples$sample))
check("3 samples in every group for each protein",
      all(count(per_sample, short_name, group)$n == 3))
check("Excel abundance equals the derived abundance (same values, same samples)",
      isTRUE(all.equal(per_sample$abundance_excel, per_sample$normalized_abundance, tolerance = 1e-6)))
check("NA abundance <=> NA peptide count",
      all(is.na(per_sample$normalized_abundance) == is.na(per_sample$n_peptides)))

per_sample <- per_sample |> select(-abundance_excel)
write_csv(per_sample, file.path(out_dir, "enzyme_per_sample_mval.csv"))

# Per group: detections out of 3, peptides, mean of the DETECTED values only.
per_group <- per_sample |>
  group_by(short_name, genotype, medium, group) |>
  summarise(n_detected        = sum(!is.na(normalized_abundance)),
            n_with_2plus_pep  = sum(n_peptides >= 2, na.rm = TRUE),
            peptides          = paste(replace_na(as.character(n_peptides), "-"), collapse = "/"),
            mean_log2_detected = round(mean(normalized_abundance, na.rm = TRUE), 2),
            .groups = "drop") |>
  mutate(mean_log2_detected = if_else(is.nan(mean_log2_detected), NA_real_, mean_log2_detected))
write_csv(per_group, file.path(out_dir, "enzyme_per_group_mval.csv"))
print(per_group, n = 18)

# Prediction check: Stk1 absent in stk (allow at most one stray detection).
stk_in_stk <- per_group |> filter(short_name == "Stk1", genotype == "stk")
check("PREDICTION Stk1 reduced in stk: <= 1 detection per medium, none with >= 2 peptides",
      all(stk_in_stk$n_detected <= 1) && all(stk_in_stk$n_with_2plus_pep == 0),
      paste0("(detections TSB/SNpulse/SN: ", paste(stk_in_stk$n_detected, collapse = "/"), ")"))

# Values that rest on a single peptide: these are the "unexpected" wild-type-like values.
single_pep <- per_sample |> filter(genotype != "anc", n_peptides == 1)
message("Mutant values resting on ONE peptide:")
print(single_pep |> select(short_name, sample, normalized_abundance, n_peptides))

# ---- 3. Figures 1a + 1b -----------------------------------------------------
anc_range <- per_sample |> filter(genotype == "anc") |>
  group_by(short_name) |> summarise(n_pep_anc = paste0(min(n_peptides), "-", max(n_peptides)))
facet_lab <- setNames(paste0(anc_range$short_name, " (anc: ", anc_range$n_pep_anc, " peptides/sample)"),
                      anc_range$short_name)

p_abund <- per_sample |>
  filter(!is.na(normalized_abundance)) |>
  ggplot(aes(group, normalized_abundance, colour = genotype)) +
  geom_point(aes(shape = n_peptides >= 2), size = 2.4,
             position = position_jitter(width = 0.15, seed = 1)) +
  geom_text(aes(label = n_peptides), size = 2.6, colour = "grey20",
            position = position_jitter(width = 0.15, seed = 1), vjust = -0.9) +
  # missing replicates written under the axis, so they are not invisible
  geom_text(data = per_group |> filter(n_detected < 3),
            aes(group, -Inf, label = paste0(n_detected, "/3")),
            inherit.aes = FALSE, vjust = -0.5, size = 2.8, colour = "firebrick") +
  facet_wrap(~ short_name, labeller = as_labeller(facet_lab)) +
  scale_colour_manual(values = GENOTYPE_COLOURS) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = ">= 2 peptides", `FALSE` = "1 peptide"), name = NULL) +
  labs(x = NULL, y = "normalised log2 abundance",
       title = "Stp1 and Stk1 per sample (number = peptides quantified in that sample)",
       subtitle = "red n/3 = samples in which the protein was detected; open circles rest on a single peptide") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p_abund, "fig_enzyme_abundance_mval.png", w = 10, h = 4.8)

p_pep <- per_sample |>
  mutate(detected = if_else(is.na(n_peptides), "not detected", "detected"),
         n_plot = replace_na(n_peptides, 0)) |>      # count scale only, never logged
  ggplot(aes(group, n_plot, colour = genotype, shape = detected)) +
  geom_point(size = 2.4, position = position_jitter(width = 0.15, height = 0, seed = 1)) +
  facet_wrap(~ short_name, scales = "free_y") +
  scale_colour_manual(values = GENOTYPE_COLOURS) +
  scale_shape_manual(values = c(detected = 16, `not detected` = 4), name = NULL) +
  labs(x = NULL, y = "peptides quantified in the sample",
       title = "Peptide evidence behind each protein value",
       subtitle = "x = protein not detected (no peptide quantified); data set total: Stp1 20, Stk1 37 peptides") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p_pep, "fig_enzyme_peptides_mval.png", w = 10, h = 4.5)

# ---- 4. Contrasts: supplied screen and our formal test, looked up only --------
dea <- read_dea("total") |> keep_je2()
treat <- read_tsv("data/derived/limma_treat_total_mval.tsv.gz", show_col_types = FALSE)
want <- c(paste0(rep(c("stk", "stp"), each = 3), "_vs_anc_at_", MEDIA),
          paste0("stp_vs_stk_at_", MEDIA))

contrasts_tbl <- dea |>
  inner_join(enzymes, by = "locus_tag") |>
  filter(contrast %in% want) |>
  transmute(short_name, contrast, supplied_diff = round(diff, 2),
            supplied_se = round(std.error, 2), supplied_FDR = signif(FDR, 2),
            supplied_estimate_type = estimate_type, supplied_model = modelName) |>
  left_join(treat |> inner_join(enzymes, by = "locus_tag") |>
              transmute(short_name, contrast, treat_diff = round(diff, 2),
                        treat_se = round(std.error, 2), treat_FDR = signif(FDR, 2),
                        treat_tested = tested, treat_min_obs = min_obs_in_groups),
            by = c("short_name", "contrast")) |>
  # direction in words (CLAUDE.md rule 3): diff < 0 = higher in the SECOND group
  mutate(direction = case_when(supplied_diff < 0 ~ paste("lower in", str_extract(contrast, "^[a-z]+")),
                               supplied_diff > 0 ~ paste("higher in", str_extract(contrast, "^[a-z]+")),
                               TRUE ~ "no difference")) |>
  arrange(short_name, contrast)
check("contrast table: 2 proteins x 9 contrasts", nrow(contrasts_tbl) == 18)
write_csv(contrasts_tbl, file.path(out_dir, "enzyme_contrasts_mval.csv"))
print(contrasts_tbl |> select(short_name, contrast, supplied_diff, supplied_FDR,
                              supplied_estimate_type, treat_diff, treat_FDR, treat_tested), n = 18)
message(explain_contrast("stp_vs_anc_at_TSB"))

# ---- 5. Stk1 phosphosite map (Stp1 has no phosphosites in our tables) ---------
sites <- bind_rows(read_sites("STYH"), read_sites("STY")) |>
  filter(locus_tag %in% enzymes$locus_tag)
check("no phosphosites on Stp1 in either search",
      !any(sites$locus_tag == "B7H15_06350"))

site_tbl <- sites |>
  group_by(search, site = paste0(residue, position), position, best_loc_prob, genotype) |>
  summarise(n = sum(!is.na(intensity)), .groups = "drop") |>
  pivot_wider(names_from = genotype, values_from = n, names_prefix = "detected_of9_") |>
  arrange(search, position)
write_csv(site_tbl, file.path(out_dir, "stk1_sites_mval.csv"))
print(site_tbl, n = 30)

p_map <- sites |>
  filter(!is.na(intensity)) |>
  mutate(sample = factor(sample, levels = rev(samples$sample[order(samples$group, samples$replicate)]))) |>
  ggplot(aes(position, sample, fill = log2_intensity)) +
  geom_tile(width = 4, height = 0.85) +
  geom_vline(xintercept = insertions$codon[insertions$short_name == "Stk1"],
             colour = "firebrick", linetype = 2) +
  annotate("text", x = 46, y = Inf, label = "stk insertion\n(codon 46)", vjust = 1.1,
           hjust = -0.05, size = 2.8, colour = "firebrick") +
  facet_wrap(~ search) +
  scale_y_discrete(drop = FALSE) +   # keep samples with NO Stk1 site visible as empty rows
  scale_x_continuous(limits = c(0, 664), breaks = c(0, 46, 100, 200, 300, 400, 500, 600, 664)) +
  scale_fill_viridis_c(name = "log2 intensity") +
  labs(x = "position in Stk1 (664 aa)", y = NULL,
       title = "Stk1 phosphosites per sample (blank = not detected)",
       subtitle = "Stp1 has no phosphosites in either search; the stp insertion (Stp1 codon 189) is not on this axis") +
  theme_bw(base_size = 9)
save_fig(p_map, "fig_stk1_phosphosite_map_mval.png", w = 10, h = 7)

# ---- 6. The HIQDSVPNVTtDVR peptide: is it real evidence for a truncated Stk1? --
fa <- readLines("data/raw/resources/JE2_proteins.fasta")
is_hdr <- str_starts(fa, ">")
seqs <- tapply(fa[!is_hdr], cumsum(is_hdr)[!is_hdr], paste, collapse = "")
names(seqs) <- str_extract(fa[is_hdr], "^>\\S+")
n_fasta_hits <- function(p) sum(str_detect(seqs, fixed(p)))
check("peptide HIQDSVPNVTTDVR occurs in exactly one FASTA entry (Stk1)",
      n_fasta_hits("HIQDSVPNVTTDVR") == 1 &&
        str_detect(names(seqs)[str_detect(seqs, fixed("HIQDSVPNVTTDVR"))], "ARG45703"))

hiq <- sites |>
  filter(position %in% c(231, 240, 242), !is.na(intensity)) |>
  select(search, site_id, position, sample, genotype, medium, log2_intensity, loc_prob, reinjected) |>
  arrange(search, position, sample)
write_csv(hiq, file.path(out_dir, "hiq_peptide_check_mval.csv"))

# Compare: the activation-loop site T164 (clear wild-type signal) vs T231
cmp <- sites |> filter(search == "STYH", position %in% c(164, 231), !is.na(intensity)) |>
  group_by(site = paste0(residue, position), genotype) |>
  summarise(n_detected_of9 = n(), median_log2 = round(median(log2_intensity), 1),
            median_loc_prob = round(median(loc_prob, na.rm = TRUE), 2), .groups = "drop")
print(cmp)
hiq_stk <- hiq |> filter(search == "STYH", position == 231, genotype == "stk")
info <- sprintf("[INFO] T231 (HIQ peptide, STYH) in stk: %d of 9 samples (%s)",
                nrow(hiq_stk), paste(hiq_stk$sample, collapse = ", "))
checks <- c(checks, info); message(info)

# ---- 7. Evidence table: one rule-based call per enzyme x strain --------------
# Rule (fixed here, before reading the result):
#   reference: anc has the protein in 3/3 samples with >= 2 peptides in every medium.
#   "supported reduction"  : mutant has <= 1 sample with >= 2 peptides in EVERY medium
#   "no reduction"         : mutant has 3/3 samples with >= 2 peptides in every medium
#                            and |supplied diff| < 1 in every medium
#   "inconsistent pattern" : anything else
# Single-peptide values are reported, not removed.
calls <- per_group |>
  filter(genotype != "anc") |>
  left_join(contrasts_tbl |> filter(str_detect(contrast, "_vs_anc_")) |>
              mutate(genotype = str_extract(contrast, "^[a-z]+"),
                     medium = str_remove(contrast, ".*_at_")) |>
              select(short_name, genotype, medium, supplied_diff, supplied_estimate_type),
            by = c("short_name", "genotype", "medium")) |>
  group_by(short_name, genotype) |>
  summarise(detections_TSB_SNpulse_SN = paste(n_detected, collapse = "/"),
            samples_2plus_peptides     = paste(n_with_2plus_pep, collapse = "/"),
            supplied_diff_TSB_SNpulse_SN = paste(supplied_diff, collapse = " / "),
            any_lod_imputed = any(supplied_estimate_type == "lod_imputed"),
            call = case_when(all(n_with_2plus_pep <= 1) ~ "supported reduction",
                             all(n_with_2plus_pep == 3) & all(abs(supplied_diff) < 1) ~ "no reduction",
                             TRUE ~ "inconsistent pattern"),
            .groups = "drop")
anc_ok <- per_group |> filter(genotype == "anc") |> pull(n_with_2plus_pep) |> (\(x) all(x == 3))()
check("reference: anc has both enzymes in 3/3 samples with >= 2 peptides in every medium", anc_ok)
write_csv(calls, file.path(out_dir, "evidence_table_mval.csv"))
print(calls)

writeLines(checks, file.path(out_dir, "validation_checks_mval.txt"))
message("Done. Outputs in ", out_dir)
