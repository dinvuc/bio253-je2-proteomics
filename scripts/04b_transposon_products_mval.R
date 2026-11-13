# Step 4b (PROJECT_WORKFLOW_mval.md): transposon products as a positive control.
#
# Run from the repository root:   Rscript scripts/04b_transposon_products_mval.R   (~30 s)
#
# ---- Plan (coding rules, section 4) -----------------------------------------
# Aim:      Is the bursa aurealis transposon present AND expressed in `stk` and
#           `stp`, and absent from `anc`? Prediction written in reference/README.md
#           before the data arrived: ErmB detected in the 18 mutant samples,
#           absent from the 9 anc samples.
# Input:    data/raw/FP24_transposonDB/ (local, git-ignored). A new FragPipe 24.0
#           search of the 27 total-proteome runs by Jonas (12 Nov 2026) against
#           p41881_db2_SA_JE2_wTransposon_20261111 = JE2 proteins + the proteins
#           of pBursa (ermB, gfp, cat, repC, orf) and of the helper plasmid pFA545
#           (tnp, bla, repC, tetB, xylR). NOT part of the checksummed course bundle.
#             combined_protein.tsv   per-sample spectral counts, MaxLFQ intensity
#             combined_peptide.tsv   per-sample peptide spectral counts, Start/End
#             combined.prot.xml      ProteinProphet probability for EVERY database
#                                    protein, incl. those that did not pass
#           read_samples()           sample sheet (join by sample name)
# Analysis: descriptive only: counts per sample. No test: the question is
#           presence/absence (0 vs > 0 spectra), not a fold change.
# Output:   mval_medi/results/step04b_transposon/
#             transposon_db_proteins_mval.csv  all 10 Tn/plasmid proteins: probability, peptides
#             ermb_per_sample_mval.csv         ErmB per sample: spectra, peptides, intensity
#             ermb_peptides_mval.csv           the ErmB peptides: position, samples seen in
#             fig_ermb_per_sample_mval.png     Fig 1e: ErmB spectra + intensity per sample
#             fig_ermb_peptide_map_mval.png    ErmB peptides along the sequence x sample
#             figure_guide_mval.md, transposon_checks_mval.txt
# Check:    27 samples, 3 per group after the join; ErmB 0 spectra in every anc
#           sample and > 0 in every stk/stp sample (= the prediction).

suppressMessages(source("R/setup.R"))
src <- "data/raw/FP24_transposonDB"
out <- "mval_medi/results/step04b_transposon"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
checks <- character()
check <- function(label, ok) {
  line <- sprintf("[%s] %s", if (isTRUE(ok)) "PASS" else "FAIL", label)
  message(line); checks <<- c(checks, line)
  invisible(ok)
}
if (!dir.exists(src)) stop("Missing ", src, ": the transposon search from Jonas is local only.")

samples <- read_samples()

# FragPipe column names look like "<raw file> Spectral Count". Split them into
# sample (canonical name via parse_sample_name()) and measure, then join BY NAME.
to_long <- function(tbl, id_cols) {
  tbl |>
    select(all_of(id_cols), matches("_input (Spectral Count|Intensity|MaxLFQ Intensity|Unique Spectral Count)$")) |>
    pivot_longer(-all_of(id_cols), names_to = "column", values_to = "value") |>
    mutate(raw_file = str_remove(column, " (Spectral Count|Intensity|MaxLFQ Intensity|Unique Spectral Count)$"),
           measure  = str_extract(column, "(Unique Spectral Count|MaxLFQ Intensity|Spectral Count|Intensity)$"),
           sample   = parse_sample_name(raw_file)$sample) |>
    select(-column)
}

# ---- 1. every transposon/plasmid protein that was in the database ----------
# combined_protein.tsv lists only proteins that passed the protein FDR, so read
# ProteinProphet's view of all of them from the prot.xml (one <protein ...> line each).
xml  <- readLines(file.path(src, "combined.prot.xml"))
prot <- str_subset(xml, '<protein protein_name="tr\\|Trps_')
# A forward protein with no peptide match at all is not listed, but its reversed
# decoy may be. Every decoy is made from a forward database entry, so the union
# of forward and decoy names is the list of searched Tn/plasmid proteins.
searched <- unique(str_match(str_subset(xml, 'protein_name="(rev_)?tr\\|Trps_'),
                             'protein_name="(?:rev_)?tr\\|(Trps_[^|]+)')[, 2])
db_proteins <- tibble(
  entry       = str_match(prot, 'protein_name="tr\\|(Trps_[^|]+)')[, 2],
  probability = as.numeric(str_match(prot, ' probability="([0-9.]+)"')[, 2]),
  n_peptides_total = as.integer(str_match(prot, 'total_number_peptides="([0-9]+)"')[, 2]),
  peptide_sequences = str_match(prot, 'unique_stripped_peptides="([^"]+)"')[, 2])
# searched but not a single peptide match: probability 0, 0 peptides
db_proteins <- db_proteins |>
  bind_rows(tibble(entry = setdiff(searched, db_proteins$entry), probability = 0, n_peptides_total = 0L)) |>
  mutate(plasmid = str_match(entry, "Trps_(pBursa|pFA545)_")[, 2],
         gene = str_remove(entry, "Trps_(pBursa|pFA545)_"),
         role = if_else(plasmid == "pBursa", "transposon delivery plasmid (bursa aurealis)",
                        "helper plasmid (transposase source)"))

protein_tbl <- read_tsv(file.path(src, "combined_protein.tsv"), show_col_types = FALSE, guess_max = 1e5)
passed <- protein_tbl |> filter(str_detect(`Protein ID`, "^Trps_")) |> pull(`Protein ID`)
db_proteins <- db_proteins |>
  mutate(in_combined_protein_table = entry %in% passed) |>
  select(plasmid, gene, entry, role, probability, n_peptides_total, in_combined_protein_table, peptide_sequences) |>
  arrange(desc(probability))
write_csv(db_proteins, file.path(out, "transposon_db_proteins_mval.csv"))
print(select(db_proteins, -peptide_sequences, -role))
check("all 10 transposon/plasmid proteins were in the searched database", nrow(db_proteins) == 10)
check("GFP was searched (in the prot.xml), so 'not found' means searched and not detected",
      "gfp" %in% db_proteins$gene[db_proteins$plasmid == "pBursa"])

# ---- 2. ErmB per sample ------------------------------------------------------
ermb <- protein_tbl |> filter(`Protein ID` == "Trps_pBursa_ermB") |>
  to_long("Protein ID") |>
  select(sample, measure, value) |>
  pivot_wider(names_from = measure, values_from = value) |>
  rename(spectral_count = `Spectral Count`, unique_spectral_count = `Unique Spectral Count`,
         intensity = Intensity, maxlfq_intensity = `MaxLFQ Intensity`)

peptide_tbl <- read_tsv(file.path(src, "combined_peptide.tsv"), show_col_types = FALSE, guess_max = 1e5) |>
  filter(`Protein ID` == "Trps_pBursa_ermB")
pep_long <- peptide_tbl |>
  to_long(c("Peptide Sequence", "Start", "End")) |>
  filter(measure == "Spectral Count") |>
  rename(peptide = `Peptide Sequence`, start = Start, end = End, spectra = value)

ermb_per_sample <- samples |>
  select(sample, genotype, medium, replicate, group) |>
  left_join(ermb, by = "sample") |>
  left_join(pep_long |> group_by(sample) |> summarise(n_peptides = sum(spectra > 0)), by = "sample") |>
  # 0 intensity in FragPipe = not quantified -> NA (never 0 on a log scale)
  mutate(detected = spectral_count > 0,
         log2_maxlfq = if_else(maxlfq_intensity > 0, log2(maxlfq_intensity), NA_real_))
stopifnot(!anyDuplicated(ermb_per_sample$sample))
check("27 samples joined by name, 3 per group, every sample found in the search",
      nrow(ermb_per_sample) == 27 && all(table(ermb_per_sample$group) == 3) &&
        !anyNA(ermb_per_sample$spectral_count))
write_csv(ermb_per_sample, file.path(out, "ermb_per_sample_mval.csv"))

by_geno <- ermb_per_sample |> group_by(genotype) |>
  summarise(detected = sum(detected), n = n(), min_spectra = min(spectral_count),
            max_spectra = max(spectral_count), min_peptides = min(n_peptides), max_peptides = max(n_peptides))
print(by_geno)
check("prediction: ErmB 0 spectra in all 9 anc samples",
      all(ermb_per_sample$spectral_count[ermb_per_sample$genotype == "anc"] == 0))
check("prediction: ErmB > 0 spectra in all 18 stk and stp samples",
      all(ermb_per_sample$spectral_count[ermb_per_sample$genotype != "anc"] > 0))

ermb_peptides <- pep_long |>
  left_join(select(samples, sample, genotype), by = "sample") |>
  group_by(peptide, start, end) |>
  summarise(samples_anc = sum(spectra > 0 & genotype == "anc"),
            samples_stk = sum(spectra > 0 & genotype == "stk"),
            samples_stp = sum(spectra > 0 & genotype == "stp"),
            total_spectra = sum(spectra), .groups = "drop") |>
  arrange(start)
write_csv(ermb_peptides, file.path(out, "ermb_peptides_mval.csv"))
covered <- unique(unlist(purrr::map2(ermb_peptides$start, ermb_peptides$end, seq)))
message(sprintf("ErmB sequence coverage: %d of 245 residues (%.1f %%)", length(covered), 100 * length(covered) / 245))

# ---- 3. figures --------------------------------------------------------------
plot_df <- ermb_per_sample
p1 <- ggplot(plot_df, aes(medium, spectral_count, colour = genotype)) +
  geom_point(aes(shape = detected), position = position_jitter(width = 0.15, height = 0, seed = 1), size = 2.4) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4), labels = c(`TRUE` = "detected", `FALSE` = "not detected (0)")) +
  scale_colour_manual(values = GENOTYPE_COLOURS) +
  facet_wrap(~genotype) +
  labs(x = NULL, y = "ErmB spectral count (count scale)", shape = NULL, colour = NULL,
       title = "a  ErmB spectra per sample")
p2 <- ggplot(filter(plot_df, !is.na(log2_maxlfq)), aes(medium, log2_maxlfq, colour = genotype)) +
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 1), size = 2.4) +
  scale_colour_manual(values = GENOTYPE_COLOURS, guide = "none") +   # legend already in panel a
  scale_x_discrete(drop = FALSE) +
  facet_wrap(~genotype, drop = FALSE) +
  labs(x = NULL, y = "log2 MaxLFQ intensity", colour = NULL,
       title = "b  ErmB intensity (samples without signal are not drawn)")
fig <- patchwork::wrap_plots(p1, p2, ncol = 1, guides = "collect") +
  patchwork::plot_annotation(
    title = "Transposon marker ErmB: a positive control for the insertions",
    subtitle = "FragPipe 24 search incl. transposon proteins (Jonas, 12 Nov 2026); 27 samples, n = 3 per group")
ggsave(file.path(out, "fig_ermb_per_sample_mval.png"), fig, width = 8, height = 7, dpi = 300)

map_df <- pep_long |>
  left_join(select(samples, sample, genotype, group), by = "sample") |>
  mutate(sample = factor(sample, rev(samples$sample)))
p3 <- ggplot(filter(map_df, spectra > 0)) +
  geom_segment(aes(x = start, xend = end, y = sample, yend = sample, colour = genotype), linewidth = 2.2) +
  scale_colour_manual(values = GENOTYPE_COLOURS, drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  coord_cartesian(xlim = c(1, 245)) +
  labs(x = "ErmB position (aa, 1-245)", y = NULL, colour = NULL,
       title = "ErmB peptides seen per sample",
       subtitle = sprintf("%d peptides, %.0f %% sequence coverage; empty row = no ErmB peptide",
                          nrow(ermb_peptides), 100 * length(covered) / 245)) +
  theme(axis.text.y = element_text(size = 7))
ggsave(file.path(out, "fig_ermb_peptide_map_mval.png"), p3, width = 8, height = 6.5, dpi = 300)

writeLines(c(paste("Step 4b transposon products, run", Sys.Date()),
             paste("Source:", src, "(FragPipe 24.0, db p41881_db2_SA_JE2_wTransposon_20261111)"),
             checks), file.path(out, "transposon_checks_mval.txt"))
if (any(grepl("^\\[FAIL", checks))) stop("At least one check FAILED, see ", out)
message("done, written to ", out)
