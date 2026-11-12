# Step 4 of our workflow (mval_medi/project_outlines/PROJECT_WORKFLOW_mval.md):
# examine the possible polar effect of the stp1::Tn insertion on stk1.
#
# Run from the repository root:   Rscript scripts/04_polar_effect_mval.R   (~20 s)
# Needs the Step 3 output (enzyme_per_sample_mval.csv): run scripts/03_* first.
#
# ---- Plan (coding rules, section 4) -----------------------------------------
# Aim:      Is the loss of Stk1 in `stp` consistent with a polar insertion?
#           Are other neighbouring proteins changed too?
# Input:    JE2_genome.gff (genes, strand, coordinates around stp1/stk1)
#           read_protein_metadata()  NCTC 8325 UniProt orthologs (for AureoWiki)
#           TRANSPOSON_INSERTIONS    insertion coordinates (R/constants.R, NTML map)
#           Step 3 per-sample table  Stk1 peptides per sample
#           read_abundances("total"), read_dea("total"), limma_treat_total_mval.tsv.gz
#           read_sites("STYH")       Stk1 phosphosites detected in anc AND stp
# Analysis: descriptive only: coordinates, counts, group means of observed values.
#           Neighbour changes are read from the supplied DEA and our treat test
#           (no new test, so no new multiple-testing correction is needed).
# Output:   mval_medi/results/step04_polar_effect/
#             fig_gene_context_mval.png       Fig 1d: gene map with insertions
#             fig_neighbours_mval.png         per-sample values of the neighbours
#             fig_stk1_peptides_by_strain_mval.png  Stk1 peptides anc vs stk vs stp
#             gene_context_mval.csv           genes, strand, gaps, ortholog, measured?
#             neighbour_contrasts_mval.csv    mutant vs anc for every neighbour
#             stk1_peptides_by_strain_mval.csv
#             stk1_site_drop_stp_mval.csv     Stk1 sites: stp vs anc (observed only)
#             polar_checks_mval.txt
# Check:    stp1/stk1 overlap by 4 bp; all genes in the window on one strand as
#           in AureoWiki; neighbour proteins detected 3/3 everywhere.
#
# Source for the transcription unit (checked 12 Nov 2026, NCTC 8325, NOT JE2):
#   AureoWiki SAOUHSC_01186 (stp1) and SAOUHSC_01187 (pknB/stk1), section Operon:
#   SAOUHSC_01182 > 01183 > 01184 > 01185 > 01186 > 01187 (MicrobesOnline, and a
#   SigA-promoter unit from Mäder et al. 2016, PLoS Genet, PMID 27035918).
#   https://aureowiki.med.uni-greifswald.de/SAOUHSC_01186
#   https://aureowiki.med.uni-greifswald.de/SAOUHSC_01187
#   UniProt Q2FZ65 = SAOUHSC_01186, Q2FZ64 = SAOUHSC_01187 (rest.uniprot.org).
#   We map SAOUHSC_01182-01185 to JE2 genes by gene order and UniProt ortholog;
#   AureoWiki gives only locus tags for those four, so this mapping is an assumption.

suppressMessages(source("R/setup.R"))

out_dir <- "mval_medi/results/step04_polar_effect"
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

step3 <- "mval_medi/results/step03_validation/enzyme_per_sample_mval.csv"
if (!file.exists(step3)) stop("Run scripts/03_mutant_validation_mval.R first.")

STP1 <- "B7H15_06350"; STK1 <- "B7H15_06355"
window <- c(1226000, 1238000)                       # ~6 kb either side of stp1-stk1
# Genes that AureoWiki lists in the same transcription unit (see header).
OPERON <- sprintf("B7H15_%05d", seq(6330, 6355, 5))

# ---- 1. Gene context from the GFF ---------------------------------------------
gff <- read_tsv("data/raw/resources/JE2_genome.gff", comment = "#", col_names = FALSE,
                col_types = "cccddcccc", progress = FALSE) |>
  filter(X3 == "gene", X4 >= window[1], X5 <= window[2]) |>
  transmute(locus_tag = str_match(X9, "locus_tag=([^;]+)")[, 2],
            start = X4, end = X5, strand = X7)

meta <- read_protein_metadata() |>
  select(locus_tag, uniprot_acc, ortholog_quality, gene_symbol_best, up_protein_name)
measured <- read_abundances("total") |> keep_je2() |>
  filter(locus_tag %in% gff$locus_tag) |>
  group_by(locus_tag) |> summarise(n_samples_detected = sum(!is.na(normalized_abundance)))

context <- gff |>
  arrange(start) |>
  mutate(gap_to_next_bp = lead(start) - end - 1,    # negative = genes overlap
         in_aureowiki_unit = locus_tag %in% OPERON) |>
  left_join(read_annotation() |> distinct(locus_tag, product), by = "locus_tag") |>
  left_join(meta, by = "locus_tag") |>
  left_join(measured, by = "locus_tag") |>
  mutate(n_samples_detected = replace_na(n_samples_detected, 0L),
         label = case_when(locus_tag == STP1 ~ "stp1",
                           locus_tag == STK1 ~ "stk1",
                           !is.na(gene_symbol_best) ~ gene_symbol_best,
                           TRUE ~ str_remove(locus_tag, "B7H15_")))
check("one row per gene, no duplicates after the joins",
      !anyDuplicated(context$locus_tag))
check("stp1 and stk1 overlap by 4 bp",
      context$gap_to_next_bp[context$locus_tag == STP1] == -4)
check("all genes of the AureoWiki unit are on the + strand in JE2",
      all(context$strand[context$in_aureowiki_unit] == "+"))
check("orthologs of stp1/stk1 are Q2FZ65/Q2FZ64 (the AureoWiki entries)",
      identical(context$uniprot_acc[match(c(STP1, STK1), context$locus_tag)], c("Q2FZ65", "Q2FZ64")))
next_after_stk1 <- context |> filter(start > context$end[context$locus_tag == STK1]) |> slice(1)
check("gap between stk1 and the next gene (bp)", TRUE,
      paste0("= ", context$gap_to_next_bp[context$locus_tag == STK1], " bp to ", next_after_stk1$locus_tag,
             " (not in the AureoWiki unit)"))
write_csv(context, file.path(out_dir, "gene_context_mval.csv"))
print(context |> select(locus_tag, label, start, end, strand, gap_to_next_bp,
                        in_aureowiki_unit, uniprot_acc, n_samples_detected))

# ---- 2. Gene map with both insertion sites -------------------------------------
ins <- TRANSPOSON_INSERTIONS |>
  mutate(codon = (insertion_pos - cds_start + 1 + 2) %/% 3,
         label = paste0(ntml_id, " (", if_else(short_name == "Stp1", "stp", "stk"),
                        ")\ncodon ", codon),
         hjust = if_else(short_name == "Stp1", 1.05, -0.05))   # the two sites are 310 bp apart: push labels apart
p_map <- ggplot(context) +
  annotate("rect", xmin = min(context$start[context$in_aureowiki_unit]),
           xmax = max(context$end[context$in_aureowiki_unit]),
           ymin = 0.55, ymax = 1.45, fill = "grey92") +
  annotate("text", x = min(context$start[context$in_aureowiki_unit]), y = 1.42, hjust = 0, vjust = 1,
           size = 2.6, colour = "grey35",
           label = "transcription unit in NCTC 8325 (AureoWiki / Mäder et al. 2016)") +
  geom_segment(aes(x = if_else(strand == "+", start, end), xend = if_else(strand == "+", end, start),
                   y = 1, yend = 1,
                   colour = case_when(locus_tag %in% c(STP1, STK1) ~ "stp1 / stk1",
                                      n_samples_detected == 0 ~ "not measured",
                                      TRUE ~ "measured")),
               linewidth = 5, arrow = arrow(length = unit(2.5, "mm"), type = "closed")) +
  geom_text(aes(x = (start + end) / 2, y = 0.8, label = label), size = 2.8) +
  geom_point(data = ins, aes(insertion_pos, 1.18), shape = 25, size = 3.5,
             fill = "firebrick", colour = "firebrick") +
  geom_text(data = ins, aes(insertion_pos, 1.3, label = label, hjust = hjust), size = 2.6,
            colour = "firebrick", lineheight = 0.9) +
  scale_colour_manual(values = c("stp1 / stk1" = "#CC79A7", measured = "grey55",
                                 "not measured" = "grey85"), name = NULL) +
  scale_x_continuous(labels = scales::label_comma(), name = "JE2 genome position (CP020619.1)") +
  scale_y_continuous(limits = c(0.55, 1.45), breaks = NULL, name = NULL) +
  labs(title = "Gene context of stp1 and stk1 in JE2",
       subtitle = "stp1 and stk1 overlap by 4 bp; triangles = bursa aurealis insertions (3,237 bp)") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
save_fig(p_map, "fig_gene_context_mval.png", w = 11, h = 3.6)

# ---- 3. Neighbouring proteins: changed in the mutants? -------------------------
neighbours <- context |> filter(n_samples_detected > 0, !locus_tag %in% c(STP1, STK1)) |> pull(locus_tag)
check("every measured neighbour is detected in all 27 samples",
      all(context$n_samples_detected[context$locus_tag %in% neighbours] == 27))

abund <- read_abundances("total") |> keep_je2() |> filter(locus_tag %in% neighbours) |>
  left_join(context |> select(locus_tag, label), by = "locus_tag") |>
  mutate(tag = str_remove(locus_tag, "B7H15_"),
         gene = if_else(label == tag, tag, paste0(label, " (", tag, ")")),
         gene = factor(gene, levels = unique(gene[order(locus_tag)])))
p_nb <- ggplot(abund, aes(group, normalized_abundance, colour = genotype)) +
  geom_point(size = 1.6, position = position_jitter(width = 0.15, seed = 1)) +
  facet_wrap(~ gene, nrow = 2) +
  scale_colour_manual(values = GENOTYPE_COLOURS) +
  labs(x = NULL, y = "normalised log2 abundance",
       title = "Neighbouring proteins per sample (all detected in 27/27 samples)") +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
save_fig(p_nb, "fig_neighbours_mval.png", w = 12, h = 5.5)

want <- c(paste0(rep(c("stk", "stp"), each = 3), "_vs_anc_at_", MEDIA))
treat <- read_tsv("data/derived/limma_treat_total_mval.tsv.gz", show_col_types = FALSE)
nb_tbl <- read_dea("total") |> keep_je2() |>
  filter(locus_tag %in% neighbours, contrast %in% want) |>
  select(locus_tag, contrast, supplied_diff = diff, supplied_FDR = FDR, estimate_type) |>
  left_join(treat |> select(locus_tag, contrast, treat_FDR = FDR), by = c("locus_tag", "contrast")) |>
  left_join(context |> select(locus_tag, label, in_aureowiki_unit), by = "locus_tag") |>
  mutate(across(c(supplied_diff), \(v) round(v, 2)), across(ends_with("FDR"), \(v) signif(v, 2)),
         formal_hit = !is.na(treat_FDR) & treat_FDR < MVAL_FDR) |>
  arrange(locus_tag, contrast)
write_csv(nb_tbl, file.path(out_dir, "neighbour_contrasts_mval.csv"))
check("neighbours: largest |supplied diff| mutant vs anc", TRUE,
      paste0("= ", max(abs(nb_tbl$supplied_diff)), " (",
             nb_tbl$label[which.max(abs(nb_tbl$supplied_diff))], ", ",
             nb_tbl$contrast[which.max(abs(nb_tbl$supplied_diff))], ")"))
check("neighbours: formal hits (treat lfc = 1, BH-FDR < 0.05)", TRUE,
      paste0("= ", sum(nb_tbl$formal_hit), " of ", nrow(nb_tbl), " protein x contrast rows"))

# ---- 4. Polar effect at peptide level: Stk1 peptides in anc vs stk vs stp ------
pep <- read_csv(step3, show_col_types = FALSE) |>
  filter(short_name == "Stk1") |>
  mutate(genotype = factor(genotype, levels = GENOTYPES), medium = factor(medium, levels = MEDIA))
pep_tbl <- pep |>
  group_by(genotype, medium) |>
  summarise(samples_detected = sum(!is.na(n_peptides)),
            peptides_per_sample = paste(replace_na(as.character(n_peptides), "-"), collapse = "/"),
            median_peptides_NA_as_0 = median(replace_na(n_peptides, 0)),
            .groups = "drop")
write_csv(pep_tbl, file.path(out_dir, "stk1_peptides_by_strain_mval.csv"))
print(pep_tbl)

p_pep <- pep |>
  mutate(n_plot = replace_na(n_peptides, 0),
         detected = if_else(is.na(n_peptides), "not detected", "detected")) |>
  ggplot(aes(genotype, n_plot, colour = genotype, shape = detected)) +
  geom_point(size = 2.6, position = position_jitter(width = 0.12, height = 0, seed = 2)) +
  facet_wrap(~ medium) +
  scale_colour_manual(values = GENOTYPE_COLOURS, guide = "none") +
  scale_shape_manual(values = c(detected = 16, `not detected` = 4), name = NULL) +
  labs(x = NULL, y = "Stk1 peptides quantified per sample",
       title = "Stk1 in stp drops to the stk level (polar-effect prediction)",
       subtitle = "x = not detected (0 peptides); count scale, never logged") +
  theme_bw(base_size = 10)
save_fig(p_pep, "fig_stk1_peptides_by_strain_mval.png", w = 8, h = 4)

# ---- 5. Polar effect at phosphosite level --------------------------------------
# Only sites seen in anc AND stp in the same medium; observed values only, no imputation.
site_drop <- read_sites("STYH") |>
  filter(locus_tag == STK1, genotype %in% c("anc", "stp"), !is.na(log2_intensity)) |>
  group_by(site = paste0(residue, position), position, medium, genotype) |>
  summarise(n = n(), mean_log2 = mean(log2_intensity), .groups = "drop") |>
  pivot_wider(names_from = genotype, values_from = c(n, mean_log2)) |>
  filter(!is.na(n_anc), !is.na(n_stp)) |>
  mutate(diff_stp_minus_anc = round(mean_log2_stp - mean_log2_anc, 2)) |>
  arrange(position, medium)
write_csv(site_drop, file.path(out_dir, "stk1_site_drop_stp_mval.csv"))
print(site_drop, n = 40)
loop <- site_drop |> filter(position %in% 161:172, n_anc >= 2, n_stp >= 2)
check("activation-loop sites (161-172, >= 2 obs per group): median stp - anc (log2)", TRUE,
      paste0("= ", round(median(loop$diff_stp_minus_anc), 2), " (", nrow(loop),
             " site x medium rows; negative = lower in stp)"))

writeLines(checks, file.path(out_dir, "polar_checks_mval.txt"))
message("Done. Outputs in ", out_dir)
