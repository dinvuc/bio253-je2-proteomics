# Step 6 add-on: the Step 6 interaction hits per strain, in genome order, with
# their strand and their distance to the transposon insertions.
# Same layout as step04_polar_effect/fig_neighbours_by_genotype_mval.png, but the
# genes are the hits instead of the stp1/stk1 neighbours.
#
# Run from the repository root:   Rscript scripts/06b_hit_genome_context_mval.R   (~10 s)
# Needs scripts/05_formal_test_mval.R and scripts/06_altered_responses_mval.R first.
#
# ---- Plan (coding rules, section 4) -----------------------------------------
# Aim:      Show the proteins that are hits in fig_response_vs_anc_mval.png and
#           fig_response_stk_vs_stp_mval.png one strain at a time, ordered along
#           the chromosome, and say where they sit relative to the transposon:
#           strand and distance from the gene's START CODON to the insertion.
# Input:    data/derived/limma_treat_total_mval.tsv.gz  which proteins are hits
#           data/raw/resources/JE2_genome.gff          start, end, strand, genome length
#           TRANSPOSON_INSERTIONS (R/constants.R)        insertion coordinates (NTML map)
#           read_abundances("total")                     per-sample values (not changed)
#           read_protein_metadata()                      gene symbols (ortholog_quality good)
# Analysis: none (descriptive). Hits = formal hits (treat lfc = 1, BH-FDR < 0.05) of
#           the three TSB-vs-SN interactions. Exploratory proteins (0.05 <= FDR < 0.1)
#           are noted; they add no protein here (both are formal hits elsewhere).
#           Distance: start codon = `start` on the + strand, `end` on the - strand.
#           Genome is circular -> shortest way round, signed:
#             + = start codon at a HIGHER genome coordinate than the insertion,
#             - = at a LOWER coordinate.
# Output:   mval_medi/results/step06_responses/
#             fig_hits_by_genotype_mval.png   one panel per strain, TSB and SN
#             hits_genome_context_mval.csv    strand, coordinates, distances, hit in
# Check:    27 samples / 3 per group joined by name; every hit found once in the GFF;
#           hit set equals the Step 6 candidate table.

suppressMessages(source("R/setup.R"))
out <- "mval_medi/results/step06_responses"

INTER <- c(stk_vs_anc = "interaction_stk_vs_anc_at_TSB_vs_SN",
           stp_vs_anc = "interaction_stp_vs_anc_at_TSB_vs_SN",
           stp_vs_stk = "interaction_stp_vs_stk_at_TSB_vs_SN")
FDR_EXPLORATORY <- 0.1   # display only, chosen after results (see Step 6 figure guide)

# ---- 1. which proteins are hits, and in which comparison ---------------------
res <- read_tsv("data/derived/limma_treat_total_mval.tsv.gz", show_col_types = FALSE) |>
  filter(contrast %in% INTER, tested, FDR < FDR_EXPLORATORY) |>
  mutate(comparison = names(INTER)[match(contrast, INTER)],
         formal = FDR < MVAL_FDR)
hits <- res |> filter(formal) |> distinct(locus_tag) |> pull()
cand <- read_csv(file.path(out, "interaction_candidates_mval.csv"), show_col_types = FALSE)
stopifnot(setequal(hits, cand$locus_tag))
extra_expl <- setdiff(res$locus_tag, hits)
message("CHECK formal-hit proteins: ", length(hits),
        "; exploratory proteins that are not formal hits elsewhere: ", length(extra_expl))

hit_in <- res |> filter(locus_tag %in% hits) |>
  mutate(txt = paste0(str_replace(comparison, "_vs_", " vs "), if_else(formal, "", " (expl.)"))) |>
  arrange(locus_tag, !formal, comparison) |>
  group_by(locus_tag) |> summarise(hit_in = paste(txt, collapse = "; "), .groups = "drop")

# ---- 2. genome position, strand, distance to the insertions -------------------
gff_lines <- readLines("data/raw/resources/JE2_genome.gff", n = 20)
genome_len <- as.numeric(str_match(str_subset(gff_lines, "^##sequence-region"), "\\s(\\d+)$")[, 2])
stopifnot(genome_len == 2874399)
gff_all <- read_tsv("data/raw/resources/JE2_genome.gff", comment = "#", col_names = FALSE,
                    col_types = "cccddcccc", progress = FALSE) |>
  filter(X3 == "gene") |>
  transmute(locus_tag = str_match(X9, "locus_tag=([^;]+)")[, 2], start = X4, end = X5, strand = X7) |>
  arrange(start)
gff <- gff_all |> filter(locus_tag %in% hits)
stopifnot(nrow(gff) == length(hits), !anyDuplicated(gff$locus_tag))

ins <- setNames(TRANSPOSON_INSERTIONS$insertion_pos, TRANSPOSON_INSERTIONS$short_name)  # Stp1, Stk1
circ_dist <- function(pos, ref) {          # signed shortest distance on a circular genome
  d <- (pos - ref) %% genome_len
  if_else(d > genome_len / 2, d - genome_len, d)
}
meta <- read_protein_metadata() |> filter(ortholog_quality == "good") |>
  distinct(locus_tag, gene = gene_symbol_best)
ctx <- gff |>
  mutate(start_codon = if_else(strand == "+", start, end),
         dist_to_stk1_Tn_bp = circ_dist(start_codon, ins[["Stk1"]]),
         dist_to_stp1_Tn_bp = circ_dist(start_codon, ins[["Stp1"]]),
         # the shortest way runs through position 1 (where the numbering starts)
         # when its sign differs from the plain coordinate difference
         path_to_Tn = if_else(sign(dist_to_stk1_Tn_bp) != sign(start_codon - ins[["Stk1"]]),
                              "via origin (position 1)", "direct")) |>
  left_join(meta, by = "locus_tag") |>
  left_join(distinct(cand, locus_tag, product), by = "locus_tag") |>
  left_join(hit_in, by = "locus_tag") |>
  arrange(start) |>
  mutate(tag = str_remove(locus_tag, "B7H15_"),
         name = if_else(is.na(gene), tag, paste0(gene, " (", tag, ")")),
         strand_txt = if_else(strand == "+", "→ forward (+)", "← reverse (−)"),
         short_product = str_trunc(str_remove(product, " biosynthetic IlvA$"), 24),
         x_label = paste(name, short_product, strand_txt, str_replace_all(hit_in, "; ", "\n"), sep = "\n"),
         x_label = factor(x_label, levels = x_label))
write_csv(ctx |> select(locus_tag, gene, product, strand, start, end, start_codon,
                        dist_to_stk1_Tn_bp, dist_to_stp1_Tn_bp, path_to_Tn, hit_in),
          file.path(out, "hits_genome_context_mval.csv"))
print(ctx |> select(locus_tag, gene, strand, start_codon, dist_to_stk1_Tn_bp, dist_to_stp1_Tn_bp))

# ---- 2b. annotation context from the metadata table and the GFF ---------------
# Metadata (reference/JE2_protein_metadata.tsv.gz): NCTC 8325 UniProt ortholog and
# what UniProt says about it. Names come from this table, never from memory (pitfall 4).
meta_full <- read_protein_metadata() |>
  filter(locus_tag %in% hits) |>
  transmute(locus_tag, uniprot_acc, ortholog_quality, up_reviewed, aa_length,
            nctc8325_locus = str_match(up_kegg, "sao:(SAOUHSC_\\d+)")[, 2],
            up_protein_name, up_ec,
            # drop the evidence codes {ECO:...} for readability
            subcellular = str_squish(str_remove_all(up_subcellular, "\\{[^}]*\\}\\.?")),
            up_keywords = str_split(up_keywords, ";") |>
              purrr::map_chr(\(k) paste(setdiff(k, "Reference proteome"), collapse = "; ")),
            gff_note)
stopifnot(nrow(meta_full) == length(hits), all(meta_full$ortholog_quality == "good"))

# GFF neighbours: 2 genes on each side, with strand and gap. A run of genes on the
# same strand with short gaps (< 150 bp) is only a HINT of co-transcription, not an
# operon (that needs transcript data; AureoWiki lists operons for NCTC 8325).
nb_name <- read_protein_metadata() |>
  transmute(locus_tag, nb = coalesce(if_else(ortholog_quality == "good", gene_symbol_best, NA_character_), product))
gff_all <- gff_all |> mutate(idx = row_number()) |> left_join(nb_name, by = "locus_tag") |>
  mutate(nb = coalesce(nb, "(no protein entry)"))   # e.g. RNA genes / pseudogenes
neighbours <- purrr::map_dfr(hits, \(h) {
  i <- gff_all$idx[gff_all$locus_tag == h]
  w <- gff_all |> filter(idx %in% (i - 2):(i + 2)) |> arrange(idx)
  me <- w |> filter(idx == i)
  # same-strand run with gaps < 150 bp, walking out from the hit (whole genome, not only the window)
  run <- i
  for (dir in c(-1, 1)) {
    j <- i
    repeat {
      k <- j + dir
      if (k < 1 || k > nrow(gff_all) || gff_all$strand[k] != me$strand) break
      gap <- if (dir > 0) gff_all$start[k] - gff_all$end[j] - 1 else gff_all$start[j] - gff_all$end[k] - 1
      if (gap >= 150) break
      run <- c(run, k); j <- k
    }
  }
  run_genes <- gff_all |> filter(idx %in% run) |> arrange(idx)
  tibble(locus_tag = h,
         neighbours_left_to_right = paste(sprintf("%s %s [%s]", str_remove(w$locus_tag, "B7H15_"),
                                                  str_trunc(w$nb, 30), w$strand), collapse = " | "),
         same_strand_run_n_genes = nrow(run_genes),
         same_strand_run = paste(sprintf("%s %s", str_remove(run_genes$locus_tag, "B7H15_"),
                                         str_trunc(run_genes$nb, 30)), collapse = ", "),
         other_hits_in_run = paste(setdiff(intersect(run_genes$locus_tag, hits), h), collapse = ", "))
})
annot <- ctx |> select(locus_tag, gene, product, strand, start, end, start_codon,
                       dist_to_stk1_Tn_bp, dist_to_stp1_Tn_bp, path_to_Tn, hit_in) |>
  left_join(meta_full, by = "locus_tag") |>
  left_join(neighbours, by = "locus_tag")
write_csv(annot, file.path(out, "hits_annotation_mval.csv"))

# Markdown overview generated from the table, so the text always matches the data
md_row <- \(r) sprintf("| %s | %s | %s | %s | %s | %s | %s | %s |",
  paste0("`", r$locus_tag, "`", if_else(is.na(r$gene), "", paste0(" ", r$gene))),
  r$up_protein_name, coalesce(r$up_ec, "–"),
  paste0(r$uniprot_acc, " (", r$up_reviewed, "; ", r$nctc8325_locus, ")"),
  coalesce(r$subcellular, "not given"),
  paste0(if_else(r$strand == "+", "forward (+)", "reverse (−)"), "; ",
         fmt_kb(r$dist_to_stk1_Tn_bp), " / ", fmt_kb(r$dist_to_stp1_Tn_bp),
         if_else(r$path_to_Tn == "direct", "", " (via origin)")),
  str_replace_all(r$hit_in, "\n", "; "),
  if_else(r$same_strand_run_n_genes > 1,
          paste0(r$same_strand_run_n_genes, " genes: ", r$same_strand_run,
                 if_else(r$other_hits_in_run != "", paste0(" (**also a hit: ", r$other_hits_in_run, "**)"), "")),
          "single gene"))
fmt_kb <- \(d) paste0(if_else(d >= 0, "+", "−"), formatC(abs(d) / 1000, format = "f", digits = 1, big.mark = ","), " kb")
writeLines(c(
  "# Step 6 hits: annotation and genome context (generated by scripts/06b_hit_genome_context_mval.R)",
  "",
  "Sources: `reference/JE2_protein_metadata.tsv.gz` (UniProt entry of the NCTC 8325 ortholog;",
  "all 8 orthologs `ortholog_quality == \"good\"`) and `data/raw/resources/JE2_genome.gff` (JE2,",
  "CP020619.1). UniProt text describes the NCTC 8325 protein, not a JE2 measurement. \"reviewed\" =",
  "Swiss-Prot (curated), \"unreviewed\" = TrEMBL (automatic). Distances: from the gene's start codon to",
  "the stk1 / stp1 insertion, + = higher genome coordinate, circular genome, shortest way.",
  "Same-strand run = neighbouring genes on the same strand with gaps < 150 bp: a hint of",
  "co-transcription only, not a proven operon.",
  "",
  "| hit | UniProt protein name | EC | UniProt (status; NCTC 8325 locus) | location (UniProt) | strand; distance to stk1-Tn / stp1-Tn | hit in | same-strand run (GFF) |",
  "|---|---|---|---|---|---|---|---|",
  purrr::map_chr(seq_len(nrow(annot)), \(i) md_row(annot[i, ])),
  "",
  "UniProt keywords per hit:",
  "",
  sprintf("- `%s`: %s", annot$locus_tag, annot$up_keywords),
  "",
  "GFF neighbours (2 each side, `[strand]`):",
  "",
  sprintf("- `%s`: %s", annot$locus_tag, annot$neighbours_left_to_right)),
  file.path(out, "hits_context_mval.md"))
print(annot |> select(locus_tag, gene, same_strand_run_n_genes, other_hits_in_run, subcellular), width = 200)

# where the insertions fall in the genome order of the hits (for the dashed line)
tn_after <- sum(ctx$start < min(ins))            # hits left of the insertions
tn_x <- tn_after + 0.5

# ---- 3. per-sample values, TSB and SN only (the hits are TSB-vs-SN comparisons) --
two_media <- c("TSB", "SN")
samples <- read_samples()
abund <- read_abundances("total") |> keep_je2() |>
  filter(locus_tag %in% hits) |>
  select(locus_tag, sample, normalized_abundance) |>
  inner_join(select(samples, sample, genotype, medium), by = "sample") |>   # by NAME (pitfall 2)
  filter(medium %in% two_media) |>
  mutate(medium = factor(medium, two_media)) |>
  left_join(select(ctx, locus_tag, x_label), by = "locus_tag")
stopifnot(!anyDuplicated(abund[c("locus_tag", "sample")]), n_distinct(abund$sample) == 18,
          all(count(distinct(abund, sample, genotype, medium), genotype, medium)$n == 3))
message("CHECK: ", n_distinct(abund$locus_tag), " proteins x ", n_distinct(abund$sample), " samples")

detect <- abund |> group_by(x_label, genotype, medium) |>
  summarise(n_det = sum(!is.na(normalized_abundance)), .groups = "drop")
y_low  <- floor(min(abund$normalized_abundance, na.rm = TRUE)) - 1
y_high <- ceiling(max(abund$normalized_abundance, na.rm = TRUE)) + 1.6   # room for the distance text

# distance text: each mutant panel shows the distance to ITS OWN insertion; anc has none
fmt_kb <- \(d) paste0(if_else(d >= 0, "+", "−"), formatC(abs(d) / 1000, format = "f", digits = 1, big.mark = ","), " kb")
dist_txt <- bind_rows(
  ctx |> transmute(x_label, genotype = "stk", via = path_to_Tn != "direct", text = fmt_kb(dist_to_stk1_Tn_bp)),
  ctx |> transmute(x_label, genotype = "stp", via = path_to_Tn != "direct", text = fmt_kb(dist_to_stp1_Tn_bp))) |>
  mutate(text = if_else(via, paste0(text, "\n(via ori)"), text)) |>
  mutate(genotype = factor(genotype, GENOTYPES))
tn_label <- tibble(genotype = factor(GENOTYPES, GENOTYPES),
                   text = c("", "stk1::Tn", "stp1::Tn"))

dodge <- position_dodge(width = 0.55)
p <- ggplot(abund |> filter(!is.na(normalized_abundance)),
            aes(x_label, normalized_abundance, colour = genotype, shape = medium, group = medium)) +
  geom_point(size = 2, position = dodge) +   # first layer: makes the x-axis discrete (gene names)
  geom_vline(xintercept = tn_x, linetype = 2, colour = "grey40") +
  geom_text(data = tn_label, aes(x = tn_x, y = y_low + 0.3, label = text), inherit.aes = FALSE,
            hjust = -0.08, size = 2.6, colour = "grey30", fontface = "italic") +
  geom_text(data = detect, aes(x_label, y_low + 0.3, group = medium,
                               label = if_else(n_det < 3, paste0(n_det, "/3"), "")),
            position = dodge, size = 2.5, colour = "firebrick", inherit.aes = FALSE) +
  geom_text(data = dist_txt, aes(x_label, y_high - 0.5, label = text), inherit.aes = FALSE,
            size = 2.7, colour = "grey25") +
  facet_wrap(~ genotype, ncol = 1,
             labeller = as_labeller(c(anc = "anc (wild type, no transposon)",
                                      stk = "stk (stk1::Tn, NE217): text on top = distance from the gene's start codon to the stk1 insertion",
                                      stp = "stp (stp1::Tn, NE1919; also low Stk1): text on top = distance to the stp1 insertion"))) +
  scale_colour_manual(values = GENOTYPE_COLOURS, guide = "none") +
  scale_shape_manual(values = c(TSB = 16, SN = 15), name = "medium") +
  scale_y_continuous(limits = c(y_low, y_high)) +
  labs(x = "hit proteins in genome order (JE2, CP020619.1): name (locus tag) / strand / comparison in which it is a hit",
       y = "normalised log2 abundance",
       title = "Step 6 hits per strain, ordered along the chromosome (each dot = one sample)",
       subtitle = paste0("Hits: treat lfc = 1, BH-FDR < 0.05 in the TSB-vs-SN interactions (stk vs anc, stp vs anc, stp vs stk); ",
                         "(expl.) = only 0.05 <= FDR < 0.1, not a hit.\n",
                         "Dashed line = position of the two insertions (stp1/stk1, ~1.233 Mb; 310 bp apart). ",
                         "Distance from the gene's start codon: + = towards higher, − = towards lower coordinates; circular genome, shortest way (via ori = round through position 1).\n",
                         "Red n/3 = replicates detected when fewer than 3; missing samples are not plotted (NA, not 0). SNpulse not shown.")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", strip.text = element_text(face = "bold", size = 8.5),
        axis.text.x = element_text(size = 8, lineheight = 0.95), plot.subtitle = element_text(size = 8.5))

ggsave(file.path(out, "fig_hits_by_genotype_mval.png"), p, width = 13, height = 10, dpi = 300)
message("Done: ", file.path(out, "fig_hits_by_genotype_mval.png"))
