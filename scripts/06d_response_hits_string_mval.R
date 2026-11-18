# Step 6 add-on: what are the proteins in fig_response_vs_anc_mval.png, how do
# their STRING neighbours behave, and is the kdpD result the same in both mutants?
#
# Run from the repository root:   Rscript scripts/06d_response_hits_string_mval.R   (~1 min)
# Needs scripts/05_formal_test_mval.R first. Uses the STRING web API ONCE and caches
# the answer in the output folder (public reference data, no course data is sent:
# only NCTC 8325 locus tags such as SAOUHSC_02314).
#
# ---- Plan (coding rules, section 4) -----------------------------------------
# Aim:      (1) Annotate the 7 proteins shown in colour in fig_response_vs_anc_mval.png
#           (formal + exploratory TSB-vs-SN interaction hits, stk and stp vs anc).
#           (2) Put them in a STRING network with their best-scoring partners and colour
#           every node by OUR measured response difference, to see whether operon /
#           pathway partners move together (a real pathway effect) or the hit is alone.
#           (3) kdpD: compare stk and stp, and check the protein value at peptide level.
# Input:    data/derived/limma_treat_total_mval.tsv.gz  our limma ~0+group, treat(lfc = 1)
#           read_abundances("total"), read_peptides()     per-sample log2 values (not changed)
#           read_protein_metadata()                       names, UniProt, SAOUHSC via up_kegg
#           STRING v12 API, species 93061 (S. aureus NCTC 8325): interaction_partners
#             (score >= 0.7, 6 per hit) and network among hits + partners (score >= 0.7)
# Analysis: descriptive. response_X = log2(X_SN) - log2(X_TSB) = -diff(TSB_vs_SN_at_X).
#           Node colour = response_mutant - response_anc = -diff(interaction_X_vs_anc...),
#           positive = rises more (falls less) in SN in the mutant than in anc.
#           kdpD peptide check: each peptide seen in SN minus the mean of the SAME peptide
#           in TSB of the SAME genotype (peptide-matched, avoids mixing different peptides).
# Output:   mval_medi/results/step06_responses/hits_biology/
#             response_hits_function_mval.csv   annotation + responses of the 7 proteins
#             string_partners_mval.tsv, string_edges_mval.tsv  cached STRING answers
#             string_nodes_mval.csv              every node with JE2 locus tag and our numbers
#             fig_string_network_mval.png        own drawing, coloured by our data
#             fig_string_official_mval.png       STRING's own picture of the same nodes
#             fig_kdp_operon_raw_mval.png         kdp proteins, per sample, 9 groups
#             fig_kdpD_peptides_mval.png          peptide counts + peptide-matched SN change
#             kdp_contrasts_mval.csv, kdpD_peptide_matched_mval.csv
#             checks_mval.txt
# Check:    7 proteins; 27 samples / 3 per group joined by name; every STRING node maps to
#           exactly one JE2 locus tag (or is listed as unmapped); kdpD numbers equal Step 6.

suppressMessages(source("R/setup.R"))
suppressPackageStartupMessages(library(patchwork))
out <- "mval_medi/results/step06_responses/hits_biology"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
checks <- character()
note <- function(...) { msg <- paste0(...); message(msg); checks <<- c(checks, msg) }

STRING_SPECIES <- 93061   # S. aureus NCTC 8325 = the strain of our UniProt orthologs
STRING_SCORE   <- 700     # "high confidence" in STRING (0-1000 scale)
STRING_LIMIT   <- 6       # partners per hit: enough for operon/pathway, still readable
FDR_EXPLORATORY <- 0.1    # display only (see Step 6 figure guide), never a hit

INTER <- c(stk = "interaction_stk_vs_anc_at_TSB_vs_SN",
           stp = "interaction_stp_vs_anc_at_TSB_vs_SN")

# ---- 1. the 7 proteins of Fig 2 ------------------------------------------------
meta <- read_protein_metadata(cds_only = TRUE) |>
  filter(!is.na(locus_tag)) |>
  mutate(saouhsc = str_extract(up_kegg, "SAOUHSC_\\d+"),
         gene = coalesce(gene_symbol_best, locus_tag))
stopifnot(!anyDuplicated(meta$locus_tag))
# CDS only: drops the frameshift/pseudo fragments that share a locus tag with a CDS
res  <- read_tsv("data/derived/limma_treat_total_mval.tsv.gz", show_col_types = FALSE) |>
  semi_join(meta, by = "protein_key")

hits_long <- res |> filter(contrast %in% INTER, tested, FDR < FDR_EXPLORATORY) |>
  mutate(mutant = names(INTER)[match(contrast, INTER)],
         call = if_else(FDR < MVAL_FDR, "formal hit", "exploratory"))
hits <- unique(hits_long$locus_tag)
stopifnot(length(hits) == 7)
note("CHECK proteins in Fig 2 (formal + exploratory): ", length(hits))

# responses per genotype (observed estimates only, tested rows only)
resp <- res |>
  filter(str_detect(contrast, "^TSB_vs_SN_at_")) |>     # all proteins: needed for STRING partners too
  transmute(locus_tag, genotype = str_remove(contrast, "TSB_vs_SN_at_"),
            response = if_else(tested, -diff, NA_real_)) |>
  pivot_wider(names_from = genotype, values_from = response, names_prefix = "resp_")
inter <- res |> filter(contrast %in% INTER) |>
  transmute(locus_tag, mutant = names(INTER)[match(contrast, INTER)],
            delta = if_else(tested, -diff, NA_real_),          # resp_mutant - resp_anc
            FDR = if_else(tested, FDR, NA_real_))

hit_table <- meta |> filter(locus_tag %in% hits) |>
  select(locus_tag, gene, product, uniprot_acc, ortholog_quality, saouhsc,
         up_protein_name, up_ec, up_keywords) |>
  left_join(resp, by = "locus_tag") |>
  left_join(hits_long |> group_by(locus_tag) |>
              summarise(hit_in = paste0(mutant, " vs anc (", call, ", delta ",
                                        sprintf("%+.2f", -diff), ", FDR ", signif(FDR, 2), ")",
                                        collapse = "; "), .groups = "drop"),
            by = "locus_tag") |>
  arrange(hit_in, locus_tag)
stopifnot(nrow(hit_table) == 7, all(hit_table$ortholog_quality == "good"), !anyNA(hit_table$saouhsc))
write_csv(hit_table, file.path(out, "response_hits_function_mval.csv"))

# ---- 2. STRING: partners and edges (cached) -------------------------------------
STRING_URL <- "https://version-12-5.string-db.org/api/"   # stable address = same answer next week
string_get <- function(endpoint, ids, ...) {
  url <- paste0(STRING_URL, "tsv/", endpoint)
  r <- httr::POST(url, body = list(identifiers = paste(ids, collapse = "\r"),
                                   species = STRING_SPECIES, caller_identity = "bio253_course", ...),
                  encode = "form", httr::timeout(60))
  httr::stop_for_status(r)
  read_tsv(I(httr::content(r, "text", encoding = "UTF-8")), show_col_types = FALSE)
}
f_part <- file.path(out, "string_partners_mval.tsv")
f_edge <- file.path(out, "string_edges_mval.tsv")
if (!file.exists(f_part) || !file.exists(f_edge)) {
  message("querying STRING (first run only) ...")
  partners <- string_get("interaction_partners", hit_table$saouhsc,
                         required_score = STRING_SCORE, limit = STRING_LIMIT)
  node_ids <- unique(c(paste0(STRING_SPECIES, ".", hit_table$saouhsc), partners$stringId_B))
  edges <- string_get("network", node_ids, required_score = STRING_SCORE)
  version <- jsonlite::fromJSON(httr::content(
    httr::GET(paste0(STRING_URL, "json/version"), httr::timeout(30)), "text", encoding = "UTF-8"))
  write_tsv(partners, f_part); write_tsv(edges, f_edge)
  writeLines(paste("STRING", version$string_version[1], "species", STRING_SPECIES, "queried", Sys.Date()),
             file.path(out, "string_version_mval.txt"))
  # STRING's own picture of the same node set (for comparison with our drawing)
  img <- httr::POST(paste0(STRING_URL, "highres_image/network"),
                    body = list(identifiers = paste(node_ids, collapse = "\r"),
                                species = STRING_SPECIES, required_score = STRING_SCORE,
                                network_flavor = "evidence", caller_identity = "bio253_course"),
                    encode = "form", httr::timeout(60))
  httr::stop_for_status(img)
  writeBin(httr::content(img, "raw"), file.path(out, "fig_string_official_mval.png"))
}
partners <- read_tsv(f_part, show_col_types = FALSE)
edges <- read_tsv(f_edge, show_col_types = FALSE) |>
  transmute(a = str_remove(stringId_A, "^\\d+\\."), b = str_remove(stringId_B, "^\\d+\\."),
            score, name_a = preferredName_A, name_b = preferredName_B) |>
  distinct(pair = paste(pmin(a, b), pmax(a, b)), .keep_all = TRUE)   # one row per pair
note("CHECK STRING: ", nrow(partners), " partner rows, ", nrow(edges),
     " edges (score >= ", STRING_SCORE / 1000, ") among hits + partners; ",
     readLines(file.path(out, "string_version_mval.txt")))

# ---- 3. nodes: map STRING (SAOUHSC) back to JE2 and add our numbers ---------------
string_names <- bind_rows(select(edges, id = a, sname = name_a), select(edges, id = b, sname = name_b))
nodes <- tibble(saouhsc = unique(c(hit_table$saouhsc, str_remove(partners$stringId_B, "^\\d+\\."))))  |>
  left_join(distinct(string_names, id, .keep_all = TRUE), by = c("saouhsc" = "id")) |>
  left_join(meta |> filter(ortholog_quality == "good", !is.na(saouhsc)) |>
              select(saouhsc, locus_tag, gene, product, blastortho_evalue),
            by = "saouhsc", relationship = "many-to-many")
# JE2 has a second copy of some arc genes (B7H15_00330-00350) that maps to the same
# NCTC 8325 protein. Keep ONE JE2 protein per STRING node: the Fig 2 protein if it is
# one, otherwise the best ortholog (lowest e-value), on a tie the one detected in more
# of the 27 samples (more data to colour the node with). The other copies are listed.
ab <- read_abundances("total")
n_det <- ab |> group_by(locus_tag) |> summarise(n_detected = sum(!is.na(normalized_abundance)))
nodes <- nodes |> left_join(n_det, by = "locus_tag") |>
  arrange(saouhsc, !(locus_tag %in% hits), blastortho_evalue, desc(n_detected))
dup <- nodes |> group_by(saouhsc) |> filter(n() > 1) |>
  summarise(txt = paste0(saouhsc[1], " -> kept ", locus_tag[1], " (", n_detected[1], "/27), not ",
                         paste0(locus_tag[-1], " (", n_detected[-1], "/27)", collapse = "/")),
            .groups = "drop")
note("CHECK nodes: ", n_distinct(nodes$saouhsc), " STRING nodes; unmapped to JE2: ",
     sum(is.na(nodes$locus_tag)), "; mapped to >1 JE2 locus: ", nrow(dup),
     if (nrow(dup)) paste0(" (", paste(dup$txt, collapse = "; "), ")") else "")
nodes <- nodes |> distinct(saouhsc, .keep_all = TRUE) |>
  mutate(is_hit = locus_tag %in% hits,
         # gene symbol, or a short product name if there is none, plus the locus number
         # (two different JE2 proteins are both called argF, so the number is needed)
         name = if_else(gene != locus_tag, gene, str_trunc(product, 18)),
         label = paste0(name, "\n", str_remove(locus_tag, "B7H15_"))) |>
  left_join(resp, by = "locus_tag")
write_csv(nodes |> left_join(inter |> pivot_wider(names_from = mutant, values_from = c(delta, FDR)),
                             by = "locus_tag"),
          file.path(out, "string_nodes_mval.csv"))

# ---- 4. layout: one small circle per STRING module ----------------------------------
# Module = connected group of nodes. Hits sit in the middle, partners on a ring. This
# is simpler to read (and to explain) than a force layout, and it is deterministic.
module_of <- function(ids, e) {
  m <- setNames(seq_along(ids), ids)
  repeat {                                  # merge until no edge joins two modules
    changed <- FALSE
    for (i in seq_len(nrow(e))) {
      ma <- m[e$a[i]]; mb <- m[e$b[i]]
      if (!is.na(ma) && !is.na(mb) && ma != mb) { m[m == max(ma, mb)] <- min(ma, mb); changed <- TRUE }
    }
    if (!changed) return(m)
  }
}
nodes$module <- module_of(nodes$saouhsc, edges)[nodes$saouhsc]
mods <- nodes |> group_by(module) |> summarise(first_hit = min(label[is_hit]), n = n(), .groups = "drop") |>
  arrange(-n, first_hit) |> mutate(slot = row_number(), cx = ((slot - 1) %% 3) * 3.2,
                                    cy = -((slot - 1) %/% 3) * 3.2)
nodes <- nodes |> left_join(mods, by = "module") |> group_by(module) |>
  mutate(k_hit = cumsum(is_hit), n_hit = sum(is_hit), k_par = cumsum(!is_hit), n_par = sum(!is_hit),
         x = if_else(is_hit, cx + (k_hit - (n_hit + 1) / 2) * 0.95, cx + 1.15 * cos(2 * pi * k_par / n_par)),
         y = if_else(is_hit, cy, cy + 1.15 * sin(2 * pi * k_par / n_par))) |>
  ungroup()
seg <- edges |> left_join(select(nodes, saouhsc, x, y), by = c("a" = "saouhsc")) |>
  left_join(select(nodes, saouhsc, xend = x, yend = y), by = c("b" = "saouhsc"))

# one panel per mutant, same layout, colour = our measured response difference
plot_nodes <- inter |> filter(locus_tag %in% nodes$locus_tag) |>
  right_join(tidyr::crossing(select(nodes, saouhsc, locus_tag, label, is_hit, x, y),
                             mutant = names(INTER)),
             by = c("locus_tag", "mutant")) |>
  mutate(panel = if_else(mutant == "stk", "stk vs anc",
                         "stp vs anc  (stp also lacks most Stk1)"),
         sig = !is.na(FDR) & FDR < MVAL_FDR)
p_net <- ggplot() +
  geom_segment(data = seg, aes(x, y, xend = xend, yend = yend, linewidth = score),
               colour = "grey65") +
  geom_point(data = plot_nodes, aes(x, y, fill = pmax(pmin(delta, 5), -5), shape = is_hit,
                                    size = is_hit, stroke = if_else(sig, 1.4, 0.4))) +
  geom_text(data = plot_nodes, aes(x, y + if_else(is_hit, 0.2, 0.13), label = label,
                                   fontface = if_else(is_hit, "bold", "plain")),
            size = 2.3, vjust = 0, lineheight = 0.8) +
  facet_wrap(~panel) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                       limits = c(-5, 5), na.value = "grey80",
                       name = "response in mutant\nminus response in anc\n(log2, SN - TSB)") +
  scale_shape_manual(values = c(`TRUE` = 22, `FALSE` = 21), labels = c("STRING partner", "Fig 2 protein"),
                     name = NULL) +
  scale_size_manual(values = c(`TRUE` = 5, `FALSE` = 3.2), guide = "none") +
  scale_linewidth(range = c(0.2, 1.2), name = "STRING score") +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  coord_equal() + theme_void(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(size = 10, margin = margin(4, 0, 4, 0))) +
  labs(title = "STRING neighbourhood of the proteins in fig_response_vs_anc_mval.png",
       subtitle = paste0("STRING v12, S. aureus NCTC 8325 orthologs, score >= 0.7, up to ", STRING_LIMIT,
                         " partners per protein. Red = rises more in SN (PASN) in the mutant than in anc;",
                         "\nblue = rises less / falls more. Grey = not testable (< 2/3 replicates in a group).",
                         " Thick outline = formal interaction hit (treat lfc = 1, BH-FDR < 0.05)."),
       caption = "B7H15_06000 has no STRING partner at score >= 0.7. Colours are observed limma estimates; n = 3 per group.")
ggsave(file.path(out, "fig_string_network_mval.png"), p_net, width = 14, height = 6.3, dpi = 200, bg = "white")

# ---- 5. kdpD in both mutants -------------------------------------------------------
KDP <- meta |> filter(str_detect(coalesce(gene_symbol_best, ""), "^kdp[A-E]$") |
                        str_detect(product, "potassium-transporting ATPase subunit F")) |>
  select(locus_tag, gene, product)
note("CHECK kdp proteins in the annotation: ", paste(KDP$gene, KDP$locus_tag, collapse = ", "))
samples <- read_samples()
stopifnot(setequal(unique(ab$sample), samples$sample), nrow(samples) == 27,
          all(count(samples, genotype, medium)$n == 3))   # pitfall 2: join by name

p_raw <- plot_protein(KDP$locus_tag, abundances = ab) +
  labs(title = "kdp operon, total proteome, per sample",
       subtitle = "kdpFABC = K+ transporter, kdpD/kdpE = sensor kinase / response regulator")
ggsave(file.path(out, "fig_kdp_operon_raw_mval.png"), p_raw, width = 10, height = 6, dpi = 200, bg = "white")

kdp_contrasts <- res |> filter(locus_tag %in% KDP$locus_tag,
                               contrast %in% c(paste0("TSB_vs_SN_at_", GENOTYPES), INTER,
                                               "stk_vs_anc_at_TSB", "stp_vs_anc_at_TSB")) |>
  left_join(select(KDP, locus_tag, gene), by = "locus_tag") |>
  select(gene, locus_tag, contrast, diff, std.error, FDR, tested, min_obs_in_groups, estimate_type) |>
  arrange(gene, contrast)
write_csv(kdp_contrasts, file.path(out, "kdp_contrasts_mval.csv"))

# peptide evidence: how many kdpD/kdpE peptides per sample, and how much does each
# peptide seen in SN change against the same peptide in TSB (same genotype)?
pep <- read_peptides() |> filter(!is.na(intensity), medium %in% c("TSB", "SN"))
n_pep <- pep |> filter(locus_tag %in% KDP$locus_tag) |>
  count(locus_tag, sample, genotype, medium, replicate, name = "n_peptides") |>
  complete(locus_tag, nesting(sample, genotype, medium, replicate), fill = list(n_peptides = 0)) |>
  left_join(select(KDP, locus_tag, gene), by = "locus_tag")
# genome-wide reference: share of TSB peptide evidence that is still seen in SN
frac <- pep |> count(locus_tag, medium, name = "n") |>
  pivot_wider(names_from = medium, values_from = n, values_fill = 0) |>
  filter(TSB >= 90) |>                      # >= 10 peptides per TSB sample on average
  mutate(frac_SN = SN / TSB, rank = rank(frac_SN, ties.method = "min"))
for (g in c("kdpD", "kdpE")) {
  r <- frac |> filter(locus_tag == KDP$locus_tag[KDP$gene == g])
  note("CHECK ", g, ": peptide observations SN/TSB = ", round(r$frac_SN, 3), " (rank ", r$rank,
       " of ", nrow(frac), " well-covered proteins; median protein ", round(median(frac$frac_SN), 2), ")")
}
note("CHECK total peptides per sample, median TSB / SN: ",
     paste(pep |> count(sample, medium) |> group_by(medium) |> summarise(m = median(n)) |>
             with(paste(medium, m)), collapse = " / "))

kdpD_tag <- KDP$locus_tag[KDP$gene == "kdpD"]
pm <- pep |> filter(locus_tag %in% KDP$locus_tag[KDP$gene %in% c("kdpD", "kdpE")]) |>
  mutate(l2 = log2(intensity)) |>
  left_join(select(KDP, locus_tag, gene), by = "locus_tag")
pm_tsb <- pm |> filter(medium == "TSB") |> group_by(gene, genotype, peptide) |>
  summarise(tsb_mean = mean(l2), n_tsb = n(), .groups = "drop")
matched <- pm |> filter(medium == "SN") |>
  select(gene, sample, genotype, replicate, peptide, protein_start, l2) |>
  left_join(pm_tsb, by = c("gene", "genotype", "peptide")) |>
  mutate(change_SN_minus_TSB = l2 - tsb_mean) |> arrange(gene, genotype, sample)
write_csv(matched, file.path(out, "kdpD_peptide_matched_mval.csv"))
note("CHECK peptide-matched SN - TSB, median per gene x genotype: ",
     paste(matched |> filter(!is.na(change_SN_minus_TSB)) |> group_by(gene, genotype) |>
             summarise(t = sprintf("%s %s %+.2f (n=%d)", first(gene), first(genotype),
                                   median(change_SN_minus_TSB), n()), .groups = "drop") |> pull(t),
           collapse = "; "))

# protein-level response for comparison (the numbers behind Fig 2)
kd <- kdp_contrasts |> filter(locus_tag == kdpD_tag, str_detect(contrast, "^TSB_vs_SN_at_")) |>
  mutate(genotype = str_remove(contrast, "TSB_vs_SN_at_"), resp = if_else(tested, -diff, NA_real_))
note("CHECK kdpD protein response (SN - TSB): ",
     paste(sprintf("%s %s", kd$genotype, ifelse(is.na(kd$resp), "not testable", sprintf("%+.2f", kd$resp))),
           collapse = "; "))
stopifnot(abs(kd$resp[kd$genotype == "stk"] - (-4.09)) < 0.01)   # equals Step 6 / notebook

p_cnt <- n_pep |> filter(gene %in% c("kdpD", "kdpE", "kdpB")) |>
  ggplot(aes(medium, n_peptides, colour = genotype)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  facet_wrap(~gene) + scale_colour_manual(values = GENOTYPE_COLOURS) +
  labs(x = NULL, y = "peptides detected per sample",
       title = "A. Peptide evidence: kdp proteins almost vanish in SN in ALL three strains")
p_pm <- matched |> filter(!is.na(change_SN_minus_TSB)) |>
  ggplot(aes(genotype, change_SN_minus_TSB, colour = genotype)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_point(position = position_jitter(width = 0.12, seed = 1), size = 2.2) +
  geom_point(data = kd |> filter(!is.na(resp)) |> mutate(gene = "kdpD"),
             aes(genotype, resp), inherit.aes = FALSE, shape = 4, size = 4, stroke = 1.3) +
  facet_wrap(~gene) + scale_colour_manual(values = GENOTYPE_COLOURS, guide = "none") +
  labs(x = NULL, y = "log2(peptide in SN) - mean log2(same peptide in TSB)",
       title = "B. Each peptide still seen in SN, against the same peptide in TSB (same genotype)",
       caption = paste("x = protein-level response used in Fig 2 (kdpD; stp not testable).",
                       "Peptides seen in SN but never in TSB of that genotype have no point."))
ggsave(file.path(out, "fig_kdpD_peptides_mval.png"), p_cnt / p_pm, width = 10, height = 8.5,
       dpi = 200, bg = "white")

writeLines(checks, file.path(out, "checks_mval.txt"))
message("done: ", out)
