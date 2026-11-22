# Step 11d (PROJECT_WORKFLOW_mval.md, after Step 11c): STRING network of the R6 and R7 proteins
#
# Run from the repository root:   Rscript C_mutant_validation_DinMete/scripts/11d_R6_R7_string_mval.R   (~15 s; the first
# run downloads the STRING links once, later runs read the cached file)
#
# Aim:      Step 11c showed WHERE R7 (up in stk) and R6 (down in stp) sit on the genome. Here we
#           ask WHAT the proteins of each block do together: which functional groups (clusters)
#           does STRING see among them, and does the whole cluster move with the block?
# Input:    R6_R7_genes_mval.csv (Step 11c; in_region == TRUE = the Step 11 block members, log2
#           diffs of the limma treat fit, observed only); STRING v12 links for the course genome
#           STRG0B09NBK (JE2_Bio253_HS26) from the STRING API, medium confidence (score >= 0.4);
#           data/raw/resources/string_enrichment_terms.txt + string_aliases.txt (term annotation);
#           read_protein_metadata() for gene symbols (only ortholog_quality == "good").
# Analysis: one network per region (nodes = block proteins). Clusters = Louvain communities on
#           the score-weighted links (seed fixed). Each cluster of >= 2 proteins gets the term
#           (KEGG, UniProt keyword, GO process/function, Reactome) shared
#           by most of its members that is most specific against the genome (Fisher test, cluster
#           vs all other annotated proteins). This label is DESCRIPTIVE: the clusters were picked
#           from a region we already selected, so the p-value is not a test result.
# Output:   fig_R6_R7_string_network_mval.png, R6_R7_string_edges_mval.tsv (cached API links),
#           R6_R7_string_clusters_mval.csv (protein, cluster, label, diffs),
#           R6_R7_string_cluster_terms_mval.csv (best terms per cluster), R6_R7_string_checks_mval.txt
# Check:    74 block proteins (49 R6 + 25 R7, as Step 11c); every protein maps to one STRING id;
#           kdpABC and rsbUVW-sigB fall into one cluster each in R7 (known operons = sanity check).

suppressMessages(source("R/setup.R"))
library(patchwork)
library(igraph, warn.conflicts = FALSE)
out <- "C_mutant_validation_DinMete/results/step11_robustness"
checks <- character()
check <- function(label, ok) {
  line <- sprintf("[%s] %s", if (isTRUE(ok)) "PASS" else "FAIL", label)
  message(line); checks <<- c(checks, line); invisible(ok)
}
note <- function(...) { line <- paste0("[NOTE] ", ...); message(line); checks <<- c(checks, line) }

STRING_ORG <- "STRG0B09NBK"   # course genome on string-db.org (see write_string_upload())
MIN_SCORE  <- 0.4             # STRING "medium confidence"
MAX_TERM_SIZE <- 300          # larger terms are too generic to name a cluster
FOCUS <- c(R7 = "stk_vs_anc_at_TSB", R6 = "stp_vs_anc_at_TSB")   # the contrast that flagged each block

# ---- 1. block proteins + annotation ---------------------------------------------------------
genes <- read_csv(file.path(out, "R6_R7_genes_mval.csv"), show_col_types = FALSE) |>
  filter(in_region) |> select(region = view, locus_tag, start, strand, product, starts_with("diff_"))
check("74 block proteins (49 R6 + 25 R7), as Step 11c",
      nrow(genes) == 74 && sum(genes$region == "R6") == 49 && sum(genes$region == "R7") == 25)

meta <- read_protein_metadata(cds_only = TRUE) |>
  transmute(locus_tag, sym = if_else(ortholog_quality %in% "good", gene_symbol_best, NA_character_))
# Node label, in this order: (1) symbol of a GOOD ortholog; (2) a gene name written at the end of
# the GenBank product ("... sigma factor SigB" -> sigB); (3) a phage module word from the product
# (capsid, portal, tail ...); (4) the B7H15 tag number. Weak ortholog symbols are NOT used: in
# R6 they put bacterial names (ureG, dnaB, ddl) on phage proteins.
PHAGE_WORDS <- "capsid|portal|terminase|tape measure|tail|head-tail|holin|integrase|recombinase|endonuclease|helicase|single-stranded DNA-binding"
genes <- genes |> left_join(meta, by = "locus_tag") |>
  mutate(prod_gene = str_match(product, "\\b([A-Z][a-z]{2}[A-Z])$")[, 2],
         prod_gene = if_else(is.na(prod_gene), NA, paste0(tolower(str_sub(prod_gene, 1, 1)), str_sub(prod_gene, 2))),
         prod_word = str_extract(product, regex(PHAGE_WORDS, ignore_case = TRUE)),
         label = coalesce(sym, prod_gene, prod_word, str_remove(locus_tag, "B7H15_")),
         label = if_else(label == str_remove(locus_tag, "B7H15_"), label,
                         paste0(label, "\n", str_remove(locus_tag, "B7H15_"))))
stopifnot(!anyDuplicated(genes$locus_tag))

# ---- 2. STRING links (API once, then cached) --------------------------------------------------
# Why the API and not data/raw/resources/string_physical_links.txt: the course file holds only
# 4 links inside R6, so almost all phage proteins would be isolated. The API gives the full
# functional network with every evidence channel, so we can say which links rest only on gene
# neighbourhood (nscore) -- such a link mostly restates "these genes sit next to each other".
edge_file <- file.path(out, "R6_R7_string_edges_mval.tsv")
if (!file.exists(edge_file)) {
  url <- paste0("https://string-db.org/api/tsv/network?identifiers=",
                paste(genes$locus_tag, collapse = "%0d"),
                "&species=", STRING_ORG, "&required_score=", MIN_SCORE * 1000,
                "&caller_identity=bio253_course")
  raw <- read_tsv(url, show_col_types = FALSE)
  stopifnot("score" %in% names(raw))
  write_tsv(raw, edge_file)
}
edges <- read_tsv(edge_file, show_col_types = FALSE) |>
  transmute(from = preferredName_A, to = preferredName_B, score,
            neighbourhood_only = score - nscore < 0.05 & pmax(fscore, pscore, ascore, escore, dscore, tscore) < 0.15) |>
  filter(score >= MIN_SCORE) |>
  mutate(a = pmin(from, to), b = pmax(from, to)) |> distinct(a, b, .keep_all = TRUE) |>
  select(from = a, to = b, score, neighbourhood_only)
check("every STRING link joins two block proteins (names = locus tags)",
      all(c(edges$from, edges$to) %in% genes$locus_tag))
edges <- edges |>
  left_join(genes |> select(from = locus_tag, r1 = region), by = "from") |>
  left_join(genes |> select(to = locus_tag, r2 = region), by = "to")
note(sprintf("%d STRING links (score >= %.1f): R7-R7 %d, R6-R6 %d, R6-R7 %d; %d rest only on gene neighbourhood",
             nrow(edges), MIN_SCORE, sum(edges$r1 == "R7" & edges$r2 == "R7"),
             sum(edges$r1 == "R6" & edges$r2 == "R6"), sum(edges$r1 != edges$r2), sum(edges$neighbourhood_only)))

# ---- 3. term annotation (course STRING files; ids map to locus tags via the alias table) -------
alias <- read_tsv("data/raw/resources/string_aliases.txt", show_col_types = FALSE) |>
  rename(sid = 1) |> filter(str_detect(alias, "^B7H15_")) |> distinct(sid, locus_tag = alias)
check("every block protein has exactly one STRING id in the course alias table",
      all(genes$locus_tag %in% alias$locus_tag) && !anyDuplicated(alias$locus_tag[alias$locus_tag %in% genes$locus_tag]))
TERM_CATS <- c("KEGG (Kyoto Encyclopedia of Genes and Genomes)", "Annotated Keywords (UniProt)",
               "Biological Process (Gene Ontology)", "Molecular Function (Gene Ontology)",
               "Reactome Pathways")   # not "Local Network Cluster": its names ("Mixed, incl. ...") say little
terms <- read_tsv("data/raw/resources/string_enrichment_terms.txt", show_col_types = FALSE) |>
  rename(sid = 1) |> filter(category %in% TERM_CATS) |>
  inner_join(alias, by = "sid") |> distinct(locus_tag, category, term, description)
n_genome <- n_distinct(terms$locus_tag)          # background = all annotated proteins
term_size <- terms |> count(term, name = "n_term")

# ---- 4. one network per region: clusters, labels, layout -------------------------------------
best_term <- function(members) {
  n_in <- sum(members %in% terms$locus_tag)
  terms |> filter(locus_tag %in% members) |> count(category, term, description, name = "k") |>
    filter(k >= 2) |> left_join(term_size, by = "term") |>
    # generic terms (> 300 proteins in the genome) and the structural keyword "Coiled coil"
    # do not describe a function, so they cannot name a cluster
    filter(n_term <= MAX_TERM_SIZE, description != "Coiled coil") |>
    mutate(n_cluster = length(members),
           p_fisher = purrr::map2_dbl(k, n_term, \(k, n) fisher.test(matrix(
             c(k, n_in - k, n - k, n_genome - n_in - n + k), 2), alternative = "greater")$p.value)) |>
    arrange(p_fisher, desc(k))
}

net_region <- function(reg) {
  nd <- genes |> filter(region == reg)
  ed <- edges |> filter(r1 == reg, r2 == reg)
  g <- graph_from_data_frame(ed |> select(from, to, weight = score), directed = FALSE,
                             vertices = nd |> select(name = locus_tag))
  set.seed(253)
  cl <- membership(cluster_louvain(g, weights = E(g)$weight))
  nd <- nd |> mutate(cluster_raw = as.integer(cl[locus_tag]))
  sizes <- table(nd$cluster_raw)
  # number clusters by size (big first); single proteins are "not linked" (cluster 0)
  ord <- names(sort(sizes[sizes >= 2], decreasing = TRUE))
  nd <- nd |> mutate(cluster = match(as.character(cluster_raw), ord), cluster = replace_na(cluster, 0L))
  tt <- nd |> filter(cluster > 0) |> group_by(cluster) |>
    group_modify(\(d, k) best_term(d$locus_tag) |> head(3)) |> ungroup() |> mutate(region = reg)
  lab <- tt |> group_by(cluster) |> slice(1) |> ungroup() |>
    transmute(cluster, term_lab = sprintf("%s (%d/%d)", str_trunc(description, 40), k, n_cluster))
  # second label line: the members' GenBank product names (from the annotation table, not memory);
  # phage proteins carry almost no STRING terms, but their products say "capsid", "portal", ...
  prod <- nd |> filter(cluster > 0, product != "hypothetical protein") |>
    mutate(p = str_replace_all(product, regex("^phage |,? ?\\b(family|protein|subunit|domain-containing)\\b", ignore_case = TRUE), " ") |>
             str_squish()) |>
    group_by(cluster) |> summarise(prod_lab = str_trunc(paste(unique(p), collapse = ", "), 60), .groups = "drop")
  nd <- nd |> left_join(lab, by = "cluster") |> left_join(prod, by = "cluster") |>
    mutate(fun = case_when(cluster == 0 ~ "not linked (score < 0.4)",
                           TRUE ~ paste0(coalesce(term_lab, "no specific STRING term"), "\n  ",
                                         coalesce(prod_lab, "hypothetical proteins only"))),
           cluster_name = if_else(cluster == 0, "-", paste0(reg, "-", cluster)))
  set.seed(253)
  xy <- layout_with_fr(g, weights = E(g)$weight, niter = 2000)
  nd <- nd |> mutate(x = xy[match(locus_tag, V(g)$name), 1], y = xy[match(locus_tag, V(g)$name), 2])
  # isolated proteins: park them in a grid (genome order) under the network so they do not clutter it
  iso <- nd$locus_tag %in% names(which(degree(g) == 0))
  if (any(iso)) {
    rng <- range(nd$x[!iso]); per_row <- 9
    i <- seq_len(sum(iso)) - 1
    step <- diff(rng) / (per_row - 1)
    nd$x[iso] <- rng[1] + (i %% per_row) * step
    nd$y[iso] <- min(nd$y[!iso]) - step * (1.3 + 1.1 * (i %/% per_row))
  }
  list(nodes = nd, edges = ed, terms = tt)
}
nets <- purrr::map(c(R7 = "R7", R6 = "R6"), net_region)
nodes <- bind_rows(purrr::map(nets, "nodes"))
cl_terms <- bind_rows(purrr::map(nets, "terms"))

# sanity: well-known operons should land in one cluster
same_cl <- function(tags)n_distinct(nodes$cluster[nodes$locus_tag %in% tags]) == 1 && all(tags %in% nodes$locus_tag)
check("R7: kdpC/kdpB/kdpA (B7H15_11530-11540) in one cluster", same_cl(sprintf("B7H15_%d", c(11530, 11535, 11540))))
check("R7: sigB/rsbW/rsbV/rsbU (B7H15_11470-11485) in one cluster", same_cl(sprintf("B7H15_%d", c(11470, 11475, 11480, 11485))))

# per-cluster effect in the contrast that flagged the block (observed values only, Step 11c)
nodes <- nodes |> mutate(focus_diff = if_else(region == "R7", diff_stk_vs_anc_at_TSB, diff_stp_vs_anc_at_TSB),
                         other_diff = if_else(region == "R7", diff_stp_vs_anc_at_TSB, diff_stk_vs_anc_at_TSB))
cl_eff <- nodes |> filter(cluster > 0) |> group_by(region, cluster_name, fun) |>
  summarise(n = n(), n_tested = sum(!is.na(focus_diff)), mean_focus = mean(focus_diff, na.rm = TRUE),
            mean_other = mean(other_diff, na.rm = TRUE), .groups = "drop")
for (i in seq_len(nrow(cl_eff))) with(cl_eff[i, ], note(sprintf(
  "%s  %s  n=%d (tested %d): mean log2 diff %s %+.2f, other mutant %+.2f",
  cluster_name, str_replace(fun, "
 *", " | "), n, n_tested, FOCUS[region], mean_focus, mean_other)))

write_csv(nodes |> select(region, cluster = cluster_name, cluster_function = fun, locus_tag, gene = sym,
                          product, strand, start, starts_with("diff_")),
          file.path(out, "R6_R7_string_clusters_mval.csv"))
write_csv(cl_terms |> mutate(cluster = paste0(region, "-", cluster), p_fisher = signif(p_fisher, 2)) |>
            select(region, cluster, category, term, description, k, n_cluster, n_term, p_fisher),
          file.path(out, "R6_R7_string_cluster_terms_mval.csv"))

# ---- 5. figure ------------------------------------------------------------------------------
hull_df <- function(nd) {      # convex hull per cluster, pushed out a bit so it frames the nodes
  nd |> filter(cluster > 0) |> group_by(cluster_name, fun) |> group_modify(\(d, k) {
    cx <- mean(d$x); cy <- mean(d$y); r <- 0.35
    pts <- purrr::map_dfr(seq(0, 2 * pi, length.out = 13)[-13], \(a)
                          tibble(x = d$x + r * cos(a), y = d$y + r * sin(a)))
    pts[chull(pts$x, pts$y), ]
  }) |> ungroup()
}
plot_net <- function(reg, title) {
  nd <- nodes |> filter(region == reg)
  ed <- nets[[reg]]$edges |>
    left_join(nd |> select(from = locus_tag, x, y), by = "from") |>
    left_join(nd |> select(to = locus_tag, xend = x, yend = y), by = "to")
  hl <- hull_df(nd)
  funs <- unique(nd$fun[nd$cluster > 0])
  pal <- setNames(rep(c("#4E79A7", "#F28E2B", "#59A14F", "#B07AA1", "#76B7B2", "#EDC948",
                        "#FF9DA7", "#9C755F", "#E15759", "#BAB0AC"), length.out = length(funs)), funs)
  ggplot() +
    geom_polygon(data = hl, aes(x, y, group = cluster_name, fill = fun), alpha = 0.22,
                 colour = NA) +
    scale_fill_manual(values = pal, guide = "none") +   # the function is written ON each cluster
    geom_segment(data = ed, aes(x, y, xend = xend, yend = yend, linewidth = score,
                                linetype = neighbourhood_only), colour = "grey45", alpha = 0.7) +
    scale_linewidth(range = c(0.2, 1.6), limits = c(0.4, 1), name = "STRING score") +
    scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"),
                          labels = c(`FALSE` = "other evidence", `TRUE` = "gene neighbourhood only"),
                          name = "link evidence") +
    geom_point(data = nd, aes(x, y), size = 6.2, colour = "grey20") +
    geom_point(data = nd, aes(x, y, colour = focus_diff), size = 5.2) +
    scale_colour_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-3, 3), oob = scales::squish, na.value = "grey80",
                         name = sprintf("%s\nlog2 diff (grey =\nnot tested)", FOCUS[reg])) +
    ggrepel::geom_text_repel(data = nd, aes(x, y, label = label), size = 2.5, lineheight = 0.8, fontface = "italic",
                             point.padding = 9, box.padding = 0.25, min.segment.length = 0.3,
                             max.overlaps = Inf, seed = 1) +
    # cluster name + function (STRING term, then member products), placed above each shaded area
    ggrepel::geom_label_repel(
      data = hl |> group_by(cluster_name, fun) |> summarise(x = mean(x), y = max(y), .groups = "drop") |>
        mutate(lab = paste0(cluster_name, ": ", fun)),
      aes(x, y, label = lab), fill = "white", size = 2.7, lineheight = 0.9, hjust = 0, alpha = 0.9,
      linewidth = 0.2, nudge_y = 0.4, direction = "both", min.segment.length = 0, seed = 2,
      max.overlaps = Inf) +
    annotate("text", x = min(nd$x[nd$cluster == 0 & !nd$locus_tag %in% c(ed$from, ed$to)], na.rm = TRUE),
             y = max(nd$y[nd$cluster == 0], na.rm = TRUE) + 0.55, hjust = 0, size = 3, fontface = "bold",
             label = sprintf("not linked in STRING (score < 0.4): %d proteins", sum(nd$cluster == 0))) +
    labs(title = title, x = NULL, y = NULL) +
    coord_equal(clip = "off") + theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11), legend.key.width = unit(0.9, "cm"))
}

eff_txt <- \(reg) cl_eff |> filter(region == reg) |>
  mutate(t = sprintf("%s %+.1f", cluster_name, mean_focus)) |> pull(t) |> paste(collapse = ", ")
pA <- plot_net("R7", sprintf("A  R7 (2.19 Mb, up in stk at TSB): %d proteins, %d links", sum(nodes$region == "R7"),
                             nrow(nets$R7$edges)))
pB <- plot_net("R6", sprintf("B  R6 (2.09 Mb, hlb prophage, down in stp at TSB): %d proteins, %d links",
                             sum(nodes$region == "R6"), nrow(nets$R6$edges)))
fig <- (pA | pB) + plot_layout(widths = c(1, 1.15), guides = "collect") +
  plot_annotation(
    title = "Step 11d: what do the R7 and R6 proteins do together? STRING network with functional clusters",
    subtitle = paste0("Node colour = log2 diff in the contrast that flagged the block (R7: stk vs anc, R6: stp vs anc, TSB; > 0 = higher in the mutant). ",
                      "Shaded areas = Louvain clusters; each label gives the most specific shared STRING term (k/n members) and, below it, the members' GenBank products.\n",
                      "Mean log2 diff per cluster: R7 ", eff_txt("R7"), "; R6 ", eff_txt("R6"), "."),
    caption = paste0("STRING v12, course genome ", STRING_ORG, ", links with combined score >= 0.4 inside each block (", sum(edges$r1 != edges$r2),
                     " R6-R7 links not drawn). Dashed = link rests only on gene neighbourhood, i.e. it mainly restates that the genes are adjacent.\n",
                     "Cluster labels: term shared by >= 2 members with the lowest Fisher p vs the genome (descriptive only, the region was pre-selected). ",
                     "Node labels: symbol of a good ortholog, else a gene name or phage-module word from the GenBank product, plus the B7H15 tag number. All values observed (no LOD imputation), n = 3 per group.\n",
                     "Reminder: the stp strain also lacks most Stk1 (polar effect in the stp1-stk1 operon)."),
    theme = theme(plot.title = element_text(face = "bold", size = 14), plot.subtitle = element_text(size = 9),
                  plot.caption = element_text(hjust = 0, size = 8)))
ggsave(file.path(out, "fig_R6_R7_string_network_mval.png"), fig, width = 20, height = 11, dpi = 160, bg = "white")

writeLines(c(sprintf("Step 11d checks, %s", Sys.Date()), checks), file.path(out, "R6_R7_string_checks_mval.txt"))
message("done")
