#!/usr/bin/env Rscript
# make_protein_metadata.R -- build a metadata table for every entry of the
# JE2 search database (db1 FASTA), combining three sources:
#
#   1. the FASTA description lines  (locus tag, product, protein id, location,
#      pseudo flag, and the BLASTORTHO UniProt ortholog + its e-value)
#   2. the NCBI GFF                 (gene biotype, gene symbol, Note, inference,
#      Dbxref, translation table -- joined on the locus tag)
#   3. the STRING v12 export        (GO / KEGG / UniProt keywords / COMPARTMENTS /
#      Reactome / local network clusters -- joined on the FASTA protein key)
#
# Writes two files next to this script:
#   JE2_protein_metadata.tsv   one row per FASTA entry (wide)
#   JE2_protein_terms.tsv      one row per protein x annotation term (long)
#
# Run:  Rscript make_protein_metadata.R

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr); library(purrr)
})

here    <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))[1]))
if (is.na(here) || !nzchar(here)) here <- "."
fasta   <- file.path(here, "p41881_db1_SA_JE2_modNmoreTranslations_20260630.fasta")
gff     <- file.path(here, "SE_SE2_GCA_002085525.1_ASM208552v1_genomic.gff")
str_dir <- file.path(here, "stringDBresources")
stopifnot(file.exists(fasta), file.exists(gff), dir.exists(str_dir))

# ---- 1. FASTA ---------------------------------------------------------------
lines   <- read_lines(fasta)
is_hdr  <- str_starts(lines, ">")
headers <- str_sub(lines[is_hdr], 2)
# amino-acid length: sequence lines belonging to each header
grp     <- cumsum(is_hdr)
aa_len  <- tapply(nchar(lines[!is_hdr]), grp[!is_hdr], sum)
aa_len  <- as.integer(aa_len[as.character(seq_along(headers))])
seqs    <- vapply(split(lines[!is_hdr], grp[!is_hdr]), paste, character(1), collapse = "")
seqs    <- unname(seqs[as.character(seq_along(headers))])

# number of fully tryptic peptides of 6-30 aa (cut after K/R, not before P).
# Useful downstream for length-normalised abundances (iBAQ-like).
n_tryptic <- function(s) {
  if (is.na(s) || !nzchar(s)) return(NA_integer_)
  pep <- str_split(str_replace_all(s, "([KR])(?!P)", "\\1\n"), "\n")[[1]]
  sum(nchar(pep) >= 6 & nchar(pep) <= 30)
}

brack <- function(x, key) str_match(x, paste0("\\[", key, "=([^\\]]*)\\]"))[, 2]

# GFF/NCBI attributes percent-encode the characters that would otherwise break
# the format (";" -> %3B, "," -> %2C, "=" -> %3D, " " -> %20, "%" -> %25).
url_decode <- function(x) {
  x |> str_replace_all("%3B", ";") |> str_replace_all("%2C", ",") |>
    str_replace_all("%3D", "=") |> str_replace_all("%20", " ") |>
    str_replace_all("%25", "%")
}
# "0" and real numbers parse; the two "BLASTORTHO None" entries carry "evalue=None".
as_evalue <- function(x) suppressWarnings(as.numeric(na_if(x, "None")))

fa <- tibble(
  fasta_header     = headers,
  fasta_id         = str_split_fixed(headers, " ", 2)[, 1],
  description      = str_split_fixed(headers, " ", 2)[, 2],
  aa_length        = aa_len,
  sequence         = seqs
) |>
  mutate(
    protein_key = str_split_fixed(fasta_id, "\\|", 3)[, 2],
    entry_class = case_when(
      str_starts(protein_key, "lcl_")                 ~ "cds",
      str_starts(protein_key, "pseudoOrframeshift")   ~ "pseudo_frameshift_translation",
      str_starts(protein_key, "zz_")                  ~ "contaminant",
      TRUE                                            ~ "database_info"),
    # --- bracketed keys (CDS entries)
    locus_tag     = brack(description, "locus_tag"),
    product       = brack(description, "protein"),
    protein_id    = brack(description, "protein_id"),
    gene_fasta    = brack(description, "gene"),
    location_raw  = brack(description, "location"),
    gbkey         = brack(description, "gbkey"),
    is_pseudo     = !is.na(brack(description, "pseudo")),
    strand        = case_when(is.na(location_raw) ~ NA_character_,
                              str_detect(location_raw, "^complement") ~ "-", TRUE ~ "+"),
    start         = as.integer(str_match(location_raw, "(\\d+)\\.\\.")[, 2]),
    end           = as.integer(str_match(location_raw, "\\.\\.(\\d+)")[, 2]),
    cds_length_nt = end - start + 1L,
    # --- BLASTORTHO: the UniProt ortholog in the NCTC 8325 reference proteome
    blastortho_raw     = na_if(str_match(description, "BLASTORTHO (\\S+)")[, 2], "None"),
    blastortho_evalue  = as_evalue(str_match(description, "evalue=(\\S+)")[, 2]),
    uniprot_db         = str_match(blastortho_raw, "^(sp|tr)_")[, 2],
    uniprot_acc        = str_match(blastortho_raw, "^(?:sp|tr)_([A-Z0-9]+)_")[, 2],
    uniprot_entry_name = str_match(blastortho_raw, "^(?:sp|tr)_[A-Z0-9]+_(.+)$")[, 2],
    # a gene symbol can be read off a reviewed (sp_) entry name, e.g. DNAA_STAA8
    uniprot_gene_guess = if_else(uniprot_db == "sp",
                                 str_remove(uniprot_entry_name, "_STAA8$"), NA_character_),
    # --- 6-frame translations of pseudogenes / frameshifted ORFs
    fs_index     = as.integer(str_match(protein_key, "^pseudoOrframeshift(\\d+)_")[, 2]),
    fs_direction = recode(str_match(protein_key, "_(fw|rv)_")[, 2],
                          fw = "forward", rv = "reverse", .default = NA_character_),
    fs_frame     = as.integer(str_match(protein_key, "_frm(\\d+)$")[, 2]),
    # Two header shapes occur: "ID=cds-LOCUS;...;Note=..." and "LOCUS; ID=gene-LOCUS;Name=LOCUS"
    fs_locus_tag = coalesce(str_match(description, "ID=(?:cds|gene)-([A-Za-z0-9_]+)")[, 2],
                            str_match(description, "^([A-Za-z0-9]+_[0-9]+);")[, 2]),
    fs_note      = url_decode(str_match(description, "Note=(.+)$")[, 2]),
    # --- contaminants: "zz_FGCZCont0000_P61626_LYSC_HUMAN blastpHomologue_5.0e-107"
    cont_uniprot_acc = str_match(description, "zz_FGCZCont\\d+_([A-Z0-9]+)_[A-Z0-9]+_")[, 2],
    cont_entry_name  = str_match(description, "zz_FGCZCont\\d+_[A-Z0-9]+_([A-Z0-9]+_[A-Z]+)")[, 2],
    cont_evalue      = as.numeric(str_match(description, "blastpHomologue_(\\S+)")[, 2]),
    # frameshift translations inherit the locus tag of the gene they come from
    locus_tag = coalesce(locus_tag, fs_locus_tag)
  ) |>
  mutate(n_tryptic_peptides = map_int(sequence, n_tryptic))

message(sprintf("FASTA: %d entries (%s)", nrow(fa),
                paste(sprintf("%s=%d", names(table(fa$entry_class)), table(fa$entry_class)),
                      collapse = ", ")))

# ---- 2. GFF -----------------------------------------------------------------
g <- read_tsv(gff, comment = "#", col_names = c("seqid", "source", "type", "gstart",
                                                "gend", "score", "gstrand", "phase", "attr"),
              col_types = cols(.default = col_character(),
                               gstart = col_integer(), gend = col_integer()))
attr_of <- function(a, key) str_match(a, paste0("(?:^|;)", key, "=([^;]*)"))[, 2]

gff_cds <- g |> filter(type == "CDS") |>
  transmute(locus_tag        = attr_of(attr, "locus_tag"),
            gff_product      = url_decode(attr_of(attr, "product")),
            gff_gene         = attr_of(attr, "gene"),
            gff_note         = url_decode(attr_of(attr, "Note")),
            gff_inference    = attr_of(attr, "inference"),
            gff_dbxref       = attr_of(attr, "Dbxref"),
            gff_transl_table = attr_of(attr, "transl_table"),
            gff_start = gstart, gff_end = gend, gff_strand = gstrand) |>
  distinct(locus_tag, .keep_all = TRUE)

gff_gene <- g |> filter(type %in% c("gene", "pseudogene")) |>
  transmute(locus_tag     = attr_of(attr, "locus_tag"),
            gene_biotype  = attr_of(attr, "gene_biotype"),
            gene_symbol   = attr_of(attr, "gene"),
            gene_name_gff = attr_of(attr, "Name"),
            gene_feature  = type) |>
  distinct(locus_tag, .keep_all = TRUE)

message(sprintf("GFF: %d CDS, %d gene/pseudogene records", nrow(gff_cds), nrow(gff_gene)))

# ---- 3. STRING --------------------------------------------------------------
strip <- function(x) str_remove(x, "^STRG0A35SHG\\.")

terms <- read_tsv(file.path(str_dir, "STRG0A35SHG.protein.enrichment.terms.v12.0.txt"),
                  col_types = cols(.default = col_character())) |>
  rename(string_protein_id = 1, category = 2, term = 3, term_description = 4) |>
  mutate(protein_key = strip(string_protein_id))

clusters <- read_tsv(file.path(str_dir, "STRG0A35SHG.clusters.proteins.v12.0.txt"),
                     col_types = cols(.default = col_character())) |>
  rename(cluster_id = 2, string_protein_id = 3) |>
  mutate(protein_key = strip(string_protein_id)) |>
  left_join(read_tsv(file.path(str_dir, "STRG0A35SHG.clusters.info.v12.0.txt"),
                     col_types = cols(.default = col_character())) |>
              rename(cluster_id = 2, cluster_size = 3, cluster_best_described_by = 4) |>
              select(cluster_id, cluster_size, cluster_best_described_by),
            by = "cluster_id")

message(sprintf("STRING: %d enrichment terms, %d cluster memberships, %d proteins covered",
                nrow(terms), nrow(clusters), n_distinct(terms$protein_key)))

# wide: collapse each category into one cell as "term=description | term=description"
cat_col <- c("Biological Process (Gene Ontology)"          = "go_bp",
             "Molecular Function (Gene Ontology)"          = "go_mf",
             "Cellular Component (Gene Ontology)"          = "go_cc",
             "KEGG (Kyoto Encyclopedia of Genes and Genomes)" = "kegg",
             "Annotated Keywords (UniProt)"                = "uniprot_keywords",
             "Subcellular localization (COMPARTMENTS)"     = "compartments",
             "Reactome Pathways"                           = "reactome",
             "Local Network Cluster (STRING)"              = "string_local_cluster")
terms_wide <- terms |>
  mutate(col = unname(cat_col[category])) |>
  filter(!is.na(col)) |>
  group_by(protein_key, col) |>
  summarise(value = paste(paste0(term, "=", term_description), collapse = " | "),
            n = n(), .groups = "drop") |>
  pivot_wider(names_from = col, values_from = c(value, n), names_glue = "{col}{ifelse(.value=='n','_n','')}")

clusters_wide <- clusters |>
  group_by(protein_key) |>
  summarise(string_clusters = paste(cluster_id, collapse = " | "),
            string_cluster_descriptions = paste(cluster_best_described_by, collapse = " | "),
            string_n_clusters = n(), .groups = "drop")

# ---- 3b. UniProt reference proteome (UP000008816, S. aureus NCTC 8325) ------
# The FASTA's BLASTORTHO field points into this proteome, so joining on the
# accession adds what the NCBI annotation lacks: real gene symbols, protein
# names, EC numbers, keywords, GO ids and subcellular locations.
# Cached locally; delete the file to refresh it.
up_file <- file.path(here, "uniprot_UP000008816.tsv.gz")
if (!file.exists(up_file)) {
  message("downloading UniProt proteome UP000008816 ...")
  url <- paste0("https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000008816",
                "&format=tsv&compressed=true&fields=accession,id,reviewed,protein_name,",
                "gene_names,gene_primary,ec,length,keyword,go_id,cc_subcellular_location,",
                "xref_kegg,protein_existence")
  utils::download.file(url, up_file, mode = "wb", quiet = TRUE)
}
up <- read_tsv(up_file, col_types = cols(.default = col_character(), Length = col_integer())) |>
  transmute(uniprot_acc        = Entry,
            up_entry_name      = `Entry Name`,
            up_reviewed        = Reviewed,
            up_protein_name    = `Protein names`,
            up_gene_primary    = `Gene Names (primary)`,
            up_gene_names      = `Gene Names`,
            up_ec              = `EC number`,
            up_length          = Length,
            up_keywords        = Keywords,
            up_go_ids          = `Gene Ontology IDs`,
            up_subcellular     = str_remove(`Subcellular location [CC]`, "^SUBCELLULAR LOCATION: "),
            up_kegg            = KEGG,
            up_protein_existence = `Protein existence`) |>
  distinct(uniprot_acc, .keep_all = TRUE)
message(sprintf("UniProt: %d entries from UP000008816", nrow(up)))

# ---- 4. Join and write ------------------------------------------------------
meta <- fa |>
  left_join(gff_cds, by = "locus_tag") |>
  left_join(gff_gene, by = "locus_tag") |>
  left_join(up, by = "uniprot_acc") |>
  left_join(terms_wide, by = "protein_key") |>
  left_join(clusters_wide, by = "protein_key") |>
  mutate(
    in_string          = protein_key %in% terms$protein_key | protein_key %in% clusters$protein_key,
    gff_location_agrees = !is.na(start) & !is.na(gff_start) & start == gff_start & end == gff_end,
    # one convenient best-guess gene symbol: NCBI first, then the UniProt ortholog
    gene_symbol_best   = coalesce(gene_fasta, gff_gene, gene_symbol,
                                  up_gene_primary, uniprot_gene_guess),
    # sanity check on the ortholog call: do the two proteins have a similar length?
    up_length_ratio    = if_else(!is.na(up_length) & up_length > 0,
                                 round(aa_length / up_length, 3), NA_real_),
    # Convenience flag so weak ortholog calls can be filtered in one step.
    # "good"  = strong BLAST hit AND the two proteins are of comparable length
    # "check" = hit exists but is weak and/or the lengths disagree
    ortholog_quality = case_when(
      is.na(uniprot_acc)                                   ~ NA_character_,
      blastortho_evalue <= 1e-5 &
        !is.na(up_length_ratio) &
        up_length_ratio >= 0.9 & up_length_ratio <= 1.1    ~ "good",
      TRUE                                                 ~ "check")
  ) |>
  select(
    protein_key, fasta_id, entry_class, locus_tag, gene_symbol_best,
    product, protein_id, aa_length, n_tryptic_peptides,
    uniprot_acc, uniprot_entry_name, uniprot_db, blastortho_evalue, ortholog_quality,
    uniprot_gene_guess, blastortho_raw,
    up_gene_primary, up_gene_names, up_protein_name, up_ec, up_reviewed,
    up_length, up_length_ratio, up_keywords, up_go_ids, up_subcellular, up_kegg,
    up_protein_existence, up_entry_name,
    location_raw, start, end, strand, cds_length_nt, gbkey, is_pseudo,
    gene_fasta,
    gff_product, gff_gene, gene_symbol, gene_name_gff, gene_biotype, gene_feature,
    gff_note, gff_inference, gff_dbxref, gff_transl_table,
    gff_start, gff_end, gff_strand, gff_location_agrees,
    fs_index, fs_direction, fs_frame, fs_locus_tag, fs_note,
    cont_uniprot_acc, cont_entry_name, cont_evalue,
    in_string, starts_with("go_"), starts_with("kegg"), starts_with("uniprot_keywords"),
    starts_with("compartments"), starts_with("reactome"), starts_with("string_"),
    description
  ) |>
  arrange(entry_class, locus_tag, protein_key)

out_wide <- file.path(here, "JE2_protein_metadata.tsv")
out_long <- file.path(here, "JE2_protein_terms.tsv")
write_tsv(meta, out_wide, na = "")
write_tsv(select(terms, protein_key, category, term, term_description) |> arrange(protein_key),
          out_long, na = "")

# ---- 5. Report --------------------------------------------------------------
message("\n--- coverage -------------------------------------------------------")
meta |> group_by(entry_class) |>
  summarise(n = n(),
            with_locus_tag = sum(!is.na(locus_tag)),
            with_uniprot   = sum(!is.na(uniprot_acc)),
            ortholog_good  = sum(ortholog_quality == "good", na.rm = TRUE),
            with_gene_sym  = sum(!is.na(gene_symbol_best)),
            with_gff       = sum(!is.na(gff_product)),
            with_string    = sum(in_string), .groups = "drop") |>
  as.data.frame() |> print(row.names = FALSE)
message(sprintf("\nwrote %s  (%d rows x %d columns)", basename(out_wide), nrow(meta), ncol(meta)))
message(sprintf("wrote %s  (%d rows)", basename(out_long), nrow(terms)))
