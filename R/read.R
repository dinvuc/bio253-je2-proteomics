# read.R -- load the tidy tables written by scripts/00_fetch_data.R.
#
# All tables are "long": one row per measurement. Every read_* function
# attaches the protein annotation (locus_tag, product, protein_id, STAA8
# ortholog), so you can always see WHICH protein a row is about.

#' Path to a file in data/derived/
derived_path <- function(file) file.path(REPO_ROOT, "data", "derived", file)

# Fail early with a helpful message if the data were not fetched yet.
check_derived <- function(file) {
  path <- derived_path(file)
  if (!file.exists(path)) {
    stop("Missing ", path, "\n",
         "Run the data set-up first (see data/README.md):\n",
         "  Rscript scripts/00_fetch_data.R /path/to/SAJE2", call. = FALSE)
  }
  path
}

#' Protein annotation: one row per protein in the search database.
read_annotation <- function() {
  read_tsv(check_derived("protein_annotation.tsv.gz"), show_col_types = FALSE)
}

DATASETS <- c("total", "phospho_STYH", "phospho_STY")

# The columns of the annotation we attach to every table.
add_annotation <- function(x) {
  annot <- read_annotation() |>
    select(protein_key, locus_tag, product, protein_id, staa8_uniprot, staa8_gene)
  left_join(x, annot, by = "protein_key") |>
    relocate(locus_tag, product, protein_id, .after = protein_key)
}

#' Differential expression results (prolfquapp), long format.
#'
#' @param dataset "total" (protein level), "phospho_STYH" or "phospho_STY"
#'   (site level).
#' @param contrasts optional character vector to keep only some contrasts.
#' @param drop_contaminants remove contaminant proteins (keratins, trypsin...).
#'
#' IMPORTANT -- direction of the fold change:
#'   contrast "A_vs_B" -> diff = log2(A) - log2(B)
#'   e.g. TSB_vs_SN_at_anc > 0 means HIGHER in TSB, i.e. LOWER in PASN (SN).
#'   Use explain_contrast("TSB_vs_SN_at_anc") if in doubt.
read_dea <- function(dataset = "total", contrasts = NULL, drop_contaminants = TRUE) {
  dataset <- match.arg(dataset, DATASETS)
  x <- read_tsv(check_derived(paste0("dea_", dataset, ".tsv.gz")), show_col_types = FALSE)
  if (!is.null(contrasts)) {
    unknown <- setdiff(contrasts, unique(x$contrast))
    if (length(unknown) > 0) stop("Unknown contrast(s): ", paste(unknown, collapse = ", "))
    x <- filter(x, contrast %in% contrasts)
  }
  if (drop_contaminants) x <- filter(x, !is_contaminant)
  add_annotation(x)
}

# Turn the sample columns into ordered factors (tsv files forget factor levels).
as_design_factors <- function(x) {
  mutate(x,
         genotype = factor(genotype, levels = GENOTYPES),
         medium   = factor(medium, levels = MEDIA),
         group    = factor(group, levels = GROUPS))
}

#' Per-sample abundances, long format (one row per protein/site x sample).
#'
#' Column `normalized_abundance` is the log2-scale, vsn-normalised value that
#' the DEA used. NA = not detected in that sample (NOT zero!).
read_abundances <- function(dataset = "total", drop_contaminants = TRUE) {
  dataset <- match.arg(dataset, DATASETS)
  x <- read_tsv(check_derived(paste0("abundances_", dataset, ".tsv.gz")),
                show_col_types = FALSE)
  if (drop_contaminants) x <- filter(x, !is_contaminant)
  add_annotation(as_design_factors(x))
}

#' Phospho sites straight from FragPipe (combined_site_*.tsv), long format.
#'
#' @param search "STYH" (phospho allowed on S, T, Y and H) or "STY".
#' Columns: residue (S/T/Y/H), position in protein, peptide, best_loc_prob
#' (best localisation probability over all samples), per sample: intensity
#' (NA = not detected), log2_intensity, loc_prob.
read_sites <- function(search = "STYH", drop_contaminants = TRUE) {
  search <- match.arg(search, c("STYH", "STY"))
  x <- read_tsv(check_derived(paste0("sites_", search, ".tsv.gz")), show_col_types = FALSE)
  if (drop_contaminants) x <- filter(x, !is_contaminant)
  add_annotation(as_design_factors(x))
}

#' Peptide-level evidence for the total proteome (one row per peptide x sample).
#'
#' This is the table the DEA was built from, BEFORE peptides were rolled up
#' into protein values. Use it when you need to know WHICH PART of a protein
#' was actually seen -- for example on either side of a transposon insertion
#' (notebook 07).
#'
#' Columns: peptide, protein_start / protein_end (position of the peptide in
#' the protein), intensity (NA = not detected in that sample), plus the usual
#' sample / genotype / medium / replicate columns.
#' Note that `peptide` still carries modification tags like "M[15.9949]"; the
#' column `peptide_plain` has them stripped.
read_peptides <- function(drop_contaminants = TRUE) {
  x <- read_tsv(check_derived("peptides_total.tsv.gz"), show_col_types = FALSE)
  if (drop_contaminants) x <- filter(x, !is_contaminant)
  add_annotation(as_design_factors(x))
}

#' Metadata for every entry of the search database (4,155 rows x 80 columns).
#'
#' Unlike the other read_* functions this one does NOT need the course data:
#' it ships with the repository in reference/, because it is built entirely
#' from public sources (the NCBI genome annotation, UniProt UP000008816 and
#' STRING v12) plus the FASTA description lines.
#'
#' Useful columns: `uniprot_acc` + `blastortho_evalue` + `ortholog_quality`
#' (for mapping to UniProt-based tools), `gene_symbol_best`, `up_protein_name`,
#' `up_ec`, `start`/`end`/`strand`, and the collapsed GO / KEGG / STRING
#' annotation. See reference/README_protein_metadata.md for all 80.
#'
#' @param cds_only keep only real coding sequences (drops the 6-frame
#'   pseudogene translations and the contaminant entries).
read_protein_metadata <- function(cds_only = FALSE) {
  x <- read_tsv(file.path(REPO_ROOT, "reference", "JE2_protein_metadata.tsv.gz"),
                show_col_types = FALSE)
  if (cds_only) x <- filter(x, entry_class == "cds")
  x
}

#' The STRING / GO / KEGG annotation in long format: one row per protein x term.
#' Easier to filter than the collapsed columns of read_protein_metadata().
read_protein_terms <- function() {
  read_tsv(file.path(REPO_ROOT, "reference", "JE2_protein_terms.tsv.gz"),
           show_col_types = FALSE)
}

#' Sample sheet: one row per sample with genotype, medium, replicate.
read_samples <- function() {
  as_design_factors(read_tsv(check_derived("samples.tsv"), show_col_types = FALSE))
}
