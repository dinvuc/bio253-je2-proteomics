#!/usr/bin/env Rscript
# make_transposon_fasta.R -- write the protein sequences carried by the
# Nebraska-library transposon (and by the two delivery plasmids) so they can be
# appended to the JE2 search database.
#
# Sources, fetched live from NCBI:
#   AY672109  pBursa  (7,383 bp) -- carries the bursa aurealis transposon
#   AY672108  pFA545  (10,079 bp) -- the mariner transposase helper plasmid
#
# GenBank annotates the transposon itself as mobile_element 4147..7383 of
# pBursa. Only TWO coding sequences lie inside it -- ermB and gfp. Everything
# else sits on the plasmid backbones, which are cured at 45 C, so those
# proteins should be ABSENT; they are included as a control for plasmid loss.
#
# Run:  Rscript make_transposon_fasta.R
# Out:  transposon_proteins.fasta

suppressPackageStartupMessages({ library(stringr) })

here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))[1]))
if (is.na(here) || !nzchar(here)) here <- "."
out  <- file.path(here, "transposon_proteins.fasta")

efetch <- function(acc) {
  f <- file.path(here, paste0(acc, ".gb"))          # cached
  if (!file.exists(f)) {
    url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                  "?db=nuccore&id=", acc, "&rettype=gbwithparts&retmode=text")
    utils::download.file(url, f, quiet = TRUE)
  }
  readLines(f, warn = FALSE)
}

# Pull every CDS (location, /gene, /product, /translation) out of a GenBank file.
parse_cds <- function(lines) {
  txt  <- paste(lines, collapse = "\n")
  feat <- str_match(txt, "(?s)\nFEATURES(.*?)\nORIGIN")[, 2]
  blocks <- str_split(feat, "\n(?=     \\S)")[[1]]
  res <- lapply(blocks, function(b) {
    m <- str_match(b, "^\\s+CDS\\s+(\\S+)")
    if (is.na(m[1])) return(NULL)
    flat <- str_replace_all(b, "\n {21}", "")
    g <- function(k) { v <- str_match(flat, paste0("/", k, '="([^"]*)"'))[, 2]; if (is.na(v)) NA_character_ else v }
    data.frame(location = m[, 2], gene = g("gene"), product = g("product"),
               translation = str_remove_all(g("translation"), "\\s"), stringsAsFactors = FALSE)
  })
  do.call(rbind, Filter(Negate(is.null), res))
}

pbursa <- parse_cds(efetch("AY672109"))
pfa545 <- parse_cds(efetch("AY672108"))

# the transposon region of pBursa, as GenBank annotates it
TN <- c(4147, 7383)
coord <- function(loc) as.integer(str_match(loc, "(\\d+)\\.\\.(\\d+)")[, 2:3])
inside <- vapply(pbursa$location, function(l) { p <- coord(l); !any(is.na(p)) && p[1] >= TN[1] && p[2] <= TN[2] }, logical(1))

entries <- rbind(
  data.frame(pbursa[inside, ],  origin = "pBursa AY672109", region = "ON THE TRANSPOSON"),
  data.frame(pbursa[!inside, ], origin = "pBursa AY672109", region = "plasmid backbone (should be cured)"),
  data.frame(pfa545,            origin = "pFA545 AY672108", region = "helper plasmid (should be cured)")
)
entries$key <- paste0(ifelse(grepl("pBursa", entries$origin), "pBursa_", "pFA545_"),
                      ifelse(is.na(entries$gene) | entries$gene == "", "orf", entries$gene))
entries$key <- make.unique(entries$key, sep = "_")

# fully tryptic 6-30 aa peptides: a rough measure of how detectable each is
n_tryptic <- function(s) {
  pep <- str_split(str_replace_all(s, "([KR])(?!P)", "\\1\n"), "\n")[[1]]
  sum(nchar(pep) >= 6 & nchar(pep) <= 30)
}

con <- file(out, "w")
for (i in seq_len(nrow(entries))) {
  e <- entries[i, ]
  writeLines(sprintf(
    ">tr|%s|%s FGCZ [protein=%s] [gene=%s] [location=%s] [gbkey=CDS] [note=%s; %s] [source=GenBank %s]",
    e$key, e$key, e$product, ifelse(is.na(e$gene), "-", e$gene), e$location, e$region,
    "Nebraska transposon mutant library construct", e$origin), con)
  writeLines(str_extract_all(e$translation, ".{1,60}")[[1]], con)
}
close(con)

entries$aa <- nchar(entries$translation)
entries$tryptic_6_30 <- vapply(entries$translation, n_tryptic, integer(1))
print(entries[, c("key", "gene", "product", "aa", "tryptic_6_30", "region")], row.names = FALSE)
message(sprintf("\nwrote %s  (%d sequences)", basename(out), nrow(entries)))
