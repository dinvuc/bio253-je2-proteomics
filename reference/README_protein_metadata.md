# Protein metadata for the JE2 search database

Built by [`make_protein_metadata.R`](make_protein_metadata.R) from
`p41881_db1_SA_JE2_modNmoreTranslations_20260630.fasta` (db1, 4,155 entries),
the NCBI GFF, the STRING v12 export and the UniProt reference proteome.

Re-run with `Rscript make_protein_metadata.R` (about 1 minute; needs network
only on the first run, to cache the UniProt table).

## Files

| file | rows | what |
|---|---|---|
| `JE2_protein_metadata.tsv` | 4,155 (one per FASTA entry) | the wide table, 80 columns |
| `JE2_protein_terms.tsv` | 121,019 | long format: `protein_key`, `category`, `term`, `term_description` — easier to join/filter than the collapsed columns |
| `uniprot_UP000008816.tsv.gz` | 2,889 | cached UniProt download; delete it to refresh |

`protein_key` is the join key used everywhere else in this project: it is the
middle field of the FASTA id (`tr|KEY|KEY`) and matches the `protein_Id` column
of the prolfquapp Excel files after stripping `tr|`.

## The four entry classes

| `entry_class` | n | what it is |
|---|---|---|
| `cds` | 2,938 | real coding sequences, fully annotated |
| `pseudo_frameshift_translation` | 720 | 6-frame translations of pseudogenes / frameshifted ORFs (one gene gives up to 6 entries) |
| `contaminant` | 496 | FGCZ contaminant list (`zz_...`) |
| `database_info` | 1 | the database description line itself |

## Column groups

**Identity** — `protein_key`, `fasta_id`, `entry_class`, `locus_tag`,
`gene_symbol_best`, `product`, `protein_id`, `aa_length`, `n_tryptic_peptides`
(fully tryptic 6–30 aa peptides, for length-normalised abundances).

**UniProt (what you asked for)** — `uniprot_acc`, `uniprot_entry_name`,
`uniprot_db` (`sp` = reviewed, `tr` = unreviewed), **`blastortho_evalue`**,
`ortholog_quality`, `blastortho_raw`. These come from the `BLASTORTHO` field of
the FASTA description, i.e. the best BLAST hit in the *S. aureus* NCTC 8325
reference proteome (UP000008816).

**UniProt annotation** (joined on `uniprot_acc`) — `up_gene_primary`,
`up_gene_names`, `up_protein_name`, `up_ec`, `up_reviewed`, `up_length`,
`up_length_ratio`, `up_keywords`, `up_go_ids`, `up_subcellular`, `up_kegg`,
`up_protein_existence`, `up_entry_name`.

**Position / NCBI** — `location_raw`, `start`, `end`, `strand`,
`cds_length_nt`, `gbkey`, `is_pseudo`, and from the GFF: `gff_product`,
`gff_gene`, `gene_symbol`, `gene_biotype`, `gene_feature`, `gff_note`,
`gff_inference`, `gff_dbxref`, `gff_transl_table`, `gff_start`, `gff_end`,
`gff_strand`, `gff_location_agrees`.

**Frameshift translations** — `fs_index`, `fs_direction` (forward/reverse),
`fs_frame` (1–3), `fs_locus_tag`, `fs_note` (e.g. "frameshifted; internal stop").
All 720 inherit the locus tag of their gene, so they join to the GFF.

**Contaminants** — `cont_uniprot_acc`, `cont_entry_name`, `cont_evalue`.

**STRING v12** — `in_string`, then per category the collapsed terms as
`id=description | id=description` plus a count: `go_bp`/`go_bp_n`,
`go_cc`, `go_mf`, `kegg`, `uniprot_keywords`, `compartments`, `reactome`,
`string_local_cluster`, and `string_clusters` / `string_cluster_descriptions` /
`string_n_clusters`.

## Coverage

| entry class | n | locus tag | UniProt acc | ortholog "good" | gene symbol | GFF | STRING |
|---|---|---|---|---|---|---|---|
| cds | 2,938 | 2,938 | 2,886 | 2,644 | 987 | 2,938 | 2,664 |
| pseudo_frameshift | 720 | 720 | – | – | – | 720 | 119 |
| contaminant | 496 | – | – | – | – | – | 453 |

## Things to watch

* **Filter on the e-value.** 102 CDS orthologs have `blastortho_evalue > 1e-5`,
  and the worst is ~1 — those are not orthologs in any useful sense. 52 CDS have
  no `BLASTORTHO` field at all (mostly "hypothetical protein"). Use
  `ortholog_quality == "good"` for a quick, defensible subset: it requires
  e-value ≤ 1e-5 **and** a length ratio between 0.9 and 1.1.
* **`up_length_ratio`** compares the JE2 protein with its UniProt ortholog.
  91 CDS differ by more than two-fold — worth a look before trusting the
  transferred annotation.
* **Gene symbols are genuinely scarce.** The FASTA has 13, the GFF 24; joining
  UniProt raises `gene_symbol_best` to 987 of 2,938. For the rest, the locus tag
  is the only identifier. Do not expect `gene_name` to be populated.
* **Frameshift entries multiply genes.** 720 translations come from far fewer
  genes (up to 6 frames each), so counting rows is not counting proteins.
* **STRING ids** are the `protein_key` prefixed with `STRG0A35SHG.`. STRING
  covers the contaminants too (453), which is usually not what you want in an
  enrichment background — filter on `entry_class == "cds"` first.
* **Joining to the DEA tables covers every real protein.** The Oct-26 total
  proteome DEA has 2,324 entries: 2,300 `tr|...` proteins, **all** of which
  match this table, plus 24 `contam_sp|...` contaminants that DIA-NN adds at
  search time and that are *not* in db1 (db1 brings its own `zz_...`
  contaminant set instead). Those 24 are flagged `CON = TRUE` in the Excel, so
  drop them before joining and nothing is lost.
* Percent-encoded GFF text (`%3B`, `%2C`) is decoded; `product` is kept as NCBI
  wrote it, while `up_protein_name` is usually the more informative name
  (e.g. "Serine/threonine-protein kinase PrkC (EC 2.7.11.1)" vs
  "serine/threonine protein kinase").

## Quick use

```r
library(readr); library(dplyr)
meta <- read_tsv("JE2_protein_metadata.tsv", show_col_types = FALSE)

# DEA results -> UniProt accession for downstream tools
dea |>
  mutate(protein_key = sub("^tr\\|([^|]+)\\|.*$", "\\1", protein_Id)) |>
  left_join(select(meta, protein_key, uniprot_acc, blastortho_evalue,
                   ortholog_quality, gene_symbol_best, up_protein_name),
            by = "protein_key")

# a clean enrichment background
bg <- meta |> filter(entry_class == "cds", ortholog_quality == "good")

# all GO BP terms of one protein, long format
read_tsv("JE2_protein_terms.tsv", show_col_types = FALSE) |>
  filter(protein_key == "lcl_CP020619.1_prot_ARG45703.1_1243",
         category == "Biological Process (Gene Ontology)")
```
