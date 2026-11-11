# Tests on the real data. Skipped if scripts/00_fetch_data.R was not run yet.

has_data <- function() file.exists(derived_path("dea_total.tsv.gz"))

test_that("the sign convention holds: diff of A_vs_B = mean(A) - mean(B)", {
  skip_if_not(has_data(), "data not fetched")
  tcyA <- "B7H15_13375"
  means <- read_abundances("total") |>
    filter(locus_tag == tcyA, group %in% c("anc_TSB", "anc_SN")) |>
    group_by(group) |>
    summarise(m = mean(normalized_abundance), n = sum(!is.na(normalized_abundance)))
  expect_equal(means$n, c(3L, 3L))                    # all observed, no imputation
  own_diff <- means$m[means$group == "anc_TSB"] - means$m[means$group == "anc_SN"]

  excel_diff <- read_dea("total", contrasts = "TSB_vs_SN_at_anc") |>
    filter(locus_tag == tcyA) |>
    pull(diff)
  expect_equal(own_diff, excel_diff, tolerance = 0.01)
  expect_lt(excel_diff, 0)                            # TcyA is higher in SN in the ancestor
})

test_that("abundance tables have 27 samples with the expected design", {
  skip_if_not(has_data(), "data not fetched")
  samples <- read_samples()
  expect_equal(nrow(samples), 27)
  expect_equal(as.vector(table(samples$group)), rep(3L, 9))
  ab <- read_abundances("total")
  expect_false(any(str_detect(ab$sample, "^2anc|^2stk|^2stp")))   # truncated names repaired
  expect_setequal(unique(ab$sample), samples$sample)
})

test_that("site tables: zeros became NA, residues and positions parsed", {
  skip_if_not(has_data(), "data not fetched")
  styh <- read_sites("STYH")
  expect_false(any(styh$intensity == 0, na.rm = TRUE))
  expect_setequal(unique(styh$residue), c("S", "T", "Y", "H"))
  expect_false(anyNA(styh$position))
  expect_true(all(is.na(styh$loc_prob[is.na(styh$intensity)])))
  sty <- read_sites("STY")
  expect_false("H" %in% sty$residue)
  expect_equal(n_distinct(styh$sample), 27)
})

test_that("find_protein finds the Tcy proteins by locus tag and protein id", {
  skip_if_not(has_data(), "data not fetched")
  expect_equal(nrow(find_protein("B7H15_13375")), 1)
  expect_equal(find_protein("ARG46976")$locus_tag, "B7H15_13375")
  expect_equal(nrow(find_protein(c("B7H15_13375", "B7H15_13370", "B7H15_13365"))), 3)
})

test_that("transposon insertion coordinates land where the NTML map says", {
  ins <- TRANSPOSON_INSERTIONS |>
    mutate(cds_length  = cds_end - cds_start + 1,
           protein_aa  = cds_length %/% 3 - 1,
           nt_into_cds = insertion_pos - cds_start + 1,
           codon       = (nt_into_cds + 2) %/% 3)
  # NE1919 (stp1): 567 nt into a 744 nt CDS -> codon 189 of a 247 aa protein
  stp <- filter(ins, short_name == "Stp1")
  expect_equal(stp$nt_into_cds, 567)
  expect_equal(stp$codon, 189)
  expect_equal(stp$protein_aa, 247)
  # NE217 (stk1): 137 nt into a 1995 nt CDS -> codon 46 of a 664 aa protein
  stk <- filter(ins, short_name == "Stk1")
  expect_equal(stk$nt_into_cds, 137)
  expect_equal(stk$codon, 46)
  expect_equal(stk$protein_aa, 664)
  # stp1 and stk1 overlap by 4 bp: one operon, hence the polar effect
  expect_equal(stk$cds_start - stp$cds_end - 1, -4)
})

test_that("peptide table resolves the transposon truncation of Stp1", {
  skip_if_not(file.exists(derived_path("peptides_total.tsv.gz")), "peptides not fetched")
  pep <- read_peptides()
  expect_equal(n_distinct(pep$sample), 27)
  stp1 <- filter(pep, locus_tag == "B7H15_06350")
  after <- filter(stp1, protein_start > 189, !is.na(intensity))
  # Peptides C-terminal of the insertion are seen all over anc and stk, but
  # essentially never in stp. Not exactly zero: one peptide shows up in a
  # single stp run, which is why the notebook discusses run-to-run spill-over
  # instead of claiming a clean absence.
  expect_gt(sum(after$genotype == "anc"), 10)
  expect_gt(sum(after$genotype == "stk"), 10)
  expect_lte(sum(after$genotype == "stp"), 2)
  # The two best-covered downstream peptides are absent from stp altogether.
  far <- filter(stp1, !is.na(intensity),
                peptide_plain %in% c("EGTIEDHGDQLMQLALDNHSK", "DNVTFILAAIEGDKV"))
  expect_gt(sum(far$genotype == "anc"), 10)
  expect_equal(sum(far$genotype == "stp"), 0)
})

test_that("protein metadata ships with the repo and maps to UniProt", {
  meta <- read_protein_metadata()
  expect_equal(nrow(meta), 4155L)
  expect_false(any(duplicated(meta$protein_key)))
  expect_setequal(unique(meta$entry_class),
                  c("cds", "pseudo_frameshift_translation", "contaminant", "database_info"))
  cds <- read_protein_metadata(cds_only = TRUE)
  expect_equal(nrow(cds), 2938L)
  expect_gt(sum(!is.na(cds$uniprot_acc)), 2800)
  # Stk1: the BLASTORTHO ortholog of B7H15_06355 in NCTC 8325
  stk1 <- filter(cds, locus_tag == "B7H15_06355")
  expect_equal(stk1$uniprot_acc, "Q2FZ64")
  expect_equal(stk1$blastortho_evalue, 0)
  expect_equal(stk1$ortholog_quality, "good")
})
