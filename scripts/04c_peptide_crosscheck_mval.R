# Step 4c: cross-check of steps 3-4 against the peptide layer used in notebook 07
# (analyses/07_transposon_mutants.qmd).
#
# Run from the repository root:   Rscript scripts/04c_peptide_crosscheck_mval.R   (~30 s)
# Needs data/derived/peptides_total.tsv.gz and the Step 3 output enzyme_per_sample_mval.csv.
#
# ---- Plan (coding rules, section 4) -----------------------------------------
# Aim:      Does the peptide-level evidence of notebook 07 agree with our Step 3/4 numbers?
#           Which peptides carry the "wild-type-like" Stp1 values in stp?
# Input:    read_peptides()                 one row per peptide x sample (NA = not detected)
#           enzyme_per_sample_mval.csv      Step 3: peptides per sample (n_peptides), abundance
#           TRANSPOSON_INSERTIONS           insertion codon per gene (R/constants.R)
# Analysis: counting only. A peptide is "before" / "after" the insertion if it ends before /
#           starts after the hit codon of THAT gene's own transposon. For Stk1 in stp this split
#           has no biological meaning (stk1 itself is intact in stp), we report it anyway
#           because notebook 07 uses it.
# Output:   mval_medi/results/step04_polar_effect/peptide_crosscheck_mval.csv
#           mval_medi/results/step04_polar_effect/peptide_crosscheck_checks_mval.txt
# Check:    our per-sample peptide counts equal Step 3 n_peptides in every mutant sample.

suppressMessages(source("R/setup.R"))

out_dir <- "mval_medi/results/step04_polar_effect"
step3 <- "mval_medi/results/step03_validation/enzyme_per_sample_mval.csv"
if (!file.exists(step3)) stop("Run scripts/03_mutant_validation_mval.R first.")

hit_codon <- TRANSPOSON_INSERTIONS |>
  transmute(locus_tag, codon = (insertion_pos - cds_start + 1 + 2) %/% 3)
stopifnot(identical(sort(hit_codon$codon), c(46, 189)))      # Stk1 46/664, Stp1 189/247

pep <- read_peptides() |>
  filter(locus_tag %in% hit_codon$locus_tag, !is.na(intensity)) |>   # detected only; NA stays out
  left_join(hit_codon, by = "locus_tag") |>
  mutate(side = case_when(protein_end < codon ~ "before",
                          protein_start > codon ~ "after",
                          TRUE ~ "spans"))

per_sample <- pep |>
  group_by(locus_tag, sample, genotype, medium) |>
  # Count peptide FORMS (modified sequence, e.g. oxidised M counts separately), because that is
  # what Step 3's n_peptides (nr_children_protein_Id) counts. Unique sequences are reported too.
  summarise(n_peptides_pep_layer = n_distinct(peptide),
            n_sequences = n_distinct(peptide_plain),
            n_after_insertion = n_distinct(peptide_plain[side == "after"]),
            peptides = paste(unique(paste0(protein_start, "-", protein_end, " (", side, ")")),
                             collapse = "; "),
            .groups = "drop")

s3 <- read_csv(step3, show_col_types = FALSE) |>
  select(short_name, locus_tag, sample, normalized_abundance, n_peptides_step3 = n_peptides)
cmp <- s3 |>
  left_join(per_sample, by = c("locus_tag", "sample")) |>      # join by name, not position
  mutate(n_peptides_pep_layer = replace_na(n_peptides_pep_layer, 0L),
         n_peptides_step3 = replace_na(n_peptides_step3, 0)) |>
  arrange(short_name, sample)
stopifnot(nrow(cmp) == 54, !anyDuplicated(cmp[c("locus_tag", "sample")]))

mut <- cmp |> filter(!str_detect(sample, "anc"))
agree <- all(mut$n_peptides_pep_layer == mut$n_peptides_step3)
checks <- c(
  sprintf("[%s] all 54 rows (anc too): peptide-layer counts equal Step 3 n_peptides",
          if (all(cmp$n_peptides_pep_layer == cmp$n_peptides_step3)) "PASS" else "FLAG"),
  sprintf("[%s] mutant samples: peptide-layer counts equal Step 3 n_peptides (%d rows)",
          if (agree) "PASS" else "FLAG", nrow(mut)),
  sprintf("[INFO] Stp1 in stp: %d of 9 samples detected; samples with a peptide after codon 189: %s",
          sum(mut$short_name == "Stp1" & str_detect(mut$sample, "stp") & mut$n_peptides_pep_layer > 0),
          paste(mut$sample[mut$short_name == "Stp1" & str_detect(mut$sample, "stp") &
                             replace_na(mut$n_after_insertion, 0) > 0], collapse = ", ")),
  sprintf("[INFO] Stk1 in stk: detections %s",
          paste(mut$sample[mut$short_name == "Stk1" & str_detect(mut$sample, "stk") &
                             mut$n_peptides_pep_layer > 0], collapse = ", "))
)
writeLines(checks, file.path(out_dir, "peptide_crosscheck_checks_mval.txt"))
message(paste(checks, collapse = "\n"))

write_csv(cmp, file.path(out_dir, "peptide_crosscheck_mval.csv"))
print(as.data.frame(mut |> filter(n_peptides_pep_layer > 0, !(short_name == "Stp1" & str_detect(sample, "stk"))) |>
                      select(short_name, sample, normalized_abundance, n_peptides_step3,
                             n_peptides_pep_layer, peptides)))
