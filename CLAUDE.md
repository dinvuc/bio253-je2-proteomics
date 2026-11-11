# Instructions for AI coding agents (Bio253 JE2 starter pack)

You are helping a small group of biology students (little programming or
statistics background) during a one-week block course. Teach while you help:
use small steps and short code, and explain *why*.

> **Before every task, read our group rules in
> [mval_medi/codingrules_mval.md](mval_medi/codingrules_mval.md) and follow
> them** (e.g. every new file name ends in `_mval`). They also link to our
> project plan and workflow.

## Hard rules

1. **Work on the real data; publish code and reports, never the raw data.** The
   lecturer allows any AI tool on the course data: read, analyse and summarise the
   files in `data/` directly (prefer running scripts over printing large tables).
   At the end the group pushes its **code and rendered reports** to its own GitHub
   repository (a clone of the starter repo) — help with that when asked. Never
   `git add` anything under `data/` (raw Excel/TSV/parquet/FASTA files); keep
   `.gitignore` intact.
2. **Use R** (tidyverse, ggplot2, limma). Do not switch to Python. Start scripts
   and notebooks with `source("R/setup.R")` (or `source("../R/setup.R")` from
   `analyses/`) and reuse the helpers in `R/`.
3. **Contrast direction:** a contrast `A_vs_B` has `diff = log2(A) - log2(B)`.
   `TSB_vs_SN_at_anc > 0` means **higher in TSB = lower in PASN (SN)**. Every
   time you report a fold change, say in words which group is higher. Use
   `explain_contrast()`.
4. **Plot raw per-sample values before trusting a fold change**
   (`plot_protein()`, `plot_sites()`). Point out missing replicates and outliers.
5. **Imputed vs observed:** in DEA tables `estimate_type == "lod_imputed"` means
   at least one group mean was not measured but imputed from the limit of
   detection. Say so explicitly. Treat it as presence/absence evidence, not as a
   precise fold change.
6. **Missing = NA, never 0.** Do not replace NA by 0 in log space.
7. **Cite sources.** For biological claims (protein function, pathway, literature)
   give a checkable source (UniProt accession, PubMed ID/DOI, the course papers).
   If you are not sure, say so. Do not invent references.
8. **Lab notebook:** after each working session, help the students add a dated
   entry in `mval_medi/notebook/` (template in `mval_medi/notebook/notebook_template_mval.md`): question, what was run,
   results (numbers), problems, next steps, which AI help was used.

## Experimental design

- *S. aureus* JE2. Genotypes: `anc` (wild type), `stk` (stk1::Tn, kinase
  Stk1/PknB), `stp` (stp1::Tn, phosphatase Stp1).
- Media: `TSB` (control); `SN` = ~15 h in 70 % TSB + 30 % cell-free,
  sterile-filtered spent medium of *P. aeruginosa* PAO1 ("PASN"); `SNpulse` = grown
  in TSB, then 15 min PASN. *S. aureus* never met living *P. aeruginosa*. PASN
  contains only secreted molecules (HQNO, PQS, pyocyanin, siderophores, proteases …),
  and it is partly nutrient-depleted. Do not call it co-culture, and do not
  say "in the presence of Pseudomonas".
- 3 replicates × 9 groups = 27 samples; total proteome and phospho-enriched (Fe-IMAC).
- Raw file names say `pulseSN`; tables/contrasts say `SNpulse`. Truncated sample names
  like `2anc_pulseSN_1` are already repaired to `JE2anc_pulseSN_1` in `data/derived/`.
- Phospho samples `stk_pulseSN_1/2` were re-injected on a different day (possible batch effect;
  column `reinjected`).
- The `stp` strain also has strongly reduced Stk1 (a polar effect). Remind
  students of this when they interpret `stp_vs_anc`, `stp_vs_stk` or any
  `interaction_stp_*`: that strain lacks the phosphatase AND most of the kinase.
- Both mutants are NTML strains carrying the 3,237 bp `bursa aurealis`
  transposon (markers `ermB` + `gfp`). NE217 hits *stk1* at codon 46 of 664;
  NE1919 hits *stp1* at codon 189 of 247 (76 % into the gene). *stp1* and *stk1*
  overlap by 4 bp, i.e. one operon. See `TRANSPOSON_INSERTIONS` in
  `R/constants.R` and notebook 07.
- `ermB` is NOT in the search FASTA, so "ErmB not found" means "never looked
  for", not "not there". See `reference/README.md`.

## Data layout (local only, git-ignored)

```
data/raw/        copies of the course files (Excel DEA results, FragPipe site tables, FASTA, GFF, STRING)
data/derived/    tidy long tables written by scripts/00_fetch_data.R
data/reference/  lecturer reference tables (e.g. pHis tiers), if provided
```

Helpers (in `R/`; all return long tibbles with locus_tag, product, protein_id
and the NCTC8325 UniProt ortholog `staa8_uniprot`):

| function | returns |
|---|---|
| `read_dea("total" / "phospho_STYH" / "phospho_STY")` | prolfquapp DEA, one row per protein/site × contrast (33 contrasts) |
| `read_abundances(...)` | normalised log2 abundances per sample (`normalized_abundance`, NA = not detected) |
| `read_sites("STYH" / "STY")` | FragPipe combined_site table, long; `intensity` NA = not detected; `best_loc_prob`, `residue`, `position`, `peptide` |
| `read_peptides()` | peptide-level evidence behind the total proteome: one row per peptide x sample, with `protein_start` / `protein_end`. Use it to ask WHICH PART of a protein was seen (notebook 07) |
| `read_protein_metadata(cds_only=)` | the big annotation table in `reference/` (ships with the repo, no course data needed): `uniprot_acc` + `blastortho_evalue` + `ortholog_quality`, `gene_symbol_best`, `up_protein_name`, `up_ec`, coordinates, GO/KEGG/STRING |
| `read_protein_terms()` | the same annotation long: protein x term |
| `read_samples()`, `read_annotation()` | sample sheet; protein annotation from the FASTA |
| `find_protein(query)` | search locus tag / product / protein_id / UniProt |
| `plot_protein()`, `plot_sites()`, `plot_peptide_coverage()`, `volcano()`, `count_hits()`, `explain_contrast()` | plotting and summaries |

Proteins are identified by **locus tag** (`B7H15_xxxxx`, JE2 genome CP020619.1);
`gene_name` is empty in the DEA tables. For a UniProt accession, a gene symbol
or GO/KEGG terms use `read_protein_metadata()`; filter weak ortholog calls with
`ortholog_quality == "good"`. Key proteins are listed in `KEY_PROTEINS` (`R/constants.R`).
The STY phospho DEA is an older run (moderated linear model, no imputation). The
STYH DEA uses `lm_impute`, like the total proteome.

## Coding style

- tidyverse, native pipe `|>`, short functions, comments that explain *why*.
- Quarto notebooks in `analyses/` must render on their own in < 1 min
  (`quarto render analyses/<file>.qmd`).
- Run the tests after you change `R/`: `Rscript tests/run_tests.R`.
- Prefer a simple, transparent method (t-test, limma, counting) over a black box.
  Mention multiple testing and n = 3.
