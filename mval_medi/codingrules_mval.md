# Coding rules — group mval

These rules apply to everyone in our group (and to any AI tool helping us).
They come on top of the course rules in `CLAUDE.md` / `AGENTS.md`.

> **Most important:** the four pitfalls in [section 11](#11-four-pitfalls-to-guard-against).
> Check every analysis against them.

## 1. File names end in `_mval`

Every file we create and save ends in `_mval` just before the file extension.
This applies to every file type: Markdown, R scripts, Quarto notebooks,
figures, tables and reports.

| type | example |
|---|---|
| Markdown | `results_summary_mval.md` |
| R script | `01_stk1_levels_mval.R` |
| Quarto notebook | `stp_vs_anc_mval.qmd` |
| rendered report | `stp_vs_anc_mval.html` |
| figure | `volcano_stp_vs_anc_mval.png` |
| table | `top_hits_stp_vs_anc_mval.csv` |

**Why:** our files stay easy to tell apart from the starter-pack files and
from other groups' files, and we can find all of them with one search
(`*_mval*`).

Exceptions: files that other tools need under a fixed name, for example
`README.md`, `.gitignore` or `R/setup.R`. We do not rename starter-pack
files.

## 2. Read the plan before you start

Before working on a task, read our project documents so the work fits our
question and the agreed order of steps:

- [Project plan](project_outlines/PROJECT_PLAN_mval.md) — Topic C question,
  hypotheses, methods and limitations
- [Project workflow](project_outlines/PROJECT_WORKFLOW_mval.md) — step-by-step
  analysis sequence, inputs and purpose of each step

## 3. Where our files live

Group documents (rules, plans, workflows, docs) go in `mval_medi/`; background
documents such as our literature review are in `mval_medi/docs/`.
The project plan and workflow are in `mval_medi/project_outlines/`.
Our lab notebook (one dated entry per session) is in `mval_medi/notebook/`,
named `YYYY-MM-DD_topic_mval.md` (template: `mval_medi/notebook/notebook_template_mval.md`).

## 4. Plan each step before coding

Before writing code for a workflow step, write down five things (in the
notebook entry or at the top of the script/notebook):

1. **Aim** — which question does this step answer?
2. **Input** — which tables/helpers (e.g. `read_dea("total")`), which columns, which scale (log2 or not)?
3. **Analysis** — which method and parameters, and why?
4. **Output** — which file, which columns, which plot?
5. **Check** — how will we know it is right (a printed count, a plot, a known answer such as Stp1 being absent in `stp`)?

One step = one script or notebook = one commit. A step uses only the inputs
it lists and writes only the outputs it lists. If the code ends up doing
something the plan did not say, explain why.

## 5. Statistics

- **Use adjusted p-values, never raw p-values**, to decide what is a hit. Our
  screen is the course default: **FDR < 0.1 and |log2 effect| > 1**, always
  reported together with the estimate and its standard error. Check which
  column a filter actually uses.
- **Fix thresholds before looking at results.** Write them in the plan or
  notebook first; do not change them afterwards to get a nicer list.
- **Name the test.** For every result say which model it came from
  (`modelName`, e.g. `lm_impute` or the moderated linear model), which
  contrast, and the direction in words (see `CLAUDE.md` rule 3).
- **Our own tests need their own correction.** If we run a new test (t-test,
  limma) on many proteins/sites, adjust with `p.adjust(method = "BH")` across
  all features tested, and state that this FDR is ours, not the supplied one.
  Do not reuse supplied p-values/FDR for a new quantity (e.g. protein-adjusted
  phosphorylation).
- **Keep the full table.** Save results for all tested features, not only the
  hits, so a filter can be checked and redone.
- **Do not re-normalise or re-log** values that are already normalised log2.
- **Batch/covariates:** check known technical factors (e.g. the `reinjected`
  phospho samples) before reading a difference as biology.
- Always mention n = 3 per group and how much of a result rests on imputed
  values or a single replicate.

## 6. Data and identifier safety

- **Join by name, never by position.** Join samples by sample name and
  proteins by locus tag. After a join or reshape, check it, e.g.
  `stopifnot(!anyDuplicated(...))`, the expected row count, 27 (or 18)
  samples, 3 per group.
- **Gene names and functions come from tables, not memory.** Annotate locus
  tags with `read_protein_metadata()` / `read_annotation()`; never type a gene
  symbol or protein function from memory into code or text.
- **Look at big files through code.** Inspect a large table with `head()`,
  `dplyr::glimpse()` or a summary, not by printing the whole file.

## 7. Run it and check the output

- Every script or notebook is run (or rendered) after it is changed. Print the
  number of samples and features after each main step.
- **Running without an error is not the same as being correct.** Compare the
  printed CHECK values with the plan and look at the plots.
- Fix the code, never the data: do not edit files in `data/` or the tests to
  make something pass, and do not hide errors (e.g. an empty `tryCatch`). If a
  step fails and the reason is unclear, stop and report it.
- Figures must use the same thresholds as the analysis they show (e.g.
  volcano lines at FDR 0.1 and |log2| 1).

## 8. Reports

- Results go into Quarto notebooks that render from scratch. Put
  `sessionInfo()` at the end and say which AI tool helped with what.
- A summary written by an AI is a draft: a group member checks every number
  and claim in it before it counts.

## 9. Git

- Work on your own branch. Pull the latest `main` into it before starting a
  session, to avoid conflicts later.
- Make small commits with a clear message saying what changed, e.g.
  `feat: add Stp1/Stk1 validation plot`, `fix: wrong contrast direction`,
  `docs: notebook entry 12 Nov`.
- Run `git status` before every commit and push. Stage only the files you
  meant to change: no `data/`, no `.RData`/`.Rhistory`, no large intermediate
  tables.
- Changes reach `main` through a pull request that another group member has
  read. An AI tool may commit and push a feature branch; it merges or pushes
  to `main` only when a group member explicitly asks.

## 10. Who decides what

- **The group decides** the question, the hypotheses, the thresholds and which
  samples are in or out.
- **The AI tool does the work in between**: writing and running code,
  checking data, making plots, reporting numbers neutrally, and pointing out
  problems (missing values, imputed estimates, odd samples, wrong direction).
- **The AI interprets the biology only when we ask it to.** When it does, it
  marks clearly what is interpretation, gives checkable sources (`CLAUDE.md`
  rule 7) and says how sure it is.

## 11. Four pitfalls to guard against

From the lecture "Agentic Coding" (slide 17, Bio-QC Guardrails), adapted to
our proteomics data. These errors are easy to make, especially by AI tools,
and they produce wrong results without any error message.

### Pitfall 1: the unadjusted p-value trap

**Risk:** filtering by raw `p.value < 0.05` across all tested proteins or
phosphosites gives many false positives.
**Guard:** decide hits on the adjusted value (FDR) only, using our screen
FDR < 0.1 and |log2 effect| > 1 (section 5). In every filter, check which
column is used: it must be the FDR column, not the raw p-value.

### Pitfall 2: silent sample-order mismatch

**Risk:** if the abundance columns and the sample sheet are in a different
order, genotype and medium labels get swapped silently, and the results
look plausible but are wrong.
**Guard:** in every script that combines measurements with sample
information, join by sample name, never by position, and assert the match
before any statistics, e.g.
`stopifnot(setequal(unique(abund$sample), samples$sample))`, plus the
expected count of 3 samples per group.

### Pitfall 3: opening huge files in the AI chat

**Risk:** pasting or printing a large table (site tables, peptide tables,
full DEA results) into an AI chat overloads it and can freeze the session.
**Guard:** large files stay in `data/` (git-ignored). Inspect them through R
with `head()`, `dplyr::glimpse()`, `nrow()` or summaries, and show the AI
only those.

### Pitfall 4: invented gene or protein names

**Risk:** asked what a locus tag (`B7H15_xxxxx`) or UniProt accession is, an
AI may guess a wrong or outdated gene name or function from memory.
**Guard:** never let an AI name proteins from memory. Get names and functions
by a join in code against `read_protein_metadata()` / `read_annotation()`
(keep `ortholog_quality == "good"`), and give a checkable source for any
functional claim (`CLAUDE.md` rule 7).
