# Coding rules — group mval

These rules apply to everyone in our group (and to any AI tool helping us).
They come on top of the course rules in `CLAUDE.md` / `AGENTS.md`.

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

Group documents (rules, plans, workflows) go in `mval_medi/`.
The project plan and workflow are in `mval_medi/project_outlines/`.
Our lab notebook (one dated entry per session) is in `mval_medi/notebook/`,
named `YYYY-MM-DD_topic_mval.md` (template: `mval_medi/notebook/README.md`).
