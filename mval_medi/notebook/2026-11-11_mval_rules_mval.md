# Group rules and file organisation — 11 November 2026

## Question
Give our group one place for its own files and a naming rule that makes them
easy to tell apart from the starter pack.

## Work performed
- Created `mval_medi/codingrules_mval.md`: every file we save ends in `_mval`
  before the extension (exceptions: files that need a fixed name, such as
  README.md). The rules link to the project plan and workflow.
- Added a note to CLAUDE.md and AGENTS.md telling every AI agent to read the
  group rules before each task.
- Moved `project_outlines/` and `notebook/` into `mval_medi/`.
- Renamed our own files to end in `_mval`: project plan, project workflow,
  `mval_medi/docs/WORKING_ON_TWO_PCS_mval.md`, `mval_medi/docs/stk_stp_rescue_evidence_mval.md`
  and the four earlier notebook entries. Updated every link and path reference.

## Results and checks
6 group documents and 4 notebook entries renamed, 2 folders moved. Searched the
repository for the old names and paths; none remain outside git history. No
analysis code or course data changed, so R tests were not needed.

## Problems and next steps
- Starter-pack files and Quarto renders in `analyses/_output/` keep their
  names (a re-render would bring back the original names).
- Continue with Step 1 of the project workflow.

## AI help
Claude Code (Opus 5.5) created the rules file, moved/renamed files, updated
links and wrote this entry.

## Later correction (same day)
- Moved the lecturer's `notebook/README.md` back to `notebook/` unchanged, and
  restored `analyses/01_explore_data.qmd` to the lecturer's version.
- Created our own template: `mval_medi/notebook/notebook_template_mval.md`.
- Moved our two background documents from `docs/` into `mval_medi/docs/`.
