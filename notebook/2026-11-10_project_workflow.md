# 2026-11-10 - Project C workflow planning

## Question for today
What must we do, in what order, to answer our mutant-validation and supernatant-response questions?

## What we did
- Created and switched to branch `docs/project-workflow` from the current `main` revision, preserving existing untracked work.
- Reviewed project instructions, PROJECT_PLAN.md, earlier notebook records, data documentation and the existing R/notebook structure.
- Created PROJECT_WORKFLOW.md with 12 ordered steps, required inputs and their purposes, outputs, decision points and completion checks.
- Included observed/imputed separation, interaction direction, phosphosite coverage/localization, protein adjustment, functional background, literature interpretation and final reporting.

## Results (key numbers)
- One new workflow document with 12 steps covering H1-H3 and conditional supporting analyses.
- Primary design: 18 samples; all 27 considered for sample QC and mutant-profile checks, with SNpulse optional.
- Rscript and Quarto are discoverable on this PC. No package installation, data recheck or biological analysis was performed in this planning session.

## Problems / limitations
- Previous documentation records intact course inputs but does not establish usable quantitative site/protein coverage.
- Predicted combined protein reduction remains a hypothesis to verify; the workflow does not treat it as an established result.
- Existing untracked documents from earlier sessions were preserved. No data were staged or published.

## AI use
Codex drafted the workflow from local project records and adapted analysis guidance, checked comparison names against the agreed plan, and reviewed the documentation changes. No new literature search or biological result was generated.

## Next steps
Follow Step 1 of PROJECT_WORKFLOW.md, then audit samples/QC and plot Stp1/Stk1 measurements before downstream interpretation.
