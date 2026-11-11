# Upstream synchronization — 11 November 2026

## Question
Is the fork up to date with the original jjGG repository?

## Work performed
Codex compared the GitHub main branches. The fork lacked three upstream commits, ending at 175730a05d2e4c6bdc6afa9c2330890043c80ed7. Merged them through pull request #2 using a merge commit to preserve both histories, including the fork's project_outlines folder and project documentation.

## Results
Imported notebook 07, transposon and peptide-evidence helpers, public protein metadata, accompanying tests and documentation. The upstream changes also update Git exclusions for worked exercise outputs and remove the previously tracked notebook 01 HTML. No local course measurement files were uploaded.

## Validation and problems
The merge had no conflicts. Attempted Rscript tests/run_tests.R on an isolated archive of the incoming revision. Tests could not start because dplyr is not installed in the available R environment; no test pass is claimed.

## Next steps
Refresh local checkouts, install the required R packages, and rerun tests/run_tests.R before using the new helpers for analysis.
