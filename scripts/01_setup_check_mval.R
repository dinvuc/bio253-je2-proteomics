# Step 1 of our workflow (mval_medi/project_outlines/PROJECT_WORKFLOW_mval.md):
# establish a reproducible starting point.
#
# What it does, and why:
#   1. records branch + git revision      -> results can be traced to a code version
#   2. runs the course checksum check     -> we only analyse intact course files
#   3. records R, package, Quarto versions -> someone else can reproduce the run
#   4. checks that data/ stays git-ignored -> we never publish the raw data
#   5. runs the helper tests + renders notebook 01 -> the helpers work on our PC
#
# Run from the repository root:
#   Rscript scripts/01_setup_check_mval.R            (all checks, ~2-4 min)
#   Rscript scripts/01_setup_check_mval.R --quick    (skip tests and rendering)
#
# Output: mval_medi/results/step01_setup/setup_record_mval.txt  (small text file, committed)

suppressMessages(source("R/setup.R"))

quick   <- "--quick" %in% commandArgs(trailingOnly = TRUE)
out_dir <- "mval_medi/results/step01_setup"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Small helper: run a shell command and keep its text output (and exit status)
run_cmd <- function(cmd, args = character()) {
  out <- suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(out = out, ok = is.null(status) || status == 0)
}

record <- c("# Step 1 setup record (group mval)",
            paste("Run at:", format(Sys.time(), "%Y-%m-%d %H:%M")), "")

# ---- 1. git provenance -------------------------------------------------------
record <- c(record, "## Git",
            paste("branch:  ", run_cmd("git", c("rev-parse", "--abbrev-ref", "HEAD"))$out),
            paste("revision:", run_cmd("git", c("rev-parse", "--short", "HEAD"))$out),
            paste("upstream:", run_cmd("git", c("remote", "get-url", "upstream"))$out), "")

# ---- 2. checksums of the course data ------------------------------------------
chk <- run_cmd("Rscript", "data/check_data.R")
summary_line <- grep("OK, .* missing", chk$out, value = TRUE)
not_ok <- grep("^(MISSING|BAD|CORRUPT)", chk$out, value = TRUE)
record <- c(record, "## Data checksums (data/check_data.R)", summary_line,
            if (length(not_ok)) paste("  ", not_ok) else "  all files OK", "")

# ---- 3. versions -------------------------------------------------------------
pkgs <- c("tidyverse", "dplyr", "ggplot2", "limma", "readxl", "testthat", "quarto")
pkg_versions <- vapply(pkgs, function(p) {
  if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "NOT INSTALLED"
}, character(1))
quarto_v <- run_cmd("quarto", "--version")
record <- c(record, "## Software versions",
            R.version.string,
            paste("Quarto CLI:", if (quarto_v$ok) quarto_v$out[1] else "NOT FOUND"),
            paste0("  ", pkgs, ": ", pkg_versions), "")

# ---- 4. git exclusions: raw data must stay ignored ---------------------------
# `git check-ignore` exits 0 when the path IS ignored, which is what we want.
probe <- c("data/derived/abundances_total.tsv.gz", "data/raw/total/DE_total.xlsx",
           "data/raw/resources/JE2_proteins.fasta", "analyses/_output/01_explore_data.html")
ignored <- vapply(probe, function(p) run_cmd("git", c("check-ignore", "-q", p))$ok, logical(1))
# Untracked large files outside data/ would be pushed by a careless `git add .`
untracked <- run_cmd("git", c("ls-files", "--others", "--exclude-standard"))$out
big_untracked <- untracked[file.exists(untracked) & file.size(untracked) > 5e6]
record <- c(record, "## Git exclusions",
            paste0("  ", probe, ": ", ifelse(ignored, "ignored (good)", "NOT IGNORED - fix .gitignore!")),
            "Untracked files > 5 MB that `git add .` would pick up:",
            if (length(big_untracked)) paste("  ", big_untracked) else "  none", "")

# ---- 5. helper tests and one rendered notebook --------------------------------
if (!quick) {
  tests <- run_cmd("Rscript", "tests/run_tests.R")
  record <- c(record, "## Helper tests (tests/run_tests.R)",
              paste("exit status OK:", tests$ok),
              paste("  ", tail(tests$out, 6)), "")

  t0 <- Sys.time()
  rend <- run_cmd("quarto", c("render", "analyses/01_explore_data.qmd"))
  secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")))
  record <- c(record, "## Render check (analyses/01_explore_data.qmd)",
              paste0("rendered OK: ", rend$ok, " (", secs, " s); output in analyses/_output/ (git-ignored)"),
              if (!rend$ok) paste("  ", tail(rend$out, 10)), "")
}

writeLines(record, file.path(out_dir, "setup_record_mval.txt"))
cat(record, sep = "\n")
