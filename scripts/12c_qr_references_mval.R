# Step 12c: QR code from the poster to our list of references.
#
# Run from the repository root:   Rscript scripts/12c_qr_references_mval.R   (a few seconds)
# Needs once:                     install.packages("qrcode")
#
# ---- Plan (coding rules, section 4) -----------------------------------------
# Aim:      Poster readers can open the references with a phone.
# Input:    none (only the address of mval_medi/docs/references_mval.md on GitHub).
# Analysis: none. The link uses main, because feature branches get deleted and the printed
#           code would then break. The list can be edited later without a new code.
#           DO NOT rename or move the references file after printing.
# Output:   mval_medi/results/step12_poster/poster/qr_references_mval.svg
# Check:    the linked file exists; the white margin is 4 modules on every side
#           (QR standard; needed on the grey poster background). Scan a printout before printing A0.

suppressMessages(source("R/setup.R"))

page <- "mval_medi/docs/references_mval.md"
url <- paste0("https://github.com/dinvuc/bio253-je2-proteomics/blob/main/", page)
out <- "mval_medi/results/step12_poster/poster/qr_references_mval.svg"
margin <- 4
stopifnot("references page is missing" = file.exists(page))

# level "M" still scans if a small part of the print is damaged
code <- qrcode::qr_code(url, ecl = "M")

# qrcode's own margin is too narrow: cut it off, then add a margin of 4 white modules
dark <- which(code, arr.ind = TRUE)
symbol <- unclass(code)[min(dark[, 1]):max(dark[, 1]), min(dark[, 2]):max(dark[, 2])]
n <- nrow(symbol)
padded <- matrix(FALSE, n + 2 * margin, n + 2 * margin)
padded[margin + seq_len(n), margin + seq_len(n)] <- symbol
attr(padded, "string") <- url
class(padded) <- class(code)

dark <- which(padded, arr.ind = TRUE)
stopifnot(
  "margin is not 4 modules" = all(range(dark[, 1]) == c(margin + 1, margin + n)),
  "margin is not 4 modules" = all(range(dark[, 2]) == c(margin + 1, margin + n))
)

qrcode::generate_svg(padded, filename = out, show = FALSE)
cat("CHECK symbol:", n, "modules, margin:", margin, "modules\n")
cat("CHECK link:  ", url, "\n")
cat("Written:     ", out, "\n")
