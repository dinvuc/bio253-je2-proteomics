# Step 12d: QR code from the poster to our list of references.
#
# Run from the repository root:   Rscript C_mutant_validation_DinMete/12d_qr_references_mval.R   (a few seconds)
# Needs once:                     install.packages("qrcode")
#
# ---- Plan (coding rules, section 4) -----------------------------------------
# Aim:      Poster readers can open the references with a phone.
# Input:    none (only the address of C_mutant_validation_DinMete/docs/references_mval.md on GitHub).
# Analysis: none. The link uses main, because feature branches get deleted and the printed
#           code would then break. The list can be edited later without a new code.
#           DO NOT rename or move the references file after printing.
# Output:   C_mutant_validation_DinMete/results/step12_poster/poster/qr_references_mval.{svg,png}
# Check:    the linked file exists; the white margin is 4 modules on every side
#           (QR standard; needed on the grey poster background). Scan a printout before printing A0.

suppressMessages(source("R/setup.R"))

page <- "C_mutant_validation_DinMete/docs/references_mval.md"
url <- paste0("https://github.com/dinvuc/bio253-je2-proteomics/blob/main/", page)
out <- "C_mutant_validation_DinMete/results/step12_poster/poster/qr_references_mval.svg"
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

# PNG copy for programs that do not take SVG: 25 whole pixels per module, so the edges stay sharp
out_png <- sub("[.]svg$", ".png", out)
px <- 25 * nrow(padded)
png(out_png, width = px, height = px, res = px / (5 / 2.54), type = "cairo", antialias = "none") # opens 5 cm wide
par(mar = rep(0, 4), xaxs = "i", yaxs = "i")
image(t(unclass(padded)[nrow(padded):1, ]) * 1, col = c("white", "black"), axes = FALSE, useRaster = TRUE)
invisible(dev.off())
cat("CHECK symbol:", n, "modules, margin:", margin, "modules\n")
cat("CHECK link:  ", url, "\n")
cat("Written:     ", out, "\n")
