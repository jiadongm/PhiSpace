# Run from the repository root. Append to a staged site or create a standalone copy.
# The Stored models vignette uses the verified CosMx inputs (scripts/cosmx-inputs.tsv)
# and always computes its results; it has no cached mode.
source("scripts/prepare-cosmx.R", local = TRUE)
site_dir <- Sys.getenv("PHISPACE_SITE_DIR", "Test/site-stored-models")
if (!dir.exists(site_dir)) {
  dir.create(site_dir, recursive = TRUE)
  files <- list.files("docs", all.files = TRUE, no.. = TRUE, full.names = TRUE)
  if (!all(file.copy(files, site_dir, recursive = TRUE))) stop("Could not stage the existing website.")
} else if (!file.exists(file.path(site_dir, "pkgdown.yml"))) {
  stop("Existing destination is not a staged pkgdown website: ", site_dir)
}
site_dir <- normalizePath(site_dir, winslash = "/")
if (site_dir == normalizePath("docs", winslash = "/")) stop("Use a separate staging directory, not docs/.")

# Stage only the inputs that the vignette reads, in a fresh folder.
run_dir <- tempfile("stored-models-run-")
dir.create(run_dir)
staged <- grep("^data/(LungFibrosis/|CosMx_lung/Lung5_Rep1_SCE)", manifest$path, value = TRUE)
if (length(staged) != 5) stop("Expected 4 references and 1 query in scripts/cosmx-inputs.tsv.")
for (rel in staged) {
  target <- file.path(run_dir, rel)
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(file.path(input_dir, rel), target)) stop("Could not stage input: ", rel)
}
Sys.setenv(PHISPACE_COSMX_DATA = normalizePath(run_dir, winslash = "/"))
pkgload::load_all("pkg", quiet = TRUE)
pkgdown::build_article(
  "StoredModels", pkg = "pkg", lazy = FALSE, new_process = FALSE,
  override = list(destination = site_dir), quiet = FALSE
)
article <- file.path(site_dir, "articles", "StoredModels.html")
html <- readLines(article, warn = FALSE)
if (any(grepl("Dropbox/Research_projects", html, fixed = TRUE)) ||
    !any(grepl("may be preprocessed differently", html, fixed = TRUE))) {
  stop("Rendered Stored models article lacks the expected preprocessing warning.")
}
source("scripts/check-results.R", local = TRUE)
check_results("StoredModels", file.path(run_dir, "output"), slug = "stored-models")
writeLines(capture.output(sessionInfo()), file.path(site_dir, "StoredModels-sessionInfo.txt"))
file.create(file.path(site_dir, ".nojekyll"))
message("Built Stored models: ", article)
