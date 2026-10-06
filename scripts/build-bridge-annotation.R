# Run from the repository root. Append to a staged site or create a standalone copy.
source("scripts/prepare-bridge-annotation.R", local = TRUE)
site_dir <- Sys.getenv("PHISPACE_SITE_DIR", "Test/site-bridge-annotation")
if (!dir.exists(site_dir)) {
  dir.create(site_dir, recursive = TRUE)
  files <- list.files("docs", all.files = TRUE, no.. = TRUE, full.names = TRUE)
  if (!all(file.copy(files, site_dir, recursive = TRUE))) stop("Could not stage the existing website.")
} else if (!file.exists(file.path(site_dir, "pkgdown.yml"))) {
  stop("Existing destination is not a staged pkgdown website: ", site_dir)
}
site_dir <- normalizePath(site_dir, winslash = "/")
if (site_dir == normalizePath("docs", winslash = "/")) stop("Use a separate staging directory, not docs/.")

# Fresh outputs prevent the result-loading branch from hiding computation errors.
run_dir <- tempfile("bridge-annotation-run-")
dir.create(run_dir)
for (rel in manifest$path) {
  target <- file.path(run_dir, rel)
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(file.path(input_dir, rel), target)) stop("Could not stage input: ", rel)
}
Sys.setenv(PHISPACE_BRIDGE_ANNOTATION_DATA = normalizePath(run_dir, winslash = "/"))
pkgload::load_all("pkg", quiet = TRUE)
pkgdown::build_article(
  "BridgeAnnotation", pkg = "pkg", lazy = FALSE, new_process = FALSE,
  override = list(destination = site_dir), quiet = FALSE
)
article <- file.path(site_dir, "articles", "BridgeAnnotation.html")
html <- readLines(article, warn = FALSE)
if (any(grepl("Dropbox/Research_projects", html, fixed = TRUE)) ||
    !any(grepl("cellTypeTable.rds", html, fixed = TRUE))) {
  stop("Rendered BridgeAnnotation article did not contain the expected path fixes.")
}
source("scripts/check-results.R", local = TRUE)
check_results("BridgeAnnotation", file.path(run_dir, "output"), slug = "bridge-annotation")
writeLines(capture.output(sessionInfo()), file.path(site_dir, "BridgeAnnotation-sessionInfo.txt"))
file.create(file.path(site_dir, ".nojekyll"))
message("Built BridgeAnnotation with fresh results: ", article)
