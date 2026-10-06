# Run from the repository root. Append to a staged site or create a standalone copy.
# PHISPACE_BRIDGE_ANNOTATION_MODE = "cached" (default) renders from the verified
# results in the shared folder; a fresh build peaks at 17-24 GB RAM, above the
# GitHub runner limit. Use "fresh" locally to recompute them.
mode <- Sys.getenv("PHISPACE_BRIDGE_ANNOTATION_MODE", "cached")
if (!mode %in% c("cached", "fresh")) stop("PHISPACE_BRIDGE_ANNOTATION_MODE must be 'cached' or 'fresh'.")
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

# A fresh build stages inputs only, so the result-loading branches cannot hide
# computation errors. A cached build also stages the verified results.
run_dir <- tempfile("bridge-annotation-run-")
dir.create(run_dir)
staged <- manifest$path[manifest$role == "input" | mode == "cached"]
for (rel in staged) {
  target <- file.path(run_dir, rel)
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(file.path(input_dir, rel), target)) stop("Could not stage input: ", rel)
}
Sys.setenv(PHISPACE_BRIDGE_ANNOTATION_DATA = normalizePath(run_dir, winslash = "/"))
cached_md5 <- tools::md5sum(file.path(run_dir, manifest$path[manifest$role == "cached"]))
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
if (mode == "cached" && !identical(cached_md5, tools::md5sum(names(cached_md5)))) {
  stop("The cached build rewrote a verified result; a result-loading branch did not run.")
}
check_results("BridgeAnnotation", file.path(run_dir, "output"), slug = "bridge-annotation")
writeLines(capture.output(sessionInfo()), file.path(site_dir, "BridgeAnnotation-sessionInfo.txt"))
file.create(file.path(site_dir, ".nojekyll"))
message("Built BridgeAnnotation with ", mode, " results: ", article)
