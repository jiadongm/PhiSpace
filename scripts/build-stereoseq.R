# Run from the repository root after installing the workflow's R dependencies.
source("scripts/prepare-stereoseq.R", local = TRUE)
site_dir <- Sys.getenv("PHISPACE_SITE_DIR", "Test/site-stereoseq")
if (file.exists(site_dir)) stop("Use a new PHISPACE_SITE_DIR; destination already exists: ", site_dir)
dir.create(site_dir, recursive = TRUE)
site_dir <- normalizePath(site_dir, winslash = "/")
site_files <- list.files("docs", all.files = TRUE, no.. = TRUE, full.names = TRUE)
if (!all(file.copy(site_files, site_dir, recursive = TRUE))) stop("Could not stage the existing website.")

# Each build gets a new output directory; supplied analysis caches are excluded.
run_dir <- tempfile("stereoseq-run-")
dir.create(run_dir)
for (rel in manifest$path) {
  target <- file.path(run_dir, rel)
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(file.path(input_dir, rel), target)) stop("Could not stage input: ", rel)
}
Sys.setenv(PHISPACE_STEREOSEQ_DATA = normalizePath(run_dir, winslash = "/"))
pkgload::load_all("pkg", quiet = TRUE)
pkgdown::build_article(
  "StereoSeq", pkg = "pkg", lazy = FALSE, new_process = FALSE,
  override = list(destination = site_dir), quiet = FALSE
)
article <- file.path(site_dir, "articles", "StereoSeq.html")
html <- readLines(article, warn = FALSE)
if (any(grepl("PhiSpaceAssay", html, fixed = TRUE)) ||
    !any(grepl("refAssay", html, fixed = TRUE))) stop("Rendered article did not contain the expected API fix.")
source("scripts/check-results.R", local = TRUE)
check_results("StereoSeq", file.path(run_dir, "output"))
writeLines(capture.output(sessionInfo()), file.path(site_dir, "StereoSeq-sessionInfo.txt"))
file.create(file.path(site_dir, ".nojekyll"))
message("Built Stereo-seq with fresh results: ", article)
