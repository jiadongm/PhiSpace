# Run from the repository root. Rebuild the function reference from pkg/man into
# a staged site, or create a standalone copy.
site_dir <- Sys.getenv("PHISPACE_SITE_DIR", "Test/site-reference")
if (!dir.exists(site_dir)) {
  dir.create(site_dir, recursive = TRUE)
  files <- list.files("docs", all.files = TRUE, no.. = TRUE, full.names = TRUE)
  if (!all(file.copy(files, site_dir, recursive = TRUE))) stop("Could not stage the existing website.")
} else if (!file.exists(file.path(site_dir, "pkgdown.yml"))) {
  stop("Existing destination is not a staged pkgdown website: ", site_dir)
}
site_dir <- normalizePath(site_dir, winslash = "/")
if (site_dir == normalizePath("docs", winslash = "/")) stop("Use a separate staging directory, not docs/.")

# Remove the staged pages so that topics deleted from pkg/man do not survive.
unlink(file.path(site_dir, "reference"), recursive = TRUE)
pkgdown::build_reference(
  pkg = "pkg", lazy = FALSE, examples = FALSE,
  override = list(destination = site_dir)
)

# Every Rd topic must have a page, and the index must link to it.
rd <- tools::file_path_sans_ext(list.files("pkg/man", pattern = "\\.Rd$"))
pages <- file.path(site_dir, "reference", paste0(rd, ".html"))
if (!all(file.exists(pages))) {
  stop("Missing reference pages: ", paste(rd[!file.exists(pages)], collapse = ", "))
}
index <- readLines(file.path(site_dir, "reference", "index.html"), warn = FALSE)
if (length(grep("href=\"[^\"]+\\.html\"", index)) == 0) stop("Reference index has no topic links.")

# The staged sitemap comes from docs/; list the rebuilt reference pages instead.
sitemap <- file.path(site_dir, "sitemap.xml")
sm <- readLines(sitemap, warn = FALSE)
is_ref <- grepl("/reference/", sm, fixed = TRUE)
if (!any(is_ref)) stop("Staged sitemap has no reference entries.")
base <- sub("/reference/.*", "/reference/", sm[which(is_ref)[1]])
html <- sort(list.files(file.path(site_dir, "reference"), pattern = "\\.html$"), method = "radix")
sm <- append(sm[!is_ref], paste0(base, html, "</loc></url>"), after = which(is_ref)[1] - 1)
writeLines(sm, sitemap)
invisible(file.create(file.path(site_dir, ".nojekyll")))
message("Built ", length(rd), " reference pages: ", file.path(site_dir, "reference"))
