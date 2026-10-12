# Run from the repository root. Rebuild the home page (index.html), the authors
# and citation page (authors.html) and 404.html from pkg/README.md,
# pkg/DESCRIPTION and pkg/inst/CITATION, and the articles index
# (articles/index.html) from pkg/_pkgdown.yml, into a staged site. Run it after
# the articles.
site_dir <- Sys.getenv("PHISPACE_SITE_DIR", "Test/site-home")
if (!dir.exists(site_dir)) {
  dir.create(site_dir, recursive = TRUE)
  files <- list.files("docs", all.files = TRUE, no.. = TRUE, full.names = TRUE)
  if (!all(file.copy(files, site_dir, recursive = TRUE))) stop("Could not stage the existing website.")
} else if (!file.exists(file.path(site_dir, "pkgdown.yml"))) {
  stop("Existing destination is not a staged pkgdown website: ", site_dir)
}
site_dir <- normalizePath(site_dir, winslash = "/")
if (site_dir == normalizePath("docs", winslash = "/")) stop("Use a separate staging directory, not docs/.")

pkgdown::build_home(pkg = "pkg", override = list(destination = site_dir), preview = FALSE)

# The citation page must list every DOI in pkg/inst/CITATION, and the home page
# must link to each of them.
cit <- utils::readCitationFile("pkg/inst/CITATION", meta = list(Encoding = "UTF-8"))
dois <- unlist(lapply(cit, function(x) x$doi))
if (length(dois) == 0) stop("pkg/inst/CITATION has no DOIs.")
for (page in c("authors.html", "index.html")) {
  html <- paste(readLines(file.path(site_dir, page), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  missing <- dois[!vapply(dois, grepl, logical(1), x = html, fixed = TRUE)]
  if (length(missing)) stop(page, " does not cite: ", paste(missing, collapse = ", "))
}

# The staged articles index and sitemap come from docs/; rebuild the index and
# list every article of pkg/_pkgdown.yml, each of which must have been built.
pkgdown::build_articles_index(pkg = "pkg", override = list(destination = site_dir))
cfg <- yaml::read_yaml("pkg/_pkgdown.yml")
articles <- unlist(lapply(cfg$articles, `[[`, "contents"))
pages <- paste0(articles, ".html")
built <- file.exists(file.path(site_dir, "articles", pages))
if (!all(built)) stop("Articles not built: ", paste(articles[!built], collapse = ", "))
index <- paste(readLines(file.path(site_dir, "articles", "index.html"), warn = FALSE), collapse = "\n")
unlinked <- pages[!vapply(paste0("href=\"", pages, "\""), grepl, logical(1), x = index, fixed = TRUE)]
if (length(unlinked)) stop("Articles index does not link to: ", paste(unlinked, collapse = ", "))
sitemap <- file.path(site_dir, "sitemap.xml")
sm <- readLines(sitemap, warn = FALSE)
is_art <- grepl("/articles/", sm, fixed = TRUE)
if (!any(is_art)) stop("Staged sitemap has no article entries.")
base <- sub("/articles/.*", "/articles/", sm[which(is_art)[1]])
html <- sort(c(pages, "index.html"), method = "radix")
sm <- append(sm[!is_art], paste0(base, html, "</loc></url>"), after = which(is_art)[1] - 1)
writeLines(sm, sitemap)
message("Built home and citation pages citing ", length(dois), " DOIs, and the index of ",
        length(articles), " articles: ", site_dir)
