# Run from the repository root. Write redirect pages for the website, which moved
# to the package repository github.com/mao-lab-org/PhiSpace on 2026-10-12. Each
# page listed in the old sitemap gets a page that redirects to the same path on
# the new website; 404.html redirects any other path, for example figures.
old_url <- "https://jiadongm.github.io/PhiSpace/"
new_url <- "https://mao-lab-org.github.io/PhiSpace/"
out_dir <- "redirects"
sitemap <- Sys.getenv("PHISPACE_OLD_SITEMAP", paste0(old_url, "sitemap.xml"))

sm <- readLines(sitemap, warn = FALSE)
loc <- regmatches(sm, regexpr("<loc>[^<]+</loc>", sm))
paths <- sub(old_url, "", gsub("</?loc>", "", loc), fixed = TRUE)
if (length(paths) == 0 || any(grepl("^https?:", paths))) stop("Unexpected sitemap entries in ", sitemap)
# The LICENSE pages are listed but do not exist; home.html is an old copy of index.html.
paths <- setdiff(paths, c("LICENSE.html", "LICENSE-text.html", "404.html"))
target <- ifelse(paths == "home.html", "index.html", paths)

page <- function(url, script = "") {
  c("<!DOCTYPE html>", "<html lang=\"en\">", "<head>", "<meta charset=\"utf-8\">",
    "<title>PhiSpace has moved</title>",
    paste0("<link rel=\"canonical\" href=\"", url, "\">"),
    paste0("<meta http-equiv=\"refresh\" content=\"0; url=", url, "\">"),
    "<meta name=\"robots\" content=\"noindex\">", script, "</head>", "<body>",
    paste0("<p>The PhiSpace website has moved to <a href=\"", url, "\">", url, "</a>.</p>"),
    "</body>", "</html>")
}

unlink(out_dir, recursive = TRUE)
for (i in seq_along(paths)) {
  file <- file.path(out_dir, paths[i])
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  writeLines(page(paste0(new_url, target[i])), file)
}
# The script runs before the refresh and keeps the path, query and fragment.
js <- paste0("<script>location.replace(\"", new_url,
             "\" + location.pathname.replace(/^\\/PhiSpace\\/?/, \"\") + location.search + location.hash);</script>")
writeLines(page(new_url, js), file.path(out_dir, "404.html"))
invisible(file.create(file.path(out_dir, ".nojekyll")))
message("Wrote ", length(paths) + 1, " redirect pages to ", out_dir)
