# Run from the repository root. Cache only these verified inputs, not results.
manifest <- read.delim("scripts/cosmx-inputs.tsv", stringsAsFactors = FALSE)
input_dir <- Sys.getenv("PHISPACE_COSMX_INPUTS", "Test/cosmx-inputs")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
paths <- file.path(input_dir, manifest$path)
if (!all(file.exists(paths))) {
  archive <- Sys.getenv("PHISPACE_COSMX_ARCHIVE")
  if (!nzchar(archive)) {
    archive <- tempfile(fileext = ".zip")
    options(timeout = max(1200, getOption("timeout")))
    download.file(
      "https://www.dropbox.com/scl/fo/z01qvst71kzxqog4bf5po/AKZzi-PgkZCDYZdUwLqfBWg?rlkey=65kg8zqatk2zdk5mnirgn1vmc&dl=1",
      archive, mode = "wb", method = "libcurl"
    )
  }
  entries <- unzip(archive, list = TRUE)$Name
  if (!all(manifest$path %in% entries)) stop("CosMx archive is missing required files.")
  unzip(archive, files = manifest$path, exdir = input_dir)
}
hashes <- vapply(paths, digest::digest, character(1), algo = "sha256", file = TRUE)
if (any(hashes != manifest$sha256)) {
  stop("Input checksum mismatch: ", paste(manifest$path[hashes != manifest$sha256], collapse = ", "),
       ". Review the data changes before updating scripts/cosmx-inputs.tsv.")
}
message("Verified ", length(paths), " CosMx files in ", normalizePath(input_dir))
