# Run from the repository root. Cache only these verified inputs, not results.
manifest <- read.delim("scripts/getting-started-inputs.tsv", stringsAsFactors = FALSE)
input_dir <- Sys.getenv("PHISPACE_GETTING_STARTED_INPUTS", "Test/getting-started-inputs")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
paths <- file.path(input_dir, manifest$path)
if (!all(file.exists(paths))) {
  archive <- Sys.getenv("PHISPACE_GETTING_STARTED_ARCHIVE")
  if (!nzchar(archive)) {
    archive <- tempfile(fileext = ".zip")
    options(timeout = max(1200, getOption("timeout")))
    download.file(
      "https://www.dropbox.com/scl/fo/vgwzp8jwmc7gu3c92je2i/ANWpygVEh0R2iygu_bT7nDI?rlkey=ouiefbmtoy8i9zettv7whxuwe&dl=1",
      archive, mode = "wb", method = "libcurl"
    )
  }
  entries <- unzip(archive, list = TRUE)$Name
  if (!all(manifest$path %in% entries)) stop("Getting Started archive is missing required inputs.")
  unzip(archive, files = manifest$path, exdir = input_dir)
}
hashes <- vapply(paths, digest::digest, character(1), algo = "sha256", file = TRUE)
if (any(hashes != manifest$sha256)) {
  stop("Input checksum mismatch: ", paste(manifest$path[hashes != manifest$sha256], collapse = ", "),
       ". Review the data changes before updating scripts/getting-started-inputs.tsv.")
}
message("Verified ", length(paths), " Getting Started inputs in ", normalizePath(input_dir))
