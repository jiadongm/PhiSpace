# Run from the repository root. Cache only these verified inputs, not results.
manifest <- read.delim("scripts/visium-inputs.tsv", stringsAsFactors = FALSE)
input_dir <- Sys.getenv("PHISPACE_VISIUM_INPUTS", "Test/visium-inputs")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
paths <- file.path(input_dir, manifest$path)
if (!all(file.exists(paths))) {
  archive <- Sys.getenv("PHISPACE_VISIUM_ARCHIVE")
  if (!nzchar(archive)) {
    archive <- tempfile(fileext = ".zip")
    options(timeout = max(1200, getOption("timeout")))
    download.file(
      "https://www.dropbox.com/scl/fo/gsaxu5jex7d8ftrwu3pf0/ANN2z_e4abgrON62BD2TEHk?rlkey=fg94wev8y6cn096khmauyaass&dl=1",
      archive, mode = "wb", method = "libcurl"
    )
  }
  entries <- unzip(archive, list = TRUE)$Name
  if (!all(manifest$path %in% entries)) stop("Visium archive is missing required inputs.")
  unzip(archive, files = manifest$path, exdir = input_dir)
}
hashes <- vapply(paths, digest::digest, character(1), algo = "sha256", file = TRUE)
if (any(hashes != manifest$sha256)) {
  stop("Input checksum mismatch: ", paste(manifest$path[hashes != manifest$sha256], collapse = ", "),
       ". Review the data changes before updating scripts/visium-inputs.tsv.")
}
message("Verified ", length(paths), " Visium inputs in ", normalizePath(input_dir))
