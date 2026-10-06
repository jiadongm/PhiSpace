# Run from the repository root. Cache only these verified inputs, not results.
manifest <- read.delim("scripts/bridge-annotation-inputs.tsv", stringsAsFactors = FALSE)
input_dir <- Sys.getenv("PHISPACE_BRIDGE_ANNOTATION_INPUTS", "Test/bridge-annotation-inputs")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
paths <- file.path(input_dir, manifest$path)
if (!all(file.exists(paths))) {
  archive <- Sys.getenv("PHISPACE_BRIDGE_ANNOTATION_ARCHIVE")
  if (!nzchar(archive)) {
    archive <- tempfile(fileext = ".zip")
    options(timeout = max(1200, getOption("timeout")))
    download.file(
      "https://www.dropbox.com/scl/fo/jeuqzjfyyr2j7doa922ve/AI-2U_wtZBpPGOswxMzReGQ?rlkey=n328yyr2llf81gz3chynjg0r6&dl=1",
      archive, mode = "wb", method = "libcurl"
    )
  }
  entries <- unzip(archive, list = TRUE)$Name
  if (!all(manifest$path %in% entries)) stop("BridgeAnnotation archive is missing required inputs.")
  unzip(archive, files = manifest$path, exdir = input_dir)
}
hashes <- vapply(paths, digest::digest, character(1), algo = "sha256", file = TRUE)
if (any(hashes != manifest$sha256)) {
  stop("Input checksum mismatch: ", paste(manifest$path[hashes != manifest$sha256], collapse = ", "),
       ". Review the data changes before updating scripts/bridge-annotation-inputs.tsv.")
}
message("Verified ", length(paths), " BridgeAnnotation inputs in ", normalizePath(input_dir))
