# Run from the repository root. Cache only these verified inputs, not results.
manifest <- read.delim("scripts/stereoseq-inputs.tsv", stringsAsFactors = FALSE)
input_dir <- Sys.getenv("PHISPACE_STEREOSEQ_INPUTS", "Test/stereoseq-inputs")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
paths <- file.path(input_dir, manifest$path)
if (!all(file.exists(paths))) {
  archive <- Sys.getenv("PHISPACE_STEREOSEQ_ARCHIVE")
  if (!nzchar(archive)) {
    archive <- tempfile(fileext = ".zip")
    options(timeout = max(1200, getOption("timeout")))
    download.file(
      "https://www.dropbox.com/scl/fo/4w5vweo2ky2vuf7g591fe/AMDVr2OAUL5W1wO7ATzM754?rlkey=aggfds07sjymsd2aomyepbylv&dl=1",
      archive, mode = "wb", method = "libcurl"
    )
  }
  entries <- unzip(archive, list = TRUE)$Name
  if (!all(manifest$path %in% entries)) stop("Stereo-seq archive is missing required inputs.")
  unzip(archive, files = manifest$path, exdir = input_dir)
}
hashes <- vapply(paths, digest::digest, character(1), algo = "sha256", file = TRUE)
if (any(hashes != manifest$sha256)) {
  stop("Input checksum mismatch: ", paste(manifest$path[hashes != manifest$sha256], collapse = ", "),
       ". Review the data changes before updating scripts/stereoseq-inputs.tsv.")
}
message("Verified ", length(paths), " Stereo-seq inputs in ", normalizePath(input_dir))
