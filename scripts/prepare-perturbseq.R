# Run from the repository root. Cache only these verified inputs, not results.
manifest <- read.delim("scripts/perturbseq-inputs.tsv", stringsAsFactors = FALSE)
input_dir <- Sys.getenv("PHISPACE_PERTURBSEQ_INPUTS", "Test/perturbseq-inputs")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
paths <- file.path(input_dir, manifest$path)
if (!all(file.exists(paths))) {
  archive <- Sys.getenv("PHISPACE_PERTURBSEQ_ARCHIVE")
  if (!nzchar(archive)) {
    archive <- tempfile(fileext = ".zip")
    options(timeout = max(1200, getOption("timeout")))
    download.file(
      "https://www.dropbox.com/scl/fo/4gm1ef27yl4wgb0f6zio7/AG0bvE9LXKgv61AAKmD5ho8?rlkey=bkd7nlwo6qlltkl21645xdw90&dl=1",
      archive, mode = "wb", method = "libcurl"
    )
  }
  entries <- unzip(archive, list = TRUE)$Name
  if (!all(manifest$path %in% entries)) stop("PerturbSeq archive is missing required files.")
  unzip(archive, files = manifest$path, exdir = input_dir)
}
hashes <- vapply(paths, digest::digest, character(1), algo = "sha256", file = TRUE)
if (any(hashes != manifest$sha256)) {
  stop("Input checksum mismatch: ", paste(manifest$path[hashes != manifest$sha256], collapse = ", "),
       ". Review the data changes before updating scripts/perturbseq-inputs.tsv.")
}
message("Verified ", length(paths), " PerturbSeq files in ", normalizePath(input_dir))
