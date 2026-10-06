# Run from the repository root. Cache only these verified inputs, not results.
manifest <- read.delim("scripts/cite-seq-inputs.tsv", stringsAsFactors = FALSE)
input_dir <- Sys.getenv("PHISPACE_CITE_SEQ_INPUTS", "Test/cite-seq-inputs")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
paths <- file.path(input_dir, manifest$path)
if (!all(file.exists(paths))) {
  archive <- Sys.getenv("PHISPACE_CITE_SEQ_ARCHIVE")
  if (!nzchar(archive)) {
    archive <- tempfile(fileext = ".zip")
    options(timeout = max(1200, getOption("timeout")))
    download.file(
      "https://www.dropbox.com/scl/fo/it8uwxd2v4k2lyoo936at/AKws5m_xeNjOOs8Vkho-wP8?rlkey=s4d5a7cfkl5sj9qiogc7mng5t&dl=1",
      archive, mode = "wb", method = "libcurl"
    )
  }
  entries <- unzip(archive, list = TRUE)$Name
  if (!all(manifest$path %in% entries)) stop("CITE-seq archive is missing required files.")
  unzip(archive, files = manifest$path, exdir = input_dir)
}
hashes <- vapply(paths, digest::digest, character(1), algo = "sha256", file = TRUE)
if (any(hashes != manifest$sha256)) {
  stop("Input checksum mismatch: ", paste(manifest$path[hashes != manifest$sha256], collapse = ", "),
       ". Review the data changes before updating scripts/cite-seq-inputs.tsv.")
}
message("Verified ", length(paths), " CITE-seq files in ", normalizePath(input_dir))
