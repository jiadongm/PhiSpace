# PhiSpace ST path loader (R).
#
# Source this from any analysis script:
#   source("scripts/00_setup/paths_loader.R")
# It walks up from the current working directory to find
# `config/paths_template.yml`, reads `paths_local.yml` if present (else the
# template, with a warning), and exposes a `paths` named list.
#
# Convenience derivations:
#   paths$data_root   = file.path(paths$phispace_data_root, "data")
#   paths$output_root = file.path(paths$phispace_data_root, "output")

local({
  find_proj_root <- function(start = getwd()) {
    d <- normalizePath(start, mustWork = TRUE)
    repeat {
      if (file.exists(file.path(d, "config", "paths_template.yml"))) return(d)
      parent <- dirname(d)
      if (parent == d) {
        stop("paths_loader: could not find sourceCode_PhiSpaceST/config/paths_template.yml ",
             "by walking up from '", start, "'. Set the working directory to the ",
             "sourceCode_PhiSpaceST/ root (or anywhere inside it) before sourcing this file.")
      }
      d <- parent
    }
  }

  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("paths_loader: the 'yaml' package is required. Install with install.packages('yaml').")
  }

  proj <- find_proj_root()
  local_yml    <- file.path(proj, "config", "paths_local.yml")
  template_yml <- file.path(proj, "config", "paths_template.yml")
  cfg_path <- if (file.exists(local_yml)) local_yml else template_yml
  if (cfg_path == template_yml) {
    message("[paths_loader] Using paths_template.yml (placeholders). ",
            "Copy it to paths_local.yml and edit before running real analyses.")
  }

  cfg <- yaml::read_yaml(cfg_path)
  cfg$project_root <- proj
  if (!is.null(cfg$phispace_data_root)) {
    if (is.null(cfg$data_root))   cfg$data_root   <- file.path(cfg$phispace_data_root, "data")
    if (is.null(cfg$output_root)) cfg$output_root <- file.path(cfg$phispace_data_root, "output")
  }
  assign("paths", cfg, envir = parent.frame(2))
})
