# Run this on the same machine and R version used to produce the paper figures,
# after loading every package the analyses depend on. Output overwrites
# environment/sessionInfo.txt.
#
# Recommended invocation from the project root:
#   Rscript environment/capture_sessionInfo.R

deps <- c(
  # PhiSpace + core SC
  "PhiSpace", "SingleCellExperiment", "SummarizedExperiment", "scran", "scuttle",
  "SpaNorm", "DropletUtils", "Seurat", "SeuratDisk", "Matrix", "sctransform",
  # Visualisation / utilities
  "ggplot2", "ggpubr", "ComplexHeatmap", "circlize", "seriation",
  "dplyr", "tidyr", "magrittr", "qs", "yaml", "devtools",
  # Spatial / deconvolution / DWD / signatures
  "Giotto", "spacexr", "SPOTlight", "kerndwd", "UCell", "msigdbr"
)

loaded <- character()
for (pkg in deps) {
  ok <- suppressPackageStartupMessages(
    tryCatch({ requireNamespace(pkg, quietly = TRUE); library(pkg, character.only = TRUE) },
             error = function(e) FALSE)
  )
  if (isTRUE(ok)) loaded <- c(loaded, pkg) else message("[skip] ", pkg, " not installed")
}

sink(file.path("environment", "sessionInfo.txt"))
cat("# R sessionInfo() captured by environment/capture_sessionInfo.R\n")
cat("# Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("# Loaded packages:", paste(loaded, collapse = ", "), "\n\n")
print(sessionInfo())
sink()
cat("Wrote environment/sessionInfo.txt\n")
