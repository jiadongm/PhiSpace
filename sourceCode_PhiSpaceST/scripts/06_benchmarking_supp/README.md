# 06_benchmarking_supp — 7-method deconvolution benchmark (Supp Fig S1)

Performance + runtime benchmark of PhiSpace against 6 spatial transcriptomics
deconvolution methods (RCTD, Cell2location, TACCO, Seurat, SPOTlight, DSTG),
on the 32 simulated datasets from Li et al. 2022 (Nature Methods, BLAST
benchmark). All methods except DSTG were re-run by us on the same prepared
CSV inputs for a fair comparison.

## Methods compared

| Method | Source | Runner | Median runtime (32 datasets) | Median composite AS |
|---|---|---|---|---|
| RCTD | spacexr | `02_run_rctd.R` | 311 s | 0.812 |
| Cell2location (GPU) | pip | `02_run_cell2location.py` | 992 s (CPU much slower) | 0.802 |
| TACCO | conda-forge | `02_run_tacco.py` | 11 s | 0.758 |
| **PhiSpace** | github.com/jiadongm/PhiSpace | `02_run_phispace.R` | 11 s | 0.676 |
| Seurat (label transfer) | CRAN 5.4.0 | `02_run_seurat.R` | 39 s | 0.664 |
| SPOTlight | Bioconductor 1.14.0 | `02_run_spotlight.R` | 32 s | 0.644 |
| DSTG | (paper values only) | — | n/a | 0.601 |

## Inputs (set via `config/paths_local.yml`)

- The 32 simulated `*.h5ad` datasets from the BLAST repository:
  `${data_root}/SpatialBenchmarking/SimualtedSpatalData/dataset{1-16}.h5ad`
  + `dataset{1-16}_r.h5ad`. See `dataset_manifest.tsv` row
  `blast_benchmark_simulated_h5ad`. Download from
  <https://github.com/QuKunLab/SpatialBenchmarking>.
- For DSTG: `${data_root}/SpatialBenchmarking/FigureData/Figure4/32SimulationData/{PCC,SSIM,RMSE,JSD}.csv`
  (from the BLAST repo; only DSTG values are kept from upstream).
- `cell2loc_python` and `tacco_python` set in `paths_local.yml`.

## Pipeline (numbered scripts run in order within each method)

```
01_prepare_data.py                       h5ad -> CSV; produces ${data_root}/SpatialBenchmarking/
                                         prepared_data/dataset{1-32}/ (NOT shipped — regenerable
                                         from the h5ad sources, ~71 GB).
02_run_<method>.{R,py}                   Run a single method on a single dataset.
03_calculate_metrics.py                  Compute PCC, SSIM, RMSE, JSD per cell type.
04_run_<method>_benchmark.py             Orchestrate one method across all 32 datasets.
04_run_full_benchmark.py                 Orchestrate PhiSpace across all 32 datasets.
05_integrate_results.py                  Merge PhiSpace + competitor metrics.
06_visualize_results.py                  (Earlier draft figures.)
07_runtime_benchmark.py                  Aggregate per-method runtime CSVs.
08_datasize.R                            Cache n_cells / n_spots / n_genes / n_celltypes per dataset.
09_visualise_results.R                   FINAL figures: Supp Fig S1A (performance) + S1B (runtime).
10_analyse_performance.R                 Per-dataset comparison + delta-from-best supplementary.
run_benchmark.sh        Shell driver for the full prepare -> PhiSpace -> integrate -> visualise loop.
run_all_datasets.sh     Per-dataset loop (for one method at a time).
```

## Outputs

Per method, each dataset:

- `${output_root}/SpatialBenchmarking/results_<method>/dataset{1-32}/{Method}_result.txt`
- `${output_root}/SpatialBenchmarking/results_<method>/all_metrics/dataset{1-32}/summary.csv`
- `${output_root}/SpatialBenchmarking/results_<method>/<dataset>/{Method}_runtime.txt`

Aggregated:

- `figures/boxplot_all_methods.{pdf,png}` — Supp Fig S1A (7-method performance)
- `figures/runtime_boxplot.{pdf,png}` — Supp Fig S1B (6-method runtime; DSTG has no runtime)
- `figures/performance_by_dataset.{pdf,png}`, `figures/performance_delta.{pdf,png}` — extra supp panels

## PhiSpace benchmark variants (in `02_run_phispace.R`)

The default PhiSpace pipeline is `normPhiScores() -> Score2Prob()` (softmax),
which is the variant reported in the paper. The script also produces:

- `PhiSpace_unnorm_result.txt` — raw scores -> Score2Prob() (no normPhiScores)
- `PhiSpace_minmax_result.txt` — normPhiScores() -> per-row min-max [0,1] -> proportions

Both alternatives are inferior; see `CLAUDE.md` of the source folder for the
delta-AS comparison. They are kept for reproducibility of the supplementary
methods discussion.

## Notes on filtering

Both RCTD (`02_run_rctd.R`) and PhiSpace (`02_run_phispace.R`) filter rare
cell types with < 25 cells before deconvolution. This matches `spacexr::Reference()`
internals. Seurat and SPOTlight were re-run with the same filter to ensure
fair comparison; their AS values therefore differ slightly from the original
BLAST paper values.

## Notes on Seurat v5 / SPOTlight

- Seurat v5 (5.4.0) requires sparse matrices (`Matrix::Matrix(x, sparse=TRUE)`)
  before `CreateSeuratObject()`; dense input triggers a `LogMap` "Rownames
  cannot be empty strings" bug. Also requires
  `options(future.globals.maxSize = 3 * 1024^3)` for SCTransform on large datasets.
- SPOTlight 1.14.0 uses `SPOTlight()` (not the older `spotlight_deconvolution()`).
  Pass spatial data as a sparse matrix, not a Seurat v5 object.
- Both runners filter out empty gene names from the spatial CSVs (occasional
  artifact from CSV index parsing).
