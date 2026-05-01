# 01_visium_case — Visium NSCLC (Section 2.2 / Fig 2)

PhiSpace ST applied to 18 10x Visium NSCLC samples (4 healthy donors + 14
NSCLC samples covering LUAD and LUSC), using the Azimuth Human Lung 2.0
(HLCA subsample) reference (61 cell types).

## Inputs (set via `config/paths_local.yml`)

- `${data_root}/Visium_NSCLC/` — per-sample Visium `filtered_feature_bc_matrix.h5`
  and `spatial/` outputs (see `dataset_manifest.tsv` row `visium_nsclc_18_samples`).
- `${data_root}/LungRef/AzimuthLung2.0_sce_0.1sub.qs` — HLCA Azimuth Lung 2.0
  reference (10% subsample, 58,759 cells, 61 cell types).
- `${data_root}/LungRef/AzimuthLung2.0_sce_0.1sub_selectFeat.h5ad` — same
  reference exported as h5ad with selected features (used by competitor
  notebooks).

## Run order

```
00_prepareData.R          Build per-sample SCE list with log1p assay (caches to ${output_root}).
01_azimuthLungMarkers.R   Extract Azimuth marker list for plotting.
02_runPhiSpace.R          Run PhiSpace (Azimuth Lung 2.0, 500 features) across all 18 samples.
03_runtimeBenchmark.R     Time PhiSpace vs RCTD on the same samples (Fig 2 timing panel).
04_DEGAnalysis.R          Differential expression analysis (supplementary).
05_visiumNSCLC_main.Rmd   Master notebook: assembles Fig 2 + Supp Figs S2–S3.
competitors/runRCTD.R                       RCTD competitor.
competitors/runTACCO.ipynb                  TACCO competitor (Python; needs tacco_env).
competitors/runCell2location_train.ipynb    Cell2location reference training (Python; needs cell2location env).
competitors/runCell2location_apply.ipynb    Cell2location query application (Python).
```

## Outputs

Written to `${output_root}/Case3/...` and to `figs/` relative to the working
directory. Major outputs (regenerated; not shipped):

- `*PhiSpace*.qs` — per-sample PhiSpace score matrices
- `figs/spatialHeatmap/` — per-sample cell type spatial heatmaps
- `figs/copresence/` — co-presence PCA, loadings, heatmaps
- `figs/specialCases/` — Fig 2 zoomed panels

Manuscript items supported: see `outputs_manifest/expected_outputs.tsv`
(rows whose `case_study=visium_nsclc`).

The pDC-fibroblast supplementary Visium panels are generated from
`scripts/05_pdc_fibroblast_supp/` and appear as Supp Figs S4-S5.

## Notes

- The competitor notebooks in `competitors/` were originally named
  `CaseCosMx_*.ipynb` despite being Visium analyses; they have been renamed
  during packaging.
- `_utils.R` provides plotting and analysis helpers shared with other case
  studies (the same file content is also referenced from CosMx legacy code).
