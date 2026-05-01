# 02_cosmx_case — CosMx NSCLC (Section 2.3 / Fig 3)

PhiSpace ST applied to 8 Nanostring CosMx WTX NSCLC samples (5 patients),
using 4 lineage-specific scRNA-seq references from Natri et al. 2024
(immune, epithelial, endothelial, mesenchymal — total 43 cell types).

## Inputs (set via `config/paths_local.yml`)

- `${data_root}/CosMx_lung/SMI_Giotto_Object.RData` — original CosMx
  Giotto object (8 samples; see `dataset_manifest.tsv` row `cosmx_nsclc`).
- `${data_root}/LungFibrosis/scRNA-seq/{immune,epithelial,endothelial,mesenchymal}_sce.qs`
  — Natri et al. 2024 lineage references (built by `00_prepareSCE.R`).

## Run order

```
00_prepareSCE.R              Process Natri Seurat objects -> per-lineage SCE; convert
                             CosMx Giotto object -> per-sample SCE (.qs).
01_runPhiSpace4Refs.R        Run PhiSpace separately per lineage; concatenate scores.
02_cosmx4Refs_main.Rmd       Master notebook: assembles Fig 3 + Supp Fig S6 / Table S1 source panels.
03_regenFig3C.R              Regenerate Fig 3C cluster spatial maps with the updated
                             niche-5 colour (gold). Run only if you need to refresh
                             the figure colours without rerunning the full notebook.
```

`_utils.R` defines `tissueNames_CosMx` (the 8 sample IDs in canonical order),
the cell type name simplifications, and plotting helpers.

## Outputs

- `${output_root}/Case3/CosMxAllLungsPhiRes4Refs_log1p.qs` — full per-sample
  per-lineage PhiSpace scores
- `${output_root}/Case3/Lung5_Rep1_PhiClusts4Refs.qs` and
  `Lung5_Rep1_OmicsClusts4Refs.qs` — clustering caches used by `03_regenFig3C.R`
- `figs/clusts/`, `figs/niches/`, `figs/spatialHeatmap/`, `figs/multiSample/`,
  `figs/signatures/`, `figs/specialCases/` (regenerated; not shipped)

Manuscript items supported: Fig 3 panels A-F, Supp Fig S6, and Supp Table S1
(see `outputs_manifest/expected_outputs.tsv` rows whose `case_study=cosmx_nsclc`).
