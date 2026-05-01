# 04_stereoseq_case — AML mouse spleen Stereo-seq (Section 2.5 / Fig 5)

PhiSpace ST applied to a BGI Stereo-seq dataset of an AML mouse spleen with
SPLINTR cellular lineage barcodes (15,454 bin50 spatial bins), using a
two-step bridging strategy through an intermediate scRNA-seq dataset from
the same mouse and 4 mouse references (Spleen, BM, Neutro, CITE — total 70
cell type scores).

> **Restricted data.** The Stereo-seq, intermediate scRNA-seq, and SPLINTR
> barcode files are unpublished collaborator data (Dawson lab, Peter
> MacCallum Cancer Centre). They are **not redistributable**. Set
> `dawson_root` in `config/paths_local.yml` only if you have authorised
> access. Without it, this case study cannot be reproduced; all other case
> studies can.

## Inputs (set via `config/paths_local.yml`)

Restricted (under `${dawson_root}`):
- `${dawson_root}/Mouse2/SS200000334TR_A2.gef`
- `${dawson_root}/Mouse2/SS200000334TR_A2.lasso.mouse2_whole_tissue.gem.gz`
- `${dawson_root}/Mouse2/DP8400028335BL_L01_read_1.unmapped_reads.counts.tsv`
- `${dawson_root}/mouse4_scRNAseq_sce.qs` (intermediate scRNA-seq, 5,595 cells)
- `${dawson_root}/scRNA-seq/refSpleen_processed.qs` (mouse spleen ref, 18 cell types)
- `${dawson_root}/scRNA-seq/refBM_processed.qs` (mouse BM ref, 13 cell types)
- `${dawson_root}/output/barcode-seq_clustering.rds` (pre-computed meta-clone clustering)

Public (under `${data_root}`):
- `${data_root}/Neutrophils/GSE244536/scRNA-seq/reference.qs` (neutrophil ref, 11 cell types)
- `${data_root}/Stereo-CITE/CITE-seq/refComboRNA.qs` (CITE-seq spleen ref, 28 cell types) —
  derived from `${data_root}/Stereo-CITE/CITE-seq/spleen_lymph_{206,111}.h5ad`

## Run order

```
00_prepareData.ipynb                Python preprocessing: GEF/GEM -> bin50 SCE; barcode counts.
01_prepareSCE.R                     Assemble Stereo-seq SCE; build CITE-seq combined reference
                                    (refComboRNA.qs from spleen_lymph_206.h5ad + spleen_lymph_111.h5ad).
02_runPhiSpace_intermediate.R       Annotate intermediate scRNA-seq with the 4 references
                                    (produces 70-column score matrix).
03_bridge_main.Rmd                  Master notebook: bridge to Stereo-seq via the intermediate
                                    scRNA-seq; clustering; enrichment; clone composition.
                                    Generates Fig 5 (B,C,D panels).
04_benchXyRange.R                   Sensitivity analysis on the spatial-enhancement xyRange parameter.
05_multiRefConsistency.R            REVISION (Reviewer 2 C4): cross-reference consistency
                                    of overlapping cell types. Supp Figs S14-S15.
06_metaCloneStates.R                REVISION (Reviewer 2 C5): quantitative cell type proportions
                                    per meta-clone + UCell scoring of 11 MSigDB Hallmark
                                    signatures. Supp Figs S11-S13.
```

## Bridging strategy

Unlike Visium / CosMx / Xenium (which annotate the spatial query directly
from references), this case study uses **two-step bridging**:

1. Annotate intermediate scRNA-seq (5,595 cells) with the 4 references →
   70-column PhiSpace score matrix.
2. Use intermediate scRNA-seq as `response` to annotate the Stereo-seq query
   via `PhiSpaceR_1ref(..., response = PhiSpace_scores, refAssay = "log1p")`.

This is necessary because the references come from different mouse tissue
contexts; the intermediate scRNA-seq from the same mouse provides a
gene-space-matched bridge.

## Outputs

- `${dawson_root}/output/stereo_bin50_PhiSpace.qs` — Stereo-seq with 70-col PhiSpace scores
- `${output_root}/barAssay_xyRange500.rds` — cached spatially enhanced barcode-seq assay
- `${output_root}/ucell_hallmark_scores.rds` — UCell signature scores (15454 × 11)
- `figs/clusts/`, `figs/cTypeEnrichBoxplots/`, `figs/divergence/`,
  `figs/heatmaps/`, `figs/spatialHeatmap/`, `figs/multiRef/`, `figs/metaClone/`

Manuscript items supported: Fig 5 (B clustering, C DWD, D clone composition
+ enriched cell types); Supp Figs S9-S15; Supp Table S2.
