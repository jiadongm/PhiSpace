# PhiSpace ST — paper analysis code

Source code for the paper:

> **PhiSpace ST: a platform-agnostic method to identify cell states in spatial transcriptomics studies**
> Jiadong Mao, Jarny Choi, Kim-Anh Lê Cao
> *Cell Reports Methods* (in press; Zenodo DOI to be inserted post-archival)

This directory contains the analysis pipelines that reproduce every main and
supplemental figure of the paper. The PhiSpace R package itself lives in
`../pkg/` of the same GitHub repository
(<https://github.com/jiadongm/PhiSpace>).

## What this reproduces

Four spatial transcriptomics platform case studies plus two supplementary
analyses:

| Section | Manuscript figure | Folder |
|---|---|---|
| 2.2 | Fig 2 — 10x Visium NSCLC (18 samples) | `scripts/01_visium_case/` |
| 2.3 | Fig 3 — Nanostring CosMx NSCLC (8 samples) | `scripts/02_cosmx_case/` |
| 2.4 | Fig 4 — 10x Xenium FFPE NSCLC (1 sample) | `scripts/03_xenium_case/` |
| 2.5 | Fig 5 — BGI Stereo-seq AML mouse spleen | `scripts/04_stereoseq_case/` |
| Supp | Reviewer 2 C3: pDC–fibroblast co-presence (Visium 40 + CosMx) | `scripts/05_pdc_fibroblast_supp/` |
| Supp | Supp Fig S1 — 7-method deconvolution benchmark on 32 simulated datasets | `scripts/06_benchmarking_supp/` |

See `outputs_manifest/expected_outputs.tsv` for the full figure-to-script map.

## Environment

- HPC cluster (developed on Univ. Melbourne *Spartan*, Linux). All R analyses
  also run on a workstation with the same package versions.
- **R 4.4.0** (foss/2022a toolchain on Spartan). Benchmark scripts also tested
  on R 4.5.0.
- **Python 3.10** for Cell2location / TACCO competitor runs (each in its own
  conda env — see `environment/`). Python 3.8 for the Stereo-seq prep notebook.
- Bioconductor 3.18 for the R 4.4.0 builds.

Module loads used on Spartan are in `hpc/module_loads.sh`. Generic SLURM
templates are in `hpc/slurm_templates/`.

## Installing dependencies

R packages required by the analyses are listed in
`environment/package_versions.tsv`. The PhiSpace package is installed from
GitHub:

```r
install.packages("devtools")
devtools::install_github("jiadongm/PhiSpace/pkg")
```

Bioconductor + CRAN dependencies can be installed in one go:

```r
install.packages(c("dplyr","tidyr","magrittr","ggplot2","ggpubr","qs","yaml",
                   "Matrix","Seurat","sctransform","circlize","seriation",
                   "kerndwd","msigdbr","devtools"))
BiocManager::install(c("SingleCellExperiment","SummarizedExperiment","scran",
                       "scuttle","SpaNorm","DropletUtils","ComplexHeatmap","UCell"))
remotes::install_github("mojaveazure/seurat-disk")
```

For competitor benchmarking (`scripts/06_benchmarking_supp/`), each method
needs its own conda environment; YAMLs are in `environment/`:
- `environment/cell2location_env.yml`
- `environment/tacco_env.yml`

## Setup before running anything

```bash
# 1. Clone the repository and cd into this directory.
git clone https://github.com/jiadongm/PhiSpace.git
cd PhiSpace/sourceCode_PhiSpaceST

# 2. Copy the path config template and edit it.
cp config/paths_template.yml config/paths_local.yml
# Edit phispace_data_root (required), dawson_root + cell2loc_python +
# tacco_python (optional, only if running those parts).

# 3. (Optional) Capture your environment for reproducibility records.
Rscript environment/capture_sessionInfo.R
```

`config/paths_local.yml` is git-ignored, so per-environment values stay local.
All scripts assume the working directory is this `sourceCode_PhiSpaceST/` root
and source `scripts/00_setup/paths_loader.R` to pick up your paths.

## Running the analyses

See `run_order.md` for the recommended end-to-end order. Each case-study
folder also has its own `README.md` with the script order and expected outputs
for that platform.

**Quick start** for one case study (Visium):

```bash
Rscript scripts/01_visium_case/00_prepareData.R
Rscript scripts/01_visium_case/01_azimuthLungMarkers.R
Rscript scripts/01_visium_case/02_runPhiSpace.R
# Then knit scripts/01_visium_case/05_visiumNSCLC_main.Rmd to regenerate Fig 2.
```

## Data availability

All datasets are listed in `config/dataset_manifest.tsv` with their
licence and redistribution status. Summary:

- **Public datasets** (HLCA / Azimuth Lung 2.0; Natri et al. 2024 GSE227136;
  10x Visium NSCLC; Nanostring CosMx NSCLC; 10x Xenium FFPE NSCLC; mouse
  neutrophil GSE244536; BLAST 32-dataset benchmark): downloaded from the
  source on first run, **not redistributed** in this repository.
- **Restricted (collaborator-only) AML mouse spleen Stereo-seq data**: not
  redistributable. The Stereo-seq scripts in `scripts/04_stereoseq_case/`
  need `dawson_root` set in `paths_local.yml`. Without access, this case
  study cannot be reproduced; all other case studies can.

## What is regenerated vs precomputed

Everything under `figs/`, `output/`, and per-script result files (`*.qs`,
`*.rds`) is **regenerated** when you run the scripts. None of these are
shipped with the repository (they are either large or trivially
reproducible from the public inputs). Only the source scripts, configs,
manifests, and metadata files are in the GitHub release.

## Compute / runtime notes

| Pipeline | Approx. runtime | Memory |
|---|---|---|
| Visium 18 samples (PhiSpace) | < 5 min total | 8–16 GB |
| CosMx 8 samples (PhiSpace, 4 refs) | ~20–40 min | 32–64 GB |
| Xenium 1 sample (4 refs, 161k cells) | ~16 s for PhiSpace; ~10 min total incl. clustering | 32 GB |
| Stereo-seq bridged annotation | ~5–10 min | 32 GB |
| pDC analysis (40 Visium samples) | ~30 min | 32 GB |
| 7-method deconvolution benchmark (32 datasets) | hours–days, dominated by Cell2location (~17 min/dataset on GPU) and RCTD (~5 min/dataset) | 64 GB; GPU recommended for Cell2location |

## Repository layout

```
sourceCode_PhiSpaceST/
├── README.md                         (this file)
├── run_order.md                      (recommended end-to-end run order)
├── .gitignore
├── config/
│   ├── paths_template.yml            (copy to paths_local.yml; git-ignored)
│   └── dataset_manifest.tsv          (every reference + query dataset)
├── environment/
│   ├── capture_sessionInfo.R         (run to populate sessionInfo.txt)
│   ├── sessionInfo.txt               (placeholder until you run the above)
│   ├── package_versions.tsv          (R + Python deps with versions)
│   ├── cell2location_env.yml
│   └── tacco_env.yml
├── hpc/
│   ├── module_loads.sh               (Spartan module loads)
│   └── slurm_templates/
│       ├── template_R.slurm
│       └── template_python_gpu.slurm
├── outputs_manifest/
│   └── expected_outputs.tsv          (every figure/table → generating script)
└── scripts/
    ├── 00_setup/                     (paths_loader.R, paths_loader.py)
    ├── 01_visium_case/               (Fig 2)
    ├── 02_cosmx_case/                (Fig 3)
    ├── 03_xenium_case/               (Fig 4)
    ├── 04_stereoseq_case/            (Fig 5)
    ├── 05_pdc_fibroblast_supp/       (Reviewer 2 C3 supplementary)
    └── 06_benchmarking_supp/         (Supp Fig S1)
```

## Citing

If you use this code, please cite both the paper and this archived release:

- **Paper**: Mao J., Choi J., Lê Cao K.-A. (2026). *PhiSpace ST: …*. Cell Reports Methods. DOI: TBD.
- **Software (this archive)**: see `CITATION.cff` at the GitHub repository
  root, and the Zenodo DOI minted from the GitHub release tag `v1.1.0`.

## License

The PhiSpace R package and this analysis code are licensed under AGPL-3.0
(see `../LICENSE` at the GitHub repository root).
