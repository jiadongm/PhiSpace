# 00_setup

Shared infrastructure used by every case-study script. You should not need to
edit anything here for a typical run.

## Files

| File | Purpose |
|---|---|
| `paths_loader.R` | Sources at the top of every R / Rmd script. Walks up from the working directory to find `config/paths_template.yml`, reads `config/paths_local.yml` if it exists (else the template), and exposes a `paths` named list. |
| `paths_loader.py` | Python equivalent of the above. Importable as `from paths_loader import paths`. |

## First-time setup

From the `sourceCode_PhiSpaceST/` root:

```bash
cp config/paths_template.yml config/paths_local.yml
# Edit config/paths_local.yml — set phispace_data_root (and dawson_root /
# cell2loc_python / tacco_python if you plan to run those parts).
```

`config/paths_local.yml` is git-ignored, so per-environment values stay local.

## Path keys exposed

| Key | Set by | Used by |
|---|---|---|
| `paths$phispace_data_root` | required (in YAML) | almost every script |
| `paths$data_root` | derived (`${phispace_data_root}/data`) | reference + query data loading |
| `paths$output_root` | derived (`${phispace_data_root}/output`) | analysis outputs / caches |
| `paths$dawson_root` | required for `scripts/04_stereoseq_case/` | AML mouse spleen Stereo-seq + intermediate scRNA-seq |
| `paths$cell2loc_python` | required for `scripts/06_benchmarking_supp/` Cell2location runs | `04_run_cell2location_benchmark.py`, `02_run_cell2location.py` |
| `paths$tacco_python` | required for `scripts/06_benchmarking_supp/` TACCO runs | `04_run_tacco_benchmark.py`, `02_run_tacco.py` |
| `paths$project_root` | injected by loader | rarely used directly |
