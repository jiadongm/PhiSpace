# Release candidate checklist — PhiSpace ST `v1.1.0`

Run through this checklist before tagging the GitHub release that Zenodo
will archive. Per the spec in `../ZENODO_HPC_PACKAGING_INSTRUCTIONS.md`.

## Identity

- **Release tag**: `v1.1.0`
- **PhiSpace package version (`pkg/DESCRIPTION`)**: confirm matches the tag
  (or document the compatible commit). Update in the GitHub clone, not on Spartan.
- **Final commit SHA**: _(fill in at tag time)_
- **Date**: _(fill in at tag time; replace placeholder `2026-05-01` in
  `release_metadata/CITATION.cff`)_

## Contents — included in the release

- [x] `pkg/` — PhiSpace R package source (DESCRIPTION, NAMESPACE, R/, man/, vignettes/ if present)
      → lives in the GitHub repo, NOT in the Spartan working tree.
- [x] `sourceCode_PhiSpaceST/` — paper analysis code (this directory)
- [x] `LICENSE` — AGPL-3.0 (already at GitHub repo root)
- [x] `README.md` (top-level) — already at GitHub repo root
- [ ] `CITATION.cff` — move from `sourceCode_PhiSpaceST/release_metadata/` to repo root
- [ ] `.zenodo.json` — move from `sourceCode_PhiSpaceST/release_metadata/` to repo root

## Contents — sourceCode_PhiSpaceST/ self-check

- [x] `README.md` — top-level, paper title, env, install, run order
- [x] `run_order.md` — recommended end-to-end run order
- [x] `.gitignore` — ignores `paths_local.yml`, `.Rhistory`, caches, etc.
- [x] `config/paths_template.yml` — placeholders only (no real paths)
- [x] `config/dataset_manifest.tsv` — every reference + query dataset
- [x] `environment/package_versions.tsv` — R + Python deps with versions
- [x] `environment/cell2location_env.yml`, `environment/tacco_env.yml`
- [x] `environment/capture_sessionInfo.R` + `sessionInfo.txt` (placeholder)
- [ ] **Run `Rscript environment/capture_sessionInfo.R`** on Spartan to populate
      `sessionInfo.txt` with the actual paper R session.
- [x] `hpc/module_loads.sh` + `hpc/slurm_templates/{template_R,template_python_gpu}.slurm`
      (sanitised; no `--account`, no `--mail-user`)
- [x] `outputs_manifest/expected_outputs.tsv` — figure/table → script map
- [x] `scripts/00_setup/{paths_loader.R, paths_loader.py, README.md}`
- [x] `scripts/01_visium_case/` (12 files) + `README.md`
- [x] `scripts/02_cosmx_case/` (5 files) + `README.md`
- [x] `scripts/03_xenium_case/` (7 files) + `README.md`
- [x] `scripts/04_stereoseq_case/` (9 files) + `README.md`
- [x] `scripts/05_pdc_fibroblast_supp/` (4 files) + `README.md`
- [x] `scripts/06_benchmarking_supp/` (22 files) + `README.md`

## Contents — explicitly EXCLUDED from the release

- [x] AML mouse spleen Stereo-seq raw / processed / barcode files (Dawson lab,
      RESTRICTED — see `dataset_manifest.tsv` rows tagged `redistribution_allowed=no`
      and `dawson_root` entries).
- [x] Manuscript LaTeX, KRT, response letter, cover letter, journal PDFs
      (`PhiSpaceST_FinalSubmission/` — separate manuscript artefact, NOT code).
- [x] Per-case `figs/`, `output/`, `data/` payloads (regenerable; gitignored).
- [x] BLAST upstream tree: `Codes/`, `Benchmarking/`, `Extenrnal/`,
      `ExampleData/`, `FigureData/`, upstream `.git/` (cite Li et al. 2022).
- [x] Benchmark `prepared_data/` (~71 GB; regenerable from h5ad sources).
- [x] All `.qs`/`.rds` intermediate caches (regeneratable; `.gitignore`d).
- [x] `.Rhistory`, `.RData`, `.ipynb_checkpoints/`, scratch dirs (`.gitignore`d).
- [x] Five legacy SLURM templates + their unshipped R scripts (per author
      decision — `runDEG.R`, `runRCTD4Refs_Lung5Rep1.R`, `timeRCTD.R`,
      `runDiProPerm.R`, `run_pseudobulking.R`).

## Sanitisation audit (last run during Step 6)

- [x] 0 hard-coded `/data/gpfs/...` or `/data/projects/...` paths in code
- [x] 0 `~/PhiSpaceR/...` HOME refs (1 commented-out doc note retained)
- [x] 0 `setwd()` calls
- [x] 0 `punim*` project account IDs
- [x] 0 personal email addresses (only Seurat `@meta.data` accessor matches)
- [x] 0 tokens / passwords / API keys / secrets
- [x] No SLURM personal info (account ID, email) in templates
- [x] `jiadongm/PhiSpace/pkg` GitHub install ref is intentional (author's GitHub handle)

## Integrity checks (last run during Step 11)

- [x] All 35 R scripts parse without errors
- [x] All 3 Rmd files purl + parse without errors
- [x] All 15 Python scripts pass `python3 -m py_compile`
- [x] `paths_loader.R` and `paths_loader.py` self-test (load template, expose
      7 keys: `phispace_data_root`, `dawson_root`, `cell2loc_python`,
      `tacco_python`, `project_root`, `data_root`, `output_root`)
- [ ] Spot-check one or two end-to-end script runs on Spartan with real
      `paths_local.yml` (recommended but not required for tagging).

## Open items / follow-up

- **Other authors' ORCIDs**: Choi and Lê Cao are name-only in
  `CITATION.cff` / `.zenodo.json`. Add `orcid:` lines if/when available.
- **Manuscript abstract**: replaced with the current manuscript Summary on
  2026-05-01. Re-check only if the Summary changes before release.
- **Manuscript DOI**: the `preferred-citation.doi` field in `CITATION.cff` is
  empty — fill in once Cell Reports Methods assigns the DOI.
- **Dataset accessions**: all `TBD` cells were removed from
  `config/dataset_manifest.tsv`. Some mouse reference accessions remain
  described as "not specified in manuscript"; verify external repositories only
  if more detailed provenance is desired.
- **Supplement numbering**: `outputs_manifest/expected_outputs.tsv` was updated
  against the current `LaTeX_source/supplemental.tex` order (S1-S15, Tables
  S1-S2).
- **Run `capture_sessionInfo.R`** on Spartan with all packages loaded so
  `environment/sessionInfo.txt` reflects the real paper environment.

## Zenodo release procedure (post-tag)

1. Confirm `release_metadata/CITATION.cff` and `.zenodo.json` are at the
   GitHub repo root (not under `sourceCode_PhiSpaceST/`).
2. Enable Zenodo GitHub integration for `jiadongm/PhiSpace`.
3. Create the GitHub release from the target commit (tag `v1.1.0`).
4. Wait for Zenodo to archive the release; copy the DOI from the Zenodo record.
5. Confirm Zenodo record title, author list, license, and repository URL.
6. Insert the DOI into:
   - `LaTeX_source/main.tex` (Data and Code Availability statement)
   - `Key resources table.docx` (PhiSpace row)
   - `release_metadata/CITATION.cff` (add a `doi:` field; commit; re-tag a
     patch release if you want Zenodo to mint a second version with the
     self-referencing DOI).
7. Rebuild / spot-check the manuscript PDF and KRT after DOI insertion.

> **Important.** Do *not* reuse the older Zenodo DOI for the prior PhiSpace
> multi-omics paper. Wait for the new Zenodo record to mint its own DOI for
> the PhiSpace ST software release.
