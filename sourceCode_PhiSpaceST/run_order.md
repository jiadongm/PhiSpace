# Recommended run order

Each case study is self-contained — you can run them in any order or skip
some entirely. Within a case study, follow the numeric prefix of each script
(`00_…`, `01_…`, etc.). The figure-to-script mapping is in
`outputs_manifest/expected_outputs.tsv`.

## End-to-end (everything)

```text
0.  One-time setup
    cp config/paths_template.yml config/paths_local.yml   # then edit
    Rscript environment/capture_sessionInfo.R              # optional but recommended

1.  Visium NSCLC (Section 2.2 / Fig 2)
    Rscript scripts/01_visium_case/00_prepareData.R
    Rscript scripts/01_visium_case/01_azimuthLungMarkers.R
    Rscript scripts/01_visium_case/02_runPhiSpace.R
    Rscript scripts/01_visium_case/03_runtimeBenchmark.R
    Rscript scripts/01_visium_case/04_DEGAnalysis.R
    Rscript -e 'rmarkdown::render("scripts/01_visium_case/05_visiumNSCLC_main.Rmd")'
    # Optional competitor runs: scripts/01_visium_case/competitors/*

2.  CosMx NSCLC (Section 2.3 / Fig 3)
    Rscript scripts/02_cosmx_case/00_prepareSCE.R
    Rscript scripts/02_cosmx_case/01_runPhiSpace4Refs.R
    Rscript -e 'rmarkdown::render("scripts/02_cosmx_case/02_cosmx4Refs_main.Rmd")'
    Rscript scripts/02_cosmx_case/03_regenFig3C.R   # only if updating Fig 3C colours

3.  Xenium FFPE NSCLC (Section 2.4 / Fig 4)
    Rscript scripts/03_xenium_case/00_prepareData.R
    "$CELL2LOC_PYTHON" scripts/03_xenium_case/01_assignPathology.py   # any python with shapely
    Rscript scripts/03_xenium_case/02_runPhiSpace4Refs.R
    Rscript scripts/03_xenium_case/03_nicheAnalysis.R
    Rscript scripts/03_xenium_case/04_clusterEnrichment.R
    Rscript scripts/03_xenium_case/05_coPresenceAnalysis_k8.R
    Rscript scripts/03_xenium_case/06_coPresenceAnalysis_multiK.R   # k=4..10 sensitivity

4.  Stereo-seq AML mouse spleen (Section 2.5 / Fig 5)
    # Requires `dawson_root` set in paths_local.yml (collaborator-restricted data).
    jupyter nbconvert --execute scripts/04_stereoseq_case/00_prepareData.ipynb
    Rscript scripts/04_stereoseq_case/01_prepareSCE.R
    Rscript scripts/04_stereoseq_case/02_runPhiSpace_intermediate.R
    Rscript -e 'rmarkdown::render("scripts/04_stereoseq_case/03_bridge_main.Rmd")'
    Rscript scripts/04_stereoseq_case/04_benchXyRange.R
    Rscript scripts/04_stereoseq_case/05_multiRefConsistency.R    # Reviewer 2 C4
    Rscript scripts/04_stereoseq_case/06_metaCloneStates.R        # Reviewer 2 C5

5.  pDC–fibroblast supplementary (Reviewer 2 C3)
    Rscript scripts/05_pdc_fibroblast_supp/00_prepareData_visium40.R
    Rscript scripts/05_pdc_fibroblast_supp/01_runPhiSpace_visium.R
    Rscript scripts/05_pdc_fibroblast_supp/02_visium_pDCFibroblast.R
    Rscript scripts/05_pdc_fibroblast_supp/03_cosmx_pDCFibroblast.R

6.  Deconvolution benchmark (Supp Fig S1)
    bash scripts/06_benchmarking_supp/run_benchmark.sh
    # or per-method orchestration:
    python scripts/06_benchmarking_supp/04_run_full_benchmark.py        --datasets all
    python scripts/06_benchmarking_supp/04_run_rctd_benchmark.py        --datasets all
    python scripts/06_benchmarking_supp/04_run_seurat_benchmark.py      --datasets all
    python scripts/06_benchmarking_supp/04_run_spotlight_benchmark.py   --datasets all
    "$TACCO_PYTHON"     scripts/06_benchmarking_supp/04_run_tacco_benchmark.py         --datasets all
    "$CELL2LOC_PYTHON"  scripts/06_benchmarking_supp/04_run_cell2location_benchmark.py --datasets all
    Rscript scripts/06_benchmarking_supp/09_visualise_results.R
```

## Five canonical steps (per the Cell Press packaging instructions)

1. **Prepare reference datasets.** HLCA/Azimuth + Natri et al. lineage
   references are processed by `scripts/02_cosmx_case/00_prepareSCE.R` (the
   four lineage `.qs` files are reused by the Xenium case study). Mouse
   refs are loaded directly from `${data_root}` and `${dawson_root}`.
2. **Prepare each spatial transcriptomics query dataset.** Visium / CosMx /
   Xenium / Stereo-seq queries are built by their respective `00_*` (or
   `01_*` for Stereo-seq R + `00_*` for Stereo-seq Python prep) scripts.
3. **Train/apply PhiSpace models.** `runPhiSpace*.R` per case study.
4. **Run downstream analyses.** Niche analysis, cluster enrichment,
   co-presence, divergence, and meta-clone state scoring per case study.
5. **Generate figures and tables.** Master `*_main.Rmd` per case study, plus
   `regenFig3C.R` (CosMx Fig 3C colours) and
   `06_benchmarking_supp/09_visualise_results.R` (Supp Fig S1).
