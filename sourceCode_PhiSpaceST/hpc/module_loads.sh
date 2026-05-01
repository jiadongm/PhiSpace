#!/bin/bash
# HPC module loads used for the PhiSpace ST analyses on Spartan (Univ. Melbourne).
# Adapt as needed for other clusters. Source this script before running R-based
# scripts in scripts/01_..04_, scripts/05_pdc_fibroblast_supp/, or the R
# portions of scripts/06_benchmarking_supp/.

module purge
module load foss/2022a
module load R/4.4.0

# For Python steps (Cell2location, TACCO, Stereo-seq prepareData), activate the
# corresponding conda environment instead — see environment/cell2location_env.yml
# and environment/tacco_env.yml. Their python paths must be set in
# config/paths_local.yml as `cell2loc_python` and `tacco_python`.
