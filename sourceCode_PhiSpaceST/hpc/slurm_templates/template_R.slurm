#!/bin/bash
# Generic R analysis SLURM template for the PhiSpace ST pipeline.
#
# To use:
#   1. Copy this file and rename to describe the job (e.g. run_visium_phispace.slurm).
#   2. Set --account, --mail-user, --time, --mem, and --cpus-per-task for your environment.
#   3. Set SCRIPT to the analysis script (path is relative to the
#      sourceCode_PhiSpaceST/ root, which is set as the working directory).
#   4. Submit with `sbatch hpc/slurm_templates/<your_copy>.slurm` from the
#      sourceCode_PhiSpaceST/ root.

#SBATCH --job-name=phispace_st
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=64GB
#SBATCH --output=slurm-%j.out

# REQUIRED: replace with your project allocation and (optionally) email
#SBATCH --account=YOUR_PROJECT_ALLOCATION
##SBATCH --mail-user=you@example.org
##SBATCH --mail-type=END,FAIL

# Resolve project root (the directory containing this template's parent).
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$PROJECT_ROOT"

# Load HPC modules (Spartan defaults — edit hpc/module_loads.sh for your cluster).
source hpc/module_loads.sh

# REQUIRED: set this to the analysis script you want to run.
SCRIPT="scripts/01_visium_case/02_runPhiSpace.R"

echo "PROJECT_ROOT=$PROJECT_ROOT"
echo "Running: $SCRIPT"
Rscript "$SCRIPT"
