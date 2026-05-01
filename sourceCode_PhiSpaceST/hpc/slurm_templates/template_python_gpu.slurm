#!/bin/bash
# Generic GPU SLURM template for Python competitor methods (Cell2location).
# Cell2location's variational inference runs much faster on GPU; TACCO and
# the prepareData notebooks are CPU-only and can use template_R.slurm with
# the relevant python interpreter substituted.

#SBATCH --job-name=cell2loc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --partition=gpgpu
#SBATCH --time=12:00:00
#SBATCH --mem=64GB
#SBATCH --output=slurm-%j.out

# REQUIRED: replace with your project allocation and (optionally) email
#SBATCH --account=YOUR_PROJECT_ALLOCATION
##SBATCH --mail-user=you@example.org
##SBATCH --mail-type=END,FAIL

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$PROJECT_ROOT"

# Load CUDA / CUDNN modules required by your conda env (edit per cluster).
module purge
# module load cuda/12.2 cudnn/8.9
# module load Anaconda3

# CELL2LOC_PYTHON should point to the python inside the cell2location conda env.
# Either export it here, or set `cell2loc_python` in config/paths_local.yml — the
# benchmark scripts will read it from there via paths_loader.py.
: "${CELL2LOC_PYTHON:?Set CELL2LOC_PYTHON or define cell2loc_python in config/paths_local.yml}"

# Example: run all 32 benchmark datasets through cell2location.
"$CELL2LOC_PYTHON" scripts/06_benchmarking_supp/04_run_cell2location_benchmark.py --datasets all
