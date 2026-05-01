#!/bin/bash

# =============================================================================
# PhiSpace Benchmarking Pipeline
# =============================================================================
# This script runs the full benchmarking pipeline for PhiSpace against
# other deconvolution methods (RCTD, Cell2location, SPOTlight, Tangram, Seurat)
# using the 32 simulated spatial transcriptomics datasets.
#
# Usage:
#   ./run_benchmark.sh [OPTIONS]
#
# Options:
#   --datasets NUMS   Comma-separated dataset numbers or "all" (default: all)
#   --skip-prepare    Skip data preparation step
#   --skip-phispace   Skip PhiSpace execution step
#   --only-metrics    Only calculate metrics (skip prepare and phispace)
#   --only-visualize  Only create visualizations
#   --help            Show this help message
# =============================================================================

set -e

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$( dirname "$SCRIPT_DIR" )"

# Default parameters
DATASETS="all"
SKIP_PREPARE=false
SKIP_PHISPACE=false
ONLY_METRICS=false
ONLY_VISUALIZE=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --datasets)
            DATASETS="$2"
            shift 2
            ;;
        --skip-prepare)
            SKIP_PREPARE=true
            shift
            ;;
        --skip-phispace)
            SKIP_PHISPACE=true
            shift
            ;;
        --only-metrics)
            ONLY_METRICS=true
            shift
            ;;
        --only-visualize)
            ONLY_VISUALIZE=true
            shift
            ;;
        --help)
            head -30 "$0" | tail -n +2
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

echo "============================================================"
echo "PhiSpace Benchmarking Pipeline"
echo "============================================================"
echo "Script directory: $SCRIPT_DIR"
echo "Base directory: $BASE_DIR"
echo "Datasets: $DATASETS"
echo "============================================================"

# Function to check if a command exists
check_command() {
    if ! command -v $1 &> /dev/null; then
        echo "ERROR: $1 is not installed or not in PATH"
        exit 1
    fi
}

# Check required dependencies
echo ""
echo "Checking dependencies..."
check_command python
check_command Rscript

# Verify Python packages
python -c "import pandas; import numpy; import scanpy" 2>/dev/null || {
    echo "ERROR: Required Python packages not found (pandas, numpy, scanpy)"
    exit 1
}
echo "  Python dependencies: OK"

# Verify R packages
Rscript -e "library(PhiSpace); library(SingleCellExperiment)" 2>/dev/null || {
    echo "ERROR: Required R packages not found (PhiSpace, SingleCellExperiment)"
    exit 1
}
echo "  R dependencies: OK"

# Create output directories
mkdir -p "$SCRIPT_DIR/prepared_data"
mkdir -p "$SCRIPT_DIR/results"
mkdir -p "$SCRIPT_DIR/integrated_results"
mkdir -p "$SCRIPT_DIR/figures"
mkdir -p "$SCRIPT_DIR/runtime_analysis"

if [ "$ONLY_VISUALIZE" = true ]; then
    echo ""
    echo "============================================================"
    echo "Running visualization only..."
    echo "============================================================"

    python "$SCRIPT_DIR/06_visualize_results.py" \
        --results_dir "PhiSpace_Benchmarking/integrated_results" \
        --output_dir "PhiSpace_Benchmarking/figures"

    python "$SCRIPT_DIR/07_runtime_benchmark.py" \
        --results_dir "PhiSpace_Benchmarking/results" \
        --output_dir "PhiSpace_Benchmarking/runtime_analysis"

    echo ""
    echo "Visualization completed!"
    exit 0
fi

# Step 1: Prepare data
if [ "$SKIP_PREPARE" = false ] && [ "$ONLY_METRICS" = false ]; then
    echo ""
    echo "============================================================"
    echo "STEP 1: Preparing data"
    echo "============================================================"

    python "$SCRIPT_DIR/01_prepare_data.py" \
        --data_dir "ExampleData/SimualtedSpatalData" \
        --output_dir "PhiSpace_Benchmarking/prepared_data" \
        --datasets "$DATASETS"
fi

# Step 2: Run PhiSpace
if [ "$SKIP_PHISPACE" = false ] && [ "$ONLY_METRICS" = false ]; then
    echo ""
    echo "============================================================"
    echo "STEP 2: Running PhiSpace"
    echo "============================================================"

    python "$SCRIPT_DIR/04_run_full_benchmark.py" \
        --data_dir "ExampleData/SimualtedSpatalData" \
        --output_dir "PhiSpace_Benchmarking/results" \
        --datasets "$DATASETS" \
        --skip_prepare
fi

# Step 3: Integrate results
echo ""
echo "============================================================"
echo "STEP 3: Integrating results"
echo "============================================================"

python "$SCRIPT_DIR/05_integrate_results.py" \
    --figure_data_dir "FigureData/Figure4/32SimulationData" \
    --phispace_results_dir "PhiSpace_Benchmarking/results" \
    --output_dir "PhiSpace_Benchmarking/integrated_results" \
    --datasets "$DATASETS" \
    --methods "RCTD,Cell2location,SPOTlight,Tangram,Seurat"

# Step 4: Create visualizations
echo ""
echo "============================================================"
echo "STEP 4: Creating visualizations"
echo "============================================================"

python "$SCRIPT_DIR/06_visualize_results.py" \
    --results_dir "PhiSpace_Benchmarking/integrated_results" \
    --output_dir "PhiSpace_Benchmarking/figures"

# Step 5: Runtime analysis
echo ""
echo "============================================================"
echo "STEP 5: Runtime analysis"
echo "============================================================"

python "$SCRIPT_DIR/07_runtime_benchmark.py" \
    --results_dir "PhiSpace_Benchmarking/results" \
    --output_dir "PhiSpace_Benchmarking/runtime_analysis"

echo ""
echo "============================================================"
echo "BENCHMARKING PIPELINE COMPLETED"
echo "============================================================"
echo ""
echo "Results saved to:"
echo "  - PhiSpace results: $SCRIPT_DIR/results/"
echo "  - Integrated results: $SCRIPT_DIR/integrated_results/"
echo "  - Figures: $SCRIPT_DIR/figures/"
echo "  - Runtime analysis: $SCRIPT_DIR/runtime_analysis/"
echo ""
echo "Key output files:"
echo "  - Metric comparison: integrated_results/summary_mean.csv"
echo "  - Method rankings: integrated_results/method_rankings.csv"
echo "  - Boxplot comparison: figures/boxplot_comparison.pdf"
echo "  - Runtime summary: runtime_analysis/runtime_summary.csv"
echo ""
