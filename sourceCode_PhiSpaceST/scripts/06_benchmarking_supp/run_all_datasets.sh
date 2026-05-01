#!/bin/bash

# Run PhiSpace benchmarking on all 32 datasets
# Configuration: sctransform (reference) + spanorm (query)

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$( dirname "$SCRIPT_DIR" )"

# Python with scanpy (cell2location conda env)
# Cell2location conda env Python — set CELL2LOC_PYTHON in your environment.
PYTHON="${CELL2LOC_PYTHON:?Set CELL2LOC_PYTHON to the cell2location conda env python}"

REF_NORM="sctransform"
QUERY_NORM="spanorm"

echo "============================================================"
echo "PhiSpace Full Benchmark: $REF_NORM + $QUERY_NORM"
echo "============================================================"
echo "Start time: $(date)"
echo ""

# Create results directory
RESULTS_DIR="$SCRIPT_DIR/results_${REF_NORM}_${QUERY_NORM}"
mkdir -p "$RESULTS_DIR"

# Track timing
TOTAL_START=$(date +%s)

# Process all 32 datasets
for DS_NUM in $(seq 1 32); do
    echo ""
    echo "============================================================"
    echo "Processing dataset $DS_NUM / 32"
    echo "============================================================"

    DS_START=$(date +%s)

    # Check if prepared data exists
    PREPARED_DIR="$SCRIPT_DIR/prepared_data/dataset${DS_NUM}"
    if [ ! -d "$PREPARED_DIR" ]; then
        echo "Preparing dataset $DS_NUM..."
        cd "$BASE_DIR"
        $PYTHON "$SCRIPT_DIR/01_prepare_data.py" \
            --data_dir "ExampleData/SimualtedSpatalData" \
            --output_dir "PhiSpace_Benchmarking/prepared_data" \
            --datasets "$DS_NUM"
    fi

    # Run PhiSpace
    OUTPUT_DIR="$RESULTS_DIR/dataset${DS_NUM}"
    mkdir -p "$OUTPUT_DIR"

    # Skip if result already exists
    if [ -f "$OUTPUT_DIR/PhiSpace_result.txt" ]; then
        echo "Skipping dataset $DS_NUM - result already exists"
        continue
    fi

    echo "Running PhiSpace on dataset $DS_NUM..."
    Rscript "$SCRIPT_DIR/02_run_phispace.R" \
        "$PREPARED_DIR" \
        "$OUTPUT_DIR" \
        "$REF_NORM" \
        "$QUERY_NORM" 2>&1 | tee "$OUTPUT_DIR/log.txt"

    DS_END=$(date +%s)
    DS_TIME=$((DS_END - DS_START))
    echo "Dataset $DS_NUM completed in ${DS_TIME}s"

    # Check if result file exists
    if [ -f "$OUTPUT_DIR/PhiSpace_result.txt" ]; then
        echo "SUCCESS: Result file created"
    else
        echo "WARNING: Result file not found"
    fi
done

TOTAL_END=$(date +%s)
TOTAL_TIME=$((TOTAL_END - TOTAL_START))

echo ""
echo "============================================================"
echo "Full Benchmark Completed"
echo "============================================================"
echo "Total time: ${TOTAL_TIME}s ($(echo "scale=1; $TOTAL_TIME/60" | bc) minutes)"
echo "Results saved to: $RESULTS_DIR"
echo "End time: $(date)"
