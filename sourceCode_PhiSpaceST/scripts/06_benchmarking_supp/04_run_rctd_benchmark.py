#!/usr/bin/env python3
"""
Orchestration script for RCTD benchmarking.
Runs RCTD on all 32 simulated datasets and computes metrics.
"""

import os
import sys
import argparse
import subprocess
import time
import pandas as pd
import numpy as np
from pathlib import Path


def run_command(cmd, description, timeout=3600):
    """Run a shell command with timeout and logging."""
    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"Command: {cmd}")
    print('='*60)

    start_time = time.time()
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            timeout=timeout
        )
        elapsed = time.time() - start_time

        if result.returncode != 0:
            print(f"ERROR: {result.stderr}")
            return False, elapsed

        print(result.stdout)
        print(f"Completed in {elapsed:.2f} seconds")
        return True, elapsed

    except subprocess.TimeoutExpired:
        print(f"TIMEOUT after {timeout} seconds")
        return False, timeout
    except Exception as e:
        print(f"ERROR: {str(e)}")
        return False, 0


def main():
    parser = argparse.ArgumentParser(description='Run full RCTD benchmarking pipeline')
    parser.add_argument('--output_dir', type=str,
                        default='PhiSpace_Benchmarking/results_rctd',
                        help='Output directory for results')
    parser.add_argument('--datasets', type=str, default='all',
                        help='Comma-separated dataset numbers or "all" (1-32)')
    parser.add_argument('--skip_rctd', action='store_true',
                        help='Skip RCTD execution step (only compute metrics)')

    args = parser.parse_args()

    # Get base directory
    script_dir = Path(__file__).parent.absolute()
    base_dir = script_dir.parent

    # Set up paths
    output_dir = base_dir / args.output_dir
    prepared_dir = script_dir / 'prepared_data'

    output_dir.mkdir(parents=True, exist_ok=True)

    # Determine datasets to process
    if args.datasets == 'all':
        datasets = list(range(1, 33))
    else:
        datasets = [int(x.strip()) for x in args.datasets.split(',')]

    print("="*60)
    print("RCTD Full Benchmarking Pipeline")
    print("="*60)
    print(f"Output directory: {output_dir}")
    print(f"Prepared data directory: {prepared_dir}")
    print(f"Datasets to process: {len(datasets)} datasets")
    print(f"Dataset numbers: {datasets}")

    # Track timing and results
    timing_results = []

    # Step 1: Run RCTD on all datasets
    if not args.skip_rctd:
        print("\n" + "#"*60)
        print("# STEP 1: Running RCTD")
        print("#"*60)

        for ds_num in datasets:
            ds_prepared_dir = prepared_dir / f"dataset{ds_num}"
            ds_output_dir = output_dir / f"dataset{ds_num}"

            if not ds_prepared_dir.exists():
                print(f"WARNING: Prepared data for dataset {ds_num} not found at {ds_prepared_dir}")
                continue

            ds_output_dir.mkdir(parents=True, exist_ok=True)

            cmd = f"Rscript {script_dir / '02_run_rctd.R'} " \
                  f"{ds_prepared_dir} {ds_output_dir}"

            success, elapsed = run_command(cmd, f"Running RCTD on dataset {ds_num}")

            timing_results.append({
                'dataset': ds_num,
                'method': 'RCTD',
                'runtime_seconds': elapsed,
                'success': success
            })

    # Step 2: Calculate metrics for all datasets
    print("\n" + "#"*60)
    print("# STEP 2: Calculating metrics")
    print("#"*60)

    all_metrics = {'pcc': [], 'ssim': [], 'rmse': [], 'jsd': []}

    for ds_num in datasets:
        ds_prepared_dir = prepared_dir / f"dataset{ds_num}"
        ds_output_dir = output_dir / f"dataset{ds_num}"
        ds_metrics_dir = output_dir / 'all_metrics' / f"dataset{ds_num}"

        gt_path = ds_prepared_dir / "ground_truth.csv"
        rctd_result = ds_output_dir / "RCTD_result.txt"

        if not gt_path.exists():
            print(f"WARNING: Ground truth for dataset {ds_num} not found")
            continue

        if not rctd_result.exists():
            print(f"WARNING: RCTD result for dataset {ds_num} not found")
            continue

        ds_metrics_dir.mkdir(parents=True, exist_ok=True)

        cmd = f"python {script_dir / '03_calculate_metrics.py'} " \
              f"--dataset_dir {ds_prepared_dir} " \
              f"--results_dir {ds_output_dir} " \
              f"--output_dir {ds_metrics_dir} " \
              f"--methods RCTD"

        success, _ = run_command(cmd, f"Calculating metrics for dataset {ds_num}")

        # Load and aggregate metrics
        if success:
            for metric in all_metrics.keys():
                metric_file = ds_metrics_dir / f"{metric.upper()}.csv"
                if metric_file.exists():
                    df = pd.read_csv(metric_file, index_col=0)
                    df['dataset'] = ds_num
                    all_metrics[metric].append(df)

    # Step 3: Combine all metrics
    print("\n" + "#"*60)
    print("# STEP 3: Combining results")
    print("#"*60)

    summary_rows = []

    for metric, dfs in all_metrics.items():
        if dfs:
            combined = pd.concat(dfs, axis=0)
            output_file = output_dir / f"RCTD_{metric.upper()}.csv"
            combined.to_csv(output_file)
            print(f"Saved combined {metric.upper()} metrics to {output_file}")

            # Aggregate per-dataset mean for RCTD
            for df in dfs:
                ds_num = df['dataset'].iloc[0]
                mean_val = df.drop(columns=['dataset']).mean(axis=0)
                if 'RCTD' in mean_val.index:
                    summary_rows.append({
                        'dataset': ds_num,
                        'metric': metric.upper(),
                        'mean': mean_val['RCTD']
                    })

    # Save all-metrics summary
    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        # Pivot to wide format: dataset x metric
        summary_wide = summary_df.pivot(index='dataset', columns='metric', values='mean')
        summary_wide.to_csv(output_dir / 'RCTD_all_metrics.csv')
        print(f"\nAll metrics summary saved to {output_dir / 'RCTD_all_metrics.csv'}")
        print(f"\nSummary across all datasets:")
        print(summary_wide.describe().loc[['mean', 'std', 'min', 'max']])

    # Save timing results
    if timing_results:
        timing_df = pd.DataFrame(timing_results)
        timing_file = output_dir / "rctd_timing.csv"
        timing_df.to_csv(timing_file, index=False)
        print(f"\nTiming results saved to {timing_file}")

        successful = timing_df[timing_df['success'] == True]
        if len(successful) > 0:
            print(f"\nTiming Summary:")
            print(f"  Mean runtime: {successful['runtime_seconds'].mean():.2f}s")
            print(f"  Median runtime: {successful['runtime_seconds'].median():.2f}s")
            print(f"  Total runtime: {successful['runtime_seconds'].sum():.2f}s")

    print("\n" + "="*60)
    print("RCTD Benchmarking Pipeline Completed!")
    print("="*60)


if __name__ == '__main__':
    main()
