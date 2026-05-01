#!/usr/bin/env python3
"""
Main orchestration script for PhiSpace benchmarking.
Runs the full pipeline across all 32 simulated datasets.
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


def aggregate_metrics(output_dir, datasets, method_name, metrics_subdir):
    """Aggregate per-dataset metrics into summary CSV files."""
    all_metrics = {'pcc': [], 'ssim': [], 'rmse': [], 'jsd': []}

    for ds_num in datasets:
        ds_metrics_dir = output_dir / f"dataset{ds_num}" / metrics_subdir
        for metric in all_metrics.keys():
            metric_file = ds_metrics_dir / f"{metric.upper()}.csv"
            if metric_file.exists():
                df = pd.read_csv(metric_file, index_col=0)
                df['dataset'] = ds_num
                all_metrics[metric].append(df)

    summary_rows = []
    for metric, dfs in all_metrics.items():
        if dfs:
            combined = pd.concat(dfs, axis=0)
            output_file = output_dir / f"{method_name}_{metric.upper()}.csv"
            combined.to_csv(output_file)
            print(f"  Saved combined {metric.upper()} to {output_file}")

            for df in dfs:
                ds_num = df['dataset'].iloc[0]
                mean_val = df.drop(columns=['dataset']).mean(axis=0)
                if method_name in mean_val.index:
                    summary_rows.append({
                        'dataset': ds_num,
                        'metric': metric.upper(),
                        'mean': mean_val[method_name]
                    })

    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_wide = summary_df.pivot(index='dataset', columns='metric', values='mean')
        summary_file = output_dir / f"{method_name}_all_metrics.csv"
        summary_wide.to_csv(summary_file)
        print(f"  Summary saved to {summary_file}")
        print(f"\n  {method_name} summary across datasets:")
        print(summary_wide.describe().loc[['mean', 'std', 'min', 'max']])


def main():
    parser = argparse.ArgumentParser(description='Run full PhiSpace benchmarking pipeline')
    parser.add_argument('--data_dir', type=str,
                        default='ExampleData/SimualtedSpatalData',
                        help='Path to SimualtedSpatalData directory')
    parser.add_argument('--output_dir', type=str,
                        default='PhiSpace_Benchmarking/results',
                        help='Output directory for results')
    parser.add_argument('--datasets', type=str, default='all',
                        help='Comma-separated dataset numbers or "all" (1-32)')
    parser.add_argument('--ref_norm', type=str, default='sctransform',
                        help='Reference normalization method (default: sctransform)')
    parser.add_argument('--query_norm', type=str, default='log1p',
                        help='Query normalization method (default: log1p)')
    parser.add_argument('--skip_prepare', action='store_true',
                        help='Skip data preparation step')
    parser.add_argument('--skip_phispace', action='store_true',
                        help='Skip PhiSpace execution step')
    parser.add_argument('--skip_deconv', action='store_true',
                        help='Skip deconvolution step (only compute metrics)')
    parser.add_argument('--n_hvg', type=int, default=0,
                        help='Number of HVGs to select (0 = use all genes)')

    args = parser.parse_args()

    # Get base directory
    script_dir = Path(__file__).parent.absolute()
    base_dir = script_dir.parent

    # Set up paths
    data_dir = base_dir / args.data_dir
    output_dir = base_dir / args.output_dir
    prepared_dir = script_dir / 'prepared_data'

    output_dir.mkdir(parents=True, exist_ok=True)

    # Determine datasets to process
    if args.datasets == 'all':
        datasets = list(range(1, 33))
    else:
        datasets = [int(x.strip()) for x in args.datasets.split(',')]

    print("="*60)
    print("PhiSpace Full Benchmarking Pipeline")
    print("="*60)
    print(f"Data directory: {data_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Datasets to process: {len(datasets)} datasets")
    print(f"Dataset numbers: {datasets}")
    print(f"Reference normalization: {args.ref_norm}")
    print(f"Query normalization: {args.query_norm}")

    # Track timing and results
    timing_results = []

    # Step 1: Prepare data for all datasets
    if not args.skip_prepare:
        print("\n" + "#"*60)
        print("# STEP 1: Preparing data")
        print("#"*60)

        for ds_num in datasets:
            ds_path = data_dir / f"dataset{ds_num}"

            if not ds_path.exists():
                print(f"WARNING: Dataset {ds_num} not found at {ds_path}")
                continue

            cmd = f"python {script_dir / '01_prepare_data.py'} " \
                  f"--data_dir {data_dir} " \
                  f"--output_dir {prepared_dir} " \
                  f"--datasets {ds_num}"

            success, elapsed = run_command(cmd, f"Preparing dataset {ds_num}")
            if not success:
                print(f"WARNING: Failed to prepare dataset {ds_num}")

    # Step 2: Run PhiSpace on all datasets
    if not args.skip_phispace and not args.skip_deconv:
        print("\n" + "#"*60)
        print("# STEP 2: Running PhiSpace")
        print("#"*60)

        for ds_num in datasets:
            ds_prepared_dir = prepared_dir / f"dataset{ds_num}"
            ds_output_dir = output_dir / f"dataset{ds_num}"

            if not ds_prepared_dir.exists():
                print(f"WARNING: Prepared data for dataset {ds_num} not found")
                continue

            ds_output_dir.mkdir(parents=True, exist_ok=True)

            cmd = f"Rscript {script_dir / '02_run_phispace.R'} " \
                  f"{ds_prepared_dir} {ds_output_dir} " \
                  f"{args.ref_norm} {args.query_norm} {args.n_hvg}"

            success, elapsed = run_command(cmd, f"Running PhiSpace on dataset {ds_num}")

            timing_results.append({
                'dataset': ds_num,
                'method': 'PhiSpace',
                'ref_norm': args.ref_norm,
                'query_norm': args.query_norm,
                'runtime_seconds': elapsed,
                'success': success
            })

    # Step 3: Calculate metrics for both normalised and unnormalised results
    print("\n" + "#"*60)
    print("# STEP 3: Calculating metrics")
    print("#"*60)

    for variant, method_name, metrics_subdir in [
        ("PhiSpace_result.txt", "PhiSpace", "metrics"),
        ("PhiSpace_unnorm_result.txt", "PhiSpace_unnorm", "metrics_unnorm"),
    ]:
        print(f"\n--- Computing metrics for {method_name} ---")

        for ds_num in datasets:
            ds_prepared_dir = prepared_dir / f"dataset{ds_num}"
            ds_output_dir = output_dir / f"dataset{ds_num}"
            ds_metrics_dir = ds_output_dir / metrics_subdir

            gt_path = ds_prepared_dir / "ground_truth.csv"
            result_file = ds_output_dir / variant

            if not gt_path.exists():
                print(f"WARNING: Ground truth for dataset {ds_num} not found")
                continue

            if not result_file.exists():
                print(f"WARNING: {variant} for dataset {ds_num} not found")
                continue

            ds_metrics_dir.mkdir(parents=True, exist_ok=True)

            cmd = f"python {script_dir / '03_calculate_metrics.py'} " \
                  f"--dataset_dir {ds_prepared_dir} " \
                  f"--results_dir {ds_output_dir} " \
                  f"--output_dir {ds_metrics_dir} " \
                  f"--methods {method_name}"

            success, _ = run_command(cmd, f"Calculating {method_name} metrics for dataset {ds_num}")

    # Step 4: Combine all metrics
    print("\n" + "#"*60)
    print("# STEP 4: Combining results")
    print("#"*60)

    for method_name, metrics_subdir in [
        ("PhiSpace", "metrics"),
        ("PhiSpace_unnorm", "metrics_unnorm"),
    ]:
        print(f"\n--- Aggregating {method_name} ---")
        aggregate_metrics(output_dir, datasets, method_name, metrics_subdir)

    # Save timing results
    if timing_results:
        timing_df = pd.DataFrame(timing_results)
        timing_file = output_dir / "phispace_timing.csv"
        timing_df.to_csv(timing_file, index=False)
        print(f"\nTiming results saved to {timing_file}")

        successful = timing_df[timing_df['success'] == True]
        if len(successful) > 0:
            print(f"\nTiming Summary:")
            print(f"  Mean runtime: {successful['runtime_seconds'].mean():.2f}s")
            print(f"  Median runtime: {successful['runtime_seconds'].median():.2f}s")
            print(f"  Total runtime: {successful['runtime_seconds'].sum():.2f}s")

    print("\n" + "="*60)
    print("PhiSpace Benchmarking Pipeline Completed!")
    print("="*60)


if __name__ == '__main__':
    main()
