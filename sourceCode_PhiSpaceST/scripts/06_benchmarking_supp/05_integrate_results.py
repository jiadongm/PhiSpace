#!/usr/bin/env python3
"""
Integrate PhiSpace results with existing benchmark results.
Adds PhiSpace column to the existing PCC, SSIM, RMSE, JS CSV files.
"""

import os
import argparse
import pandas as pd
import numpy as np
from pathlib import Path


def load_existing_results(figure_data_dir, metric):
    """Load existing benchmark results for a metric."""
    metric_file = figure_data_dir / f"{metric}.csv"
    if not metric_file.exists():
        print(f"WARNING: Existing results not found: {metric_file}")
        return None

    df = pd.read_csv(metric_file, index_col=0)
    print(f"Loaded existing {metric} results: {df.shape}")
    return df


def load_phispace_results(results_dir, metric, datasets):
    """Load PhiSpace results for all datasets."""
    all_results = []

    for ds_num in datasets:
        metric_file = results_dir / f"dataset{ds_num}" / "metrics" / f"{metric}.csv"

        if not metric_file.exists():
            print(f"  Dataset {ds_num}: PhiSpace results not found")
            continue

        df = pd.read_csv(metric_file, index_col=0)
        df['dataset'] = ds_num

        # Keep only PhiSpace column
        if 'PhiSpace' in df.columns:
            result = df[['PhiSpace']].copy()
            result['dataset'] = ds_num
            all_results.append(result)
            print(f"  Dataset {ds_num}: Loaded {len(result)} cell types")
        else:
            print(f"  Dataset {ds_num}: PhiSpace column not found")

    if not all_results:
        return None

    combined = pd.concat(all_results, axis=0)
    return combined


def integrate_results(existing_df, phispace_df, methods_to_keep):
    """
    Integrate PhiSpace results with existing benchmark results.

    Parameters:
    -----------
    existing_df : pd.DataFrame
        Existing benchmark results (cell types x methods, with dataset column)
    phispace_df : pd.DataFrame
        PhiSpace results (cell types x PhiSpace, with dataset column)
    methods_to_keep : list
        List of methods to include in final comparison

    Returns:
    --------
    pd.DataFrame : Combined results
    """
    if existing_df is None or phispace_df is None:
        return None

    # Filter existing results to keep only specified methods
    cols_to_keep = [c for c in methods_to_keep if c in existing_df.columns]
    cols_to_keep.append('dataset')

    filtered_existing = existing_df[cols_to_keep].copy()

    # Create a combined index from cell type and dataset
    # The existing data has cell type as index and dataset as column
    filtered_existing = filtered_existing.reset_index()
    filtered_existing.columns = ['celltype'] + cols_to_keep

    phispace_df = phispace_df.reset_index()
    phispace_df.columns = ['celltype', 'PhiSpace', 'dataset']

    # Merge on celltype and dataset
    merged = pd.merge(
        filtered_existing,
        phispace_df[['celltype', 'dataset', 'PhiSpace']],
        on=['celltype', 'dataset'],
        how='left'
    )

    # Set celltype back as index
    merged = merged.set_index('celltype')

    # Reorder columns to put PhiSpace first for visibility
    method_cols = ['PhiSpace'] + [c for c in methods_to_keep if c in merged.columns]
    merged = merged[method_cols + ['dataset']]

    return merged


def create_summary_table(integrated_results, output_dir):
    """
    Create summary statistics table across all datasets.

    Parameters:
    -----------
    integrated_results : dict
        Dictionary of metric -> DataFrame
    output_dir : Path
        Output directory for summary
    """
    summary_data = []

    for metric, df in integrated_results.items():
        if df is None:
            continue

        method_cols = [c for c in df.columns if c != 'dataset']

        for method in method_cols:
            values = df[method].dropna()

            summary_data.append({
                'Metric': metric,
                'Method': method,
                'Mean': values.mean(),
                'Std': values.std(),
                'Median': values.median(),
                'Q1': values.quantile(0.25),
                'Q3': values.quantile(0.75),
                'N_celltypes': len(values)
            })

    summary_df = pd.DataFrame(summary_data)

    # Create pivot table for easier comparison
    pivot_mean = summary_df.pivot(index='Method', columns='Metric', values='Mean')
    pivot_std = summary_df.pivot(index='Method', columns='Metric', values='Std')

    # Save summaries
    summary_df.to_csv(output_dir / 'summary_detailed.csv', index=False)
    pivot_mean.to_csv(output_dir / 'summary_mean.csv')
    pivot_std.to_csv(output_dir / 'summary_std.csv')

    print("\nSummary (Mean values):")
    print(pivot_mean.to_string())

    return summary_df


def main():
    parser = argparse.ArgumentParser(description='Integrate PhiSpace results with existing benchmarks')
    parser.add_argument('--figure_data_dir', type=str,
                        default='FigureData/Figure4/32SimulationData',
                        help='Directory containing existing benchmark results')
    parser.add_argument('--phispace_results_dir', type=str,
                        default='PhiSpace_Benchmarking/results',
                        help='Directory containing PhiSpace results')
    parser.add_argument('--output_dir', type=str,
                        default='PhiSpace_Benchmarking/integrated_results',
                        help='Output directory for integrated results')
    parser.add_argument('--datasets', type=str, default='all',
                        help='Comma-separated dataset numbers or "all"')
    parser.add_argument('--methods', type=str,
                        default='RCTD,Cell2location,SPOTlight,Tangram,Seurat',
                        help='Comma-separated list of methods to compare against')

    args = parser.parse_args()

    # Get base directory
    script_dir = Path(__file__).parent.absolute()
    base_dir = script_dir.parent

    # Set up paths
    figure_data_dir = base_dir / args.figure_data_dir
    phispace_results_dir = base_dir / args.phispace_results_dir
    output_dir = base_dir / args.output_dir

    output_dir.mkdir(parents=True, exist_ok=True)

    # Parse methods
    methods_to_keep = [m.strip() for m in args.methods.split(',')]

    # Determine datasets
    if args.datasets == 'all':
        datasets = list(range(1, 33))
    else:
        datasets = [int(x.strip()) for x in args.datasets.split(',')]

    print("="*60)
    print("Integrating PhiSpace Results with Existing Benchmarks")
    print("="*60)
    print(f"Figure data directory: {figure_data_dir}")
    print(f"PhiSpace results directory: {phispace_results_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Methods to compare: {methods_to_keep}")
    print(f"Datasets: {len(datasets)}")

    metrics = ['PCC', 'SSIM', 'RMSE', 'JS']
    integrated_results = {}

    for metric in metrics:
        print(f"\n{'='*40}")
        print(f"Processing {metric}...")
        print('='*40)

        # Load existing results
        existing_df = load_existing_results(figure_data_dir, metric)

        # Load PhiSpace results
        print(f"Loading PhiSpace {metric} results...")
        phispace_df = load_phispace_results(phispace_results_dir, metric, datasets)

        if phispace_df is not None:
            print(f"PhiSpace results shape: {phispace_df.shape}")
        else:
            print("No PhiSpace results found")

        # Integrate results
        integrated = integrate_results(existing_df, phispace_df, methods_to_keep)

        if integrated is not None:
            # Save integrated results
            output_file = output_dir / f"{metric}_integrated.csv"
            integrated.to_csv(output_file)
            print(f"Saved integrated {metric} results to {output_file}")
            print(f"Shape: {integrated.shape}")

            integrated_results[metric] = integrated
        else:
            print(f"WARNING: Could not integrate {metric} results")

    # Create summary tables
    print(f"\n{'='*60}")
    print("Creating summary tables...")
    print('='*60)

    create_summary_table(integrated_results, output_dir)

    print(f"\n{'='*60}")
    print("Integration completed!")
    print(f"Results saved to: {output_dir}")
    print('='*60)


if __name__ == '__main__':
    main()
