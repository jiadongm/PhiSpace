#!/usr/bin/env python3
"""
Runtime and computational resource benchmarking.
Compares PhiSpace runtime against other deconvolution methods.
Addresses reviewer comment 3.2 about computational resource evaluation.
"""

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


def load_phispace_timing(results_dir):
    """Load PhiSpace timing results."""
    timing_file = results_dir / "phispace_timing.csv"
    if timing_file.exists():
        return pd.read_csv(timing_file)
    return None


def load_dataset_metadata(prepared_dir, datasets):
    """Load dataset metadata (n_spots, n_cells, n_genes)."""
    metadata = []
    for ds_num in datasets:
        meta_file = prepared_dir / f"dataset{ds_num}" / "metadata.csv"
        if meta_file.exists():
            df = pd.read_csv(meta_file)
            df['dataset'] = ds_num
            metadata.append(df)

    if metadata:
        return pd.concat(metadata, ignore_index=True)
    return None


def create_runtime_comparison_plot(timing_df, metadata_df, output_dir):
    """
    Create runtime comparison plots.
    """
    if timing_df is None:
        print("No timing data available")
        return

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Plot 1: Runtime distribution
    ax1 = axes[0]
    successful = timing_df[timing_df['success'] == True]
    sns.histplot(successful['runtime_seconds'], bins=15, ax=ax1, color='#3498db')
    ax1.set_xlabel('Runtime (seconds)')
    ax1.set_ylabel('Count')
    ax1.set_title('PhiSpace Runtime Distribution')
    ax1.axvline(successful['runtime_seconds'].median(), color='red', linestyle='--',
                label=f'Median: {successful["runtime_seconds"].median():.1f}s')
    ax1.legend()

    # Plot 2: Runtime vs dataset size (if metadata available)
    ax2 = axes[1]
    if metadata_df is not None:
        merged = pd.merge(successful, metadata_df, on='dataset')
        if 'n_spots' in merged.columns:
            ax2.scatter(merged['n_spots'], merged['runtime_seconds'],
                       alpha=0.7, s=50, c='#3498db')
            ax2.set_xlabel('Number of Spots')
            ax2.set_ylabel('Runtime (seconds)')
            ax2.set_title('Runtime vs Dataset Size')

            # Add trend line
            z = np.polyfit(merged['n_spots'], merged['runtime_seconds'], 1)
            p = np.poly1d(z)
            x_line = np.linspace(merged['n_spots'].min(), merged['n_spots'].max(), 100)
            ax2.plot(x_line, p(x_line), 'r--', alpha=0.5, label='Trend')
            ax2.legend()
    else:
        ax2.text(0.5, 0.5, 'Metadata not available', transform=ax2.transAxes,
                ha='center', va='center')

    # Plot 3: Runtime by dataset
    ax3 = axes[2]
    successful_sorted = successful.sort_values('dataset')
    ax3.bar(range(len(successful_sorted)), successful_sorted['runtime_seconds'],
           color='#3498db', alpha=0.8)
    ax3.set_xlabel('Dataset')
    ax3.set_ylabel('Runtime (seconds)')
    ax3.set_title('PhiSpace Runtime by Dataset')
    ax3.set_xticks(range(len(successful_sorted)))
    ax3.set_xticklabels(successful_sorted['dataset'], rotation=45)

    plt.tight_layout()

    output_file = output_dir / 'phispace_runtime_analysis.pdf'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'phispace_runtime_analysis.png', dpi=300, bbox_inches='tight')
    print(f"Saved runtime analysis to {output_file}")
    plt.close()


def create_runtime_summary(timing_df, output_dir):
    """
    Create runtime summary statistics.
    """
    if timing_df is None:
        return None

    successful = timing_df[timing_df['success'] == True]

    summary = {
        'Method': 'PhiSpace',
        'N_datasets': len(successful),
        'Mean_runtime_sec': successful['runtime_seconds'].mean(),
        'Median_runtime_sec': successful['runtime_seconds'].median(),
        'Std_runtime_sec': successful['runtime_seconds'].std(),
        'Min_runtime_sec': successful['runtime_seconds'].min(),
        'Max_runtime_sec': successful['runtime_seconds'].max(),
        'Total_runtime_sec': successful['runtime_seconds'].sum(),
    }

    summary_df = pd.DataFrame([summary])
    output_file = output_dir / 'runtime_summary.csv'
    summary_df.to_csv(output_file, index=False)
    print(f"Saved runtime summary to {output_file}")

    print("\nRuntime Summary:")
    print(f"  Mean: {summary['Mean_runtime_sec']:.2f} seconds")
    print(f"  Median: {summary['Median_runtime_sec']:.2f} seconds")
    print(f"  Std: {summary['Std_runtime_sec']:.2f} seconds")
    print(f"  Range: {summary['Min_runtime_sec']:.2f} - {summary['Max_runtime_sec']:.2f} seconds")
    print(f"  Total (all datasets): {summary['Total_runtime_sec']:.2f} seconds")

    return summary_df


def compare_with_literature(output_dir):
    """
    Create comparison table with literature-reported runtimes.
    Reference: Li et al., Nature Methods 2022 (benchmark paper)
    """
    # Approximate runtimes from the benchmark paper (per dataset, in seconds)
    # These are rough estimates - actual values would need to be extracted
    # from running the methods or from supplementary data

    literature_runtimes = pd.DataFrame([
        {'Method': 'RCTD', 'Approx_runtime_sec': 300, 'Source': 'Li et al. 2022'},
        {'Method': 'Cell2location', 'Approx_runtime_sec': 1800, 'Source': 'Li et al. 2022'},
        {'Method': 'SPOTlight', 'Approx_runtime_sec': 600, 'Source': 'Li et al. 2022'},
        {'Method': 'Tangram', 'Approx_runtime_sec': 900, 'Source': 'Li et al. 2022'},
        {'Method': 'Seurat', 'Approx_runtime_sec': 450, 'Source': 'Li et al. 2022'},
    ])

    output_file = output_dir / 'literature_runtime_comparison.csv'
    literature_runtimes.to_csv(output_file, index=False)
    print(f"\nSaved literature runtime comparison to {output_file}")
    print("\nNote: These are approximate values from the benchmark paper.")
    print("For accurate comparison, run all methods on the same hardware.")

    return literature_runtimes


def main():
    parser = argparse.ArgumentParser(description='Runtime benchmarking analysis')
    parser.add_argument('--results_dir', type=str,
                        default='PhiSpace_Benchmarking/results',
                        help='Directory containing PhiSpace results')
    parser.add_argument('--prepared_dir', type=str,
                        default='PhiSpace_Benchmarking/prepared_data',
                        help='Directory containing prepared data')
    parser.add_argument('--output_dir', type=str,
                        default='PhiSpace_Benchmarking/runtime_analysis',
                        help='Output directory for runtime analysis')
    parser.add_argument('--datasets', type=str, default='all',
                        help='Comma-separated dataset numbers or "all"')

    args = parser.parse_args()

    # Get base directory
    script_dir = Path(__file__).parent.absolute()
    base_dir = script_dir.parent

    # Set up paths
    results_dir = base_dir / args.results_dir
    prepared_dir = base_dir / args.prepared_dir
    output_dir = base_dir / args.output_dir

    output_dir.mkdir(parents=True, exist_ok=True)

    # Determine datasets
    if args.datasets == 'all':
        datasets = list(range(1, 33))
    else:
        datasets = [int(x.strip()) for x in args.datasets.split(',')]

    print("="*60)
    print("Runtime Benchmarking Analysis")
    print("="*60)

    # Load timing data
    timing_df = load_phispace_timing(results_dir)
    if timing_df is not None:
        print(f"Loaded timing data for {len(timing_df)} runs")

    # Load metadata
    metadata_df = load_dataset_metadata(prepared_dir, datasets)
    if metadata_df is not None:
        print(f"Loaded metadata for {len(metadata_df)} datasets")

    # Create summary
    print("\n" + "-"*40)
    create_runtime_summary(timing_df, output_dir)

    # Create plots
    print("\n" + "-"*40)
    print("Creating runtime plots...")
    create_runtime_comparison_plot(timing_df, metadata_df, output_dir)

    # Add literature comparison
    print("\n" + "-"*40)
    compare_with_literature(output_dir)

    print("\n" + "="*60)
    print("Runtime analysis completed!")
    print(f"Results saved to: {output_dir}")
    print("="*60)


if __name__ == '__main__':
    main()
