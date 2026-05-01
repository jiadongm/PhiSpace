#!/usr/bin/env python3
"""
Visualize benchmarking results comparing PhiSpace with other methods.
Creates publication-quality figures for the benchmarking comparison.
"""

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set plotting style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")


def load_integrated_results(results_dir):
    """Load all integrated metric results."""
    metrics = ['PCC', 'SSIM', 'RMSE', 'JS']
    results = {}

    for metric in metrics:
        metric_file = results_dir / f"{metric}_integrated.csv"
        if not metric_file.exists():
            # Try alternate name (e.g., JSD vs JS)
            metric_file = results_dir / f"{metric}D_integrated.csv"
        if metric_file.exists():
            df = pd.read_csv(metric_file, index_col=0)
            results[metric] = df
            print(f"Loaded {metric}: {df.shape}")
        else:
            print(f"WARNING: {metric} results not found")

    return results


def plot_boxplot_comparison(results, output_dir, figsize=(12, 10)):
    """
    Create boxplot comparison of all methods across all metrics.
    """
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    axes = axes.flatten()

    metric_order = ['PCC', 'SSIM', 'RMSE', 'JS']
    metric_titles = {
        'PCC': 'Pearson Correlation Coefficient (PCC)',
        'SSIM': 'Structural Similarity Index (SSIM)',
        'RMSE': 'Root Mean Square Error (RMSE)',
        'JS': 'Jensen-Shannon Divergence (JSD)'
    }

    # Higher is better for PCC, SSIM; Lower is better for RMSE, JSD
    higher_better = {'PCC': True, 'SSIM': True, 'RMSE': False, 'JS': False}

    for idx, metric in enumerate(metric_order):
        ax = axes[idx]

        if metric not in results:
            ax.set_visible(False)
            continue

        df = results[metric]
        method_cols = [c for c in df.columns if c != 'dataset']

        # Melt data for seaborn
        melted = df[method_cols].melt(var_name='Method', value_name='Value')
        melted = melted.dropna()

        # Order methods by median performance
        method_order = melted.groupby('Method')['Value'].median()
        if higher_better[metric]:
            method_order = method_order.sort_values(ascending=False)
        else:
            method_order = method_order.sort_values(ascending=True)

        # Create boxplot
        colors = ['#e74c3c' if m == 'PhiSpace' else '#3498db' for m in method_order.index]
        palette = dict(zip(method_order.index, colors))

        sns.boxplot(
            data=melted,
            x='Method',
            y='Value',
            order=method_order.index,
            palette=palette,
            ax=ax,
            width=0.6
        )

        ax.set_title(metric_titles.get(metric, metric), fontsize=12, fontweight='bold')
        ax.set_xlabel('')
        ax.set_ylabel(metric)
        ax.tick_params(axis='x', rotation=45)

        # Add statistical annotation for PhiSpace
        if 'PhiSpace' in method_cols:
            phispace_median = df['PhiSpace'].median()
            ax.axhline(y=phispace_median, color='red', linestyle='--', alpha=0.5, linewidth=1)

    plt.tight_layout()

    output_file = output_dir / 'boxplot_comparison.pdf'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'boxplot_comparison.png', dpi=300, bbox_inches='tight')
    print(f"Saved boxplot comparison to {output_file}")
    plt.close()


def plot_heatmap_by_dataset(results, output_dir, figsize=(14, 8)):
    """
    Create heatmap showing method performance by dataset.
    """
    for metric, df in results.items():
        method_cols = [c for c in df.columns if c != 'dataset']

        # Calculate mean per dataset
        dataset_means = df.groupby('dataset')[method_cols].mean()

        fig, ax = plt.subplots(figsize=figsize)

        # Create heatmap
        sns.heatmap(
            dataset_means.T,
            annot=False,
            cmap='RdYlBu_r' if metric in ['RMSE', 'JS'] else 'RdYlBu',
            ax=ax,
            linewidths=0.5
        )

        ax.set_title(f'{metric} by Dataset and Method', fontsize=14, fontweight='bold')
        ax.set_xlabel('Dataset')
        ax.set_ylabel('Method')

        plt.tight_layout()

        output_file = output_dir / f'heatmap_{metric.lower()}.pdf'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.savefig(output_dir / f'heatmap_{metric.lower()}.png', dpi=300, bbox_inches='tight')
        print(f"Saved {metric} heatmap to {output_file}")
        plt.close()


def plot_summary_barplot(results, output_dir, figsize=(10, 8)):
    """
    Create barplot of mean performance with error bars.
    """
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    axes = axes.flatten()

    metric_order = ['PCC', 'SSIM', 'RMSE', 'JS']
    higher_better = {'PCC': True, 'SSIM': True, 'RMSE': False, 'JS': False}

    for idx, metric in enumerate(metric_order):
        ax = axes[idx]

        if metric not in results:
            ax.set_visible(False)
            continue

        df = results[metric]
        method_cols = [c for c in df.columns if c != 'dataset']

        # Calculate mean and std
        means = df[method_cols].mean()
        stds = df[method_cols].std()

        # Sort by performance
        if higher_better[metric]:
            order = means.sort_values(ascending=False).index
        else:
            order = means.sort_values(ascending=True).index

        means = means[order]
        stds = stds[order]

        # Create barplot
        colors = ['#e74c3c' if m == 'PhiSpace' else '#3498db' for m in order]
        bars = ax.bar(range(len(means)), means.values, yerr=stds.values,
                      color=colors, capsize=3, alpha=0.8)

        ax.set_xticks(range(len(means)))
        ax.set_xticklabels(order, rotation=45, ha='right')
        ax.set_ylabel(metric)
        ax.set_title(metric, fontsize=12, fontweight='bold')

        # Add value labels on bars
        for i, (bar, mean) in enumerate(zip(bars, means.values)):
            ax.annotate(f'{mean:.3f}',
                       xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
                       xytext=(0, 3),
                       textcoords="offset points",
                       ha='center', va='bottom', fontsize=8)

    plt.tight_layout()

    output_file = output_dir / 'summary_barplot.pdf'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'summary_barplot.png', dpi=300, bbox_inches='tight')
    print(f"Saved summary barplot to {output_file}")
    plt.close()


def plot_ranking_table(results, output_dir):
    """
    Create ranking table showing relative performance.
    """
    rankings = {}
    higher_better = {'PCC': True, 'SSIM': True, 'RMSE': False, 'JS': False}

    for metric, df in results.items():
        method_cols = [c for c in df.columns if c != 'dataset']
        means = df[method_cols].mean()

        # Rank methods
        if higher_better[metric]:
            rank = means.rank(ascending=False)
        else:
            rank = means.rank(ascending=True)

        rankings[metric] = rank

    ranking_df = pd.DataFrame(rankings)
    ranking_df['Mean_Rank'] = ranking_df.mean(axis=1)
    ranking_df = ranking_df.sort_values('Mean_Rank')

    # Save ranking table
    output_file = output_dir / 'method_rankings.csv'
    ranking_df.to_csv(output_file)
    print(f"Saved ranking table to {output_file}")

    # Print ranking
    print("\nMethod Rankings (1 = best):")
    print(ranking_df.to_string())

    return ranking_df


def plot_violin_comparison(results, output_dir, figsize=(14, 5)):
    """
    Create violin plot comparison focusing on key metrics.
    """
    fig, axes = plt.subplots(1, 4, figsize=figsize)

    metric_order = ['PCC', 'SSIM', 'RMSE', 'JS']

    for idx, metric in enumerate(metric_order):
        ax = axes[idx]

        if metric not in results:
            ax.set_visible(False)
            continue

        df = results[metric]
        method_cols = [c for c in df.columns if c != 'dataset']

        # Melt data
        melted = df[method_cols].melt(var_name='Method', value_name='Value')
        melted = melted.dropna()

        # Create violin plot
        colors = ['#e74c3c' if m == 'PhiSpace' else '#3498db' for m in method_cols]
        palette = dict(zip(method_cols, colors))

        sns.violinplot(
            data=melted,
            x='Method',
            y='Value',
            palette=palette,
            ax=ax,
            inner='box',
            scale='width'
        )

        ax.set_title(metric, fontsize=12, fontweight='bold')
        ax.set_xlabel('')
        ax.tick_params(axis='x', rotation=45)

    plt.tight_layout()

    output_file = output_dir / 'violin_comparison.pdf'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'violin_comparison.png', dpi=300, bbox_inches='tight')
    print(f"Saved violin comparison to {output_file}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Visualize benchmarking results')
    parser.add_argument('--results_dir', type=str,
                        default='PhiSpace_Benchmarking/integrated_results',
                        help='Directory containing integrated results')
    parser.add_argument('--output_dir', type=str,
                        default='PhiSpace_Benchmarking/figures',
                        help='Output directory for figures')

    args = parser.parse_args()

    # Get base directory
    script_dir = Path(__file__).parent.absolute()
    base_dir = script_dir.parent

    # Set up paths
    results_dir = base_dir / args.results_dir
    output_dir = base_dir / args.output_dir

    output_dir.mkdir(parents=True, exist_ok=True)

    print("="*60)
    print("Visualizing Benchmarking Results")
    print("="*60)
    print(f"Results directory: {results_dir}")
    print(f"Output directory: {output_dir}")

    # Load results
    results = load_integrated_results(results_dir)

    if not results:
        print("ERROR: No results found to visualize")
        return

    # Create visualizations
    print("\n" + "-"*40)
    print("Creating boxplot comparison...")
    plot_boxplot_comparison(results, output_dir)

    print("\n" + "-"*40)
    print("Creating summary barplot...")
    plot_summary_barplot(results, output_dir)

    print("\n" + "-"*40)
    print("Creating violin comparison...")
    plot_violin_comparison(results, output_dir)

    print("\n" + "-"*40)
    print("Creating heatmaps by dataset...")
    plot_heatmap_by_dataset(results, output_dir)

    print("\n" + "-"*40)
    print("Creating ranking table...")
    plot_ranking_table(results, output_dir)

    print("\n" + "="*60)
    print("Visualization completed!")
    print(f"Figures saved to: {output_dir}")
    print("="*60)


if __name__ == '__main__':
    main()
