#!/usr/bin/env python3
"""
Calculate benchmarking metrics for deconvolution results.
Computes PCC, SSIM, RMSE, and JS divergence against ground truth.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from scipy.spatial.distance import jensenshannon
from sklearn.metrics import mean_squared_error
import warnings
warnings.filterwarnings('ignore')


def ssim(im1, im2, M=1):
    """
    Compute Structural Similarity Index (SSIM).

    Parameters:
    -----------
    im1, im2 : array-like
        Two vectors to compare
    M : float
        Dynamic range (default 1 for normalized data)

    Returns:
    --------
    float : SSIM value
    """
    im1 = np.asarray(im1, dtype=float)
    im2 = np.asarray(im2, dtype=float)

    # Normalize
    if im1.max() > 0:
        im1 = im1 / im1.max()
    if im2.max() > 0:
        im2 = im2 / im2.max()

    mu1 = im1.mean()
    mu2 = im2.mean()
    sigma1 = np.sqrt(((im1 - mu1) ** 2).mean())
    sigma2 = np.sqrt(((im2 - mu2) ** 2).mean())
    sigma12 = ((im1 - mu1) * (im2 - mu2)).mean()

    k1, k2, L = 0.01, 0.03, M
    C1 = (k1 * L) ** 2
    C2 = (k2 * L) ** 2
    C3 = C2 / 2

    l12 = (2 * mu1 * mu2 + C1) / (mu1 ** 2 + mu2 ** 2 + C1)
    c12 = (2 * sigma1 * sigma2 + C2) / (sigma1 ** 2 + sigma2 ** 2 + C2)
    s12 = (sigma12 + C3) / (sigma1 * sigma2 + C3)

    return l12 * c12 * s12


def rmse(x1, x2):
    """Compute Root Mean Square Error."""
    return mean_squared_error(x1, x2, squared=False)


def str_remove_special(string):
    """Remove special characters from string (for consistent column naming)."""
    return ''.join(c if c.isalnum() else '_' for c in str(string))


def load_results(results_path, method_name):
    """
    Load deconvolution results from a CSV file.

    Parameters:
    -----------
    results_path : str
        Path to results file
    method_name : str
        Name of the method for parsing specific formats

    Returns:
    --------
    pd.DataFrame : Results with spots as rows and cell types as columns
    """
    if not os.path.exists(results_path):
        print(f"  Warning: Results file not found: {results_path}")
        return None

    results = pd.read_csv(results_path, index_col=0)

    # Handle method-specific formatting
    if method_name == 'Cell2location':
        # Cell2location has 'q05cell_abundance_w_sf_' prefix
        new_cols = []
        for c in results.columns:
            if 'q05cell_abundance_w_sf_' in c:
                new_cols.append(str_remove_special(c.split('q05cell_abundance_w_sf_')[1]))
            else:
                new_cols.append(str_remove_special(c))
        results.columns = new_cols
        # Normalize to proportions
        results = results.div(results.sum(axis=1), axis=0).fillna(0)

    elif method_name == 'Seurat':
        # Seurat has 'prediction.score.' prefix and extra columns
        if 'predicted.id' in results.columns:
            results = results.drop(columns=['predicted.id'], errors='ignore')
        if 'prediction.score.max' in results.columns:
            results = results.drop(columns=['prediction.score.max'], errors='ignore')
        new_cols = []
        for c in results.columns:
            if 'prediction.score.' in c:
                new_cols.append(str_remove_special(c.split('prediction.score.')[1]))
            else:
                new_cols.append(str_remove_special(c))
        results.columns = new_cols

    elif method_name in ['RCTD', 'Tangram', 'SPOTlight', 'PhiSpace', 'PhiSpace_unnorm', 'PhiSpace_minmax', 'PhiSpace_rawscore', 'TACCO']:
        # Standard format, just clean column names
        results.columns = [str_remove_special(c) for c in results.columns]

    # Ensure all values are numeric and clip to [0, 1]
    results = results.apply(pd.to_numeric, errors='coerce').fillna(0)
    results = results.clip(lower=0, upper=1)

    # Normalize rows to sum to 1 if they don't already
    row_sums = results.sum(axis=1)
    if not np.allclose(row_sums, 1.0, atol=0.1):
        results = results.div(row_sums, axis=0).fillna(0)

    return results


def compare_results(ground_truth, predicted, metric='pcc', axis=1):
    """
    Compare predicted results against ground truth.

    Parameters:
    -----------
    ground_truth : pd.DataFrame
        Ground truth proportions (spots x cell types)
    predicted : pd.DataFrame
        Predicted proportions (spots x cell types)
    metric : str
        Metric to compute ('pcc', 'ssim', 'rmse', 'jsd')
    axis : int
        0 for by-spot comparison, 1 for by-celltype comparison

    Returns:
    --------
    list : Metric values
    """
    # Define metric functions
    metric_funcs = {
        'pcc': lambda x, y: pearsonr(x, y)[0] if len(set(x)) > 1 and len(set(y)) > 1 else np.nan,
        'ssim': ssim,
        'rmse': rmse,
        'jsd': jensenshannon
    }

    func = metric_funcs.get(metric)
    if func is None:
        raise ValueError(f"Unknown metric: {metric}")

    results = []

    if axis == 1:  # By cell type
        for col in ground_truth.columns:
            if col in predicted.columns:
                gt_vals = ground_truth[col].values
                pred_vals = np.clip(predicted[col].values, 0, 1)
                pred_vals = np.nan_to_num(pred_vals)
                try:
                    results.append(func(gt_vals, pred_vals))
                except:
                    results.append(np.nan)
            else:
                results.append(np.nan)
    else:  # By spot (axis=0)
        for idx in ground_truth.index:
            if idx in predicted.index:
                gt_vals = ground_truth.loc[idx].values
                pred_vals = np.clip(predicted.loc[idx].values, 0, 1)
                pred_vals = np.nan_to_num(pred_vals)
                try:
                    results.append(func(gt_vals, pred_vals))
                except:
                    results.append(np.nan)
            else:
                results.append(np.nan)

    return results


def calculate_all_metrics(ground_truth, results_dict, output_path):
    """
    Calculate all metrics for all methods.

    Parameters:
    -----------
    ground_truth : pd.DataFrame
        Ground truth proportions
    results_dict : dict
        Dictionary of method_name -> results DataFrame
    output_path : str
        Path to save metric results
    """
    metrics = ['pcc', 'ssim', 'rmse', 'jsd']

    # Clean ground truth column names
    ground_truth.columns = [str_remove_special(c) for c in ground_truth.columns]

    # Normalize ground truth to proportions
    ground_truth = ground_truth.div(ground_truth.sum(axis=1), axis=0).fillna(0)

    print(f"  Ground truth: {ground_truth.shape[0]} spots x {ground_truth.shape[1]} cell types")

    all_metrics = {}

    for metric in metrics:
        print(f"\n  Computing {metric.upper()}...")
        metric_results = {}

        for method_name, predicted in results_dict.items():
            if predicted is None:
                continue

            # Align columns
            common_cols = list(set(ground_truth.columns) & set(predicted.columns))
            if len(common_cols) == 0:
                print(f"    Warning: No common cell types between ground truth and {method_name}")
                continue

            gt_aligned = ground_truth[common_cols]
            pred_aligned = predicted[common_cols]

            # Ensure same number of rows
            if len(gt_aligned) != len(pred_aligned):
                min_len = min(len(gt_aligned), len(pred_aligned))
                gt_aligned = gt_aligned.iloc[:min_len]
                pred_aligned = pred_aligned.iloc[:min_len]

            # Compute metric by cell type (axis=1)
            values = compare_results(gt_aligned, pred_aligned, metric=metric, axis=1)
            metric_results[method_name] = values

            mean_val = np.nanmean(values)
            print(f"    {method_name}: {mean_val:.4f} (mean)")

        # Create DataFrame
        if metric_results:
            df = pd.DataFrame(metric_results, index=common_cols)
            all_metrics[metric] = df

    return all_metrics


def main():
    parser = argparse.ArgumentParser(description='Calculate benchmarking metrics')
    parser.add_argument('--dataset_dir', type=str, required=True,
                        help='Directory containing prepared dataset with ground_truth.csv')
    parser.add_argument('--results_dir', type=str, required=True,
                        help='Directory containing deconvolution results')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Directory for metric results')
    parser.add_argument('--methods', type=str,
                        default='RCTD,Cell2location,SPOTlight,Tangram,Seurat,PhiSpace',
                        help='Comma-separated list of methods to compare')

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("="*60)
    print("Calculating Benchmarking Metrics")
    print("="*60)

    # Load ground truth
    gt_path = os.path.join(args.dataset_dir, "ground_truth.csv")
    if not os.path.exists(gt_path):
        print(f"Error: Ground truth not found at {gt_path}")
        sys.exit(1)

    ground_truth = pd.read_csv(gt_path, index_col=0)
    print(f"Ground truth loaded: {ground_truth.shape}")

    # Load results for each method
    methods = [m.strip() for m in args.methods.split(',')]
    results_dict = {}

    for method in methods:
        result_file = os.path.join(args.results_dir, f"{method}_result.txt")
        print(f"\nLoading {method} results from {result_file}")
        results = load_results(result_file, method)
        if results is not None:
            print(f"  Loaded: {results.shape}")
            results_dict[method] = results
        else:
            print(f"  Not found or error loading")

    if not results_dict:
        print("Error: No results loaded")
        sys.exit(1)

    # Calculate metrics
    print("\n" + "="*60)
    print("Computing metrics...")
    print("="*60)

    all_metrics = calculate_all_metrics(ground_truth, results_dict, args.output_dir)

    # Save metrics
    print("\n" + "="*60)
    print("Saving results...")
    print("="*60)

    for metric_name, df in all_metrics.items():
        output_file = os.path.join(args.output_dir, f"{metric_name.upper()}.csv")
        df.to_csv(output_file)
        print(f"  {metric_name.upper()} saved to {output_file}")

    # Create summary
    summary_data = []
    for metric_name, df in all_metrics.items():
        for method in df.columns:
            summary_data.append({
                'Metric': metric_name.upper(),
                'Method': method,
                'Mean': df[method].mean(),
                'Std': df[method].std(),
                'Median': df[method].median()
            })

    summary_df = pd.DataFrame(summary_data)
    summary_file = os.path.join(args.output_dir, "summary.csv")
    summary_df.to_csv(summary_file, index=False)
    print(f"\n  Summary saved to {summary_file}")

    print("\n" + "="*60)
    print("Metric calculation completed!")
    print("="*60)


if __name__ == '__main__':
    main()
