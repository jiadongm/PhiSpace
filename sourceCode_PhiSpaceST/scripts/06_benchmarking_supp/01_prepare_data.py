#!/usr/bin/env python3
"""
Prepare data for PhiSpace benchmarking.
Converts h5ad files to formats suitable for R/PhiSpace.
Extracts ground truth from simulated datasets.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmwrite
from scipy.sparse import csr_matrix

def prepare_dataset(dataset_num, data_dir, output_dir):
    """
    Prepare a single simulated dataset for PhiSpace benchmarking.

    Parameters:
    -----------
    dataset_num : int
        Dataset number (1-32)
        Datasets 1-16 are named "dataset1" to "dataset16"
        Datasets 17-32 are named "dataset1_r" to "dataset16_r"
    data_dir : str
        Path to SimualtedSpatalData directory
    output_dir : str
        Path to output directory
    """

    dataset_path = os.path.join(data_dir, f"dataset{dataset_num}")

    if not os.path.exists(dataset_path):
        print(f"Dataset {dataset_num} not found at {dataset_path}")
        return False

    print(f"\n{'='*50}")
    print(f"Processing dataset {dataset_num}")
    print(f"{'='*50}")

    # Create output directory for this dataset
    ds_output_dir = os.path.join(output_dir, f"dataset{dataset_num}")
    os.makedirs(ds_output_dir, exist_ok=True)

    # Load scRNA-seq reference data
    scrna_path = os.path.join(dataset_path, "scRNA.h5ad")
    print(f"Loading scRNA-seq data from {scrna_path}")
    scrna = sc.read_h5ad(scrna_path)
    print(f"  Shape: {scrna.shape}")
    print(f"  Cell types: {scrna.obs.columns.tolist()}")

    # Identify celltype column
    celltype_col = None
    for col in ['celltype', 'cell_type', 'CellType', 'cellType', 'cluster']:
        if col in scrna.obs.columns:
            celltype_col = col
            break

    if celltype_col is None:
        # Try to find any column with 'type' or 'cluster' in name
        for col in scrna.obs.columns:
            if 'type' in col.lower() or 'cluster' in col.lower():
                celltype_col = col
                break

    if celltype_col is None:
        print(f"  WARNING: Could not identify celltype column. Available: {scrna.obs.columns.tolist()}")
        celltype_col = scrna.obs.columns[0]

    print(f"  Using celltype column: {celltype_col}")
    print(f"  Cell types found: {scrna.obs[celltype_col].nunique()}")

    # Load spatial data
    spatial_path = os.path.join(dataset_path, "Spatial.h5ad")
    print(f"Loading spatial data from {spatial_path}")
    spatial = sc.read_h5ad(spatial_path)
    print(f"  Shape: {spatial.shape}")

    # Extract ground truth cell type proportions
    if 'density' in spatial.uns:
        ground_truth = spatial.uns['density']
        print(f"  Ground truth shape: {ground_truth.shape}")
    else:
        print("  WARNING: No ground truth 'density' found in spatial.uns")
        ground_truth = None

    # Save ground truth
    if ground_truth is not None:
        gt_path = os.path.join(ds_output_dir, "ground_truth.csv")
        ground_truth.to_csv(gt_path)
        print(f"  Saved ground truth to {gt_path}")

    # Save scRNA-seq count matrix (genes x cells)
    scrna_counts_path = os.path.join(ds_output_dir, "scRNA_counts.csv")
    if hasattr(scrna.X, 'toarray'):
        counts_df = pd.DataFrame(
            scrna.X.toarray().T,
            index=scrna.var_names,
            columns=scrna.obs_names
        )
    else:
        counts_df = pd.DataFrame(
            scrna.X.T,
            index=scrna.var_names,
            columns=scrna.obs_names
        )
    counts_df.to_csv(scrna_counts_path)
    print(f"  Saved scRNA counts to {scrna_counts_path}")

    # Save scRNA-seq cell type annotations
    celltype_path = os.path.join(ds_output_dir, "scRNA_celltype.csv")
    celltype_df = pd.DataFrame({
        'cell': scrna.obs_names,
        'celltype': scrna.obs[celltype_col].values
    })
    celltype_df.to_csv(celltype_path, index=False)
    print(f"  Saved cell types to {celltype_path}")

    # Save spatial count matrix (genes x spots)
    spatial_counts_path = os.path.join(ds_output_dir, "spatial_counts.csv")
    if hasattr(spatial.X, 'toarray'):
        spatial_df = pd.DataFrame(
            spatial.X.toarray().T,
            index=spatial.var_names,
            columns=spatial.obs_names
        )
    else:
        spatial_df = pd.DataFrame(
            spatial.X.T,
            index=spatial.var_names,
            columns=spatial.obs_names
        )
    spatial_df.to_csv(spatial_counts_path)
    print(f"  Saved spatial counts to {spatial_counts_path}")

    # Save metadata file
    metadata = {
        'dataset_num': dataset_num,
        'n_cells': scrna.n_obs,
        'n_spots': spatial.n_obs,
        'n_genes_scrna': scrna.n_vars,
        'n_genes_spatial': spatial.n_vars,
        'celltype_column': celltype_col,
        'n_celltypes': scrna.obs[celltype_col].nunique()
    }
    metadata_path = os.path.join(ds_output_dir, "metadata.csv")
    pd.DataFrame([metadata]).to_csv(metadata_path, index=False)
    print(f"  Saved metadata to {metadata_path}")

    return True


def main():
    parser = argparse.ArgumentParser(description='Prepare data for PhiSpace benchmarking')
    parser.add_argument('--data_dir', type=str,
                        default='ExampleData/SimualtedSpatalData',
                        help='Path to SimualtedSpatalData directory')
    parser.add_argument('--output_dir', type=str,
                        default='PhiSpace_Benchmarking/prepared_data',
                        help='Output directory for prepared data')
    parser.add_argument('--datasets', type=str, default='all',
                        help='Comma-separated dataset numbers or "all"')

    args = parser.parse_args()

    # Get absolute paths
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(base_dir, args.data_dir)
    output_dir = os.path.join(base_dir, args.output_dir)

    os.makedirs(output_dir, exist_ok=True)

    # Determine which datasets to process
    if args.datasets == 'all':
        datasets = list(range(1, 33))
    else:
        datasets = [int(x.strip()) for x in args.datasets.split(',')]

    print(f"Data directory: {data_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Datasets to process: {datasets}")

    # Process each dataset
    success_count = 0
    for ds_num in datasets:
        if prepare_dataset(ds_num, data_dir, output_dir):
            success_count += 1

    print(f"\n{'='*50}")
    print(f"Completed: {success_count}/{len(datasets)} datasets processed successfully")
    print(f"{'='*50}")


if __name__ == '__main__':
    main()
