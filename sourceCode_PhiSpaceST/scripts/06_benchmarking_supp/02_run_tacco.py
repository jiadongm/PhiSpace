#!/usr/bin/env python3
"""
TACCO Deconvolution Pipeline for Benchmarking.

Reads prepared CSV data from the benchmarking framework and runs TACCO
(optimal transport-based cell type annotation).

Usage:
    python 02_run_tacco.py <dataset_dir> <output_dir>

Must be run with the TACCO conda environment (set its absolute path
in config/paths_local.yml as `tacco_python`).
"""

import sys
import os
import time
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import tacco as tc

import warnings
warnings.filterwarnings('ignore')


def main():
    if len(sys.argv) < 3:
        print("Usage: python 02_run_tacco.py <dataset_dir> <output_dir>")
        sys.exit(1)

    dataset_dir = sys.argv[1]
    output_dir = sys.argv[2]

    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print("TACCO Benchmarking Pipeline")
    print("=" * 60)
    print(f"Dataset directory: {dataset_dir}")
    print(f"Output directory: {output_dir}")

    # ==========================================================================
    # Load data from CSVs
    # ==========================================================================
    print("\nLoading data...")

    # Load scRNA-seq counts (genes x cells CSV, first column = gene names)
    scrna_counts = pd.read_csv(os.path.join(dataset_dir, "scRNA_counts.csv"), index_col=0)
    print(f"  scRNA counts: {scrna_counts.shape[0]} genes x {scrna_counts.shape[1]} cells")

    # Load cell type annotations
    celltype_df = pd.read_csv(os.path.join(dataset_dir, "scRNA_celltype.csv"))
    celltypes = celltype_df.set_index('cell')['celltype']
    print(f"  Cell types: {celltypes.nunique()} unique types")

    # Load spatial counts (genes x spots CSV, first column = gene names)
    spatial_counts = pd.read_csv(os.path.join(dataset_dir, "spatial_counts.csv"), index_col=0)
    print(f"  Spatial counts: {spatial_counts.shape[0]} genes x {spatial_counts.shape[1]} spots")

    # ==========================================================================
    # Convert to AnnData (cells x genes)
    # ==========================================================================
    print("\nPreparing AnnData objects...")

    # Transpose: CSV is genes x cells -> AnnData needs cells x genes
    adata_ref = anndata.AnnData(
        X=scrna_counts.values.T.astype(np.float32),
        obs=pd.DataFrame(index=scrna_counts.columns),
        var=pd.DataFrame(index=scrna_counts.index)
    )
    adata_ref.obs['celltype'] = celltypes.reindex(adata_ref.obs_names).values
    print(f"  Reference AnnData: {adata_ref.shape}")

    adata_spatial = anndata.AnnData(
        X=spatial_counts.values.T.astype(np.float32),
        obs=pd.DataFrame(index=spatial_counts.columns),
        var=pd.DataFrame(index=spatial_counts.index)
    )
    print(f"  Spatial AnnData: {adata_spatial.shape}")

    # ==========================================================================
    # Filter cell types with <=1 cell
    # ==========================================================================
    ct_counts = adata_ref.obs['celltype'].value_counts()
    rare_cts = ct_counts[ct_counts <= 1].index
    if len(rare_cts) > 0:
        print(f"  Removing {len(rare_cts)} cell types with <=1 cell: {list(rare_cts)}")
        adata_ref = adata_ref[~adata_ref.obs['celltype'].isin(rare_cts)].copy()

    adata_ref.obs['celltype'] = pd.Categorical(adata_ref.obs['celltype'])

    print(f"  Reference after filtering: {adata_ref.shape}")
    print(f"  Cell types: {adata_ref.obs['celltype'].nunique()}")

    # ==========================================================================
    # Run TACCO
    # ==========================================================================
    print("\nRunning TACCO (OT method)...")
    start_time = time.time()

    result_df = tc.tl.annotate(
        adata_spatial,
        adata_ref,
        annotation_key='celltype',
        method='OT',
        result_key=None
    )

    end_time = time.time()
    runtime = end_time - start_time

    print(f"  TACCO completed in {runtime:.2f} seconds")
    print(f"  Result shape: {result_df.shape}")

    # ==========================================================================
    # Normalize rows to sum to 1
    # ==========================================================================
    row_sums = result_df.sum(axis=1)
    result_df = result_df.div(row_sums, axis=0).fillna(0)

    # Use spot names as index
    result_df.index = adata_spatial.obs_names

    # ==========================================================================
    # Save results
    # ==========================================================================
    print("\nSaving results...")

    result_file = os.path.join(output_dir, "TACCO_result.txt")
    result_df.to_csv(result_file)
    print(f"  Results saved to {result_file}")

    # Save runtime
    runtime_file = os.path.join(output_dir, "TACCO_runtime.txt")
    pd.DataFrame({
        'method': ['TACCO'],
        'runtime_seconds': [runtime]
    }).to_csv(runtime_file, index=False)
    print(f"  Runtime saved to {runtime_file}")

    print(f"\nTACCO completed in {runtime:.2f} seconds")
    print("=" * 60)


if __name__ == '__main__':
    main()
