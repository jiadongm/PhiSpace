#!/usr/bin/env python3
"""
Cell2location Deconvolution Pipeline for Benchmarking.

Reads prepared CSV data from the benchmarking framework and runs Cell2location.

Usage:
    python 02_run_cell2location.py <dataset_dir> <output_dir>

Must be run with the cell2location conda environment (set its absolute
path in config/paths_local.yml as `cell2loc_python`).
"""

import sys
import os
import time
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from scipy.sparse import csr_matrix

import cell2location
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel

import warnings
warnings.filterwarnings('ignore')


def main():
    if len(sys.argv) < 3:
        print("Usage: python 02_run_cell2location.py <dataset_dir> <output_dir>")
        sys.exit(1)

    dataset_dir = sys.argv[1]
    output_dir = sys.argv[2]

    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print("Cell2location Benchmarking Pipeline")
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
    celltype_df['cell'] = celltype_df['cell'].astype(str)
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
        X=csr_matrix(scrna_counts.values.T.astype(np.float32)),
        obs=pd.DataFrame(index=scrna_counts.columns),
        var=pd.DataFrame(index=scrna_counts.index)
    )
    # Add cell type annotations
    adata_ref.obs['celltype'] = celltypes.reindex(adata_ref.obs_names).values
    print(f"  Reference AnnData: {adata_ref.shape}")

    adata_spatial = anndata.AnnData(
        X=csr_matrix(spatial_counts.values.T.astype(np.float32)),
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
    adata_ref = adata_ref[~adata_ref.obs['celltype'].isna(), :].copy()

    # Remove cells and genes with 0 counts
    sc.pp.filter_genes(adata_ref, min_cells=1)
    sc.pp.filter_cells(adata_ref, min_genes=1)

    # ==========================================================================
    # Filter genes using cell2location's filter
    # ==========================================================================
    print("\nFiltering genes...")
    selected = filter_genes(
        adata_ref,
        cell_count_cutoff=5,
        cell_percentage_cutoff2=0.03,
        nonz_mean_cutoff=1.12
    )
    adata_ref = adata_ref[:, selected].copy()
    print(f"  After filtering: {adata_ref.shape}")

    # ==========================================================================
    # Train Regression Model
    # ==========================================================================
    print("\nTraining RegressionModel...")
    start_time = time.time()

    RegressionModel.setup_anndata(adata=adata_ref, labels_key='celltype')
    mod = RegressionModel(adata_ref)
    mod.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002)

    # Export posterior (100 samples is sufficient and reduces memory usage)
    adata_ref = mod.export_posterior(
        adata_ref,
        sample_kwargs={'num_samples': 100, 'batch_size': 2500}
    )

    # Get inferred average expression per cluster
    if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
        inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[
            f'means_per_cluster_mu_fg_{i}'
            for i in adata_ref.uns['mod']['factor_names']
        ]].copy()
    else:
        inf_aver = adata_ref.var[[
            f'means_per_cluster_mu_fg_{i}'
            for i in adata_ref.uns['mod']['factor_names']
        ]].copy()
    inf_aver.columns = adata_ref.uns['mod']['factor_names']

    print(f"  Inferred signatures: {inf_aver.shape}")

    # ==========================================================================
    # Intersect genes and train Cell2location model
    # ==========================================================================
    print("\nTraining Cell2location model...")

    intersect_genes = np.intersect1d(adata_spatial.var_names, inf_aver.index)
    print(f"  Intersecting genes: {len(intersect_genes)}")

    adata_spatial = adata_spatial[:, intersect_genes].copy()
    inf_aver = inf_aver.loc[intersect_genes, :].copy()

    cell2location.models.Cell2location.setup_anndata(adata=adata_spatial)

    mod = cell2location.models.Cell2location(
        adata_spatial,
        cell_state_df=inf_aver,
        N_cells_per_location=30,
        detection_alpha=200
    )

    mod.train(
        max_epochs=30000,
        batch_size=None,
        train_size=1
    )

    # Export posterior (100 samples is sufficient and reduces memory usage)
    adata_spatial = mod.export_posterior(
        adata_spatial,
        sample_kwargs={
            'num_samples': 100,
            'batch_size': mod.adata.n_obs
        }
    )

    end_time = time.time()
    runtime = end_time - start_time

    # ==========================================================================
    # Save results
    # ==========================================================================
    print("\nSaving results...")

    result = adata_spatial.obsm['q05_cell_abundance_w_sf']
    result_file = os.path.join(output_dir, "Cell2location_result.txt")
    result.to_csv(result_file)
    print(f"  Results saved to {result_file}")

    # Save runtime
    runtime_file = os.path.join(output_dir, "Cell2location_runtime.txt")
    pd.DataFrame({
        'method': ['Cell2location'],
        'runtime_seconds': [runtime]
    }).to_csv(runtime_file, index=False)
    print(f"  Runtime saved to {runtime_file}")

    print(f"\nCell2location completed in {runtime:.2f} seconds")
    print("=" * 60)


if __name__ == '__main__':
    main()
