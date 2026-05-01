#!/usr/bin/env python3
import numpy as np
import sys
import pandas as pd
import squidpy as sq
import matplotlib.pyplot as plt
import stereo as st


# from stereopy package but added offset from raw gef file
def parse_bin_coor(df, x_min, y_min, bin_size, x_raw_offset = 0, y_raw_offset = 0):
    """
    merge bins to a bin unit according to the bin size, also calculate the center coordinate of bin unit,
    and generate cell id of bin unit using the coordinate after merged.

    :param df: a dataframe of the bin file.
    :param bin_size: the size of bin to merge.
    :return:
    """
    df['bin_x'] = merge_bin_coor(df['x'].values, x_min, bin_size, x_raw_offset)
    df['bin_y'] = merge_bin_coor(df['y'].values, y_min, bin_size, y_raw_offset)
    df['cell_id'] = df['bin_x'].astype(str) + '_' + df['bin_y'].astype(str)
    df['x_center'] = get_bin_center(df['bin_x'], x_min, bin_size)
    df['y_center'] = get_bin_center(df['bin_y'], y_min, bin_size)
    
def merge_bin_coor(coor: np.ndarray, coor_min: int, bin_size: int, raw_offset: int):
    return np.floor((coor - coor_min - raw_offset) / bin_size).astype(np.int32)


def get_bin_center(bin_coor: np.ndarray, coor_min: int, bin_size: int):
    return bin_coor * bin_size + coor_min + int(bin_size / 2)


def bin_barcodes(bc, minX=False, minY=False, gem_tissue_path=None, gef_raw_path=None, bin_size=100):

    # in case gem file was already loaded, can directly specify minX and minY
    # this saves time
    if not minX and not minY:
        if gem_tissue_path == None:
            sys.exit("Specify either minX and minY or path to gem tissue file.")
        else:
            attr_tissue = st.io.read_gem(file_path=gem_tissue_path, bin_size=bin_size).attr
            minX = attr_tissue["minX"]
            minY = attr_tissue["minY"]

    attr_raw = st.io.read_gef_info(file_path=gef_raw_path)

    bc[['x','y']] = bc['cell'].str.split('_',expand=True)
    bc['x'] = bc['x'].astype("int")
    bc['y'] = bc['y'].astype("int")

    parse_bin_coor(
        bc, 
        minX,
        minY,
        bin_size,
        attr_raw["offsetX"],
        attr_raw["offsetY"]
    )

    bc["count_binned"] = bc.groupby(["cell_id", "barcode"])["count"].transform(sum)

    bc = bc[["barcode", "count_binned", "bin_x", "bin_y", "cell_id", "x_center", "y_center"]].drop_duplicates().reset_index(drop=True)

    return bc


# function to plot barcode
def plot_barcode(adata, barcode, axs, **kwargs):
    adata.obs["barcode"] = adata.obs["barcode"].replace('nan', np.nan)
    if "nan" not in adata.obs.columns:
        # make nan column for plotting
        adata.obs["nan"] = np.nan
        adata.obs["nan"] = adata.obs["nan"].astype("category")

    sq.pl.spatial_scatter(
        adata, color="nan", size=1, shape=None, frameon=False, legend_loc=False, ax=axs,
        na_color=(0.9, 0.9, 0.9, 0.01)
    )
    sq.pl.spatial_scatter(
        adata[adata.obs["barcode"] == barcode, :], 
        color="barcode", size=1, shape=None, frameon=False,
        title=barcode,
        palette="Set1",
        ax=axs,
        legend_loc=False, 
        **kwargs
    )

def plot_barcode_grid(adata, barcodes, figsize=(20,20), cols=3, invert=False, **kwargs):
    fig = plt.figure(figsize=figsize)
    for i, barcode in enumerate(barcodes):
        ax = fig.add_subplot(int(np.ceil(len(barcodes) / cols)), cols, i + 1)
        plot_barcode(adata, barcode, ax, **kwargs)
        if invert:
            ax.invert_yaxis()
    plt.tight_layout()