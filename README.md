# Φ-Space: Continuous phenotyping for single-cell and spatial multiomics data.

Mao, Jiadong, Deng, Yidi and Lê Cao, Kim-Anh. (2024). Φ-Space: Continuous phenotyping of single-cell multi-omics data. *bioRxiv*.

Check out our talk and bioRxiv manuscript at [mixOmics website](http://mixomics.org/2024/06/phispace/).

For more case studies on spatial transcriptomics data, see our [poster](https://github.com/jiadongm/PhiSpace/blob/main/PhiSpace_poster_CBA.pdf).

## Installation

Installing BioConductor dependencies:
``` r
if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("SingleCellExperiment", "scran", "scuttle", "ComplexHeatmap"))
```

Install the GitHub version of PhiSpace:
``` r
if(!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github('jiadongm/PhiSpace/pkg')
```



## What is Φ-Space?

Φ-Space ([Mao et al., 2024](https://www.biorxiv.org/content/10.1101/2024.06.19.599787v1)) is a computational framework for reference-based continuous annotation of single-cell and spatial multiomics data. Given a bulk or single-cell RNA-seq **reference** dataset with potentially multiple layers of phenotypes defined in the metadata (e.g. cell type and sample source), Φ-Space can phenotype on a continuum the cells and cell-like objects in the **query** dataset. The core of Φ-Space is continuous phenotyping based on partial least squares (PLS) regression. Compared to conventional cell type annotation methods, Φ-Space has the following strengths:

- Identifying continuous and out-of-reference cell states;
- Robust against batch effects;
- Utilising bulk and multiomic refereneces and queries;
- More suitable for exploring the spatial patterns of rare cell types. 

<img src="./figs/schema_core.png" width="60%" style="display: block; margin: auto;" />


We have applied Φ-Space to many different use cases, including

| Reference     |      Query    | Note   |
| :------ |    :------   | :------  |
| bulk RNA-seq  |   scRNA-seq   |   |
| scRNA-seq     |   scRNA-seq   |  |
| scRNA-seq     |   scATAC-seq  | requires a bimodal bridge dataset |
| CITE-seq (scRNA+Protein-seq)  |  CITE-seq   | using both modalities |
| scRNA-seq    |   subcellular spatial transcriptomics | e.g. Stereo-seq, CosMx, 10x Xenium |
| scRNA-seq    |   supercellular spatial transcriptomics | e.g. 10x Visium, Slide-seqV2 |


## Example: transitional identities of induced dendritic cells (DCs)

We illustrate how Φ-Space works using the first case study in our [manuscript](https://www.biorxiv.org/content/10.1101/2024.06.19.599787v1).

- Reference bulk RNA-seq: [Stemformatics DC atlas](https://journals.aai.org/jimmunol/article/209/12/2352/237295/The-Human-Dendritic-Cell-Atlas-An-Integrated)
- Query scRNA-seq: [Rosa et al. (2018)](https://www.science.org/doi/10.1126/sciimmunol.aau4292) 

Dendritic cells (DCs) are a type of immune cells. DCs are relatively rare in human blood samples. Hence it is desirable to culture *in vivo* like 
DCs using *in vitro* methods. Rosa et al. (2018) claimed that they successfully reprogrammed human esophagus fibroblasts (HEFs) into induced DCs after 9 days *in vitro* cell culturing. 

The Stemformatics DC atlas is a bulk RNA-seq atlas of different subtypes of human DC samples (FACs sorted). The DC atlas contains 


### Read data

Load packages
``` r
suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ComplexHeatmap))
```

Download the processed reference dataset [ref_dc.rds](https://unimelbcloud-my.sharepoint.com/:u:/g/personal/jiadong_mao_unimelb_edu_au/EZVy-qceLC5Ik9YQ9yiASM8BRy0eKn4KYj_fy5A8LVdifA?e=NhIuTt) and the processed and downsampled query dataset [query_Rosa_sub.rds](https://unimelbcloud-my.sharepoint.com/:u:/g/personal/jiadong_mao_unimelb_edu_au/Eep7PpTnTHJIirmK8EM6JGsBRrxRlx_Soqk5DT-8KiheNQ?e=MvFzNA). In addtion, we download the selected genes [ref_dc_feat.rds](https://unimelbcloud-my.sharepoint.com/:u:/g/personal/jiadong_mao_unimelb_edu_au/EYW4m1WMtxhNg9vTUFQdZAQB12sF0VOj3u2pmz3Uce5U6A?e=zdyv2a). See our [manuscript](https://www.biorxiv.org/content/10.1101/2024.06.19.599787v1) for a description of feature selection. 

``` r
dat_dir <- "/data/projects/punim0613/JiaDong/PhiSpace/" # replace this by your own directory where you store ref_dc.rds and query_dc.rds
query <- readRDS(paste0(dat_dir, "query_Rosa_sub.rds"))
reference <- readRDS(paste0(dat_dir,"ref_dc.rds"))
selectedFeat <- readRDS(paste0(dat_dir, "ref_dc_test.rds"))

# Rank normalise reference and query
query <- RankTransf(query, "counts")
reference <- RankTransf(reference, "data", sparse = F)

# Customised colour code
DC_cols <- c(
  `DC precursor` = "#B3B3B3",DC_prec = "#FFFFB3",
  MoDC = "#CCEBC5", cDC1 = "#1B9E77", 
  cDC2 = "#D95F02",`dendritic cell` = "#B3B3B3", 
  DC = "#B3DE69", monocyte = "#B3B3B3",
  mono = "#80B1D3",`plasmacytoid dendritic cell` = "#7570B3",
  DC1 = "#1B9E77",DC2 = "#D95F02",
  pDC = "#7570B3",HEF = "#FEE5D9", 
  Day3 = "#FCAE91",Day6 = "#FB6A4A",
  Day9_SP = "#DE2D26", Day9_DP = "#A50F15",
  Day9 = "#A50F15",other = "#B3B3B3"
)

DC_cols_source <- c(
  ex_vivo = "#B3B3B3",in_vitro = "#A50026",
  in_vivo = "#313695",`in_vivo_HuMouse` = "#80B1D3", query = "#B3B3B3"
)

DC_symbs <- c(
  query = 4, reference = 16
)
```

In the above code, we applied rank transform to both reference and query. This is because the rank transform is more appropriate for this particular reference dataset ([Elahi et al., 2022](https://journals.aai.org/jimmunol/article/209/12/2352/237295/The-Human-Dendritic-Cell-Atlas-An-Integrated)). In general, normalisation methods should be chosen to suit individual cases. The only requirement is that both reference and query are normalised in the same way. **No** additional harmonisation of reference and query is needed.


### Continuous phenotyping

Now we are ready to apply PhiSpace to continuously phenotype the query cells. 
``` r
PhiSpaceAssay <- "rank"
phenotypes <- c("Cell Type", "Sample Source")
PhiMethod <- "PLS"

PhiResPath <- paste0(dat_dir, "PhiRes.rds")

if(!file.exists(PhiResPath)){
  
  PhiRes <- PhiSpaceR(
    reference, 
    query,
    ncomp = 30,
    selectedFeat = selectedFeat,
    phenotypes = phenotypes, 
    PhiSpaceAssay = PhiSpaceAssay,
    regMethod = PhiMethod,
    scale = FALSE, 
    DRinfo = TRUE
  )
  
  saveRDS(PhiRes, PhiResPath)
} else {
  
  PhiRes <- readRDS(PhiResPath)
}
```

### Visualise annotation





