# Φ-Space: continuous cell state annotation for single-cell and spatial omics studies


## Installation

Installing BioConductor dependencies:
``` r
if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("SingleCellExperiment", "scran", "scuttle", "ComplexHeatmap", "SpatialExperiment"))
```

Install the GitHub version of PhiSpace:
``` r
if(!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github('jiadongm/PhiSpace/pkg')
```



## What is Φ-Space?

Φ-Space is a computational framework for reference-based continuous annotation 
of cell states in single-cell and spatial multiomics data. Currently it has two
modules:

- PhiSpace multiomics for single-cell multiomics data ([Mao et al., 2024](https://www.biorxiv.org/content/10.1101/2024.06.19.599787v1), and
- PhiSpace ST for spatial transcriptomics data [Mao et al., 2025](https://www.biorxiv.org/content/10.1101/2024.06.19.599787v1))

Given a bulk or single-cell RNA-seq **reference** dataset with potentially multiple 
layers of phenotypes defined in the metadata (e.g. cell type and sample source), 
Φ-Space can phenotype on a continuum the cells and cell-like objects in the **query** dataset. 
The core of Φ-Space is continuous phenotyping based on partial least squares (PLS) regression. 
Compared to conventional cell type annotation methods, Φ-Space has the following strengthss

- Identifying continuous and out-of-reference cell states;
- Robust against batch effects;
- Utilising bulk and multiomic refereneces and queries;
- More suitable for exploring the spatial patterns of rare cell types. 


We have applied Φ-Space to many different use cases, including

| Reference     |      Query    | Note   | Vignettes |
| :------ |    :------   | :------  | :------ |
| bulk RNA-seq  |   scRNA-seq   |   | [Stemformatics DC atlas](articles/getting_started.html) |
| scRNA-seq     |   scRNA-seq   |  | |
| scRNA-seq     |   scATAC-seq  | requires a bimodal bridge dataset | |
| CITE-seq (scRNA+Protein-seq)  |  CITE-seq   | using both modalities | |
| scRNA-seq    |   subcellular spatial transcriptomics | e.g. Stereo-seq, CosMx, 10x Xenium | [CosMx lung cancer microenvironment](articles/CosMx.html) |
| scRNA-seq    |   supercellular spatial transcriptomics | e.g. 10x Visium, Slide-seqV2 | [Cell type deconvolution for Visium](articles/Visium.html) |



## Cite PhiSpace

Mao, Jiadong, Choi, Jarny and Lê Cao, Kim-Anh. (2025). Φ-Space ST: a platform-agnostic method to identify cell states in spatial transcriptomics studies. [*bioRxiv*](https://www.biorxiv.org/content/10.1101/2025.02.05.636735v1).

Mao, Jiadong, Deng, Yidi and Lê Cao, Kim-Anh. (2024). Φ-Space: Continuous phenotyping of single-cell multi-omics data. [*bioRxiv*](https://www.biorxiv.org/content/10.1101/2024.06.19.599787v1.full).

Check out our talk for PhiSpace single-cell multiomics at [mixOmics website](http://mixomics.org/2024/06/phispace/).

