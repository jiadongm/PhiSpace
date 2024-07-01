# Φ-Space: Continuous phenotyping for single-cell and spatial multiomics data.

Check out our talk and bioRxiv manuscript at [mixOmics website](http://mixomics.org/2024/06/phispace/).

For more case studies on spatial transcriptomics data, see our [poster](https://github.com/jiadongm/PhiSpace/blob/main/PhiSpace_poster_CBA.pdf).

## Installation

Install the GitHub version of PhiSpace:
``` r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github('jiadongm/PhiSpace/pkg')
```

## What is Φ-Space?

Φ-Space ([Mao et al., 2024](https://www.biorxiv.org/content/10.1101/2024.06.19.599787v1)) is a computational framework for reference-based continuous annotation of single-cell and spatial multiomics data. Given a bulk or single-cell RNA-seq **reference** dataset with potentially multiple layers of phenotypes defined in the metadata (e.g. cell type and sample source), Φ-Space can phenotype on a continuum the cells and cell-like objects in the **query** dataset. The core of Φ-Space is continuous phenotyping based on partial least squares (PLS) regression. Compared to conventional cell type annotation methods, Φ-Space has the following strengths:

- Identifying continuous and out-of-reference cell states;
- Robust against batch effects;
- Utilising bulk and multiomic refereneces and queries;
- Superior 

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







