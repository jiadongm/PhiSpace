# Φ-Space: Continuous phenotyping of single-cell and spatial omics data.


<p float="left">
  <img src="./figs/PhiSpace.png" width="30%" /> 
  <img src="./figs/PhiSpaceST.png" width="30%" />
</p>


## Vignettes

[PhiSpace pkgdown site](https://jiadongm.github.io/PhiSpace/) 



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


## Other resources

Mao, Jiadong, Choi, Jarny and Lê Cao, Kim-Anh. (2025). Φ-Space ST: a platform-agnostic method to identify cell states in spatial transcriptomics studies. [*bioRxiv*](https://www.biorxiv.org/content/10.1101/2025.02.05.636735v1).

Mao, Jiadong, Deng, Yidi and Lê Cao, Kim-Anh. (2024). Φ-Space: Continuous phenotyping of single-cell multi-omics data. [*bioRxiv*](https://www.biorxiv.org/content/10.1101/2024.06.19.599787v1.full).

Check out our talk for PhiSpace single-cell multiomics at [mixOmics website](http://mixomics.org/2024/06/phispace/).

Check out our poster for PhiSpace ST at [Lorne Genome 2025](https://github.com/jiadongm/PhiSpace/blob/main/PhiSpaceST_LorneGenome2025.pdf)





