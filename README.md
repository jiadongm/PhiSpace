# Φ-Space: Continuous phenotyping of single-cell and spatial omics data.


<p float="left">
  <img src="./figs/PhiSpace.png" width="30%" /> 
  <img src="./figs/PhiSpaceST.png" width="30%" />
</p>


The PhiSpace R package moved to [github.com/mao-lab-org/PhiSpace](https://github.com/mao-lab-org/PhiSpace) on 12 October 2026. Report issues and follow development there. This repository keeps the papers' analysis code (`sourceCode_PhiSpaceMultiomics/`, `sourceCode_PhiSpaceST/`) and posters. Its `pkg/` folder is the package as of the move and is no longer updated.


## Vignettes

[PhiSpace pkgdown site](https://mao-lab-org.github.io/PhiSpace/) 

[A workflow for spatial transcriptomics featuring SpaNorm, PhiSpace and jazzPanda](https://github.com/bhuvad/CosMxSpatialAnalysisWorkshop)


## Installation


``` r
if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install('mao-lab-org/PhiSpace')
```


## Citing PhiSpace

Mao, Jiadong, Choi, Jarny and Lê Cao, Kim-Anh. (2026). Φ-Space ST: a platform-agnostic method to identify cell states in spatial transcriptomics studies. *Cell Reports Methods* 6(8), 101483. [doi:10.1016/j.crmeth.2026.101483](https://doi.org/10.1016/j.crmeth.2026.101483).

Mao, Jiadong, Deng, Yidi and Lê Cao, Kim-Anh. (2025). Φ-Space: continuous phenotyping of single-cell multi-omics data. *Genome Biology* 26, 323. [doi:10.1186/s13059-025-03755-8](https://doi.org/10.1186/s13059-025-03755-8).






