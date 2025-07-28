#' PhiSpace continuous phenotyping
#'
#' PhiSpace annotates a query dataset or a list of query datasets, given an annotated bulk or single-cell
#' RNA-seq references. PhiSpace can simultaneously model multiple layers of cell phenotypes, e.g. cell type and disease condtion.
#'
#' @param reference The references. A `SingleCellExperiment` (SCE) object or a list of SCE objects. Each must contain an assay named by `refAssay`.
#' @param query The queries. An SCE object or a list of SCE object. Each must contain an assay named by `queryAssay`.
#' @param phenotypes Which phenotypes (e.g. "cell type") to predict. If `NULL`, then have to specify `response`.
#' @param response Named matrix. Rows correpond to cells (columns) in reference; columns correspond to phenotypes. If not `NULL`, then will override `phenotypes`.
#' @param refAssay Character. Which assay in reference to use to train PhiSpace.
#' @param queryAssay Character. Which assay in query to use to predict.
#' @param regMethod Character. Regression method: one of "PLS" and "PCA".
#' @param ncomp Integer. Number of components. If `NULL`, will use the default, i.e. same as the total number of phenotypes.
#' @param nfeat Integer. Number of features to choose to predict each phenotype. See details.
#' @param selectedFeat Character. Alternatively, can provide a vector of pre-selected features.
#' @param center Logic. Whether to perform centering. See details.
#' @param scale Logic. Whether to perform scaling. See details.
#' @param DRinfo Logic. Whether to return dimension reduction information from PCA or PLS. By default disabled to save memory.
#' @param storeUnNorm Store unnormalised raw PhiSpace scores or not. Default is `FALSE`.
#' @param updateRef Update reference (store reference PhiSpace scores in reference sce object) or not.
#'
#' @return
#' - If `updateRef = FALSE` (default): An updated query SCE object with PhiSpace annotation results stored in reducedDim slot "PhiSpace";
#'
#' - If `updateRef = TRUE`: A list of updated reference and query SCE objects with PhiSpace annotation results stored in reducedDim slot "PhiSpace".
#'
#' @details
#' - **Parameter tuning** By default, PhiSpace takes `ncomp` equal to the total number of phenotypes. For example, if there are
#' 20 cell types and 5 sample sources (e.g. *in vivo* and *in vitro*) defined in the reference, `ncomp` will be set to be 25.
#' By default, PhiSpace doesn't do any feature selection. However, feature selection is very straightforward in PhiSpace. The partial least squares model automatically rank all features
#' according to their contribution to predicting each phenotype. For example, if there are 25 phenotypes then PLS will return
#' 25 rankings of all features. Based on this, customised feature selection can be done in two ways. The user either
#' specify `nfeat`, i.e. the number of features (or markers) used to predict each phenotype; or specify `selectedFeat`, i.e.
#' a subset of features used to predict all phenotypes (e.g. highly variable genes). If `nfeat` is provided, then `selectedFeat`
#' will be automatically defined as the union of the `nfeat` markers of all phenotypes.
#'
#' - **Center and scale** By default, PhiSpace only center the data but do not scale. This is because we often use more features
#' than standard single-cell pipelines (e.g. Seurat selects 2,000 highly variable genes) for a better prediction. This means
#' that we might have included genes with small variance and scaling would greatly inflate the expression levels of these genes.
#'
#'
#' @references
#' Mao J., Deng, Y. and Lê Cao, K.-A. (2024).Φ-Space: Continuous phenotyping of single-cell multi-omics data. bioRxiv.
#'
#' @export
PhiSpace <- function(
    reference,
    query,
    phenotypes = NULL,
    response = NULL,
    refAssay = "rank",
    queryAssay = NULL,
    regMethod = c("PLS", "PCA"),
    ncomp = NULL,
    nfeat = NULL,
    selectedFeat = NULL,
    center = TRUE,
    scale = FALSE,
    DRinfo = FALSE,
    storeUnNorm = FALSE,
    updateRef = FALSE
){

  # Check if multiple references are provided
  if(is.list(reference)){

    # Check if references are named
    if(is.null(names(reference))){

      warning("Names of reference datasets not provided. Will use generic names.")
      names(reference) <- paste0("Ref", 1:length(reference))
    }


    sc_list <- scUnnorm_list <- vector("list", length(reference))
    names(sc_list) <- names(reference)
    for(ii in 1:length(reference)){

      refDataName <- names(reference)[ii]
      refSingle <- reference[[refDataName]]

      PhiRes <- PhiSpaceR_1ref(
        reference = refSingle,
        query = query,
        phenotypes = phenotypes,
        response = response,
        refAssay = refAssay,
        queryAssay = queryAssay,
        regMethod = regMethod,
        ncomp = ncomp,
        nfeat = nfeat,
        selectedFeat = selectedFeat,
        center = center,
        scale = scale,
        DRinfo = DRinfo
      )

      # Rename cell type names by appending reference dataset name
      sc <- PhiRes$PhiSpaceScore

      if(is.list(sc)){

        sc <- lapply(
          sc,
          function(x){
            colnames(x) <- paste0(colnames(x), "(", refDataName, ")")
            return(x)
          }
        )

        sc_norm <- lapply(sc, normPhiScores)

      } else {

        colnames(sc) <- paste0(colnames(sc), "(", refDataName, ")")
        sc_norm <- normPhiScores(sc)
      }

      sc_list[[refDataName]] <- sc_norm
      if(storeUnNorm) scUnnorm_list[[refDataName]] <- sc
    }


    # If multiple queries were provided
    if(is.list(query)){

      for(i in 1:length(query)){

        scSingleQuery_list <- lapply(sc_list, function(x) x[[i]])
        reducedDim(query[[i]], "PhiSpace") <- Reduce(cbind, scSingleQuery_list)

        if(storeUnNorm){

          scSingleQuery_list <- lapply(scUnnorm_list, function(x) x[[i]])
          reducedDim(query[[i]], "PhiSpaceNonNorm") <- Reduce(cbind, scSingleQuery_list)
        }
      }
    } else {

      reducedDim(query, "PhiSpace") <- Reduce(cbind, sc_list)
      if(storeUnNorm) reducedDim(query, "PhiSpaceNonNorm") <-  Reduce(cbind, scUnnorm_list)
    }

    return(query)

  } else { # If single reference was provided

    PhiRes <- PhiSpaceR_1ref(
      reference = reference,
      query = query,
      phenotypes = phenotypes,
      response = response,
      refAssay = refAssay,
      queryAssay = queryAssay,
      regMethod = regMethod,
      ncomp = ncomp,
      nfeat = nfeat,
      selectedFeat = selectedFeat,
      center = center,
      scale = scale,
      DRinfo = DRinfo
    )

    if(!is.list(query)){

      if(storeUnNorm) reducedDim(query, "PhiSpaceNonNorm") <- PhiRes$PhiSpaceScore
      reducedDim(query, "PhiSpace") <- normPhiScores(PhiRes$PhiSpaceScore)
    } else {

      for(i in 1:length(query)){

        if(storeUnNorm) reducedDim(query[[i]], "PhiSpaceNonNorm") <- PhiRes$PhiSpaceScore[[i]]
        reducedDim(query[[i]], "PhiSpace") <- normPhiScores(PhiRes$PhiSpaceScore[[i]])
      }
    }

    if(updateRef){

      if(storeUnNorm) reducedDim(reference, "PhiSpaceNonNorm") <- PhiRes$YrefHat
      reducedDim(reference, "PhiSpace") <- normPhiScores(PhiRes$YrefHat)

      return(
        list(
          reference = reference,
          query = query
        )
      )
    } else {

      return(query)
    }
  }







}
