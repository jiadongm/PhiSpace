#' Tune ncomp using cross-validation.
#'
#' @param reference `SingleCellExperiment` object.
#' @param YY Matrix.
#' @param ncompLimits Vector.
#' @param gridSize Integer.
#' @param ncompGrid Vector.
#' @param selectedFeat Vector.
#' @param assayName Character.
#' @param Kfolds Integer.
#' @param seed Integer.
#' @param regMethod Character.
#' @param Ncores Integer.
#' @param labelCode Character.
#' @param phenotypes Character.
#' @param center Logic.
#' @param scale Logic.
#'
#' @return A list.
#'
#' @import ggplot2
CVTune_ncomp <- function(reference,
                         assayName = "logcounts",
                         phenotypes = NULL,
                         YY = NULL,
                         #
                         center = TRUE,
                         scale = FALSE,
                         #
                         ncompLimits = c(1,100),
                         gridSize = 20,
                         ncompGrid = NULL,
                         selectedFeat = NULL,
                         #
                         Kfolds = 5,
                         regMethod = 'PCA',
                         labelCode = "-1,1",
                         Ncores = 1,
                         seed = 5202056)
{
  ## Prepare predictor
  XX <- as.matrix(t(assay(reference, assayName)))

  ## Prepare ncompGrid
  if(is.null(ncompGrid)){

    ncompMax <- max(ncompLimits)
    ncompMin <- min(ncompLimits)
    ncompGrid <-
      round(
        exp(
          seq(log(ncompMin), log(ncompMax), length.out = gridSize)
        )
      )
    ncompGrid <- unique(ncompGrid)

  }

  ## Prepare response: dummy matrix
  if(is.null(YY)){
    if(is.null(phenotypes)) stop("Have to specify either YY or phenotypes.")
    YY <- codeY(reference, phenotypes)
  }


  ## Cross-validation
  # Create data partition
  set.seed(seed)
  Nsamples <- nrow(XX)
  allIdx <- 1:Nsamples # indices of all samples
  permuIdx <- sample(allIdx, Nsamples) # random permutation of allIdx
  foldSize <- floor(Nsamples/Kfolds) # input of split f below
  splitIdx <- split2(permuIdx, Kfolds)

  # Function defining what to do for k-th fold partition
  CVjob <- function(k){

    valIdx <- splitIdx[[k]]
    trainIdx <- setdiff(allIdx, valIdx)

    errs <- getErr(
      ncompGrid = ncompGrid,
      selectedFeat = selectedFeat,
      XXtrain = XX[trainIdx, ], YYtrain = YY[trainIdx, ],
      XXtest = XX[valIdx, ], YYtest = YY[valIdx, ],
      regMethod = regMethod,
      assayName = assayName,
      center = center,
      scale = scale
    )

    return(errs)
  }

  # Parallel computing CV
  if(Ncores > 1){
    Ncores <- min(Ncores, parallel::detectCores()-1)
    err_re <- parallel::mclapply(1:Kfolds, FUN = CVjob, mc.cores = Ncores)
  } else {
    err_re <- lapply(1:Kfolds, FUN = CVjob)
  }
  errs <- Reduce(`+`, err_re)/length(err_re)


  ### Plotting --------------------------------------------------------

  # Useful for multiple choices of errs
  # errs <- errs %>%
  #   t() %>%
  #   as.data.frame() %>%
  #   dplyr::mutate(ncomp = ncompGrid)

  plot_dat <- data.frame(
    FrobRatio = errs,
    ncomp = ncompGrid
  )

  p <- plot_dat %>%
    ggplot(mapping = aes(x = ncomp)) +
    scale_x_continuous(breaks = plot_dat$ncomp, guide = guide_axis(angle = 90, n.dodge=2)) +
    geom_line(mapping = aes(y = FrobRatio)) +
    geom_point(mapping = aes(y = FrobRatio)) +
    ylab("Frobenius ratio") +
    ggtitle("Frobenius ratio")

  print(p)

  # TODO: use gap statistic for data-driven choice of ncomp
  ncomp <- readline("What is your choice of ncomp? ")
  ncomp <- as.numeric(ncomp)

  return(
    list(
      tuneNcompErrs = plot_dat,
      regMethod = regMethod,
      ncompGrid = ncompGrid,
      selectedFeat = selectedFeat,
      assayName = assayName,
      ncomp = ncomp,
      center = center,
      scale = scale
    )
  )





}
