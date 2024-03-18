#' Tune nfeat using cross-validation.
#'
#' @param reference `SingleCellExperiment` object.
#' @param ncomp Integer.
#' @param impScores Matrix.
#' @param YY Matrix.
#' @param assayName Character.
#' @param nfeatLimits Vector.
#' @param gridSize Integer.
#' @param nfeatV Vector.
#' @param Kfolds Integer.
#' @param seed Integer.
#' @param regMethod Character.
#' @param Ncores Integer.
#' @param labelCode Character.
#' @param phenotypes Vector.
#' @param center Logic.
#' @param scale Logic.
#'
#' @return A list.
#'
#' @import ggplot2
CVTune_nfeat <- function(reference,
                         assayName = "logcounts",
                         phenotypes = NULL,
                         YY = NULL,
                         ncomp = NULL,
                         impScores = NULL,
                         #
                         center = TRUE,
                         scale = FALSE,
                         #
                         nfeatLimits = c(10, 15000),
                         gridSize = 30,
                         nfeatV = NULL,
                         #
                         Kfolds = 5,
                         regMethod = 'PCA',
                         labelCode = "-1,1",
                         seed = 5202056,
                         Ncores = 1)
{
  ## Prepare response: dummy matrix
  if(is.null(YY)){
    if(is.null(phenotypes)) stop("Have to specify either YY or phenotypes.")
    YY <- codeY(reference, phenotypes)
  }

  if(is.null(ncomp)) ncomp <- ncol(YY)

  ## Prepare predictor
  XX <- as.matrix(t(assay(reference, assayName)))


  ## impScores
  if(is.null(impScores)){
    impScores <- mvr(XX,
                     YY,
                     ncomp,
                     method = regMethod,
                     center = center,
                     scale = scale)$coefficients[,,ncomp]
  } else {
    impScores <- as.matrix(impScores)
  }
  if(nrow(impScores) != nrow(reference)) stop("Incorrect dimension of impScores.")

  ## Prepare nfeatV
  if(is.null(nfeatV)){

    nfeatMin <- min(nfeatLimits)
    if(nfeatMin > ncol(XX)) stop("nfeat has to be smaller than the total number of features.")
    nfeatMax <- min(max(nfeatLimits), ncol(XX))


    # If min nfeat results in too few features, fewer than ncomp
    actualnfeatMin <- length(
      unique(
        as.vector(
          selectFeat(impScores, nfeatMin)
        )
      )
    )
    if(ncomp > actualnfeatMin){
      nfeatMin <- ncomp
    }

    nfeatV <-
      round(
        exp(
          seq(log(nfeatMin), log(nfeatMax), length.out = gridSize)
        )
      )
    nfeatV <- unique(nfeatV)

  }




  ## Prepare response: dummy matrix
  if(is.null(YY)){
    if(is.null(phenotypes)) stop("Have to specify either phenotypes or YY.")
    YY <- codeY(reference, phenotypes)
  }

  ## CV
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

    errs <- getErr_nfeat(
      nfeatV = nfeatV,
      impScores = impScores,
      ncomp = ncomp,
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





  ### Plotting ----------------------------------------

  # errs <- errs %>%
  #   t() %>%
  #   as.data.frame() %>%
  #   dplyr::mutate(nfeat = nfeatV)

  plot_dat <- data.frame(
    FrobRatio = errs,
    nfeat = nfeatV
  )

  p <-
    ggplot(data = plot_dat, mapping = aes(x = nfeat)) +
    scale_x_continuous(
      breaks = plot_dat$nfeat,
      trans = "log",
      guide = guide_axis(angle = 90, n.dodge=2)
    ) +
    geom_line(mapping = aes(y = plot_dat$FrobRatio)) +
    geom_point(mapping = aes(y = plot_dat$FrobRatio)) +
    ylab("Loss")

  print(p)

  ## TODO: gap statistic for data-driven selection of nfeat
  nfeat <- readline("What is your choice of nfeat? ")
  nfeat <- as.numeric(nfeat)
  selectedFeat <- selectFeat(impScores, nfeat)

  return(
    list(
      ncomp = ncomp,
      impScores = impScores,
      tuneNfeatErrs = plot_dat,
      selectedFeat = selectedFeat,
      nfeat = nfeat,
      nfeatV = nfeatV,
      regMethod = regMethod,
      assayName = assayName,
      center = center,
      scale = scale
    )
  )
}


