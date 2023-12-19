#' Tune ncomp using cross-validation.
#'
#' @param reference
#' @param YY
#' @param labPerSample
#' @param labelName
#' @param ncompLimits
#' @param gridSize
#' @param ncompGrid
#' @param selectFeat
#' @param assayName
#' @param Kfolds
#' @param seed
#' @param regMethod
#' @param mode
#' @param normYY
#' @param Ncores
#' @param labelCode
#'
#' @return
#'
#' @import ggplot2
CVTuneNcomp <- function(reference,
                        YY = NULL,
                        labPerSample = NULL,
                        labelName = NULL,
                        ncompLimits = c(1,100),
                        gridSize = 20,
                        ncompGrid = NULL,
                        selectFeat = NULL,
                        assayName = "logcounts",
                        Kfolds = 5,
                        seed = 5202056,
                        regMethod = 'PCA',
                        mode = 'supervised',
                        normYY = F,
                        Ncores = 2,
                        labelCode = "-1,1")
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
  if(mode == "supervised" & is.null(YY)){

    if(is.null(labelName)) stop("Have to specify either labelName or YY.")

    Ytrain <- colData(reference)[,labelName]
    classLabels <- names(table(Ytrain))
    YY <- sapply(1:length(Ytrain),
                 function(x){

                   if(labelCode == "-1,1"){

                     out <- as.numeric(classLabels == Ytrain[x])
                     out[out == 0] <- -1

                   } else {

                     out <- as.numeric(classLabels == Ytrain[x])

                   }

                   return(out)
                 })
    YY <- t(YY)
    dimnames(YY) <- list(rownames(XX), classLabels)
  }

  if(mode != "supervised" & is.null(YY)){
    stop("YY cannot be NULL for unsupervised mode.")
  }

  ## By default each sample has one label, eg cell type
  if(is.null(labPerSample)){
    labPerSample <- 1
  }


  if(length(ncompGrid) > 1){
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

      errs <- getErr(
        ncompGrid = ncompGrid,
        selectFeat = selectFeat,
        XX = XX[trainIdx, ], YY = YY[trainIdx, ],
        XXtest = XX[valIdx, ], YYtest = YY[valIdx, ],
        labPerSample = labPerSample,
        regMethod = regMethod,
        assayName = assayName,
        normYY = normYY,
        mode = mode
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
    ## Prepare data.frame for visualisation
    if(mode == "supervised"){
      errs <- errs %>%
        t() %>%
        as.data.frame() %>%
        dplyr::mutate(ncomp = ncompGrid)

      p <- errs %>%
        ggplot(mapping = aes(x = ncomp)) +
        scale_x_continuous(breaks = errs$ncomp, guide = guide_axis(angle = 90, n.dodge=2))

      ## Errors using normalised scores
      rescale1 <- diff(range(errs$overall))/diff(range(errs$balanced))
      intercep1 <- min(errs$balanced) - min(errs$overall) / rescale1
      p1 <-
        p +
        geom_line(mapping = aes(y = overall)) +
        geom_point(mapping = aes(y = overall)) +
        geom_line(mapping = aes(y = (balanced - intercep1) * rescale1), color = "blue") +
        geom_point(
          aes(y = (balanced - intercep1) * rescale1),
          color = "blue"
        ) +
        ggrepel::geom_text_repel(
          aes(
            y = (balanced - intercep1) * rescale1,
            label = ncomp
          ),
          color = "blue"
        ) +
        scale_y_continuous(
          "Overall",
          sec.axis = sec_axis(~./rescale1 + intercep1, name = "Balanced"),
        )  +
        theme(axis.line.y.right = element_line(color = "blue"),
              axis.ticks.y.right = element_line(color = "blue"),
              axis.text.y.right = element_text(color = "blue"),
              axis.title.y.right = element_text(color = "blue")
        ) +
        ggtitle("Classif errors - Original")


      ## Error using non-normalised scores
      rescale2 <- diff(range(errs$overall_norm))/diff(range(errs$balanced_norm))
      intercep2 <- min(errs$balanced_norm) - min(errs$overall_norm) / rescale2
      p2 <-
        p +
        geom_line(mapping = aes(y = overall_norm)) +
        geom_point(mapping = aes(y = overall_norm)) +
        geom_line(mapping = aes(y = (balanced_norm - intercep2) * rescale2), color = "blue") +
        geom_point(
          aes(y = (balanced_norm - intercep2) * rescale2),
          color = "blue"
        ) +
        ggrepel::geom_text_repel(
          aes(
            y = (balanced_norm - intercep2) * rescale2,
            label = ncomp
          ),
          color = "blue"
        ) +
        scale_y_continuous(
          "Overall - normalised",
          sec.axis = sec_axis(~./rescale2 + intercep2, name = "Balanced - normalised")
        )  +
        theme(axis.line.y.right = element_line(color = "blue"),
              axis.ticks.y.right = element_line(color = "blue"),
              axis.text.y.right = element_text(color = "blue"),
              axis.title.y.right = element_text(color = "blue")
        ) +
        ggtitle("Classif errors - Normalised")

      p3 <-
        p +
        geom_line(mapping = aes(y = `Frob ratio`)) +
        geom_point(mapping = aes(y = `Frob ratio`)) +
        ylab("Scaled Frobenius") +
        ggtitle("Scaled Frobenius error")

      p4 <-
        p +
        geom_line(mapping = aes(y = `inf norm`)) +
        geom_point(mapping = aes(y = `inf norm`)) +
        ylab("Inf norm") +
        ggtitle("Inf norm")

      print(p1 + p2 + p3 + p4 + patchwork::plot_layout(ncol = 2))
    } else {
      errs <- errs %>%
        t() %>%
        as.data.frame() %>%
        dplyr::mutate(ncomp = ncompGrid)

      p <- errs %>%
        ggplot(mapping = aes(x = ncomp)) +
        scale_x_continuous(breaks = errs$ncomp, guide = guide_axis(angle = 90, n.dodge=2))

      p1 <-
        p +
        geom_line(mapping = aes(y = `Frob ratio`)) +
        geom_point(mapping = aes(y = `Frob ratio`)) +
        ylab("Frobenius ratio") +
        ggtitle("Frobenius ratio")

      p2 <-
        p +
        geom_line(mapping = aes(y = `inf norm`)) +
        geom_point(mapping = aes(y = `inf norm`)) +
        ylab("Inf norm") +
        ggtitle("Inf norm")

      print(p1 + p2 + plot_layout(ncol = 2))
    }

    ncomp <- readline("What is your choice of ncomp? ")
    ncomp <- as.numeric(ncomp)

    return(
      list(
        XX = XX,
        YY = YY,
        tuneNcompErrs = errs,
        regMethod = regMethod,
        ncompGrid = ncompGrid,
        selectFeat = selectFeat,
        assayName = assayName,
        mode = mode,
        normYY = normYY,
        normType = normType,
        ncomp = ncomp
      )
    )

  } else {

    ## when ncompGrid has only 1 element, simply calcualte impScores
    ncomp <- ncompGrid[1]


    return(
      list(
        XX = XX,
        YY = YY,
        regMethod = regMethod,
        selectFeat = selectFeat,
        assayName = assayName,
        mode = mode,
        normYY = normYY,
        ncomp = ncomp
      )
    )

  }


}
