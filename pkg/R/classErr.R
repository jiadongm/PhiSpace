#' Calculate classification errors.
#'
#' @param classQuery Vector.
#' @param classOriginal Vector.
#' @param labPerSample Integer.
#'
#' @return A list.
classErr <- function(classQuery, classOriginal, labPerSample = NULL){

  if(is.null(labPerSample)) labPerSample <- 1

  if(labPerSample == 1){

    out_errs <- rep(NA, 2)

    # Overall error
    out_errs[1] <- sum(classQuery != classOriginal)/length(classOriginal)

    # Balanced error
    classL <- unique(classOriginal)
    errs <- sapply(1:length(classL),
                   function(x){
                     cl <- classL[x]
                     idx <- (classOriginal == cl)
                     # Per class classification error
                     sum(classQuery[idx] != classOriginal[idx])/sum(idx)
                   })
    names(errs) <- classL
    out_errs[2] <- mean(errs)

    return(list(err = out_errs, byClassErrs = sort(errs)))

  } else {

    if(labPerSample <= 0) stop("labPerSample has to be positive integer.")

    # Overall errors
    ove_errs <-
      sapply(1:nrow(classQuery),
             function(x){
               out <- sum(2 - sum(classQuery[x,] %in% classOriginal[x,]))/labPerSample
             })
    ove_errs <- sum(ove_errs)/length(ove_errs)

    # Balanced error
    out_list <- vector("list", labPerSample)
    for(erIdx in 1:ncol(classOriginal)){

      classL <- unique(classOriginal[,erIdx])
      errs <- sapply(1:length(classL),
                     function(x){
                       cl <- classL[x]
                       idx <- (classOriginal[,erIdx] == cl)
                       sum(classQuery[idx, erIdx] != classOriginal[idx, erIdx])/sum(idx)
                     })
      names(errs) <- classL
      out_list[[erIdx]] <- errs

    }

    byClassErrs <- do.call(`c`, out_list)
    BER <- mean(byClassErrs)

    return(
      list(
        err = c(ove_errs, BER),
        byClassErrs = byClassErrs
      )
    )

  }




}
