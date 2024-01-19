
# PhiSpaceR <- function(references,
#                       query,
#                       phenotypes,
#                       PhiSpaceAssay = "rank",
#                       regMethod = "PLS",
#                       # arguments for CVTune_ncomp
#                       ncompLimits = c(1,100),
#                       gridSize = 20,
#                       Kfolds = 5,
#                       CVseed = 5202056,
#                       TuneCVmode = 'unsupervised',
#                       normYY = F,
#                       Ncores = 10
#                       )
# {
#   ## KeepCommonGenes happen internally
#   query_anns <-
#     lapply(references,
#            PhiSpaceR_1ref,
#            query = query,
#            phenotypes = phenotypes,
#            PhiSpaceAssay = PhiSpaceAssay,
#            regMethod = regMethod,
#            ncompLimits = ncompLimits,
#            gridSize = gridSize,
#            Kfolds = Kfolds,
#            CVseed = CVseed,
#            TuneCVmode = TuneCVmode,
#            normYY = normYY,
#            Ncores = Ncores)
# }
