#' An example data set with added DE results
#'
#' An example data set with added DE results. Data is from the airway package,
#' DE results are generated with `DESeq2`.
#' @details
#' This is a `DeeDeeExperiment` object, a derivative of the airway data set.
#' In the DEA slot are two separate differential expression analyses, both
#' generated with `DESeq2`.
#'
#' @return A DeeDeeExperiment with two DEA analyses
#'
#' @format A DeeDeeExperiment
#'
#' @name se_with_de
#' @docType data
NULL

#' An example data set with added FE results
#'
#' An example data set with added FE results. Data is from the airway package,
#' DE results are generated with `clusterProfiler`.
#' @details
#' This is a `DeeDeeExperiment` object, a derivative of the airway data set.
#' In the DEA slot are two separate functional enrichment analyses, both
#' generated with `clusterProfiler`.
#'
#' @return A DeeDeeExperiment with two FEA analyses
#'
#' @format A DeeDeeExperiment
#'
#' @name se_with_go
#' @docType data
NULL
