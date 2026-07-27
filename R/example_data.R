#' DE results for airway - baseline
#'
#' Example DE results for airway data set
#'
#' @details
#' This is a `DESEqResults` object, with the simplest model, using `~dex` as
#' the design. Data is from the airway package, DE results are generated with `DESeq2`.
#'
#' @return A DESEqResults
#'
#' @format A DESEqResults
#'
#' @name baseline
#' @docType data
NULL

#' DE results for airway - cell-controlled
#'
#' Example DE results for airway data set
#'
#' @details
#' This is a `DESEqResults` object, using `~ cell + dex` as the design.
#' Data is from the airway package, DE results are generated with `DESeq2`,
#' using the `dex untrt vs trt` contrast
#'
#' @return A DESEqResults
#'
#' @format A DESEqResults
#'
#' @name cell_controlled
#' @docType data
NULL


#' Example GO enrichments for the baseline model
#'
#' GO enrichments based on the `baseline` DE results
#' @details
#'
#' GO enrichmens are generated with `clusterProfiler`, filtered on up-regulated
#' genes only; genes that underwent DE testing were taken as the universe.
#'
#'
#' The code to create said object can be found in the folder `/inst/scripts` in
#' the exploreSE package, the file is called `creating_dataset_examples.R`.
#'
#' @return A EnrichResults object
#'
#' @format A EnrichResults
#'
#' @name baseline_up_go
#' @docType data
NULL


#' Example GO enrichments for the baseline model
#'
#' GO enrichments based on the `baseline` DE results
#' @details
#'
#' GO enrichmens are generated with `clusterProfiler`, filtered on down-regulated
#' genes only; genes that underwent DE testing were taken as the universe.
#'
#' The code to create said object can be found in the folder `/inst/scripts` in
#' the exploreSE package, the file is called `creating_dataset_examples.R`.
#'
#' @return A EnrichResults object
#'
#' @format A EnrichResults
#'
#' @name baseline_dn_go
#' @docType data
NULL


#' Example GO enrichments for the cell controlled model
#'
#' GO enrichments based on the `cell_controlled` DE results
#' @details
#'
#' GO enrichmens are generated with `clusterProfiler`, filtered on up-regulated
#' genes only; genes that underwent DE testing were taken as the universe.
#'
#' The code to create said object can be found in the folder `/inst/scripts` in
#' the exploreSE package, the file is called `creating_dataset_examples.R`.
#'
#' @return A EnrichResults object
#'
#' @format A EnrichResults
#'
#' @name controlled_up_go
#' @docType data
NULL


#' Example GO enrichments for the cell_controlled model
#'
#' GO enrichments based on the `cell_controlled` DE results
#' @details
#'
#' GO enrichmens are generated with `clusterProfiler`, filtered on down-regulated
#' genes only; genes that underwent DE testing were taken as the universe.
#'
#' The code to create said object can be found in the folder `/inst/scripts` in
#' the exploreSE package, the file is called `creating_dataset_examples.R`.
#'
#' @return A EnrichResults object
#'
#' @format A EnrichResults
#'
#' @name controlled_dn_go
#' @docType data
NULL
