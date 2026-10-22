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
#' @name baseline_down_go
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
#' @name controlled_down_go
#' @docType data
NULL

#' Example context lists
#'
#' Example gene lists to use as contexts in the Context tab of [exploreSE()].
#'
#' @details
#' A named list of two gene sets, as gene symbols:
#' \describe{
#'   \item{dex_response}{11 well-known glucocorticoid response genes,
#'   compiled from studies in human dermal fibroblasts, human airway smooth
#'   muscle cells (the `airway` data set itself) and mouse keratinocytes.}
#'   \item{aortic_valve_stenosis}{261 genetic risk factors for aortic valve
#'   stenosis. This set is a placeholder and will be revised.}
#' }
#' The same contexts are available as a file with one context per column,
#' which can be used with the "Add context" upload:
#' `system.file("extdata", "example_contexts.csv", package = "exploreSE")`.
#' Both are created by `inst/scripts/creating_example_contexts.R`.
#'
#' @source
#' dex_response: \doi{10.3390/biom11050619}, \doi{10.1371/journal.pone.0099625}
#' and \doi{10.1016/j.mce.2015.05.015}.
#'
#' aortic_valve_stenosis: \doi{10.1038/s41588-025-02417-6}.
#'
#' @format A named list of two character vectors.
#'
#' @examples
#' data("example_contexts", package = "exploreSE")
#' lengths(example_contexts)
#'
#' @name example_contexts
#' @docType data
NULL
