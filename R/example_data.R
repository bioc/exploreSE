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
#' A named list of 11 gene sets, as gene symbols:
#' \describe{
#'   \item{dex_response}{11 well-known glucocorticoid response genes,
#'   compiled from studies in human dermal fibroblasts, human airway smooth
#'   muscle cells (the `airway` data set itself) and mouse keratinocytes.}
#'   \item{aortic_valve_stenosis_high}{89 high priority genetic risk factors
#'   for aortic valve stenosis.}
#'   \item{aortic_valve_stenosis_high_medium}{261 high and medium priority
#'   genetic risk factors for aortic valve stenosis.}
#'   \item{bicuspid_aortic_valve}{46 genetic risk factors for bicuspid aortic
#'   valve.}
#'   \item{thoracic_aortic_aneurysm}{100 genetic risk factors for thoracic
#'   aortic aneurysm.}
#'   \item{abdominal_aortic_aneurysm}{130 genetic risk factors for abdominal
#'   aortic aneurysm.}
#'   \item{ascending_aorta_diameter}{82 genetic risk factors for an increased
#'   ascending aorta diameter.}
#'   \item{descending_aorta_diameter}{47 genetic risk factors for an increased
#'   descending aorta diameter.}
#'   \item{coronary_artery_disease_high}{123 high priority genetic risk
#'   factors for coronary artery disease.}
#'   \item{coronary_artery_disease_high_medium}{220 high and medium priority
#'   genetic risk factors for coronary artery disease.}
#'   \item{monogenic_aortic_disease}{44 genes associated with
#'   familial/monogenic aortic diseases.}
#' }
#' The same contexts are available as a file with one context per column,
#' which can be used with the "Add context" upload:
#' `system.file("extdata", "example_contexts.csv", package = "exploreSE")`.
#' Both are created by `inst/scripts/creating_example_contexts.R`.
#'
#' @source
#' \itemize{
#'   \item dex_response: \doi{10.3390/biom11050619},
#'   \doi{10.1371/journal.pone.0099625} and \doi{10.1016/j.mce.2015.05.015}.
#'   \item aortic_valve_stenosis_high, aortic_valve_stenosis_high_medium:
#'   Small et al. (2025), \doi{10.1038/s41588-025-02417-6}.
#'   \item bicuspid_aortic_valve: \enc{Thériault}{Theriault} et al. (2026),
#'   \doi{10.1161/CIRCULATIONAHA.125.074752}.
#'   \item thoracic_aortic_aneurysm: Klarin et al. (2023),
#'   \doi{10.1038/s41588-023-01420-z}, and Murdock et al. (2026, preprint),
#'   \doi{10.64898/2026.09.16.26363254}.
#'   \item abdominal_aortic_aneurysm: Roychowdhury et al. (2023),
#'   \doi{10.1038/s41588-023-01510-y}.
#'   \item ascending_aorta_diameter, descending_aorta_diameter: Pirruccello
#'   et al. (2022), \doi{10.1038/s41588-021-00962-4}.
#'   \item coronary_artery_disease_high, coronary_artery_disease_high_medium:
#'   Aragam et al. (2022), \doi{10.1038/s41588-022-01233-6}.
#'   \item monogenic_aortic_disease: source to be added.
#' }
#'
#' @format A named list of 11 character vectors.
#'
#' @examples
#' data("example_contexts", package = "exploreSE")
#' lengths(example_contexts)
#'
#' @name example_contexts
#' @docType data
NULL
