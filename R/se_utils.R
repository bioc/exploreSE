#' Calculate GO Enrichments for a given comparison
#'
#' A function to calculate GO term enrichments from a given summarizedExperiment
#' or DeeDee experiment. By default, it calculates separate enrichments for
#' up- downregulated genes, at the cutoffs of |logFC| > 1 and adjusted p-
#' value < 0.05. It utilises the [clusterProfiler::enrichGO()] function.
#'
#' @param NAME character, name of the DE comparison
#' @param obj DeeDeeExperiment or SummarizedExperiment, the object
#' @param species character, the species hs/mm
#' @param gene_type character, the type of gene identifier
#'
#' @returns a DeeDeeExperiment or SummarizedExperiment
#' @export
#'
#' @examples
#' data("se_with_des", package = "exploreSE")
#' se_with_go <- get.gos(obj = se_with_de,
#'                       NAME = "baseline",
#'                       gene_type = "ENSEMBL")
#'

get.gos <- function(NAME, obj, species = "hs", gene_type = "SYMBOL") {
  if (species == "hs") {
    DB <- org.Hs.eg.db::org.Hs.eg.db
  } else if (species == "mm") {
    DB <- org.Mm.eg.db::org.Mm.eg.db
  } else {
    stop("Please specify a valid species. Supported: hs, mm")
  }

  if (gene_type %notin% c("SYMBOL", "ENTREZID", "ENSEMBL")) {
    stop(
      "Please specify a valid gene identifier. Supported: SYMBOL, ENSEMBL, ENTREZID"
    )
  }

  res <- .de_result(obj, NAME)

  up_genes <- res %>%
    BiocGenerics::as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(padj < 0.05, log2FoldChange > 1) %>%
    dplyr::pull(gene)
  dn_genes <- res %>%
    BiocGenerics::as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(padj < 0.05, log2FoldChange < -1) %>%
    dplyr::pull(gene)

  universe <- res %>%
    BiocGenerics::as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::pull(gene)

  up_go <- clusterProfiler::enrichGO(
    up_genes,
    DB,
    keyType = gene_type,
    ont = "BP",
    universe = universe
  ) %>%
    BiocGenerics::as.data.frame()
  dn_go <- clusterProfiler::enrichGO(
    dn_genes,
    DB,
    keyType = gene_type,
    ont = "BP",
    universe = universe
  ) %>%
    BiocGenerics::as.data.frame()

  if (!methods::is(obj, "DeeDeeExperiment")) {
    obj <- DeeDeeExperiment::DeeDeeExperiment(obj)
  }

  if (length(up_genes) > 0) {
    up_go <- clusterProfiler::enrichGO(
      up_genes,
      DB,
      keyType = gene_type,
      ont = "BP",
      universe = universe
    )
    if (!is.null(up_go)) {
      obj <- DeeDeeExperiment::addFEA(obj, up_go, NAME)
      obj <- DeeDeeExperiment::renameFEA(obj, "up_go", paste0(NAME, "_up_go"))
    }
  }

  if (length(dn_genes) > 0) {
    dn_go <- clusterProfiler::enrichGO(
      dn_genes,
      DB,
      keyType = gene_type,
      ont = "BP",
      universe = universe
    )
    if (!is.null(dn_go)) {
      obj <- DeeDeeExperiment::addFEA(obj, dn_go, NAME)
      obj <- DeeDeeExperiment::renameFEA(obj, "dn_go", paste0(NAME, "_down_go"))
    }
  }
  return(obj)
}


#' Calculate GSEA Enrichments for a given comparison
#'
#' Calculates the GSEA enrichments with the counts within a summarizedExperiment
#' or DeeDeeExperiment and adds them in the respective slots. It follows
#' the [original signal-to-noise ratio of the method](https://www.pnas.org/doi/10.1073/pnas.0506580102); other approaches for
#' ranking of genes are available in other implementations,
#' such as the [clusterProfiler::GSEA()] implementation which calls for
#' logFC ranked genes.
#'
#' @param NAME character, name of the DE comparison
#' @param obj DeeDeeExperiment or SummarizedExperiment, the object
#' @param type character, the type of gene set, HALLMARK or REACTOME
#' @param conditions a vector of length 2, which levels of the condition variable are being compared,
#' @param species character, supported values are "hs" and "mm"
#' @param condition_var symbol or character, the variable name for this comparison
#' @param gene_type character, the type of gene identifier
#'
#' @returns a DeeDeeExperiment or SummarizedExperiment
#' @export
#'
#' @examples
#' data("se_with_des", package = "exploreSE")
#' se_with_gsea <- get.gsea(NAME = "baseline",
#'                          obj = se_with_de,
#'                          type = "HALLMARK",
#'                          condition_var = "dex",
#'                          conditions = c("trt", "untrt"),
#'                          gene_type = "ENSEMBL")

get.gsea <- function(
  NAME,
  obj,
  type = "HALLMARK",
  conditions,
  species = "hs",
  condition_var = "condition",
  gene_type = "SYMBOL"
) {
  condition_var <- rlang::ensym(condition_var)

  if (species == "hs") {
    SPECIES <- "HS"
  } else if (species == "mm") {
    SPECIES <- "MM"
  } else {
    stop("Please specify a valid species. Supported: hs, mm")
  }

  if (type == "HALLMARK") {
    gene_sets <- msigdbr::msigdbr(collection = "H", db_species = SPECIES)
  } else if (type == "REACTOME") {
    gene_sets <- msigdbr::msigdbr(
      collection = "C2",
      subcollection = "CP:REACTOME",
      db_species = SPECIES
    )
  } else {
    stop(
      "Please specify a valid gene set identifier. Supported: HALLMARK, REACTOME"
    )
  }

  if (gene_type == "SYMBOL") {
    gene_sets <- gene_sets %>%
      dplyr::select(gs_name, gene_symbol) %>%
      dplyr::distinct()
  } else if (gene_type == "ENSEMBL") {
    gene_sets <- gene_sets %>%
      dplyr::select(gs_name, ensembl_gene) %>%
      dplyr::distinct()
  } else if (type == "ENTREZID") {
    gene_sets <- gene_sets %>%
      dplyr::select(gs_name, ncbi_gene) %>%
      dplyr::distinct()
  } else {
    stop(
      "Please specify a valid gene identifier. Supported: SYMBOL, ENSEMBL, ENTREZID"
    )
  }

  cdata <- BiocGenerics::as.data.frame(SummarizedExperiment::colData(obj)) |>
    tibble::rownames_to_column("sample_ident")

  rankings <- BiocGenerics::counts(obj, normalized = TRUE) %>%
    BiocGenerics::as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene) %>%
    dplyr::left_join(
      cdata,
      dplyr::join_by(name == sample_ident)
    ) %>%
    dplyr::mutate(condition = !!condition_var) %>%
    dplyr::filter(condition %in% conditions) %>%
    dplyr::select(gene, condition, value) %>%
    dplyr::group_by(condition, gene) %>%
    dplyr::summarise(
      sd = stats::sd(value, na.rm = TRUE),
      mean_expr = mean(value, na.rm = TRUE)
    ) %>%
    dplyr::mutate(
      sd = ifelse(sd >= 0.2 * abs(mean_expr), sd, 0.2 * abs(mean_expr)),
      mean_expr = ifelse(mean_expr == 0, 1, mean_expr)
    ) %>%
    dplyr::mutate(
      condition = ifelse(condition == conditions[2], "group2", "group1")
    ) %>%
    dplyr::select(condition, mean_expr:sd, gene) %>%
    tidyr::pivot_wider(names_from = condition, values_from = mean_expr:sd) %>%
    dplyr::filter(
      dplyr::if_all(sd_group1:sd_group2, \(x) x > 0),
      dplyr::if_any(mean_expr_group1:mean_expr_group2, \(x) x > 0)
    ) %>%
    dplyr::mutate(
      ranking = (mean_expr_group1 - mean_expr_group2) / (sd_group1 + sd_group2),
      gene = forcats::fct_reorder(gene, ranking, .desc = TRUE)
    ) %>%
    dplyr::filter(!is.na(ranking)) %>%
    dplyr::arrange(dplyr::desc(ranking)) %>%
    dplyr::pull(ranking, name = gene)

  gsea <- clusterProfiler::GSEA(rankings, TERM2GENE = gene_sets)

  if (!methods::is(obj, "DeeDeeExperiment")) {
    obj <- DeeDeeExperiment::DeeDeeExperiment(obj)
  }

  obj <- DeeDeeExperiment::addFEA(obj, gsea, NAME)
  obj <- DeeDeeExperiment::renameFEA(obj, "gsea", paste0(NAME, "_gsea_", type))

  return(obj)
}
