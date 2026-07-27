## Tests for the exported orchestration functions get.gos()/get.gsea().
##
## Both functions delegate the actual enrichment computation to
## clusterProfiler (enrichGO()/GSEA()), external annotation databases and
## DeeDeeExperiment (addFEA()/renameFEA()). Rather than depending on real
## biological results (slow, non-deterministic, and would really be testing
## clusterProfiler rather than this package), we mock those calls and assert
## on the package's own logic: which genes get classified as up/down/universe,
## how the GSEA ranking statistic is computed, and how results get attached
## to the returned object.
##
## Where possible, the input object is the package's own bundled real-world
## DeeDeeExperiment (se_with_de) rather than a synthetic fixture - the
## expected gene sets/ranking direction below are computed independently
## from that real data rather than hand-picked. Two edge cases aren't
## naturally present in se_with_de (a comparison with zero significant genes
## in one direction; a colData column literally named "condition", to
## exercise get.gsea()'s default condition_var), so those two tests still
## use small hand-built fixtures.

se_de <- se_with_de

make_de_se <- function(res_df) {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      counts = matrix(
        1,
        nrow = nrow(res_df),
        ncol = 1,
        dimnames = list(rownames(res_df), "s1")
      )
    )
  )
  S4Vectors::metadata(se)$de_results <- list(comparisonA = res_df)
  se
}

test_that("get.gos classifies up/down/universe genes and forwards them to enrichGO()", {
  res <- .de_result(se_de, "baseline")
  up_expected <- rownames(res)[which(res$padj < 0.05 & res$log2FoldChange > 1)]
  dn_expected <- rownames(res)[which(res$padj < 0.05 & res$log2FoldChange < -1)]
  universe_expected <- rownames(res)[which(!is.na(res$padj))]

  enrichGO_calls <- list()
  testthat::local_mocked_bindings(
    enrichGO = function(gene, OrgDb, keyType, ont, universe) {
      enrichGO_calls[[length(enrichGO_calls) + 1]] <<- list(
        gene = gene,
        keyType = keyType,
        ont = ont,
        universe = universe
      )
      if (length(gene) == 0) {
        return(NULL)
      }
      data.frame(ID = paste0("GO:", seq_along(gene)), Description = "term")
    },
    .package = "clusterProfiler"
  )

  addFEA_calls <- list()
  renameFEA_calls <- list()
  testthat::local_mocked_bindings(
    addFEA = function(x, y, name) {
      addFEA_calls[[length(addFEA_calls) + 1]] <<- list(name = name)
      x
    },
    renameFEA = function(x, old_name, new_name) {
      renameFEA_calls[[length(renameFEA_calls) + 1]] <<- list(
        old_name = old_name,
        new_name = new_name
      )
      x
    },
    .package = "DeeDeeExperiment"
  )

  result <- get.gos(
    "baseline",
    obj = se_de,
    species = "hs",
    gene_type = "ENSEMBL"
  )

  expect_true(methods::is(result, "DeeDeeExperiment"))

  # up/down classification is based on padj < 0.05 & |log2FC| > 1, NAs excluded
  up_calls <- Filter(function(x) setequal(x$gene, up_expected), enrichGO_calls)
  dn_calls <- Filter(function(x) setequal(x$gene, dn_expected), enrichGO_calls)
  expect_true(length(up_calls) > 0)
  expect_true(length(dn_calls) > 0)

  # universe excludes genes with NA padj
  expect_true(all(vapply(
    enrichGO_calls,
    function(x) setequal(x$universe, universe_expected),
    logical(1)
  )))

  expect_true(all(vapply(
    enrichGO_calls,
    function(x) x$keyType == "ENSEMBL",
    logical(1)
  )))
  expect_true(all(vapply(
    enrichGO_calls,
    function(x) x$ont == "BP",
    logical(1)
  )))

  expect_equal(length(addFEA_calls), 2)
  expect_true(all(vapply(
    addFEA_calls,
    function(x) x$name == "baseline",
    logical(1)
  )))

  expect_equal(
    vapply(renameFEA_calls, function(x) x$old_name, character(1)),
    c("up_go", "dn_go")
  )
  expect_equal(
    vapply(renameFEA_calls, function(x) x$new_name, character(1)),
    c("baseline_up_go", "baseline_down_go")
  )
})

test_that("get.gos skips addFEA/renameFEA when a direction has no significant genes", {
  # se_de has significant genes in both directions for every comparison, so
  # this edge case needs a hand-built comparison with only one direction hit.
  res_df <- data.frame(
    padj = c(0.001, 0.9),
    log2FoldChange = c(3, 0.1),
    row.names = c("UP1", "NS1")
  )
  se <- make_de_se(res_df)

  testthat::local_mocked_bindings(
    enrichGO = function(gene, OrgDb, keyType, ont, universe) {
      if (length(gene) == 0) {
        return(NULL)
      }
      data.frame(ID = paste0("GO:", seq_along(gene)))
    },
    .package = "clusterProfiler"
  )

  addFEA_calls <- list()
  testthat::local_mocked_bindings(
    addFEA = function(x, y, name) {
      addFEA_calls[[length(addFEA_calls) + 1]] <<- list(name = name)
      x
    },
    renameFEA = function(x, old_name, new_name) x,
    .package = "DeeDeeExperiment"
  )

  get.gos("comparisonA", obj = se, species = "hs", gene_type = "SYMBOL")

  # only the "up" direction has significant genes, so addFEA should be called once
  expect_equal(length(addFEA_calls), 1)
})

test_that("get.gsea builds a descending gene ranking and forwards it to GSEA()", {
  msigdbr_calls <- list()
  testthat::local_mocked_bindings(
    msigdbr = function(collection, db_species, ...) {
      msigdbr_calls[[length(msigdbr_calls) + 1]] <<- list(
        collection = collection,
        db_species = db_species
      )
      data.frame(gs_name = "SET_A", gene_symbol = "GENE1")
    },
    .package = "msigdbr"
  )

  gsea_calls <- list()
  testthat::local_mocked_bindings(
    GSEA = function(geneList, TERM2GENE, ...) {
      gsea_calls[[length(gsea_calls) + 1]] <<- list(
        geneList = geneList,
        TERM2GENE = TERM2GENE
      )
      structure(list(), class = "mock_gsea_result")
    },
    .package = "clusterProfiler"
  )

  addFEA_calls <- list()
  renameFEA_calls <- list()
  testthat::local_mocked_bindings(
    addFEA = function(x, y, name) {
      addFEA_calls[[length(addFEA_calls) + 1]] <<- list(name = name)
      x
    },
    renameFEA = function(x, old_name, new_name) {
      renameFEA_calls[[length(renameFEA_calls) + 1]] <<- list(
        old_name = old_name,
        new_name = new_name
      )
      x
    },
    .package = "DeeDeeExperiment"
  )

  result <- get.gsea(
    NAME = "baseline",
    obj = se_de,
    type = "HALLMARK",
    conditions = c("trt", "untrt"),
    species = "hs",
    condition_var = "dex"
  )

  expect_true(methods::is(result, "DeeDeeExperiment"))

  expect_equal(length(msigdbr_calls), 1)
  expect_equal(msigdbr_calls[[1]]$collection, "H")
  expect_equal(msigdbr_calls[[1]]$db_species, "HS")

  expect_equal(length(gsea_calls), 1)
  rankings <- gsea_calls[[1]]$geneList

  expect_true(is.numeric(rankings))
  expect_true(!is.null(names(rankings)))
  # results are returned sorted from most to least "up"
  expect_true(all(diff(rankings) <= 0))

  # Cross-check the ranking direction against a plain, independently
  # computed mean-expression comparison (not the package's sd-scaled
  # formula) for the most extreme gene at each end of the ranking.
  counts_norm <- BiocGenerics::counts(se_de, normalized = TRUE)
  dex <- SummarizedExperiment::colData(se_de)$dex
  top_gene <- names(rankings)[1]
  bottom_gene <- names(rankings)[length(rankings)]
  expect_gt(
    mean(counts_norm[top_gene, dex == "trt"]),
    mean(counts_norm[top_gene, dex == "untrt"])
  )
  expect_gt(
    mean(counts_norm[bottom_gene, dex == "untrt"]),
    mean(counts_norm[bottom_gene, dex == "trt"])
  )

  expect_equal(length(addFEA_calls), 1)
  expect_equal(addFEA_calls[[1]]$name, "baseline")
  expect_equal(renameFEA_calls[[1]]$old_name, "gsea")
  expect_equal(renameFEA_calls[[1]]$new_name, "baseline_gsea_HALLMARK")
})

test_that("get.gsea requests the Reactome collection for type = 'REACTOME'", {
  msigdbr_calls <- list()
  testthat::local_mocked_bindings(
    msigdbr = function(collection, subcollection = NULL, db_species, ...) {
      msigdbr_calls[[length(msigdbr_calls) + 1]] <<- list(
        collection = collection,
        subcollection = subcollection,
        db_species = db_species
      )
      data.frame(gs_name = "SET_A", gene_symbol = "GENE1")
    },
    .package = "msigdbr"
  )
  testthat::local_mocked_bindings(
    GSEA = function(geneList, TERM2GENE, ...) {
      structure(list(), class = "mock_gsea_result")
    },
    .package = "clusterProfiler"
  )
  testthat::local_mocked_bindings(
    addFEA = function(x, y, name) x,
    renameFEA = function(x, old_name, new_name) x,
    .package = "DeeDeeExperiment"
  )

  get.gsea(
    NAME = "baseline",
    obj = se_de,
    type = "REACTOME",
    conditions = c("trt", "untrt"),
    species = "mm",
    condition_var = "dex"
  )

  expect_equal(msigdbr_calls[[1]]$collection, "C2")
  expect_equal(msigdbr_calls[[1]]$subcollection, "CP:REACTOME")
  expect_equal(msigdbr_calls[[1]]$db_species, "MM")
})

test_that("get.gsea's default condition_var works without being supplied explicitly", {
  # se_de's grouping column is called "dex", not "condition" - this test is
  # specifically about the literal default value, so it needs an object with
  # a colData column named "condition".
  counts <- matrix(
    c(
      100,
      110,
      90,
      10,
      12,
      8,
      10,
      12,
      8,
      100,
      110,
      90,
      50,
      52,
      48,
      50,
      51,
      49
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- paste0("GENE", 1:3)
  colnames(counts) <- paste0("S", 1:6)
  coldata <- data.frame(
    sample = colnames(counts),
    condition = factor(rep(c("groupA", "groupB"), each = 3)),
    row.names = colnames(counts)
  )
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, coldata, design = ~condition)
  dds <- DESeq2::estimateSizeFactors(dds)

  testthat::local_mocked_bindings(
    msigdbr = function(...) {
      data.frame(gs_name = "SET_A", gene_symbol = "GENE1")
    },
    .package = "msigdbr"
  )
  gsea_calls <- list()
  testthat::local_mocked_bindings(
    GSEA = function(geneList, TERM2GENE, ...) {
      gsea_calls[[length(gsea_calls) + 1]] <<- list(geneList = geneList)
      structure(list(), class = "mock_gsea_result")
    },
    .package = "clusterProfiler"
  )
  testthat::local_mocked_bindings(
    addFEA = function(x, y, name) x,
    renameFEA = function(x, old_name, new_name) x,
    .package = "DeeDeeExperiment"
  )

  # condition_var intentionally omitted - relies on the "condition" default
  expect_no_error(
    get.gsea(
      NAME = "comparisonA",
      obj = dds,
      type = "HALLMARK",
      conditions = c("groupA", "groupB"),
      species = "hs"
    )
  )

  expect_equal(length(gsea_calls), 1)
  expect_true(is.numeric(gsea_calls[[1]]$geneList))
})

test_that("get.gsea accepts condition_var as either a bare symbol or a string", {
  testthat::local_mocked_bindings(
    msigdbr = function(...) {
      data.frame(gs_name = "SET_A", gene_symbol = "GENE1")
    },
    .package = "msigdbr"
  )
  testthat::local_mocked_bindings(
    addFEA = function(x, y, name) x,
    renameFEA = function(x, old_name, new_name) x,
    .package = "DeeDeeExperiment"
  )
  gsea_calls <- list()
  testthat::local_mocked_bindings(
    GSEA = function(geneList, TERM2GENE, ...) {
      gsea_calls[[length(gsea_calls) + 1]] <<- list(geneList = geneList)
      structure(list(), class = "mock_gsea_result")
    },
    .package = "clusterProfiler"
  )

  expect_no_error(
    get.gsea(
      NAME = "baseline",
      obj = se_de,
      type = "HALLMARK",
      conditions = c("trt", "untrt"),
      species = "hs",
      condition_var = dex # bare symbol, not a string
    )
  )
  expect_no_error(
    get.gsea(
      NAME = "baseline",
      obj = se_de,
      type = "HALLMARK",
      conditions = c("trt", "untrt"),
      species = "hs",
      condition_var = "dex"
    )
  )

  expect_equal(length(gsea_calls), 2)
  expect_equal(gsea_calls[[1]]$geneList, gsea_calls[[2]]$geneList)
})
