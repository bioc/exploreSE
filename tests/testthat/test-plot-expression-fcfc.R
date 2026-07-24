## Tests for .plot_expression and .plot_fcfc, using the package's bundled
## real-world DeeDeeExperiment (se_with_de) rather than synthetic fixtures.

data("se_with_des", package = "exploreSE")
se_de <- se_with_de

test_that(".plot_expression returns a ggplot titled with the gene name", {
  p <- .plot_expression(
    se_de,
    GENE = "TSPAN6",
    BY = "dex",
    LEVELS = c("trt", "untrt"),
    PLOT_TYPE = "box",
    GENE_VAR = "symbol"
  )

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Expression: TSPAN6")
  expect_equal(p$labels$x, "dex")
  expect_equal(p$labels$y, "Counts")
})

test_that(".plot_expression uses boxplot vs violin geoms depending on PLOT_TYPE", {
  p_box <- .plot_expression(
    se_de,
    GENE = "TSPAN6",
    BY = "dex",
    LEVELS = c("trt", "untrt"),
    PLOT_TYPE = "box",
    GENE_VAR = "symbol"
  )
  p_violin <- .plot_expression(
    se_de,
    GENE = "TSPAN6",
    BY = "dex",
    LEVELS = c("trt", "untrt"),
    PLOT_TYPE = "violin",
    GENE_VAR = "symbol"
  )

  box_geoms <- sapply(p_box$layers, function(l) class(l$geom)[1])
  violin_geoms <- sapply(p_violin$layers, function(l) class(l$geom)[1])

  expect_true("GeomBoxplot" %in% box_geoms)
  expect_true("GeomViolin" %in% violin_geoms)
})

test_that(".plot_expression normalises counts by sizeFactor when present", {
  gene <- rownames(SummarizedExperiment::rowData(se_de)[
    which(SummarizedExperiment::rowData(se_de)[, "symbol"] == "TSPAN6"),
  ])[1]

  raw_counts <- SummarizedExperiment::assay(se_de, "counts")[gene, ]
  sf <- SummarizedExperiment::colData(se_de)[, "sizeFactor"]
  expected <- unname(raw_counts / sf)

  p <- .plot_expression(
    se_de,
    GENE = "TSPAN6",
    BY = "dex",
    LEVELS = c("trt", "untrt"),
    PLOT_TYPE = "box",
    GENE_VAR = "symbol"
  )
  built <- ggplot2::ggplot_build(p)$data[[1]]
  plotted <- built$y[match(colnames(se_de), built$text)]

  expect_equal(plotted, unname(expected), tolerance = 1e-6)
})

test_that(".plot_expression filters samples down to the requested LEVELS", {
  p <- .plot_expression(
    se_de,
    GENE = "TSPAN6",
    BY = "dex",
    LEVELS = "trt",
    PLOT_TYPE = "box",
    GENE_VAR = "symbol"
  )
  built <- ggplot2::ggplot_build(p)$data[[1]]

  expect_equal(nrow(built), sum(SummarizedExperiment::colData(se_de)$dex == "trt"))
})

test_that(".plot_fcfc returns a plotly comparing two DE comparisons", {
  p <- .plot_fcfc(
    se_de,
    NAME = "baseline",
    PARTNER_NAME = "cell_controlled",
    GENE_VAR = "symbol",
    GENE_VARS = c("symbol")
  )

  expect_s3_class(p, "plotly")
  expect_equal(p$x$layout$xaxis$title$text, "logFC baseline")
  expect_equal(p$x$layout$yaxis$title$text, "logFC cell_controlled")
})

test_that(".plot_fcfc buckets genes into ns/both/NAME/PARTNER_NAME significance groups", {
  p <- .plot_fcfc(
    se_de,
    NAME = "baseline",
    PARTNER_NAME = "cell_controlled",
    GENE_VAR = "symbol",
    GENE_VARS = c("symbol")
  )
  trace_names <- vapply(p$x$data, function(d) d$name, character(1))

  expect_setequal(trace_names, c("ns", "both", "baseline", "cell_controlled"))
})

test_that(".plot_fcfc labels hover text with the gene identifier from GENE_VAR", {
  p <- .plot_fcfc(
    se_de,
    NAME = "baseline",
    PARTNER_NAME = "cell_controlled",
    GENE_VAR = "symbol",
    GENE_VARS = c("symbol")
  )
  all_text <- unlist(lapply(p$x$data, `[[`, "text"))

  expect_true(any(grepl("Gene: TSPAN6", all_text, fixed = TRUE)))
})
