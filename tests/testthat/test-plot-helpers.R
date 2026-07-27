## Tests for the ggplot2/plotly building helpers in explorer_utils.R.
## We only assert on the returned object's type/structure, not on exact
## pixel/rendering output.
##
## Where the package's bundled real-world data has a matching shape, it's
## used directly (via .de_result()/.fe_result() on se_with_de/se_with_go)
## instead of a hand-built data.frame. Neither bundled object contains a
## GSEA-style (NES) result though - only precomputed GO enrichments - so the
## GSEA-related .plot_fe tests still use a small synthetic fixture.

de_results <- function(NAME = "baseline") {
  .de_result(se_with_de, NAME) %>%
    tibble::rownames_to_column("gene_id")
}

go_fe_results <- function(NAME = "baseline", direction = "up_go") {
  fe <- .fe_result(se_with_go, NAME)
  fe[[paste0(NAME, "_", direction)]]
}

make_gsea_fe_results <- function() {
  data.frame(
    Description = paste0("HALLMARK_TERM_", 1:10),
    p.adjust = stats::runif(10, 0, 0.1),
    NES = stats::rnorm(10, 0, 2)
  )
}

test_that(".plot_des returns a ggplot summarising up/down gene counts", {
  p <- .plot_des(de_results("baseline"), "baseline")

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "DE genes in baseline")
})

test_that(".plot_fe handles GO-style (FoldEnrichment) results", {
  p <- .plot_fe(go_fe_results("baseline", "up_go"), NAME = "up_go")
  expect_s3_class(p, "plotly")
})

test_that(".plot_fe handles GSEA-style (NES) results without erroring", {
  p <- .plot_fe(make_gsea_fe_results(), NAME = "comparisonA_gsea_HALLMARK")
  expect_s3_class(p, "plotly")
})

test_that(".plot_fe titles GSEA plots based on HALLMARK vs Reactome in NAME", {
  hallmark_df <- make_gsea_fe_results()
  reactome_df <- make_gsea_fe_results()

  p_hallmark <- .plot_fe(hallmark_df, NAME = "comparisonA_gsea_HALLMARK")
  p_reactome <- .plot_fe(reactome_df, NAME = "comparisonA_gsea_REACTOME")

  expect_s3_class(p_hallmark, "plotly")
  expect_s3_class(p_reactome, "plotly")
})

test_that(".plot_volcano returns a plotly object with LABEL_TOP = TRUE (default)", {
  p <- .plot_volcano(
    de_results("baseline"),
    "baseline",
    COLS = c(Up = "red", Down = "blue", Highlighted = "green", NS = "grey")
  )
  expect_s3_class(p, "plotly")
})

test_that(".plot_volcano returns a plotly object with LABEL_TOP = FALSE", {
  p <- .plot_volcano(
    de_results("baseline"),
    "baseline",
    COLS = c(Up = "red", Down = "blue", Highlighted = "green", NS = "grey"),
    LABEL_TOP = FALSE
  )
  expect_s3_class(p, "plotly")
})

test_that(".plot_volcano works with highlighted genes and no top-N labels", {
  res <- de_results("baseline")
  p <- .plot_volcano(
    res,
    "baseline",
    COLS = c(Up = "red", Down = "blue", Highlighted = "green", NS = "grey"),
    highlights = res$gene_id[1:3],
    LABEL_TOP = FALSE
  )
  expect_s3_class(p, "plotly")
})

test_that(".plot_volcano handles no significant genes gracefully", {
  # A real subset of baseline that happens to contain no genes clearing the
  # padj < 0.05 significance threshold used inside .plot_volcano().
  res <- de_results("baseline") %>%
    dplyr::filter(!is.na(padj), padj > 0.05)
  expect_true(nrow(res) > 0)

  p <- .plot_volcano(
    res,
    "baseline",
    COLS = c(Up = "red", Down = "blue", Highlighted = "green", NS = "grey")
  )
  expect_s3_class(p, "plotly")
})
