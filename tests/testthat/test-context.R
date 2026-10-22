## Tests for the context helpers in explorer_utils.R:
## .check_context_enrichment, .parse_gene_list and .read_context_file

make_context_res <- function() {
  data.frame(
    SYMBOL = paste0("G", 1:1000),
    padj = c(rep(0.001, 100), rep(0.5, 900)),
    log2FoldChange = c(rep(3, 50), rep(-3, 50), rep(0, 900))
  )
}

write_context_file <- function(content, ext = ".txt") {
  path <- tempfile(fileext = ext)
  writeBin(charToRaw(content), path)
  path
}

# .check_context_enrichment -----

test_that(".check_context_enrichment returns one row per direction", {
  res <- .check_context_enrichment(
    make_context_res(),
    paste0("G", 1:40),
    gene_var = "SYMBOL"
  )

  expect_s3_class(res, "data.frame")
  expect_equal(res$direction, c("up regulated genes", "down regulated genes"))
  expect_equal(
    names(res),
    c(
      "direction",
      "n_genes",
      "context_size",
      "context_size_tested",
      "overlap",
      "fold_enrichment",
      "pvalue"
    )
  )
  expect_equal(res$n_genes, c(50, 50))
  expect_equal(res$overlap, c(40, 0))
})

test_that(".check_context_enrichment computes fold enrichment and p value", {
  res <- .check_context_enrichment(
    make_context_res(),
    paste0("G", 1:40),
    gene_var = "SYMBOL"
  )
  up <- res[res$direction == "up regulated genes", ]

  expect_equal(up$fold_enrichment, round((40 / 50) / (40 / 1000), 2))
  expect_equal(
    up$pvalue,
    stats::fisher.test(
      matrix(c(40, 10, 0, 950), nrow = 2),
      alternative = "greater"
    )$p.value
  )
  expect_lt(up$pvalue, 1e-10)
})

test_that(".check_context_enrichment separates raw and tested context size", {
  res <- .check_context_enrichment(
    make_context_res(),
    c("G1", "G2", "G2", "NOT_TESTED", NA),
    gene_var = "SYMBOL"
  )

  expect_equal(unique(res$context_size), 3)
  expect_equal(unique(res$context_size_tested), 2)
})

test_that(".check_context_enrichment drops genes with missing statistics from the universe", {
  de <- make_context_res()
  de$padj[1000] <- NA

  res <- .check_context_enrichment(de, c("G1", "G1000"), gene_var = "SYMBOL")

  expect_equal(unique(res$context_size_tested), 1)
})

test_that(".check_context_enrichment drops directions without DE genes", {
  de <- make_context_res()
  de <- de[de$log2FoldChange >= 0, ]

  res <- .check_context_enrichment(de, paste0("G", 1:40), gene_var = "SYMBOL")

  expect_equal(res$direction, "up regulated genes")
})

test_that(".check_context_enrichment returns an empty data.frame without matches", {
  res <- .check_context_enrichment(
    make_context_res(),
    "NOT_TESTED",
    gene_var = "SYMBOL"
  )
  expect_equal(nrow(res), 0)

  res <- .check_context_enrichment(
    make_context_res(),
    paste0("G", 1:40),
    padj_CO = 1e-10,
    gene_var = "SYMBOL"
  )
  expect_equal(nrow(res), 0)
})

test_that(".check_context_enrichment uses strict fold change cutoffs", {
  de <- make_context_res()
  de$log2FoldChange[1:50] <- 1

  res <- .check_context_enrichment(de, paste0("G", 1:40), gene_var = "SYMBOL")

  expect_equal(res$direction, "down regulated genes")
})

# .parse_gene_list -----

test_that(".parse_gene_list splits on whitespace, commas and semicolons", {
  expect_equal(
    .parse_gene_list("TP53 MYC,EGFR;KRAS\nBRCA1\tCDK4"),
    c("TP53", "MYC", "EGFR", "KRAS", "BRCA1", "CDK4")
  )
})

test_that(".parse_gene_list strips quotes, empty entries and duplicates", {
  expect_equal(
    .parse_gene_list(c("\"TP53, MYC\"", "", "  EGFR  ", "\"TP53\"")),
    c("TP53", "MYC", "EGFR")
  )
})

test_that(".parse_gene_list returns an empty vector for empty input", {
  expect_equal(.parse_gene_list(""), character())
  expect_equal(.parse_gene_list(NULL), character())
})

# .read_context_file -----

test_that(".read_context_file reads a single column as one unnamed list", {
  path <- write_context_file("TP53\nMYC\n\nEGFR\n")
  res <- .read_context_file(path)

  expect_null(names(res))
  expect_equal(res[[1]], c("TP53", "MYC", "EGFR"))
})

test_that(".read_context_file reads a separated single line as one list", {
  path <- write_context_file("TP53 MYC, EGFR;KRAS\n")
  res <- .read_context_file(path)

  expect_null(names(res))
  expect_equal(res[[1]], c("TP53", "MYC", "EGFR", "KRAS"))
})

test_that(".read_context_file handles quotes, a BOM and Windows line endings", {
  path <- write_context_file("﻿\"TP53\"\r\n\"MYC\"\r\n", ".csv")
  res <- .read_context_file(path)

  expect_equal(res[[1]], c("TP53", "MYC"))
})

test_that(".read_context_file reads multi-column csv files per column", {
  path <- write_context_file(
    "\"up\",\"down\"\r\n\"TP53\",\"KRAS\"\r\n\"MYC\",\"\"\r\n",
    ".csv"
  )
  res <- .read_context_file(path)

  expect_equal(res, list(up = c("TP53", "MYC"), down = "KRAS"))
})

test_that(".read_context_file reads ragged tsv files", {
  path <- write_context_file(
    "setA\tsetB\tsetC\nTP53\tKRAS\tBRCA1\nMYC\tEGFR\nCDK4\n",
    ".tsv"
  )
  res <- .read_context_file(path)

  expect_equal(
    res,
    list(setA = c("TP53", "MYC", "CDK4"), setB = c("KRAS", "EGFR"), setC = "BRCA1")
  )
})

test_that(".read_context_file reads semicolon separated files", {
  path <- write_context_file("a;b\nTP53;MYC\n", ".csv")

  expect_equal(.read_context_file(path), list(a = "TP53", b = "MYC"))
})

test_that(".read_context_file names columns with an empty header", {
  path <- write_context_file(",x\nTP53,MYC\n", ".csv")

  expect_equal(names(.read_context_file(path)), c("Column 1", "x"))
})

# .validate_context -----

test_that(".validate_context returns NULL for NULL", {
  expect_null(.validate_context(NULL))
})

test_that(".validate_context cleans the gene vectors", {
  res <- .validate_context(list(
    a = c("TP53", " MYC ", NA, "", "TP53"),
    b = factor(c("EGFR", "KRAS")),
    c = c(7157, 4609)
  ))

  expect_equal(
    res,
    list(a = c("TP53", "MYC"), b = c("EGFR", "KRAS"), c = c("7157", "4609"))
  )
})

test_that(".validate_context only accepts named lists", {
  expect_error(.validate_context(c("TP53", "MYC")), "named list")
  expect_error(.validate_context(data.frame(a = "TP53")), "named list")
  expect_error(.validate_context(list("TP53", "MYC")), "must be named")
  expect_error(.validate_context(list(a = "TP53", "MYC")), "must be named")
  expect_error(
    .validate_context(list(a = "TP53", a = "MYC")),
    "must be unique"
  )
})

test_that(".validate_context rejects non-vector and empty elements", {
  expect_error(
    .validate_context(list(a = list("TP53"))),
    "vector of gene identifiers"
  )
  expect_error(
    .validate_context(list(a = "TP53", b = c(NA, ""))),
    "contain no genes: b"
  )
})

test_that("exploreSE validates the context argument", {
  expect_error(exploreSE(context = "TP53"), "named list")
})

# .get_context_overlap -----

test_that(".get_context_overlap returns the sorted overlap per direction", {
  res <- .get_context_overlap(
    make_context_res(),
    c("G60", "G3", "G1", "G500", "NOT_TESTED"),
    gene_var = "SYMBOL"
  )

  expect_equal(
    res,
    list(
      "up regulated genes" = c("G1", "G3"),
      "down regulated genes" = "G60"
    )
  )
})

test_that(".get_context_overlap matches the overlap column of the enrichment table", {
  ctx <- paste0("G", c(1:40, 95:120))
  overlap <- .get_context_overlap(make_context_res(), ctx, gene_var = "SYMBOL")
  tbl <- .check_context_enrichment(make_context_res(), ctx, gene_var = "SYMBOL")

  expect_equal(unname(lengths(overlap)), tbl$overlap)
})

test_that(".get_context_overlap drops directions without DE genes", {
  de <- make_context_res()
  de <- de[de$log2FoldChange >= 0, ]

  res <- .get_context_overlap(de, c("G1", "G2"), gene_var = "SYMBOL")

  expect_equal(names(res), "up regulated genes")
})
