## Tests for .get_split_des(), which cross-tabulates two DE comparisons
## (metadata()$de_results path) into "significant in A / B / both / ns"
## calls, filtered by a minimum log2FoldChange difference (CUTOFF).

make_se_split <- function(de_a, de_b, gene_data) {
  gene_ids <- rownames(gene_data)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      counts = matrix(
        1,
        nrow = length(gene_ids),
        ncol = 2,
        dimnames = list(gene_ids, c("s1", "s2"))
      )
    ),
    rowData = S4Vectors::DataFrame(gene_data)
  )
  S4Vectors::metadata(se)$de_results <- list(A = de_a, B = de_b)
  se
}

test_that(".get_split_des returns an empty data.frame when NAME == PARTNER_NAME", {
  gene_data <- data.frame(SYMBOL = "Ga", row.names = "g1")
  de <- data.frame(padj = 0.01, log2FoldChange = 2, row.names = "g1")
  se <- make_se_split(de, de, gene_data)

  result <- .get_split_des(se, NAME = "A", PARTNER_NAME = "A", CUTOFF = 1)

  expect_identical(result, data.frame())
})

test_that(".get_split_des classifies genes into both/NAME-only/PARTNER-only, filtered by CUTOFF", {
  gene_data <- data.frame(
    SYMBOL = c("Ga", "Gb", "Gc", "Gd", "Ge", "Gf"),
    row.names = paste0("g", 1:6)
  )

  de_a <- data.frame(
    padj = c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
    log2FoldChange = c(2, 2, 0.3, 2, 0.5, 3),
    row.names = paste0("g", 1:6)
  )
  de_b <- data.frame(
    padj = c(0.01, 0.01, 0.01, 0.9, 0.9, NA),
    log2FoldChange = c(2.1, 5, 2, 0.1, 0.1, 3),
    row.names = paste0("g", 1:6)
  )

  se <- make_se_split(de_a, de_b, gene_data)
  result <- .get_split_des(se, NAME = "A", PARTNER_NAME = "B", CUTOFF = 1)

  expect_setequal(result$gene_id, c("Gb", "Gc", "Gd"))
  expect_named(
    result,
    c(
      "gene_id",
      "log2FoldChange A",
      "log2FoldChange B",
      "Significant_In",
      "Effect_Classification"
    )
  )

  by_gene <- \(g) dplyr::filter(result, gene_id == g)

  expect_equal(
    by_gene("Gb")$Significant_In,
    factor("both", levels = levels(result$Significant_In))
  )
  expect_equal(
    as.character(by_gene("Gb")$Effect_Classification),
    "stronger effect in B"
  )

  expect_equal(as.character(by_gene("Gc")$Significant_In), "B")
  expect_equal(
    as.character(by_gene("Gc")$Effect_Classification),
    "stronger effect in B"
  )

  expect_equal(as.character(by_gene("Gd")$Significant_In), "A")
  expect_equal(
    as.character(by_gene("Gd")$Effect_Classification),
    "stronger effect in A"
  )
})

test_that(".get_split_des drops genes with NA padj, NA log2FoldChange or infinite log2FoldChange in NAME", {
  gene_data <- data.frame(
    SYMBOL = c("Keep", "BadPadj", "BadLfc", "BadInf"),
    row.names = c("keep", "bad_padj", "bad_lfc", "bad_inf")
  )
  de_a <- data.frame(
    padj = c(0.01, NA, 0.01, 0.01),
    log2FoldChange = c(2, 2, NA, Inf),
    row.names = c("keep", "bad_padj", "bad_lfc", "bad_inf")
  )
  de_b <- data.frame(
    padj = c(0.01, 0.01, 0.01, 0.01),
    log2FoldChange = c(3, 2, 2, 2),
    row.names = c("keep", "bad_padj", "bad_lfc", "bad_inf")
  )

  se <- make_se_split(de_a, de_b, gene_data)
  result <- .get_split_des(
    se,
    NAME = "A",
    PARTNER_NAME = "B",
    CUTOFF = 0,
    padj_CO = 1,
    fc_CO = 0
  )

  expect_equal(result$gene_id, "Keep")
})

test_that(".get_split_des uses GENE_VAR/GENE_VARS to pick the gene_id column", {
  gene_data <- data.frame(
    SYMBOL = "Ga",
    ENTREZ = "1234",
    row.names = "g1"
  )
  de_a <- data.frame(padj = 0.01, log2FoldChange = 2, row.names = "g1")
  de_b <- data.frame(padj = 0.01, log2FoldChange = 5, row.names = "g1")

  se <- make_se_split(de_a, de_b, gene_data)
  result <- .get_split_des(
    se,
    NAME = "A",
    PARTNER_NAME = "B",
    CUTOFF = 1,
    GENE_VAR = "ENTREZ",
    GENE_VARS = c("SYMBOL", "ENTREZ")
  )

  expect_equal(result$gene_id, "1234")
})

test_that(".get_split_des honors a custom padj_CO threshold for significance calls", {
  gene_data <- data.frame(SYMBOL = "Flex", row.names = "g1")
  de_a <- data.frame(padj = 0.03, log2FoldChange = 1.5, row.names = "g1")
  de_b <- data.frame(padj = 0.9, log2FoldChange = 0.1, row.names = "g1")
  se <- make_se_split(de_a, de_b, gene_data)

  sig <- .get_split_des(se, NAME = "A", PARTNER_NAME = "B", CUTOFF = 1)
  expect_equal(as.character(sig$Significant_In), "A")

  not_sig <- .get_split_des(
    se,
    NAME = "A",
    PARTNER_NAME = "B",
    CUTOFF = 1,
    padj_CO = 0.01
  )
  expect_equal(as.character(not_sig$Significant_In), "ns")
})

test_that(".get_split_des honors a custom fc_CO threshold for significance calls", {
  gene_data <- data.frame(SYMBOL = "Flex", row.names = "g1")
  de_a <- data.frame(padj = 0.01, log2FoldChange = 1.5, row.names = "g1")
  de_b <- data.frame(padj = 0.9, log2FoldChange = 0.1, row.names = "g1")
  se <- make_se_split(de_a, de_b, gene_data)

  sig <- .get_split_des(se, NAME = "A", PARTNER_NAME = "B", CUTOFF = 1)
  expect_equal(as.character(sig$Significant_In), "A")

  not_sig <- .get_split_des(
    se,
    NAME = "A",
    PARTNER_NAME = "B",
    CUTOFF = 1,
    fc_CO = 2
  )
  expect_equal(as.character(not_sig$Significant_In), "ns")
})
