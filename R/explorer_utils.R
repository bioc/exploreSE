#' checks for the presence of pre-computed DE results in the
#' the object. It checks the DeeDee slot if applicable, otherwise
#' it checks metadata(object)$de_results
#' @noRd
.check_precomputed_de <- function(SE) {
    if (is.null(SE)) {
        return(FALSE)
    }
    if (
        methods::is(SE, "DeeDeeExperiment") &&
            !is.null(DeeDeeExperiment::getDEANames(SE))
    ) {
        return(TRUE)
    } else {
        !is.null(S4Vectors::metadata(SE)$de_results) &&
            is.list(S4Vectors::metadata(SE)$de_results) &&
            length(S4Vectors::metadata(SE)$de_results) > 0
    }
}

#' checks for the presence of pre-computed FE results in the
#' the object. It checks the DeeDee slot if applicable, otherwise
#' it checks metadata(object)$fe_results
#' @noRd
.check_precomputed_fe <- function(SE) {
    if (is.null(SE)) {
        return(FALSE)
    }
    if (
        methods::is(SE, "DeeDeeExperiment") &&
            !is.null(DeeDeeExperiment::getFEANames(SE))
    ) {
        return(TRUE)
    } else {
        !is.null(S4Vectors::metadata(SE)$fe_results) &&
            is.list(S4Vectors::metadata(SE)$fe_results) &&
            length(S4Vectors::metadata(SE)$fe_results) > 0
    }
}

#' helper funciton to fetch the names of the DE comparisons in the object
#' @noRd
.de_results_names <- function(SE) {
    if (methods::is(SE, "DeeDeeExperiment")) {
        DeeDeeExperiment::getDEANames(SE)
    } else {
        names(S4Vectors::metadata(SE)$de_results)
    }
}

#' helper funciton to fetch the names of the FE comparisons in the object,
#' belonging to a specified DE comparison
#' @noRd
.fe_results_names <- function(SE, NAME) {
    if (methods::is(SE, "DeeDeeExperiment")) {
        return(names(DeeDeeExperiment::getFEAList(SE, NAME)))
    }
    return(names(S4Vectors::metadata(SE)$fe_results[[NAME]]))
}

# fetches a specified DE comparison by name
#' @noRd
.de_result <- function(SE, NAME) {
    if (methods::is(SE, "DeeDeeExperiment")) {
        res <- DeeDeeExperiment::getDEA(SE, NAME, type = "data.frame") %>%
            dplyr::rename_with(.fn = \(x) {
                stringr::str_remove(x, paste0(NAME, "_"))
            })
    } else {
        res <- as.data.frame(S4Vectors::metadata(SE)$de_results[[NAME]])
    }
}

# fetches a specified DE comparison by name
#' @noRd
.fe_result <- function(SE, NAME) {
    fe_list <- DeeDeeExperiment::getFEAList(SE, NAME, format = "original") %>%
        purrr::map(BiocGenerics::as.data.frame)
    return(fe_list)
}

# plotting function for DE genes, returns a bar chart
#' @noRd
.plot_des <- function(
    RES,
    NAME,
    padj_CO = 0.05,
    fc_CO = 1,
    COLS = c("darkred", "darkblue")
) {
    p <- RES %>%
        dplyr::filter(padj < padj_CO, abs(log2FoldChange) >= fc_CO) %>%
        dplyr::mutate(
            direction = ifelse(log2FoldChange > fc_CO, "Up", "Down")
        ) %>%
        dplyr::count(direction) %>%
        ggplot2::ggplot(ggplot2::aes(direction, y = n, fill = direction)) +
        ggplot2::geom_col() +
        ggplot2::geom_label(
            ggplot2::aes(label = n),
            fill = "white",
            show.legend = FALSE
        ) +
        ggplot2::scale_fill_manual(values = COLS) +
        ggplot2::labs(title = paste("DE genes in", NAME), ) +
        ggplot2::theme_light(base_size = 14) +
        ggplot2::theme(legend.position = "none")
    return(p)
}
#' takes the functional enrichments and creates a bar chart
#' recognises GSEA results  by the use of NES in the colnames
#' for ORA analyiss it plots the Fold Enrichment when available or
#' able to be calculated
#' @noRd
.plot_fe <- function(
    FE,
    NAME = "up_go",
    padj_CO = 0.05,
    N_terms = 5,
    COLS = c("Up" = "darkred", "Down" = "blue"),
    SHOW_GENES = TRUE
) {
    # this can probably go when the DeeDee devel branch hits main, but for now needs to stay
    padj_col <- ifelse("gs_p.adjust" %in% names(FE), "gs_p.adjust", "p.adjust")

    if (sum(stringr::str_detect(names(FE), "NES")) > 0) {
        gsea_df <- FE %>%
            dplyr::filter(!!rlang::sym(padj_col) < padj_CO) %>%
            dplyr::group_by(sign(NES)) %>%
            dplyr::slice_max(abs(NES), n = round(N_terms / 2)) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(
                Description = stringr::str_remove(Description, "HALLMARK_") %>%
                    stringr::str_remove("REACTOME") %>%
                    stringr::str_replace_all("_", " ") %>%
                    forcats::fct_reorder(NES),
                dir = ifelse(NES > 0, "Up", "Down")
            )

        p <- ggplot2::ggplot(gsea_df, ggplot2::aes(NES, Description)) +
            ggplot2::geom_col(ggplot2::aes(fill = dir)) +
            ggplot2::scale_y_discrete(labels = \(x) {
                stringr::str_wrap(x, width = 60)
            }) +
            ggplot2::scale_fill_manual(values = COLS) +
            ggplot2::theme_light(base_size = 14) +
            ggplot2::labs(
                x = "NES (Normalized Enrichment Score)",
                title = if (stringr::str_detect(NAME, "HALLMARK")) {
                    paste("Top", N_terms, "enriched hallmark gene sets")
                } else {
                    paste("Top", N_terms, "enriched Reactome gene sets")
                }
            ) +
            ggplot2::theme(
                axis.title.y = ggplot2::element_blank(),
                legend.position = "none"
            )
    } else {
        # this can probably go, when the DeeDee devel branch hits main, but for now needs to stay
        filter_col <- dplyr::case_when(
            "FoldEnrichment" %in% names(FE) ~ "FoldEnrichment",
            "GeneRatio" %in%
                names(FE) &
                "BgRatio" %in% names(FE) ~ "FoldEnrichment",
            "gs_p.adjust" %in% names(FE) ~ "gs_p.adjust",
            TRUE ~ "p.adjust"
        )

        if (
            filter_col == "FoldEnrichment" & !("FoldEnrichment" %in% names(FE))
        ) {
            FE <- FE %>%
                dplyr::rowwise() %>%
                dplyr::mutate(
                    GR = stringr::str_split(GeneRatio, "/") %>%
                        unlist() %>%
                        as.numeric() %>%
                        (\(x) x[[1]] / x[[2]]),
                    BR = stringr::str_split(BgRatio, "/") %>%
                        unlist() %>%
                        as.numeric() %>%
                        (\(x) x[[1]] / x[[2]])
                ) %>%
                dplyr::ungroup() %>%
                dplyr::mutate(FoldEnrichment = GR / BR) %>%
                dplyr::select(-GR, -BR)
        }

        if (filter_col == "FoldEnrichment") {
            df <- FE %>%
                dplyr::filter(!!rlang::sym(padj_col) < padj_CO) %>%
                dplyr::slice_max(FoldEnrichment, n = N_terms) %>%
                dplyr::mutate(
                    Description = forcats::fct_reorder(
                        Description,
                        FoldEnrichment
                    ),
                    dir = ifelse(
                        stringr::str_detect(NAME, "[Uu][Pp]"),
                        "Up",
                        "Down"
                    )
                )
            p <- ggplot2::ggplot(
                df,
                ggplot2::aes(
                    FoldEnrichment,
                    Description,
                    text = stringr::str_wrap(
                        stringr::str_replace_all(geneID, "\\/", ", "),
                        width = 60
                    )
                )
            ) +
                ggplot2::geom_col(ggplot2::aes(fill = dir)) +
                ggplot2::scale_y_discrete(labels = \(x) {
                    stringr::str_wrap(x, width = 60)
                }) +
                ggplot2::scale_fill_manual(values = COLS) +
                ggplot2::theme_light(base_size = 14) +
                ggplot2::labs(
                    x = "Fold Enrichment over Background",
                    title = if (stringr::str_detect(NAME, "up")) {
                        paste("Top", N_terms, "upregulated GO Terms")
                    } else {
                        paste("Top", N_terms, "downregulated GO Terms")
                    }
                ) +
                ggplot2::theme(
                    axis.title.y = ggplot2::element_blank(),
                    legend.position = "none"
                )
        } else {
            df <- FE %>%
                dplyr::filter(!!rlang::sym(padj_col) < padj_CO) %>%
                dplyr::slice_min(!!rlang::sym(filter_col), n = N_terms) %>%
                dplyr::mutate(
                    padj = -log10(!!rlang::sym(filter_col)),
                    Description = forcats::fct_reorder(Description, padj),
                    dir = ifelse(
                        stringr::str_detect(NAME, "[Uu][Pp]"),
                        "Up",
                        "Down"
                    )
                )
            p <- ggplot2::ggplot(
                df,
                ggplot2::aes(
                    padj,
                    Description,
                    text = stringr::str_wrap(
                        stringr::str_replace_all(geneID, "\\/", ", "),
                        width = 60
                    )
                )
            ) +
                ggplot2::geom_col(ggplot2::aes(fill = dir)) +
                ggplot2::scale_y_discrete(labels = \(x) {
                    stringr::str_wrap(x, width = 60)
                }) +
                ggplot2::scale_fill_manual(values = COLS) +
                ggplot2::theme_light(base_size = 14) +
                ggplot2::labs(
                    x = "-log10 adjusted p-Value",
                    title = if (stringr::str_detect(NAME, "[Uu][Pp]")) {
                        paste("Top", N_terms, "upregulated GO Terms")
                    } else {
                        paste("Top", N_terms, "downregulated GO Terms")
                    }
                ) +
                ggplot2::theme(
                    axis.title.y = ggplot2::element_blank(),
                    legend.position = "none"
                )
        }
    }

    if (!SHOW_GENES) {
        p <- p + ggplot2::aes(text = NULL)
    }
    p_p <- plotly::ggplotly(p, tooltip = c("text")) %>%
        plotly::layout(hovermode = "closest") %>%
        plotly::style(textposition = "right")
    return(p_p)
}
#' plots a volcano plot, highlights selected genes.
#' @noRd
.plot_volcano <- function(
    RES,
    NAME,
    padj_CO = 0.05,
    fc_CO = 1,
    highlights = NULL,
    COLS,
    LABEL_TOP = TRUE,
    TOPN = 5,
    gene_var = "gene_id",
    LABEL_HIGHLIGHTS = TRUE
) {
    gene_var <- rlang::sym(gene_var)
    # Prep data
    volcano_df <- RES %>%
        dplyr::filter(
            !is.na(padj) & !is.na(log2FoldChange) & !is.infinite(log2FoldChange)
        ) %>%
        dplyr::mutate(gene_id = !!gene_var) %>%
        dplyr::mutate(
            highlighted = gene_id %in% highlights,
            sig = dplyr::case_when(
                padj < padj_CO & log2FoldChange > fc_CO ~ "Up",
                padj < padj_CO & log2FoldChange < -fc_CO ~ "Down",
                TRUE ~ "NS"
            )
        ) %>%
        dplyr::mutate(
            display_category = dplyr::case_when(
                highlighted ~ "Highlighted",
                TRUE ~ sig
            ) |>
                factor(levels = c("NS", "Down", "Up", "Highlighted")),
            emphasised = highlighted & sig != "NS"
        )

    # Identify top genes
    genes_to_label <- data.frame()
    if (LABEL_TOP && TOPN > 0) {
        top_genes <- volcano_df %>%
            dplyr::filter(sig != "NS" & !highlighted) %>%
            dplyr::slice_min(padj, n = TOPN)
        genes_to_label <- rbind(genes_to_label, top_genes)
    }

    # highlighted genes
    if (LABEL_HIGHLIGHTS && length(highlights) > 0) {
        highlighted_genes <- volcano_df %>%
            dplyr::filter(highlighted)
        genes_to_label <- rbind(genes_to_label, highlighted_genes)
    }

    #  plot
    p <- ggplot2::ggplot(
        volcano_df,
        ggplot2::aes(
            x = log2FoldChange,
            y = -log10(padj),
            color = display_category,
            text = paste(
                "log2 FC:",
                round(log2FoldChange, 2),
                "\nGene:",
                gene_id,
                "\nadjusted p-value:",
                format(padj, digits = 4)
            )
        )
    ) +
        ggplot2::geom_point(ggplot2::aes(
            size = emphasised,
            alpha = ifelse(emphasised, 1, 0.6)
        )) +
        ggplot2::scale_size_manual(
            values = c("TRUE" = 3, "FALSE" = 1.5),
            guide = "none"
        ) +
        ggplot2::scale_alpha_identity() +
        ggplot2::scale_color_manual(
            values = COLS,
            name = "Category",
            breaks = c("Up", "Down", "Highlighted", "NS"),
            labels = c(
                "Up-regulated",
                "Down-regulated",
                "Highlighted",
                "Not significant"
            )
        ) +
        ggplot2::geom_vline(
            xintercept = c(-fc_CO, fc_CO),
            linetype = "dashed",
            color = "grey30",
            linewidth = 0.5
        ) +
        ggplot2::geom_hline(
            yintercept = -log10(padj_CO),
            linetype = "dashed",
            color = "grey30",
            linewidth = 0.5
        ) +
        ggplot2::labs(
            x = "Log2 Fold Change",
            y = "-Log10 Adjusted P-value",
            title = if (!is.null(NAME)) {
                paste("Volcano Plot:", NAME)
            } else {
                "Volcano Plot"
            }
        ) +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(
            legend.position = "right",
            panel.grid.minor = ggplot2::element_blank()
        )

    # Add labels
    if (nrow(genes_to_label) > 0) {
        p <- p +
            ggplot2::geom_text(
                data = genes_to_label,
                ggplot2::aes(label = gene_id),
                size = 5,
                show.legend = FALSE
            )
    }

    p_p <- plotly::ggplotly(p, tooltip = c("text")) %>%
        plotly::layout(
            hovermode = "closest",

            # threshold note at the right end of the padj line
            annotations = list(
                x = 1,
                xref = "paper",
                xanchor = "right",
                y = -log10(padj_CO),
                yref = "y",
                yanchor = "bottom",
                text = paste(
                    "FDR:",
                    format(padj_CO, scientific = FALSE, drop0trailing = TRUE)
                ),
                showarrow = FALSE,
                font = list(color = "grey30")
            )
        ) %>%
        plotly::style(textposition = "right")

    return(p_p)
}

#' creates the expression plot; still needs work re: normalisation of counts
#' @noRd
.plot_expression <- function(
    SE,
    GENE,
    BY,
    LEVELS,
    PLOT_TYPE,
    GENE_VAR,
    LOG = FALSE
) {
    # get the right row to plot
    gene <- rownames(SummarizedExperiment::rowData(SE)[
        which(
            SummarizedExperiment::rowData(SE)[, GENE_VAR] == GENE
        ),
    ])[1]

    if ("counts" %notin% SummarizedExperiment::assayNames(SE)) {
        message("counts not found in assays, defaulting to first assay")
        counts_data <- SummarizedExperiment::assay(SE)[[1]][gene, ]
    } else if ("sizeFactor" %in% names(SummarizedExperiment::colData(SE))) {
        counts_data <- SummarizedExperiment::assay(SE, "counts")
        sf <- SummarizedExperiment::colData(SE)[, "sizeFactor"]
        counts_data <- sweep(counts_data, 2, sf, "/")[gene, ]
    } else {
        counts_data <- SummarizedExperiment::assay(SE, "counts")[gene, ]
    }

    plot_df <- data.frame(
        s_a_m_p_l_e = colnames(SE),
        expression = counts_data,
        group = forcats::as_factor(SummarizedExperiment::colData(SE)[[
            BY
        ]])
    ) %>%
        dplyr::filter(group %in% LEVELS)

    p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
            x = group,
            y = expression,
            fill = group,
            text = s_a_m_p_l_e
        )
    ) +
        ggplot2::labs(
            title = paste("Expression:", GENE),
            x = BY,
            y = "Counts"
        ) +
        ggplot2::scale_x_discrete(labels = \(x) {
            stringr::str_wrap(stringr::str_replace_all(x, "_", " "), 10)
        }) +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(legend.position = "none")

    if (PLOT_TYPE == "violin" && any(table(droplevels(plot_df$group)) < 2)) {
        # group size check
        shiny::showNotification(
            "Some selected groups have fewer than two samples; showing a box plot instead of violin.",
            type = "warning"
        )
        PLOT_TYPE <- "box"
    }

    if (PLOT_TYPE == "box") {
        p <- p +
            ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
            ggplot2::geom_jitter(width = 0.2, alpha = 0.5, size = 2)
    } else {
        # removing the text geom from violin, cant handle it
        p <- p +
            ggplot2::geom_violin(ggplot2::aes(text = NULL), alpha = 0.7) +
            ggplot2::geom_jitter(width = 0.1, alpha = 0.5, size = 2)
    }

    if (LOG) {
        p <- p + ggplot2::scale_y_log10()
    }

    p
}


#' filters all DEAs for the gene at hand.
#' @noRd
.get_de_status_table <- function(
    .SE,
    GENE,
    GENE_VAR,
    padj_CO = 0.05,
    fc_CO = 1
) {
    row_data <- SummarizedExperiment::rowData(.SE) |>
        as.data.frame() |>
        dplyr::select(tidyselect::any_of(GENE_VAR)) |>
        tibble::rownames_to_column("gene_ident")

    df_signif <- .de_results_names(.SE) |>
        purrr::map(\(x) .de_result(SE = .SE, NAME = x)) |>
        purrr::map(as.data.frame) |>
        purrr::map(\(x) tibble::rownames_to_column(x, "gene_ident")) |>
        purrr::set_names(.de_results_names(.SE)) |>
        dplyr::bind_rows(.id = "comparison") |>
        dplyr::filter(padj < padj_CO, abs(log2FoldChange) > fc_CO) |>
        dplyr::left_join(row_data, by = dplyr::join_by(gene_ident)) |>
        dplyr::filter(.data[[GENE_VAR]] == GENE) |>
        dplyr::select(
            tidyselect::any_of(c(
                "comparison",
                GENE_VAR,
                "baseMean",
                "log2FoldChange",
                "pvalue",
                "padj"
            ))
        ) |>
        dplyr::mutate(
            dplyr::across(
                tidyselect::any_of(c("baseMean", "log2FoldChange")),
                \(x) round(x, 2)
            ),
            dplyr::across(tidyselect::any_of(c("pvalue", "padj")), \(x) {
                round(x, 6)
            })
        )
    return(df_signif)
}


#' Creates some random data to demonstrate the colData, PCA and expression elements.
#' @noRd
.create_demo_data <- function() {
    n_genes <- 1000
    n_samples <- 12

    # Simulate count data
    counts <- matrix(
        stats::rnbinom(n_genes * n_samples, mu = 100, size = 1 / 0.5),
        nrow = n_genes
    )
    rownames(counts) <- paste0("GENE", seq_len(n_genes))
    colnames(counts) <- paste0("Sample", seq_len(n_samples))

    # Add some differential expression
    de_genes <- seq_len(100)
    counts[de_genes, seq_len(6)] <- counts[de_genes, seq_len(6)] * 3

    # Sample metadata
    colData <- S4Vectors::DataFrame(
        condition = factor(rep(c("Control", "Treatment"), each = 6)),
        batch = factor(rep(c("A", "B"), times = 6)),
        replicate = rep(seq_len(6), times = 2)
    )
    rownames(colData) <- colnames(counts)

    # Gene metadata
    rowData <- S4Vectors::DataFrame(
        gene_id = rownames(counts),
        gene_name = paste0("Gene_", seq_len(n_genes))
    )

    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts),
        colData = colData,
        rowData = rowData
    )
}

#' creates fc-fc plots
#' @noRd
.plot_fcfc <- function(
    OBJ,
    NAME,
    PARTNER_NAME,
    GENE_VAR = "SYMBOL",
    GENE_VARS = c("SYMBOL"),
    padj_CO = 0.05,
    fc_CO = 1,
    COLS = rev(c("red", "blue", "purple", "grey")),
    LABEL_TOP = TRUE,
    TOPN = 5
) {
    COLS <- COLS |>
        purrr::set_names(c(NAME, PARTNER_NAME, "ns", "both"))
    gene_var <- rlang::sym(GENE_VAR)

    df_a <- .de_result(OBJ, NAME = NAME) |>
        tibble::rownames_to_column("gene_ident")

    df_b <- .de_result(OBJ, NAME = PARTNER_NAME) |>
        tibble::rownames_to_column("gene_ident")

    df_list <- list(df_a, df_b) |>
        purrr::set_names(c(NAME, PARTNER_NAME)) |>
        purrr::map(\(x) {
            dplyr::filter(
                x,
                !is.na(padj),
                !is.na(log2FoldChange),
                !is.infinite(log2FoldChange)
            )
        }) |>
        purrr::imap(\(DF, .NAME) {
            dplyr::mutate(
                DF,
                sig = ifelse(
                    padj < padj_CO & abs(log2FoldChange) > fc_CO,
                    .NAME,
                    "ns"
                )
            )
        })

    df_comb <- dplyr::left_join(
        df_list[[NAME]],
        df_list[[PARTNER_NAME]],
        by = dplyr::join_by(gene_ident),
        suffix = paste0("_", c("A", "B"))
    ) |>
        dplyr::mutate(
            sig = dplyr::case_when(
                sig_A == "ns" & sig_B == "ns" ~ "ns",
                sig_A == NAME & sig_B == PARTNER_NAME ~ "both",
                sig_A == "ns" ~ PARTNER_NAME,
                sig_B == "ns" ~ NAME
            ) |>
                forcats::fct_relevel("ns", "both", PARTNER_NAME)
        ) |>
        dplyr::select(-sig_A, -sig_B) |>
        dplyr::filter(!is.na(sig))

    merging_data <- SummarizedExperiment::rowData(OBJ) |>
        as.data.frame() |>
        dplyr::select(tidyselect::any_of(GENE_VARS)) |>
        tibble::rownames_to_column("gene_ident")

    df_comb <- merging_data |>
        dplyr::inner_join(df_comb, by = dplyr::join_by(gene_ident)) |>
        dplyr::mutate(gene_id = !!gene_var)

    p <- ggplot2::ggplot(
        df_comb,
        ggplot2::aes(
            log2FoldChange_A,
            log2FoldChange_B,
            color = sig,
            text = paste(
                "Gene:",
                gene_id,
                "\n",
                paste("logFC", NAME),
                round(log2FoldChange_A, 2),
                "\n",
                paste("logFC", PARTNER_NAME),
                round(log2FoldChange_B, 2)
            )
        )
    ) +
        # added line for visibility of the "ideal" overlap of a=b
        ggplot2::geom_abline(
            slope = 1,
            intercept = 0,
            alpha = 0.5,
            linetype = 2,
            show.legend = FALSE
        ) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::scale_color_manual(values = COLS) +
        ggplot2::labs(
            x = paste("logFC", NAME),
            y = paste("logFC", PARTNER_NAME),
            color = "Signficant in:"
        ) +
        ggplot2::theme_minimal(base_size = 14)

    p <- plotly::ggplotly(p, tooltip = c("text"))
    return(p)
}


#' creates fc-fc difference table
#' @noRd
.get_split_des <- function(
    OBJ,
    NAME,
    PARTNER_NAME,
    CUTOFF,
    GENE_VAR = "SYMBOL",
    GENE_VARS = c("SYMBOL"),
    padj_CO = 0.05,
    fc_CO = 1
) {
    if (NAME == PARTNER_NAME) {
        return(data.frame())
    }

    gene_var <- rlang::sym(GENE_VAR)
    df_a <- .de_result(OBJ, NAME = NAME) |>
        tibble::rownames_to_column("gene_ident")

    df_b <- .de_result(OBJ, NAME = PARTNER_NAME) |>
        tibble::rownames_to_column("gene_ident")

    df_list <- list(df_a, df_b) |>
        purrr::set_names(c(NAME, PARTNER_NAME)) |>
        purrr::map(\(x) {
            dplyr::filter(
                x,
                !is.na(padj),
                !is.na(log2FoldChange),
                !is.infinite(log2FoldChange)
            )
        }) |>
        purrr::imap(\(DF, .NAME) {
            dplyr::mutate(
                DF,
                sig = ifelse(
                    padj < padj_CO & abs(log2FoldChange) > fc_CO,
                    .NAME,
                    "ns"
                )
            )
        })

    df_comb <- dplyr::left_join(
        df_list[[NAME]],
        df_list[[PARTNER_NAME]],
        by = dplyr::join_by(gene_ident),
        suffix = paste0("_", c("A", "B"))
    ) |>
        dplyr::mutate(
            sig = dplyr::case_when(
                sig_A == "ns" & sig_B == "ns" ~ "ns",
                sig_A == NAME & sig_B == PARTNER_NAME ~ "both",
                sig_A == "ns" ~ PARTNER_NAME,
                sig_B == "ns" ~ NAME
            ) |>
                forcats::fct_relevel("ns", "both", PARTNER_NAME)
        ) |>
        dplyr::select(-sig_A, -sig_B) |>
        dplyr::filter(!is.na(sig))

    merging_data <- SummarizedExperiment::rowData(OBJ) |>
        as.data.frame() |>
        dplyr::select(tidyselect::any_of(GENE_VARS)) |>
        tibble::rownames_to_column("gene_ident")

    df_comb <- merging_data |>
        dplyr::inner_join(df_comb, by = dplyr::join_by(gene_ident)) |>
        dplyr::mutate(gene_id = !!gene_var)

    genes_with_diff <- df_comb |>
        dplyr::filter(
            sig != "ns",
            abs(log2FoldChange_A - log2FoldChange_B) > CUTOFF
        ) |>
        dplyr::mutate(
            Effect_Classification = dplyr::case_when(
                log2FoldChange_A > log2FoldChange_B ~ paste(
                    "stronger effect in",
                    NAME
                ),
                log2FoldChange_B > log2FoldChange_A ~ paste(
                    "stronger effect in",
                    PARTNER_NAME
                ),
                T ~ "No Diff"
            ) |>
                forcats::fct_relevel("No Diff", after = Inf),
            dplyr::across(tidyselect::contains("log2Fol"), \(x) round(x, 2)),
            dplyr::across(tidyselect::contains("padj"), \(x) round(x, 4))
        ) |>
        dplyr::rename(Significant_In = sig) |>
        dplyr::rename_with(\(x) {
            stringr::str_replace(
                x,
                "log2FoldChange_A",
                paste("log2FoldChange", NAME)
            )
        }) |>
        dplyr::rename_with(\(x) {
            stringr::str_replace(
                x,
                "log2FoldChange_B",
                paste("log2FoldChange", PARTNER_NAME)
            )
        }) |>
        dplyr::select(
            gene_id,
            tidyselect::contains("log2Fol"),
            Significant_In,
            Effect_Classification
        )

    return(genes_with_diff)
}


#' finds the rowData columns that can serve as gene identifiers: character (or
#' entrez) columns whose non-missing values are mostly unique. drops columns
#' like gene_biotype or seq_name, which are character but not identifiers
#' @noRd
.identifier_columns <- function(row_data, min_unique = 0.5) {
    row_data |>
        as.data.frame() |>
        dplyr::select(
            tidyselect::where(is.character),
            tidyselect::contains("entrez")
        ) |>
        dplyr::select(tidyselect::where(\(x) {
            x <- x[!is.na(x)]
            length(x) > 0 &&
                dplyr::n_distinct(x) > 1 &&
                dplyr::n_distinct(x) / length(x) >= min_unique
        })) |>
        names()
}

#' matches every context entry against all gene identifier columns and
#' translates the matched genes into the selected identifier, so contexts can
#' use any (or a mix of) identifiers. each entry is attributed to the first
#' column it is found in, so put the selected identifier first
#' @noRd
.match_context <- function(RES, CONTEXT, gene_var, GENE_VARS) {
    context <- unique(CONTEXT[!is.na(CONTEXT) & CONTEXT != ""])
    GENE_VARS <- intersect(GENE_VARS, colnames(RES))

    matched_via <- rep(NA_character_, length(context))
    for (v in rev(GENE_VARS)) {
        matched_via[context %in% RES[[v]]] <- v
    }

    # rows matched by any identifier, translated into the selected one
    rows <- purrr::map(GENE_VARS, \(v) RES[[v]] %in% context) |>
        purrr::reduce(`|`, .init = rep(FALSE, nrow(RES)))
    genes <- unique(RES[[gene_var]][rows])

    list(
        genes = genes[!is.na(genes)],
        context_size = length(context),
        matched_via = table(factor(matched_via, levels = GENE_VARS)),
        unmatched = context[is.na(matched_via)]
    )
}

#' splits the tested genes into up- and downregulated genes and restricts the
#' context to the tested genes. directions without DE genes are dropped
#' @noRd
.split_context_genes <- function(
    RES,
    CONTEXT,
    padj_CO = 0.05,
    fc_CO = 1,
    gene_var = "gene_id"
) {
    tested_genes <- RES |>
        dplyr::filter(
            !is.na(padj) & !is.na(log2FoldChange) & !is.infinite(log2FoldChange)
        ) |>
        dplyr::mutate(gene_id = .data[[gene_var]]) |>
        dplyr::filter(!is.na(gene_id))

    universe <- unique(tested_genes$gene_id)
    context_raw <- unique(CONTEXT[!is.na(CONTEXT)])

    de_genes <- list(
        "up regulated genes" = dplyr::filter(
            tested_genes,
            padj < padj_CO,
            log2FoldChange > fc_CO
        ),
        "down regulated genes" = dplyr::filter(
            tested_genes,
            padj < padj_CO,
            log2FoldChange < -fc_CO
        )
    ) |>
        purrr::map(\(x) unique(x$gene_id)) |>
        purrr::discard(\(x) length(x) == 0)

    list(
        universe = universe,
        context_raw = context_raw,
        context = intersect(context_raw, universe),
        de_genes = de_genes
    )
}

#' tests for overrepresentation of a gene set (CONTEXT) among up- and
#' down-regulated genes, using all tested genes as the universe.
#' directions without DE genes are dropped from the table. context_size
#' overrides the reported size of the context, e.g. the size of the list
#' before it was translated by .match_context()
#' @noRd
.check_context_enrichment <- function(
    RES,
    CONTEXT,
    padj_CO = 0.05,
    fc_CO = 1,
    gene_var = "gene_id",
    context_size = NULL
) {
    split <- .split_context_genes(RES, CONTEXT, padj_CO, fc_CO, gene_var)
    if (length(split$de_genes) == 0 || length(split$context) == 0) {
        return(data.frame())
    }

    n_universe <- length(split$universe)
    n_context <- length(split$context)
    purrr::imap(split$de_genes, \(genes, direction) {
        n_genes <- length(genes)
        n_overlap <- length(intersect(genes, split$context))
        data.frame(
            direction = direction,
            n_genes = n_genes,
            context_size = context_size %||% length(split$context_raw),
            context_size_tested = n_context,
            overlap = n_overlap,
            fold_enrichment = round(
                (n_overlap / n_genes) / (n_context / n_universe),
                2
            ),
            pvalue = stats::phyper(
                n_overlap - 1,
                n_context,
                n_universe - n_context,
                n_genes,
                lower.tail = FALSE
            )
        )
    }) |>
        dplyr::bind_rows()
}

#' returns the context genes found among the up- and downregulated genes,
#' one sorted vector per direction
#' @noRd
.get_context_overlap <- function(
    RES,
    CONTEXT,
    padj_CO = 0.05,
    fc_CO = 1,
    gene_var = "gene_id"
) {
    split <- .split_context_genes(RES, CONTEXT, padj_CO, fc_CO, gene_var)
    purrr::map(split$de_genes, \(genes) sort(intersect(genes, split$context)))
}

#' splits pasted text or the lines of a file into a clean gene vector.
#' splits on whitespace, commas and semicolons and strips surrounding quotes
#' @noRd
.parse_gene_list <- function(x) {
    if (length(x) == 0) {
        return(character())
    }
    genes <- x |>
        strsplit("[\\s,;]+", perl = TRUE) |>
        unlist() |>
        stringr::str_remove_all('^"+|"+$') |>
        trimws()
    unique(genes[nchar(genes) > 0])
}

#' reads an uploaded context file. A file with several lines whose first line
#' has several tab-, semicolon- or comma-separated fields is read as a table
#' with a header, giving one gene list per column named by its header.
#' Anything else is read as a single, unnamed gene list.
#' @noRd
.read_context_file <- function(path) {
    lines <- readLines(path, warn = FALSE)
    lines <- lines[nchar(trimws(lines)) > 0]

    sep <- NA_character_
    if (length(lines) >= 2) {
        sep <- dplyr::case_when(
            grepl("\t", lines[[1]]) ~ "\t",
            grepl(";", lines[[1]]) ~ ";",
            grepl(",", lines[[1]]) ~ ",",
            TRUE ~ NA_character_
        )
    }
    if (is.na(sep)) {
        return(list(.parse_gene_list(lines)))
    }

    tbl <- utils::read.table(
        text = lines,
        sep = sep,
        header = TRUE,
        quote = "\"",
        colClasses = "character",
        check.names = FALSE,
        na.strings = c("", "NA"),
        strip.white = TRUE,
        fill = TRUE,
        comment.char = ""
    )
    col_names <- names(tbl)
    col_names[col_names == ""] <- paste(
        "Column",
        which(col_names == "")
    )
    tbl |>
        as.list() |>
        purrr::set_names(col_names) |>
        purrr::map(\(x) unique(x[!is.na(x) & nchar(x) > 0])) |>
        purrr::discard(\(x) length(x) == 0)
}

#' checks the context passed to exploreSE(): a named list of gene vectors.
#' removes NA, empty and duplicated genes from each element
#' @noRd
.validate_context <- function(context) {
    if (is.null(context)) {
        return(NULL)
    }
    if (!is.list(context) || is.data.frame(context)) {
        stop("`context` must be a named list of gene vectors.", call. = FALSE)
    }
    ctx_names <- names(context)
    if (
        is.null(ctx_names) ||
            any(is.na(ctx_names) | trimws(ctx_names) == "")
    ) {
        stop("All elements of `context` must be named.", call. = FALSE)
    }
    if (anyDuplicated(ctx_names)) {
        stop("The names of `context` must be unique.", call. = FALSE)
    }
    if (!all(purrr::map_lgl(context, is.atomic))) {
        stop(
            "Every element of `context` must be a vector of gene identifiers.",
            call. = FALSE
        )
    }

    context <- purrr::map(context, \(x) {
        x <- trimws(as.character(x))
        unique(x[!is.na(x) & x != ""])
    })
    empty <- names(context)[lengths(context) == 0]
    if (length(empty) > 0) {
        stop(
            "These elements of `context` contain no genes: ",
            paste(empty, collapse = ", "),
            call. = FALSE
        )
    }
    context
}

#' the tabs presented on the About tab: the tab title (as used in
#' `main_tabs`), the thumbnail file and a short description in markdown
#' @noRd
.about_tabs <- function() {
    list(
        list(
            tab = "Data Overview",
            thumbnail = "01_overview.png",
            help = "Displays an overview of the data loaded into the app, with a focus on its metadata. Here you will find additional information on your samples, for example which patient a sample was taken from or whether it was treated."
        ),
        list(
            tab = "PCA",
            thumbnail = "02_pca.png",
            help = "Calculated on your data, the PCA shows the spread of your samples. Each sample is coloured by the \"Color/Group by:\" option set on the left."
        ),
        list(
            tab = "Gene Expression",
            thumbnail = "03_gene_expression.png",
            help = "Check the expression of a gene of interest in your data. The expression is shown for each level of the selected grouping variable, and you can select/deselect certain groups. Different gene identifiers can be selected on the left; gene symbols are the most human-readable. Below the plot, the expression values are shown for each sample, and below that, the differential expression status of the gene in the different comparisons."
        ),
        list(
            tab = "DE Results",
            thumbnail = "04_de_results.png",
            help = "A big-picture overview of the differential expression result selected on the left. A bar chart displays the number of differentially expressed genes, and below it is a table with the test statistics per gene."
        ),
        list(
            tab = "Volcano Plot",
            thumbnail = "05_volcano_plot.png",
            help = "Explore your differential expression results as a volcano plot. You can select a number of top genes (by p-value) to label, highlight genes of interest or just hover over a point to get the underlying gene, p-value and fold change."
        ),
        list(
            tab = "Enrichment Results",
            thumbnail = "06_enrichment_results.png",
            help = "Shows the enrichment results for the selected DE result. The genes driving each enrichment can be shown as well, using the checkbox at the top of the tab."
        ),
        list(
            tab = "FC FC",
            thumbnail = "07_fc_fc.png",
            help = "Compare DE results to one another, using significance and fold change. Below the plot, you get a table showing all genes that are significant in at least one of the results; you can filter this list by a minimum difference in fold change between the two results in question."
        ),
        list(
            tab = "Context",
            thumbnail = "08_context.png",
            help = "It can help to put your results into context, for example with known genetic risk factors. The genes in a context are overlaid on the selected result, and their overrepresentation is tested. You can select a context from the list or add one yourself, either by pasting a list of genes or by uploading a table where each column is a context."
        )
    )
}

#' builds the grid of tab thumbnails for the About tab. Clicking a thumbnail
#' or tab name sets `input$about_goto`, clicking the question mark sets
#' `input$about_help`, both to the tab title
#' @noRd
.about_tab_cards <- function() {
    cards <- purrr::map(.about_tabs(), \(x) {
        shiny::column(
            3,
            htmltools::div(
                class = "about-card",
                htmltools::tags$a(
                    href = "#",
                    class = "about-goto",
                    `data-tab` = x$tab,
                    htmltools::tags$img(
                        src = paste0("exploreSE/thumbnails/", x$thumbnail),
                        alt = paste("Preview of the", x$tab, "tab")
                    )
                ),
                htmltools::div(
                    class = "about-card-title",
                    htmltools::tags$a(
                        href = "#",
                        class = "about-goto",
                        `data-tab` = x$tab,
                        x$tab
                    ),
                    htmltools::tags$a(
                        href = "#",
                        class = "about-help",
                        `data-tab` = x$tab,
                        title = paste("More about the", x$tab, "tab"),
                        `aria-label` = paste("More about the", x$tab, "tab"),
                        shiny::icon("circle-question")
                    )
                )
            )
        )
    })

    htmltools::tagList(
        htmltools::tags$style(htmltools::HTML(
            "
            .about-card { text-align: center; margin-bottom: 24px; }
            .about-card img {
                width: 150px; height: 150px; max-width: 100%;
                border-radius: 16px;
                transition: box-shadow 0.15s;
            }
            .about-card a.about-goto:hover img,
            .about-card a.about-goto:focus img {
                box-shadow: 0 2px 10px rgba(0, 0, 0, 0.25);
            }
            .about-card-title { margin-top: 8px; font-weight: 600; }
            .about-help { margin-left: 6px; color: #777; }
            .about-help:hover, .about-help:focus { color: #333; }
            "
        )),
        htmltools::tags$script(htmltools::HTML(
            "
            $(document).on('click', '.about-goto, .about-help', function(e) {
              e.preventDefault();
              var id = $(this).hasClass('about-help') ? 'about_help' : 'about_goto';
              Shiny.setInputValue(id, $(this).data('tab'), {priority: 'event'});
            });
            "
        )),
        shiny::fluidRow(cards)
    )
}
