# UI
ui <- function(request) {
  shiny::fluidPage(
    shiny::titlePanel("RNA-seq SummarizedExperiment Explorer"),
    # sidebar ----------
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        width = 3,
        shiny::h4("Data Input"),
        shiny::fileInput(
          "se_file",
          "Upload SummarizedExperiment (.rds)",
          accept = ".rds"
        ),
        shiny::checkboxInput("use_demo", "Use Demo Data", value = TRUE),
        shiny::hr(),
        shiny::h4("Analysis Options"),

        shiny::conditionalPanel(
          condition = "input.main_tabs == 'PCA' | input.main_tabs == 'Gene Expression'",
          shiny::selectInput(
            "color_var",
            "Color/Group by:",
            choices = NULL
          ),
        ),

        shiny::conditionalPanel(
          condition = "input.main_tabs == 'Gene Expression' | input.main_tabs == 'Volcano Plot' | input.main_tabs == 'FC FC'",

          shiny::selectInput(
            "row_data_var",
            "Select the gene identifier",
            choices = NULL
          ),
        ),

        shiny::conditionalPanel(
          condition = "input.main_tabs == 'DE Results' | input.main_tabs == 'Volcano Plot' | input.main_tabs == 'Enrichment Results' | input.main_tabs == 'FC FC'",
          shiny::conditionalPanel(
            condition = "output.data_loaded",
            shiny::conditionalPanel(
              condition = "output.has_precomputed_de",
              shiny::selectInput(
                "de_comparison",
                "Select Comparison:",
                choices = NULL
              )
            ),
          )
        ),

        shiny::conditionalPanel(
          condition = "input.main_tabs == 'DE Results' | input.main_tabs == 'Volcano Plot' | input.main_tabs == 'Enrichment Results'",
          colourpicker::colourInput(
            "up_col",
            "Color for Upregulated",
            "#d62728"
          ),
          colourpicker::colourInput(
            "dn_col",
            "Color for Downregulated",
            "#1f77b4"
          ),
          colourpicker::colourInput(
            "high_col",
            "Color for Highlights",
            "#FFD700"
          )
        ),
        shiny::conditionalPanel(
          condition = "input.main_tabs == 'FC FC'",
          colourpicker::colourInput(
            "both_col",
            "Color for DE Overlap",
            "#ce1ae6"
          ),
          colourpicker::colourInput(
            "comp1_col",
            "Color for DE Comparison 1",
            "#0a31db"
          ),
          colourpicker::colourInput(
            "comp2_col",
            "Color for DE Comparison 2",
            "#f03405"
          ),
        )
      ),

      shiny::mainPanel(
        width = 9,
        shiny::tabsetPanel(
          id = "main_tabs",

          # tab overview ------
          shiny::tabPanel(
            "Overview",
            shiny::h3("Dataset Summary"),
            shiny::verbatimTextOutput("data_summary"),
            shiny::hr(),
            shiny::h4("Sample Metadata"),
            DT::DTOutput("metadata_table")
          ),
          # tab PCA ----
          shiny::tabPanel(
            "PCA",
            shiny::h3("Principal Component Analysis"),
            shiny::fluidRow(
              shiny::column(
                3,
                shiny::numericInput(
                  "top_genes",
                  "Top variable genes for PCA:",
                  value = 500,
                  min = 100,
                  max = 5000,
                  step = 100
                ),
                shiny::hr(),
              )
            ),
            plotly::plotlyOutput("pca_plot", height = "600px"),
            shiny::hr(),
            shiny::verbatimTextOutput("pca_variance")
          ),
          # tab gene expression ----
          shiny::tabPanel(
            "Gene Expression",
            shiny::h3("Gene Expression Plot"),
            shiny::fluidRow(
              shiny::column(
                3,
                shiny::checkboxGroupInput(
                  "groups_to_show",
                  "Include levels:"
                )
              ),
              shiny::column(
                3,
                shinyWidgets::pickerInput(
                  "gene_id",
                  "Select Gene:",
                  choices = NULL,
                  options = list(
                    `live-search` = TRUE,
                    `live-search-placeholder` = "Search genes...",
                    size = 10
                  )
                ),
                shiny::selectInput(
                  "plot_type",
                  "Plot Type:",
                  choices = c(
                    "Boxplot" = "box",
                    "Violin" = "violin"
                  )
                )
              )
            ),
            plotly::plotlyOutput("expr_plot", height = "500px"),
            shiny::hr(),
            shiny::h4("Expression Values"),
            DT::DTOutput("expr_table"),
            shiny::downloadButton(
              "download_expr",
              "Download Expression Table"
            )
          ),
          # tab de res ----------
          shiny::tabPanel(
            "DE Results",
            shiny::h3("Differential Expression Results"),
            shiny::uiOutput("de_status_message"),
            shiny::fluidRow(),
            shiny::conditionalPanel(
              condition = "!output.has_precomputed_de",
              shiny::actionButton(
                "run_de",
                "Run Basic DE Analysis",
                class = "btn-primary"
              )
            ),
            shiny::hr(),
            shiny::plotOutput("de_plot", height = "700px"),
            shiny::hr(),
            DT::DTOutput("de_table"),
            shiny::downloadButton("download_de", "Download DE Table")
          ),

          # tab volcano --------
          shiny::tabPanel(
            "Volcano Plot",
            shiny::h3("Volcano Plot"),
            shiny::fluidRow(
              shiny::column(
                3,
                shiny::numericInput(
                  "padj_cutoff_volcano",
                  "Adjusted p-value cutoff:",
                  value = 0.05,
                  min = 0,
                  max = 1,
                  step = 0.01
                )
              ),
              shiny::column(
                3,
                shiny::numericInput(
                  "lfc_cutoff_volcano",
                  "Log2 Fold Change cutoff:",
                  value = 1,
                  min = 0,
                  max = 10,
                  step = 0.5
                )
              ),
              shiny::column(
                3,
                shiny::checkboxInput(
                  "label_top",
                  "Label top genes",
                  value = TRUE
                ),
                shiny::numericInput(
                  "n_labels",
                  "Number to label:",
                  value = 10,
                  min = 0,
                  max = 50,
                  step = 5
                )
              )
            ),
            shiny::fluidRow(
              shiny::column(
                12,
                shiny::textAreaInput(
                  "highlight_genes",
                  "Highlight specific genes (one per line or comma-separated):",
                  value = "",
                  placeholder = "GENE1, GENE2, GENE3\nor\nGENE1\nGENE2\nGENE3",
                  rows = 3,
                  width = "100%"
                ),
                shiny::helpText(
                  "Enter gene IDs or gene names to highlight in yellow on the plot."
                )
              )
            ),
            shiny::hr(),
            plotly::plotlyOutput("volcano_plot", height = "700px"),
            shiny::hr(),
            shiny::h4("Summary Statistics")
          ),
          # tab enrichment ---
          shiny::tabPanel(
            "Enrichment Results",
            shiny::h3("Enrichment Results"),
            shiny::fluidRow(
              shiny::column(
                3,
                shiny::numericInput(
                  "padj_cutoff_enrichment",
                  "Adjusted p-value cutoff:",
                  value = 0.05,
                  min = 0,
                  max = 1,
                  step = 0.01
                )
              ),
              shiny::column(
                3,
                shiny::numericInput(
                  "n_terms_enrichment",
                  "Number of top terms to show:",
                  value = 10,
                  min = 3,
                  max = 20,
                  step = 1
                )
              )
            ),
            shiny::uiOutput("fe_status_message"),
            shiny::hr(),
            shiny::uiOutput("enrichment_plots", width = "700px")
          ),
          # tab fcfc ------
          shiny::tabPanel(
            "FC FC",
            shiny::h3("Fold Change comparisons"),
            shiny::fluidRow(),
            shiny::selectInput(
              "comparison_comparison",
              "Select Comparison:",
              choices = NULL
            ),
            shiny::hr(),
            plotly::plotlyOutput("fcfc_plot", height = "700px")
          )
        )
      )
    )
  )
}
# nocov start
server <- function(input, output, session) {
  # Increase upload size limit to 500MB
  options(shiny.maxRequestSize = 500 * 1024^2)

  # Reactive values --------
  rv <- shiny::reactiveValues(
    se = NULL,
    vst_data = NULL,
    pca_result = NULL,
    de_results = NULL,
    # up_col = "#d62728",
    # dn_col = "#1f77b4",
    # highlight_col = "#FFD700",
    # both_col = "#ce1ae6",
    # comp1_col = "#0a31db",
    # comp2_col = "#f03405",
    # color_var = NULL,
    # row_var = NULL,
    gene_idents = NULL,
    # second_comparison = NULL
  )

  has_precomputed_de <- shiny::reactive({
    shiny::req(rv$se)
    .check_precomputed_de(rv$se)
  })

  has_precomputed_fe <- shiny::reactive({
    shiny::req(rv$se)
    .check_precomputed_fe(rv$se)
  })

  # run observers --------------
  # .create_dir_color_observers(input, session, rv)
  .create_interest_color_observers(input, session, rv)
  .observe_demo_data(input, session, rv)
  .observe_inital_obj(input, session, rv)
  .observe_load_file(input, session, rv)
  .create_rowvar_observer(input, session, rv)
  # precomputed results ---------

  # Get precomputed DE comparisons
  de_comparisons <- shiny::reactive({
    shiny::req(has_precomputed_de())
    .de_results_names(rv$se)
  })

  fes <- shiny::reactive({
    shiny::req(has_precomputed_fe())
    .fe_results_names(rv$se, input$de_comparison)
  })

  .observe_se_load(input, session, rv, has_precomputed_de, de_comparisons)

  # VST transformation -------
  vst_data <- shiny::reactive({
    shiny::req(rv$se)

    if (is.null(rv$vst_data)) {
      shiny::withProgress(message = "Transforming data...", {
        dds <- DESeq2::DESeqDataSet(rv$se, design = ~1)
        rv$vst_data <- SummarizedExperiment::assay(DESeq2::vst(
          dds,
          blind = TRUE
        ))
      })
    }
    rv$vst_data
  })

  # PCA calculation
  pca_data <- shiny::reactive({
    shiny::req(vst_data(), input$color_var, input$top_genes)

    vst_mat <- vst_data()

    # Select top variable genes
    rv_genes <- matrixStats::rowVars(vst_mat)
    select_genes <- order(rv_genes, decreasing = TRUE)[
      seq_len(min(input$top_genes, nrow(vst_mat)))
    ]

    # Run PCA
    pca <- stats::prcomp(t(vst_mat[select_genes, ]), scale. = FALSE)
    rv$pca_result <- pca

    # Create data frame for plotting
    pca_df <- data.frame(
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2],
      sample = colnames(rv$se),
      group = SummarizedExperiment::colData(rv$se)[[input$color_var]]
    )

    # Calculate variance explained
    var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

    list(data = pca_df, var = var_explained)
  })

  # Run DE analysis ----------
  shiny::observeEvent(input$run_de, {
    shiny::req(rv$se, input$color_var)

    shiny::withProgress(message = "Running DESeq2...", {
      tryCatch(
        {
          # Create DESeq2 object
          col_var <- input$color_var
          design_formula <- stats::as.formula(paste("~", col_var))
          dds <- DESeq2::DESeqDataSet(rv$se, design = design_formula)

          # Filter low counts
          keep <- rowSums(BiocGenerics::counts(dds)) >= 10
          dds <- dds[keep, ]

          # Run DESeq2
          dds <- DESeq2::DESeq(dds)

          # Get results
          res <- DESeq2::results(dds)
          rv$de_results <- as.data.frame(res) %>%
            tibble::rownames_to_column("gene_id") %>%
            dplyr::arrange(padj)

          shiny::showNotification(
            "DE analysis complete!",
            type = "message"
          )
        },
        error = function(e) {
          shiny::showNotification(
            paste("Error in DE analysis:", e$message),
            type = "error"
          )
        }
      )
    })
  })

  # current DE results ---------
  current_de_results <- shiny::reactive({
    if (has_precomputed_de() && !is.null(input$de_comparison)) {
      # Load precomputed results

      de_res <- .de_result(rv$se, input$de_comparison) %>%
        tibble::rownames_to_column("gene_ident")

      merging_data <- SummarizedExperiment::rowData(rv$se) |>
        as.data.frame() |>
        dplyr::select(tidyselect::any_of(rv$gene_idents)) |>
        tibble::rownames_to_column("gene_ident")
      de_res <- merging_data |>
        dplyr::left_join(de_res, by = dplyr::join_by(gene_ident))
      return(de_res)
    } else if (!is.null(rv$de_results)) {
      # Use computed results
      return(rv$de_results)
    } else {
      return(NULL)
    }
  })

  # current FE results ------
  current_fe_results <- shiny::reactive({
    if (has_precomputed_fe()) {
      # Load precomputed results
      if (methods::is(rv$se, "DeeDeeExperiment")) {
        fe_res <- .fe_result(rv$se, input$de_comparison)
      } else {
        fe_res <- S4Vectors::metadata(rv$se)$fe_results[[
          input$de_comparison
        ]]
      }
      return(fe_res)
    } else {
      return(NULL)
    }
  })

  # Outputs ---------------
  output$data_loaded <- shiny::reactive({
    !is.null(rv$se)
  })
  shiny::outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)

  output$has_precomputed_de <- shiny::reactive({
    has_precomputed_de()
  })
  output$has_precomputed_fe <- shiny::reactive({
    has_precomputed_fe()
  })
  shiny::outputOptions(
    output,
    "has_precomputed_de",
    suspendWhenHidden = FALSE
  )
  shiny::outputOptions(
    output,
    "has_precomputed_fe",
    suspendWhenHidden = FALSE
  )

  output$de_status_message <- shiny::renderUI({
    if (has_precomputed_de()) {
      comparisons <- de_comparisons()
      htmltools::tagList(
        htmltools::p(
          shiny::icon("check-circle", class = "text-success"),
          htmltools::strong(paste(
            "Found",
            length(comparisons),
            "precomputed DE comparison(s):"
          )),
          htmltools::br(),
          paste(comparisons, collapse = ", ")
        )
      )
    } else {
      htmltools::p(
        shiny::icon("info-circle"),
        "No precomputed DE results found. Run basic DE analysis or upload data with results in metadata(se)$de_results."
      )
    }
  })
  output$fe_status_message <- shiny::renderUI({
    if (has_precomputed_fe()) {
      comparisons <- names(current_fe_results())
      htmltools::tagList(
        htmltools::p(
          shiny::icon("check-circle", class = "text-success"),
          htmltools::strong(paste(
            "Found",
            length(comparisons),
            "precomputed functional enrichment(s):"
          )),
          htmltools::br(),
          paste(comparisons, collapse = ", ")
        )
      )
    } else {
      htmltools::p(
        shiny::icon("info-circle"),
        "No precomputed FE results found. Run basic FE analysis or upload data with results in metadata(se)$fe_results."
      )
    }
  })

  output$data_summary <- shiny::renderText({
    shiny::req(rv$se)

    assay_names <- paste(
      names(SummarizedExperiment::assays(rv$se)),
      collapse = ", "
    )
    sample_metadata_cols <- paste(
      colnames(SummarizedExperiment::colData(rv$se)),
      collapse = ", "
    )

    summary_lines <- c(
      "SummarizedExperiment Object",
      "===========================",
      "",
      "Dimensions:",
      paste("  Genes:", nrow(rv$se)),
      paste("  Samples:", ncol(rv$se)),
      "",
      paste("Assays:", assay_names),
      "",
      "Sample Metadata Columns:",
      paste(" ", sample_metadata_cols)
    )

    if (ncol(SummarizedExperiment::rowData(rv$se)) > 0) {
      gene_metadata_cols <- paste(
        colnames(SummarizedExperiment::rowData(rv$se)),
        collapse = ", "
      )
      summary_lines <- c(
        summary_lines,
        "",
        "Gene Metadata Columns:",
        paste(" ", gene_metadata_cols)
      )
    }

    paste(summary_lines, collapse = "\n")
  })

  output$metadata_table <- DT::renderDT({
    shiny::req(rv$se)
    DT::datatable(
      as.data.frame(SummarizedExperiment::colData(rv$se)),
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = TRUE
    )
  })

  output$pca_plot <- plotly::renderPlotly({
    shiny::req(pca_data())

    pca_info <- pca_data()
    df <- pca_info$data
    var <- pca_info$var

    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = PC1, y = PC2, color = group, text = sample)
    ) +
      ggplot2::geom_point(size = 4, alpha = 0.8) +
      ggplot2::labs(
        x = paste0("PC1 (", var[1], "%)"),
        y = paste0("PC2 (", var[2], "%)"),
        color = input$color_var
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(legend.position = "right")

    plotly::ggplotly(p, tooltip = c("text", "group")) %>%
      plotly::layout(hovermode = "closest")
  })

  output$pca_variance <- shiny::renderText({
    shiny::req(rv$pca_result)
    var_explained <- round(
      100 * rv$pca_result$sdev^2 / sum(rv$pca_result$sdev^2),
      2
    )
    top_pcs <- seq_len(min(10, length(var_explained)))
    paste(
      c(
        "Variance Explained by PCs:",
        sprintf("  PC%d: %.2f%%", top_pcs, var_explained[top_pcs])
      ),
      collapse = "\n"
    )
  })
  ## expression plot-------
  output$expr_plot <- plotly::renderPlotly({
    shiny::req(rv$se, input$gene_id, input$color_var)

    p <- .plot_expression(
      rv$se,
      input$gene_id,
      input$color_var,
      input$groups_to_show,
      input$plot_type,
      input$row_data_var
    )
    plotly::ggplotly(p)
  })

  output$expr_table <- DT::renderDT({
    shiny::req(
      rv$se,
      input$gene_id,
      input$groups_to_show,
      input$row_data_var
    )

    gene <- rownames(SummarizedExperiment::rowData(rv$se)[
      which(
        SummarizedExperiment::rowData(rv$se)[, input$row_data_var] ==
          input$gene_id
      ),
    ])[1]
    expr_data <- data.frame(
      Sample = colnames(rv$se),
      SummarizedExperiment::colData(rv$se),
      Gene = input$gene_id,
      Count = SummarizedExperiment::assay(rv$se, "counts")[gene, ],
      g_r_o_u_p = SummarizedExperiment::colData(rv$se)[[input$color_var]]
    ) %>%
      dplyr::filter(g_r_o_u_p %in% input$groups_to_show) %>%
      dplyr::select(-g_r_o_u_p)

    DT::datatable(
      expr_data,
      options = list(pageLength = 12, scrollX = TRUE),
      rownames = FALSE
    ) %>%
      DT::formatRound("Count", digits = 0)
  })

  output$de_table <- DT::renderDT({
    de_data <- current_de_results()
    de_data <- dplyr::select(
      de_data,
      tidyselect::any_of(c(
        "gene_ident",
        rv$gene_idents,
        "baseMean",
        "log2FoldChange",
        "pvalue",
        "padj"
      ))
    ) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::any_of(c("baseMean", "log2FoldChange")),
          \(x) round(x, 2)
        ),
        dplyr::across(tidyselect::any_of(c("pvalue", "padj")), \(x) {
          round(x, 4)
        })
      )
    shiny::req(de_data)

    DT::datatable(
      de_data,
      options = list(pageLength = 25, scrollX = TRUE),
      rownames = FALSE,
      filter = "top"
    )
  })

  output$de_plot <- shiny::renderPlot({
    colors_acute <- c(
      "Up" = input$up_col,
      "Down" = input$dn_col,
      "NS" = "grey70",
      "Highlighted" = input$highlight_col
    )
    de_data <- current_de_results()

    .plot_des(
      de_data,
      input$de_comparison,
      input$padj_cutoff_volcano,
      input$lfc_cutoff_volcano,
      colors_acute
    )
  })

  output$download_expr <- shiny::downloadHandler(
    filename = function() {
      comparison_name <- stringr::str_replace(
        input$gene_id,
        "[^a-zA-Z0-9_-]",
        "_"
      )

      paste0(comparison_name, "_expression_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      shiny::req(rv$se, input$gene_id, input$groups_to_show)

      gene <- input$gene_id
      expr_data <- data.frame(
        Sample = colnames(rv$se),
        SummarizedExperiment::colData(rv$se),
        Gene = input$gene_id,
        Count = SummarizedExperiment::assay(rv$se, "counts")[gene, ],
        g_r_o_u_p = SummarizedExperiment::colData(rv$se)[[
          input$color_var
        ]]
      ) %>%
        dplyr::filter(g_r_o_u_p %in% input$groups_to_show) %>%
        dplyr::select(-g_r_o_u_p)
      readr::write_excel_csv2(expr_data, file)
    }
  )

  output$download_de <- shiny::downloadHandler(
    filename = function() {
      comparison_name <- if (
        has_precomputed_de() && !is.null(input$de_comparison)
      ) {
        gsub("[^a-zA-Z0-9_-]", "_", input$de_comparison)
      } else {
        "DE_results"
      }
      paste0(comparison_name, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      de_data <- current_de_results()
      shiny::req(de_data)
      readr::write_excel_csv2(de_data, file)
    }
  )

  ## Volcano plot ----
  output$volcano_plot <- plotly::renderPlotly({
    de_data <- current_de_results()
    shiny::req(
      de_data,
      # input$padj_cutoff_volcano,
      # input$lfc_cutoff_volcano,
      # input$up_col,
      # input$dn_col,
      # input$highlight_col,
      input$row_data_var
    )
    colors_acute <- c(
      "Up" = input$up_col,
      "Down" = input$dn_col,
      "NS" = "grey70",
      "Highlighted" = input$highlight_col
    )

    highlight_vec <- c()
    if (
      !is.null(input$highlight_genes) &&
        nchar(trimws(input$highlight_genes)) > 0
    ) {
      # Split by newlines and commas, trim whitespace
      highlight_vec <- input$highlight_genes %>%
        strsplit("[\n,]") %>%
        unlist() %>%
        trimws() %>%
        .[nchar(.) > 0]
    }

    .plot_volcano(
      RES = de_data,
      NAME = input$de_comparison,
      padj_CO = input$padj_cutoff_volcano,
      fc_CO = input$lfc_cutoff_volcano,
      highlights = highlight_vec,
      COLS = colors_acute,
      LABEL_TOP = input$label_top,
      TOPN = input$n_labels,
      gene_var = input$row_data_var
    )
  })

  ## enrichment plots -------
  output$enrichment_plots <- shiny::renderUI({
    current_fes <- current_fe_results()
    colors_acute <- c(
      "Up" = input$up_col,
      "Down" = input$dn_col,
      "NS" = "grey70",
      "Highlighted" = input$highlight_col
    )
    shiny::req(
      current_fes,
      input$padj_cutoff_enrichment,
      input$n_terms_enrichment,
      colors_acute
    )

    plots <- purrr::imap(current_fes, function(enrich, name) {
      id <- paste0("e_plot_", name)

      output[[id]] <- plotly::renderPlotly({
        .plot_fe(
          FE = enrich,
          NAME = name,
          padj_CO = input$padj_cutoff_enrichment,
          N_terms = input$n_terms_enrichment,
          COLS = colors_acute
        )
      })
    })
  })

  output$fcfc_plot <- plotly::renderPlotly({
    shiny::req(
      rv$se,
      input$padj_cutoff_volcano,
      input$lfc_cutoff_volcano,
      input$comp1_col,
      input$comp2_col
    )
    colors_acute <- c(
      "Comp1" = input$comp1_col,
      "Comp2" = input$comp2_col,
      "ns" = "grey70",
      "both" = input$both_col
    )

    #     highlight_vec <- c()
    #     if (
    #         !is.null(input$highlight_genes) &&
    #             nchar(trimws(input$highlight_genes)) > 0
    #     ) {
    #         # Split by newlines and commas, trim whitespace
    #         highlight_vec <- input$highlight_genes %>%
    #             strsplit("[\n,]") %>%
    #             unlist() %>%
    #             trimws() %>%
    #             .[nchar(.) > 0]
    #     }

    .plot_fcfc(
      OBJ = rv$se,
      NAME = input$de_comparison,
      PARTNER_NAME = input$comparison_comparison,
      GENE_VAR = input$row_data_var,
      GENE_VARS = rv$gene_idents,
      padj_CO = input$padj_cutoff_volcano,
      fc_CO = input$lfc_cutoff_volcano,
      COLS = colors_acute
    )
  })
}
# nocov end

# Run app
#' exploreSE
#'
#'
#' @description
#' This runs the explorer app. You can add a path to a .rds file or a the
#' relevant R object as an argument to load it on start-up.
#' The app looks for DE results in the metadata of the metadata(object)$de_results
#' and, if it is a DeeDeeExperiment, in the DEA slots of the object.
#' Similarly, the enrichment results are looked for in the metadata(object)$fe_results
#' or FEA slots (DeeDeeExperiment). As the DeeDeeExperiment structure does
#' a lot of work to clean and standardize the inputs and outputs for both DE
#' and enrichment results, using it is highly recommended.
#'
#' @param file a string leading to a file you want loaded as a default
#' @param object a summarizedExperiment object (or DeeDeeExperiment) to load in
#'
#' @returns an application
#' @export
#'
#' @examples
#' app <- exploreSE()
#' if(interactive()){
#' shiny::runApp(app, port = 1234)
#'}

exploreSE <- function(file = NULL, object = NULL) {
  shiny::shinyOptions(set_file_name = file, set_object = object)
  shiny::shinyApp(ui, server)
}
