.create_dir_color_observers <- function(input, session, rv) {
    # nocov start
    shiny::observeEvent(input$up_col, {
        rv$up_col <- input$up_col
    })
    shiny::observeEvent(input$dn_col, {
        rv$dn_col <- input$dn_col
    })
    shiny::observeEvent(input$high_col, {
        rv$highlight_col <- input$high_col
    })
    # nocov end
}

.create_interest_color_observers <- function(input, session, rv) {
    # nocov start
    shiny::observe({
        if (!is.null(rv$se)) {
            shiny::updateCheckboxGroupInput(
                session,
                "groups_to_show",
                choices = levels(SummarizedExperiment::colData(rv$se)[[
                    rv$color_var
                ]]),
                selected = levels(SummarizedExperiment::colData(rv$se)[[
                    rv$color_var
                ]])
            )
        }
    })

    # nocov end

    # nocov start
    shiny::observeEvent(input$color_var, {
        rv$color_var <- input$color_var
    })
    # nocov end
}

.create_rowvar_observer <- function(input, session, rv) {
    # nocov start
    shiny::observeEvent(input$row_data_var, {
        shiny::req(input$row_data_var)
        rv$row_var <- input$row_data_var
    })
    # nocov end

    # nocov start

    shiny::observeEvent(rv$row_var, {
        if (
            !is.null(rv$se) &&
                !is.null(SummarizedExperiment::rowData(rv$se)) &&
                rv$row_var %in% colnames(SummarizedExperiment::rowData(rv$se))
        ) {
            gene_choices <- SummarizedExperiment::rowData(rv$se)[,
                rv$row_var,
                drop = TRUE
            ]
            # if ("gene_name" %in% colnames(SummarizedExperiment::rowData(rv$se))) {
            #   names(gene_choices) <- SummarizedExperiment::rowData(rv$se)$gene_name
            # }

            shinyWidgets::updatePickerInput(
                session,
                "gene_id",
                choices = gene_choices,
                selected = gene_choices[1]
            )
        }
    })
    # nocov end
}

# .gene_ident_observer <- function(input, session, rv) {}

.observe_demo_data <- function(input, session, rv) {
    # nocov start

    shiny::observe({
        if (input$use_demo) {
            rv$se <- .create_demo_data()
        }
    })
    # nocov end
}


.observe_inital_obj <- function(input, session, rv) {
    set_file <- shiny::getShinyOption("set_file_name")
    set_object <- shiny::getShinyOption("set_object")

    # nocov start

    shiny::observe({
        shiny::req(set_file)
        rv$se <- readRDS(set_file)
        shiny::updateCheckboxInput(session, "use_demo", value = FALSE)
        shiny::showNotification("Data loaded successfully!", type = "message")
    })

    shiny::observe({
        shiny::req(set_object)
        rv$se <- set_object
        shiny::updateCheckboxInput(session, "use_demo", value = FALSE)
        shiny::showNotification("Data loaded successfully!", type = "message")
    })

    # nocov end
}


.observe_se_load <- function(
    input,
    session,
    rv,
    has_precomputed_de,
    de_comparisons
) {
    # nocov start
    shiny::observeEvent(rv$se, {
        shiny::req(rv$se)

        # Update grouping variable choices
        col_vars <- colnames(SummarizedExperiment::colData(rv$se))
        row_vars <- colnames(SummarizedExperiment::rowData(rv$se))

        # Set default to "condition" if it exists, otherwise first column
        default_col_var <- if ("condition" %in% col_vars) {
            "condition"
        } else {
            col_vars[1]
        }

        shiny::updateSelectInput(
            session,
            "color_var",
            choices = col_vars,
            selected = default_col_var
        )

        default_row_var <- if ("gene_name" %in% row_vars) {
            "gene_name"
        } else if ("SYMBOL" %in% stringr::str_to_upper(row_vars)) {
            row_vars[which(stringr::str_to_upper(row_vars) == "SYMBOL")[1]]
        } else {
            row_vars[1]
        }
        shiny::updateSelectInput(
            session,
            "row_data_var",
            choices = row_vars,
            selected = default_row_var
        )

        # Update gene choices with searchable picker
        gene_choices <- SummarizedExperiment::rowData(rv$se)[,
            default_row_var,
            drop = TRUE
        ]
        # if ("gene_name" %in% colnames(SummarizedExperiment::rowData(rv$se))) {
        #   names(gene_choices) <- SummarizedExperiment::rowData(rv$se)$gene_name
        # }

        shinyWidgets::updatePickerInput(
            session,
            "gene_id",
            choices = gene_choices,
            selected = gene_choices[1]
        )

        # Update DE comparison choices if precomputed results exist
        if (has_precomputed_de()) {
            comparisons <- de_comparisons()
            shiny::updateSelectInput(
                session,
                "de_comparison",
                choices = comparisons,
                selected = comparisons[1]
            )
            shiny::showNotification(
                paste(
                    "Found",
                    length(comparisons),
                    "precomputed DE comparison(s)"
                ),
                type = "message"
            )
        }
    })
    # nocov end
}


.observe_load_file <- function(input, session, rv) {
    # nocov start
    shiny::observeEvent(input$se_file, {
        shiny::req(input$se_file)
        tryCatch(
            {
                rv$se <- readRDS(input$se_file$datapath)
                shiny::updateCheckboxInput(session, "use_demo", value = FALSE)
                shiny::showNotification(
                    "Data loaded successfully!",
                    type = "message"
                )
            },
            error = function(e) {
                shiny::showNotification(
                    paste("Error loading file:", e$message),
                    type = "error"
                )
            }
        )
    })
    # nocov end
}
