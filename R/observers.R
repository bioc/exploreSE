.create_interest_color_observers <- function(input, session, rv) {
  # nocov start
  shiny::observe({
    if (!is.null(rv$se)) {
      shiny::updateCheckboxGroupInput(
        session,
        "groups_to_show",
        choices = levels(SummarizedExperiment::colData(rv$se)[[
          input$color_var
        ]]),
        selected = levels(SummarizedExperiment::colData(rv$se)[[
          input$color_var
        ]])
      )
    }
  })

  # nocov end
}

.create_rowvar_observer <- function(input, session, rv) {
  # nocov start

  shiny::observeEvent(input$row_data_var, {
    if (
      !is.null(rv$se) &&
        !is.null(SummarizedExperiment::rowData(rv$se)) &&
        input$row_data_var %in%
          colnames(SummarizedExperiment::rowData(rv$se))
    ) {
      gene_choices <- SummarizedExperiment::rowData(rv$se)[,
        input$row_data_var,
        drop = TRUE
      ] |>
        unname()

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

    col_vars <- colnames(SummarizedExperiment::colData(rv$se))

    row_vars <- SummarizedExperiment::rowData(rv$se) |>
      as.data.frame() |>
      dplyr::select(
        tidyselect::where(is.character),
        # tidyselect::contains("id"),
        tidyselect::contains("entrez")
      ) |>
      dplyr::select(tidyselect::where(\(x) {
        !all(is.na(x)) & dplyr::n_distinct(x, na.rm = TRUE) > 1
      })) |>
      names()
    rv$gene_idents <- row_vars

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
      # } else if ("SYMBOL" %in% stringr::str_to_upper(row_vars)) {
      #   row_vars[which(stringr::str_to_upper(row_vars) == "SYMBOL")[1]]
    } else {
      row_vars[1]
    }

    shiny::updateSelectInput(
      session,
      "row_data_var",
      choices = row_vars,
      selected = default_row_var
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
      #for fcfc
      shiny::updateSelectInput(
        session,
        "comparison_comparison",
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


.observe_context <- function(input, session, rv) {
  # nocov start
  # contexts passed to exploreSE(context = ...)
  set_context <- shiny::getShinyOption("set_context")
  if (!is.null(set_context)) {
    rv$context_lists <- set_context
    shiny::updateSelectInput(
      session,
      "context_item",
      choices = names(set_context),
      selected = names(set_context)[1]
    )
  }

  shiny::observeEvent(input$add_context, {
    shiny::showModal(shiny::modalDialog(
      title = "Add context",
      shiny::textInput(
        "context_name",
        "Name of the context list:",
        placeholder = "e.g. Interferon response"
      ),
      shiny::fileInput(
        # new id each time the pop-up opens, so an earlier upload isn't reused
        paste0("context_file_", input$add_context),
        "Upload a gene list (.txt, .csv, .tsv):",
        accept = c(".txt", ".csv", ".tsv")
      ),
      shiny::helpText(
        "A file with one column, or genes separated by spaces, commas or",
        "semicolons, is added as one context under the name above.",
        "A file with several tab-, semicolon- or comma-separated columns is",
        "read as a table: the first row must hold the column headers, and",
        "each column is added as its own context named after its header",
        "(prefixed with the name above, if given). Gene IDs must match the",
        "selected gene identifier exactly, including upper/lower case."
      ),
      shiny::textAreaInput(
        "context_genes",
        "Or paste genes (one per line, or separated by spaces or commas):",
        value = "",
        placeholder = "GENE1, GENE2, GENE3\nor\nGENE1\nGENE2\nGENE3",
        rows = 5,
        width = "100%"
      ),
      footer = htmltools::tagList(
        shiny::modalButton("Cancel"),
        shiny::actionButton(
          "confirm_context",
          "Add",
          class = "btn-primary"
        )
      ),
      easyClose = TRUE
    ))
  })

  shiny::observeEvent(input$confirm_context, {
    context_file <- input[[paste0("context_file_", input$add_context)]]
    user_name <- trimws(input$context_name)
    context_name <- if (user_name == "") {
      paste("Context", length(rv$context_lists) + 1)
    } else {
      user_name
    }

    pasted <- .parse_gene_list(input$context_genes)
    file_contexts <- if (!is.null(context_file)) {
      .read_context_file(context_file$datapath)
    } else {
      list()
    }

    new_contexts <- list()
    if (is.null(names(file_contexts))) {
      # single gene list: combine with the pasted genes
      genes <- unique(c(pasted, unlist(file_contexts)))
      if (length(genes) > 0) {
        new_contexts[[context_name]] <- genes
      }
    } else {
      # table: one context per column
      if (user_name != "") {
        names(file_contexts) <- paste0(user_name, ": ", names(file_contexts))
      }
      new_contexts <- file_contexts
      if (length(pasted) > 0) {
        new_contexts[[context_name]] <- pasted
      }
    }

    if (length(new_contexts) == 0) {
      shiny::showNotification(
        "No genes found. Upload a file or paste gene names.",
        type = "error"
      )
      return()
    }

    context_lists <- rv$context_lists
    context_lists[names(new_contexts)] <- new_contexts
    rv$context_lists <- context_lists
    # set directly, the select observer won't fire if the name is unchanged
    rv$current_context <- new_contexts[[1]]

    shiny::updateSelectInput(
      session,
      "context_item",
      choices = names(rv$context_lists),
      selected = names(new_contexts)[1]
    )
    shiny::removeModal()
    shiny::showNotification(
      paste0(
        "Added ",
        paste0(
          "'",
          names(new_contexts),
          "' (",
          lengths(new_contexts),
          " genes)",
          collapse = ", "
        )
      ),
      type = "message"
    )
  })

  shiny::observeEvent(input$context_item, {
    shiny::req(input$context_item %in% names(rv$context_lists))
    rv$current_context <- rv$context_lists[[input$context_item]]
  })
  # nocov end
}
