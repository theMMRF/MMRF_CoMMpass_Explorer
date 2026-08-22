# Server function ------
shinyServer(function(input, output, session) {
  picker_state <- reactiveValues()  # per-input cache: labels and selection
  shinyjs::useShinyjs()

  .inputs_identical <- function(old, new) {
    if (is.null(old)) return(FALSE)
    identical(old$labels, new$labels) && identical(old$selected, new$selected)
  }

  .filter_survival_for_cohort <- function(data, cohort_id, read_input = function(id) input[[id]]) {
    if (!isTRUE(read_input(paste0("enable_survival_filter_", cohort_id)))) return(data)

    threshold_type <- read_input(paste0("surv_threshold_type_", cohort_id))
    surv_var <- read_input(paste0("surv_variable_", cohort_id))
    require_event <- isTRUE(read_input(paste0("require_surv_event_", cohort_id)))
    event_var <- survival_event_column(surv_var)

    if (threshold_type == "percentile") {
      return(filter_by_survival(
        data,
        surv_var = surv_var,
        threshold_type = "percentile",
        min_percentile = read_input(paste0("surv_threshold_min_percentile_", cohort_id)),
        max_percentile = read_input(paste0("surv_threshold_max_percentile_", cohort_id)),
        require_event = require_event,
        event_var = event_var
      ))
    }

    if (threshold_type == "value") {
      return(filter_by_survival(
        data,
        surv_var = surv_var,
        threshold_type = "value",
        min_value = read_input(paste0("surv_threshold_min_value_", cohort_id)),
        max_value = read_input(paste0("surv_threshold_max_value_", cohort_id)),
        require_event = require_event,
        event_var = event_var
      ))
    }

    data
  }

  .cohort_name <- function(cohort_id) {
    default <- if (identical(cohort_id, "cohort1")) "Cohort 1" else "Cohort 2"
    value <- input[[paste0("cohort_name_", cohort_id)]]
    if (is.null(value) || !nzchar(trimws(value))) default else trimws(value)
  }

  .cohort_labels <- function() {
    if (exists("cohort_metadata", inherits = TRUE) && !is.null(cohort_metadata$labels)) {
      return(cohort_metadata$labels)
    }
    labels <- c(Cohort1 = .cohort_name("cohort1"), Cohort2 = .cohort_name("cohort2"))
    if (identical(labels[["Cohort1"]], labels[["Cohort2"]])) {
      labels[["Cohort1"]] <- paste0(labels[["Cohort1"]], " (1)")
      labels[["Cohort2"]] <- paste0(labels[["Cohort2"]], " (2)")
    }
    labels
  }

  .with_display_cohort <- function(data) {
    labels <- .cohort_labels()
    data$cohort <- factor(labels[as.character(data$cohort)], levels = unname(labels))
    data
  }

  .observe_percentile_bounds <- function(input_id) {
    observeEvent(input[[input_id]], {
      value <- suppressWarnings(as.numeric(input[[input_id]]))
      if (!is.finite(value)) return()
      bounded <- max(0, min(100, value))
      if (!identical(value, bounded)) {
        updateNumericInput(session, input_id, value = bounded)
      }
    }, ignoreNULL = TRUE)
  }

  lapply(
    c(
      "surv_threshold_min_percentile_cohort1",
      "surv_threshold_max_percentile_cohort1",
      "surv_threshold_min_percentile_cohort2",
      "surv_threshold_max_percentile_cohort2",
      "gene_expr_percentile_min_cohort1",
      "gene_expr_percentile_max_cohort1",
      "gene_expr_percentile_min_cohort2",
      "gene_expr_percentile_max_cohort2"
    ),
    .observe_percentile_bounds
  )

  .count_card <- function(num_total, num_cohort1, num_cohort2, extra = NULL) {
    labels <- .cohort_labels()
    extra_html <- if (is.null(extra)) "" else paste0("<br>", extra)
    HTML(sprintf(
      "<div style='padding:10px; border: 1px solid #ccc; border-radius: 5px; background-color:#f5f5f5;'>
     <strong>Total Samples:</strong> %d<br>
     <span style='color: #E41A1C;'>%s:</span> %d (%.2f%%)<br>
     <span style='color: #4DBBD5;'>%s:</span> %d (%.2f%%)%s
   </div>",
      num_total,
      htmltools::htmlEscape(labels[["Cohort1"]]), num_cohort1, ifelse(num_total == 0, 0, num_cohort1 / num_total * 100),
      htmltools::htmlEscape(labels[["Cohort2"]]), num_cohort2, ifelse(num_total == 0, 0, num_cohort2 / num_total * 100),
      extra_html
    ))
  }

  .fmt_values <- function(values) {
    values <- values[!is.na(values) & nzchar(as.character(values))]
    paste(as.character(values), collapse = ", ")
  }

  .cohort_filter_labels <- c(
    sex = "Sex", race = "Race", genetic_ancestry = "Genetic ancestry",
    stage = "ISS stage", risk = "IMWG risk", cyto_risk = "Cytogenetic high risk",
    cgs_risk = "CGS risk", rna_subtype = "RNA subtype", cna_subtype = "CNA subtype",
    triplet = "Triplet firstline", asct = "ASCT firstline", regimen = "Regimen firstline",
    q21_gain = "1q21 gain", q21_amp = "1q21 amplification", del13q14 = "13q14 deletion",
    del13q34 = "13q34 deletion", del17p13 = "17p13 deletion", diploidy = "Hyperdiploidy",
    chromothripsis = "Chromothripsis", t11_14 = "t(11;14)", t4_14 = "t(4;14)",
    maf = "MAF/MAFB", apobec = "APOBEC", tp53 = "TP53 functional copies",
    tp53_ns = "TP53 non-synonymous mutation count"
  )

  .describe_cohort_selection <- function(cohort_id) {
    filters <- c(
      get_cohort_filters(input, cohort_id, "clinical"),
      get_cohort_filters(input, cohort_id, "molecular")
    )
    parts <- character()

    for (key in names(.cohort_filter_labels)) {
      values <- filters[[key]]
      if (!is.null(values) && length(values) > 0) {
        parts <- c(parts, paste0(.cohort_filter_labels[[key]], ": ", .fmt_values(values)))
      }
    }

    age <- filters$age
    if (!is.null(age) && length(age) == 2) {
      min_age <- min(clinical_data$Age, na.rm = TRUE)
      max_age <- max(clinical_data$Age, na.rm = TRUE)
      if (!identical(as.numeric(age), c(min_age, max_age))) {
        parts <- c(parts, sprintf("Age: %s-%s", age[1], age[2]))
      }
    }

    mut_rows <- if (identical(cohort_id, "cohort1")) mut_row_counter$cohort1 else mut_row_counter$cohort2
    if (mut_rows > 0) {
      mut_parts <- character()
      for (i in seq_len(mut_rows)) {
        selection <- .mutation_rule_selection(input, cohort_id, i)
        state <- input[[paste0("state_mut_", i, "_", cohort_id)]]
        logic <- input[[paste0("logic_mut_", i, "_", cohort_id)]]
        if (!is.null(selection) && nzchar(selection) && !is.null(state) && nzchar(state)) {
          connector <- if (!is.null(logic) && !identical(logic, "END")) paste0(" ", logic) else ""
          mut_parts <- c(mut_parts, paste0(
            .mutation_selector_display(selection), " ", state, connector
          ))
        }
      }
      if (length(mut_parts)) parts <- c(parts, paste0("Mutation rules: ", paste(mut_parts, collapse = "; ")))
    }

    gene <- input[[paste0("gene_expr_search_", cohort_id)]]
    if (!is.null(gene) && nzchar(gene)) {
      expr_type <- input[[paste0("expr_threshold_type_", cohort_id)]]
      if (identical(expr_type, "percentile")) {
        parts <- c(parts, sprintf(
          "Gene expression: %s percentile %s-%s",
          gene,
          input[[paste0("gene_expr_percentile_min_", cohort_id)]],
          input[[paste0("gene_expr_percentile_max_", cohort_id)]]
        ))
      } else {
        parts <- c(parts, sprintf(
          "Gene expression: %s value %s-%s",
          gene,
          input[[paste0("gene_expr_min_", cohort_id)]],
          input[[paste0("gene_expr_max_", cohort_id)]]
        ))
      }
    }

    if (isTRUE(input[[paste0("enable_survival_filter_", cohort_id)]])) {
      surv_var <- input[[paste0("surv_variable_", cohort_id)]]
      stype <- input[[paste0("surv_threshold_type_", cohort_id)]]
      event_text <- if (isTRUE(input[[paste0("require_surv_event_", cohort_id)]])) ", event = 1" else ""
      if (identical(stype, "percentile")) {
        parts <- c(parts, sprintf(
          "Survival: %s percentile %s-%s%s",
          surv_var,
          input[[paste0("surv_threshold_min_percentile_", cohort_id)]],
          input[[paste0("surv_threshold_max_percentile_", cohort_id)]],
          event_text
        ))
      } else {
        parts <- c(parts, sprintf(
          "Survival: %s days %s-%s%s",
          surv_var,
          input[[paste0("surv_threshold_min_value_", cohort_id)]],
          input[[paste0("surv_threshold_max_value_", cohort_id)]],
          event_text
        ))
      }
    }

    loaded_n <- if (identical(cohort_id, "cohort1")) length(user_cohorts$c1_tsb) else length(user_cohorts$c2_tsb)
    if (loaded_n > 0) {
      parts <- c(parts, sprintf("Uploaded cohort constraint: %d matched IDs", loaded_n))
    }

    if (!length(parts)) "All eligible samples; no cohort-specific filters were applied." else paste(parts, collapse = "; ")
  }

  tour_stage <- reactiveVal("idle")

  observeEvent(input$start_tour, {
    # -------- Phase 1: Custom Cohorts side --------
    steps1 <- data.frame(
      element = c(
        "body",
        "#how_it_works",
        "#cohort1_upload",
        "#cohort2_upload",
        "#filters_sidebar",
        "#apply_filters_wrap",
        'a[data-value="Overall Summary"]'  # <- ask user to click
      ),
      intro = c(
        "Welcome to CoMMpass Explorer! This tour highlights the main parts of the app.",
        "Start here. This panel explains how to load public_ids as CSV/TXT or paste them directly.",
        "Load or paste Cohort 1 patient IDs here. Click <b>Load Cohort 1 IDs</b>.",
        "Do the same for Cohort 2. You can also download any unmatched IDs.",
        "Use these filters to define cohorts (demographics, molecular features, etc.). Counts update live as you filter.",
        "<b>Apply Filters</b> locks in your current settings across tabs.",
        "Now click <b>Overall Summary</b> to continue the tour."
      ),
      position = c("floating","right","right","left","right","right","bottom"),
      stringsAsFactors = FALSE
    )

    rintrojs::introjs(
      session,
      options = list(
        steps = steps1,
        nextLabel = "Next", prevLabel = "Back", skipLabel = "Skip", doneLabel = "Finish",
        showProgress = TRUE, showBullets = FALSE, scrollToElement = TRUE,
        disableInteraction = FALSE, overlayOpacity = 0.35
      )
    )

    # Arm phase 2
    tour_stage("await_summary")
  })

  # -------- Phase 2: start when user actually switches to Overall Summary --------
  observeEvent(input$main_tabs, {
    if (identical(tour_stage(), "await_summary") && identical(input$main_tabs, "Overall Summary")) {
      # small delay lets the tab content render before anchoring steps
      later::later(function() {
        steps2 <- data.frame(
          element = c("#counts_card", "#summary_g1", "#summary_g2"),
          intro = c(
            "This card shows how many samples are in each cohort after filtering.",
            "A quick summary of Cohort 1's clinical distributions.",
            "And the same for Cohort 2."
          ),
          position = c("left","top","top"),
          stringsAsFactors = FALSE
        )

        rintrojs::introjs(
          session,
          options = list(
            steps = steps2,
            nextLabel = "Next", prevLabel = "Back", skipLabel = "Skip", doneLabel = "Finish",
            showProgress = TRUE, showBullets = FALSE, scrollToElement = TRUE,
            disableInteraction = FALSE, overlayOpacity = 0.35
          )
        )
        tour_stage("done")
      }, delay = 0.15)
    }
  })



  # --- user-defined cohorts  --------------------------------------------------
  user_cohorts <- reactiveValues(
    c1_public = character(), c2_public = character(),
    c1_tsb = character(), c2_tsb = character(),
    c1_unmatched = character(), c2_unmatched = character()
  )

  .has_loaded_cohort_ids <- function(cohort_num) {
    if (identical(cohort_num, 1L)) {
      return(length(user_cohorts$c1_public) > 0 ||
               length(user_cohorts$c1_tsb) > 0 ||
               length(user_cohorts$c1_unmatched) > 0)
    }
    length(user_cohorts$c2_public) > 0 ||
      length(user_cohorts$c2_tsb) > 0 ||
      length(user_cohorts$c2_unmatched) > 0
  }

  observe({
    toggleState("clear_loaded_cohort1", condition = .has_loaded_cohort_ids(1L))
    toggleState("clear_loaded_cohort2", condition = .has_loaded_cohort_ids(2L))
  })

  cohort_metadata <- reactiveValues(
    labels = c(Cohort1 = "Cohort 1", Cohort2 = "Cohort 2"),
    desc1 = "All eligible samples; no cohort-specific filters were applied.",
    desc2 = "All eligible samples; no cohort-specific filters were applied."
  )

  # Load Cohort 1
  observeEvent(input$load_cohort1, {
    ids <- unique(c(.parse_public_ids_from_file(input$upload_cohort1),
                    .parse_public_ids_from_text(input$paste_cohort1)))
    map <- .map_public_ids_to_tsb(ids, clinical_data)
    user_cohorts$c1_public <- map$matched_public
    user_cohorts$c1_tsb <- map$tsb
    user_cohorts$c1_unmatched <- map$unmatched

    if (length(map$unmatched)) showNotification(sprintf("%d Cohort 1 IDs were unmatched", length(map$unmatched)), type = "warning")
  })

  observeEvent(input$clear_loaded_cohort1, {
    user_cohorts$c1_public <- character()
    user_cohorts$c1_tsb <- character()
    user_cohorts$c1_unmatched <- character()
    showNotification("Cleared loaded Cohort 1 IDs. Click Apply Filters to refresh downstream panels.", type = "message")
  })

  # Load Cohort 2
  observeEvent(input$load_cohort2, {
    ids <- unique(c(.parse_public_ids_from_file(input$upload_cohort2),
                    .parse_public_ids_from_text(input$paste_cohort2)))
    map <- .map_public_ids_to_tsb(ids, clinical_data)
    user_cohorts$c2_public <- map$matched_public
    user_cohorts$c2_tsb <- map$tsb
    user_cohorts$c2_unmatched <- map$unmatched

    if (length(map$unmatched)) showNotification(sprintf("%d Cohort 2 IDs were unmatched", length(map$unmatched)), type = "warning")
  })

  observeEvent(input$clear_loaded_cohort2, {
    user_cohorts$c2_public <- character()
    user_cohorts$c2_tsb <- character()
    user_cohorts$c2_unmatched <- character()
    showNotification("Cleared loaded Cohort 2 IDs. Click Apply Filters to refresh downstream panels.", type = "message")
  })

  # Status + previews
  output$cohort1_status <- renderText({
    if (!.has_loaded_cohort_ids(1L)) {
      return("No uploaded ID constraint active.")
    }
    sprintf("Matched: %d  |  Unmatched: %d", length(user_cohorts$c1_tsb), length(user_cohorts$c1_unmatched))
  })
  output$cohort2_status <- renderText({
    if (!.has_loaded_cohort_ids(2L)) {
      return("No uploaded ID constraint active.")
    }
    sprintf("Matched: %d  |  Unmatched: %d", length(user_cohorts$c2_tsb), length(user_cohorts$c2_unmatched))
  })

  output$cohort1_preview <- renderDT({
    if (!length(user_cohorts$c1_tsb)) return(datatable(data.frame(Message = "No Cohort 1 IDs loaded")))
    dat <- clinical_data[clinical_data$Tumor_Sample_Barcode %in% user_cohorts$c1_tsb, ]
    datatable(dat[, c("public_id","Tumor_Sample_Barcode","Age","Sex","Race")], options = list(pageLength = 5, scrollX = TRUE))
  })
  output$cohort2_preview <- renderDT({
    if (!length(user_cohorts$c2_tsb)) return(datatable(data.frame(Message = "No Cohort 2 IDs loaded")))
    dat <- clinical_data[clinical_data$Tumor_Sample_Barcode %in% user_cohorts$c2_tsb, ]
    datatable(dat[, c("public_id","Tumor_Sample_Barcode","Age","Sex","Race")], options = list(pageLength = 5, scrollX = TRUE))
  })

  output$download_unmatched_c1 <- downloadHandler(
    filename = function() paste0("unmatched_cohort1_", Sys.Date(), ".txt"),
    content = function(file) writeLines(user_cohorts$c1_unmatched, file)
  )
  output$download_unmatched_c2 <- downloadHandler(
    filename = function() paste0("unmatched_cohort2_", Sys.Date(), ".txt"),
    content = function(file) writeLines(user_cohorts$c2_unmatched, file)
  )


  # Mutational Profile. Store rule row counters
  mut_row_counter <- reactiveValues(cohort1 = 0, cohort2 = 0)
  mut_rule_cache <- reactiveValues(cohort1 = list(), cohort2 = list())
  mut_rule_restore <- reactiveValues(cohort1 = FALSE, cohort2 = FALSE)

  # Disable remove when no rules exist
  observe({
    toggleState("remove_mut_row_cohort1", condition = mut_row_counter$cohort1 > 0)
    toggleState("remove_mut_row_cohort2", condition = mut_row_counter$cohort2 > 0)
  })

  # Add new rows on click
  observeEvent(input$add_mut_row_cohort1, {
    isolate({
      mut_rule_cache$cohort1 <- lapply(seq_len(mut_row_counter$cohort1), function(i) {
        list(
          gene = input[[paste0("gene_mut_", i, "_cohort1")]],
          codon = input[[paste0("codon_mut_", i, "_cohort1")]],
          variant = input[[paste0("variant_mut_", i, "_cohort1")]],
          state = input[[paste0("state_mut_", i, "_cohort1")]],
          logic = input[[paste0("logic_mut_", i, "_cohort1")]]
        )
      })
      mut_rule_restore$cohort1 <- TRUE
    })
    mut_row_counter$cohort1 <- mut_row_counter$cohort1 + 1
  })


  # Remove new rows
  observeEvent(input$remove_mut_row_cohort1, {
    if (mut_row_counter$cohort1 > 0) {
      # Save inputs from remaining rows before decrement
      isolate({
        mut_rule_cache$cohort1 <- lapply(seq_len(mut_row_counter$cohort1 - 1), function(i) {
          list(
            gene = input[[paste0("gene_mut_", i, "_cohort1")]],
            codon = input[[paste0("codon_mut_", i, "_cohort1")]],
            variant = input[[paste0("variant_mut_", i, "_cohort1")]],
            state = input[[paste0("state_mut_", i, "_cohort1")]],
            logic = input[[paste0("logic_mut_", i, "_cohort1")]]
          )
        })
        mut_rule_restore$cohort1 <- TRUE
      })
      mut_row_counter$cohort1 <- mut_row_counter$cohort1 - 1
    }
  })


  observeEvent(input$add_mut_row_cohort2, {
    isolate({
      mut_rule_cache$cohort2 <- lapply(seq_len(mut_row_counter$cohort2), function(i) {
        list(
          gene = input[[paste0("gene_mut_", i, "_cohort2")]],
          codon = input[[paste0("codon_mut_", i, "_cohort2")]],
          variant = input[[paste0("variant_mut_", i, "_cohort2")]],
          state = input[[paste0("state_mut_", i, "_cohort2")]],
          logic = input[[paste0("logic_mut_", i, "_cohort2")]]
        )
      })
      mut_rule_restore$cohort2 <- TRUE
    })
    mut_row_counter$cohort2 <- mut_row_counter$cohort2 + 1
  })


  observeEvent(input$remove_mut_row_cohort2, {
    if (mut_row_counter$cohort2 > 0) {
      isolate({
        mut_rule_cache$cohort2 <- lapply(seq_len(mut_row_counter$cohort2 - 1), function(i) {
          list(
            gene = input[[paste0("gene_mut_", i, "_cohort2")]],
            codon = input[[paste0("codon_mut_", i, "_cohort2")]],
            variant = input[[paste0("variant_mut_", i, "_cohort2")]],
            state = input[[paste0("state_mut_", i, "_cohort2")]],
            logic = input[[paste0("logic_mut_", i, "_cohort2")]]
          )
        })
        mut_rule_restore$cohort2 <- TRUE
      })

      mut_row_counter$cohort2 <- mut_row_counter$cohort2 - 1
    }
  })

  mutation_selectize_options <- function(placeholder) list(
    maxOptions = 100,
    allowEmptyOption = TRUE,
    placeholder = placeholder,
    score = I(paste0(
      "function(search) {",
      "  var terms = search.toLowerCase().trim().split(/\\s+/).filter(Boolean);",
      "  return function(item) {",
      "    var text = ((item.label || '') + ' ' + (item.value || '')).toLowerCase();",
      "    for (var i = 0; i < terms.length; i++) {",
      "      if (text.indexOf(terms[i]) === -1) return 0;",
      "    }",
      "    var match = (item.label || '').match(/\\((\\d+) patients?\\)$/);",
      "    var prevalence = match ? parseInt(match[1], 10) : 0;",
      "    return 1 + Math.min(prevalence, 999999) / 1000000;",
      "  };",
      "}"
    ))
  )

  # Render dynamic rule UI
  output$mutation_rules_cohort1 <- renderUI({
    if (mut_row_counter$cohort1 == 0) return(NULL)

    lapply(seq_len(mut_row_counter$cohort1), function(i) {
      cached <- if (length(mut_rule_cache$cohort1) >= i) mut_rule_cache$cohort1[[i]] else NULL

      gene_val <- if (!is.null(cached) && !is.null(cached$gene)) cached$gene else ""
      state_val <- if (!is.null(cached)) cached$state else "Mutated"
      logic_val <- if (!is.null(cached)) cached$logic else "END"

      div(
        selectInput(
          paste0("logic_mut_", i, "_cohort1"),
          label = tags$div(style = "font-size: 12px;", "Logic"),
          choices = c("AND", "OR", "END"),
          selected = logic_val,
          width = "100%"
        ),
        selectInput(
          paste0("state_mut_", i, "_cohort1"),
          label = tags$div(style = "font-size: 12px;", "State"),
          choices = c("Mutated", "Not Mutated"),
          selected = state_val,
          width = "100%"
        ),
        selectizeInput(
          paste0("gene_mut_", i, "_cohort1"),
          label = tags$div(style = "font-size: 12px;", "Gene"),
          choices = NULL,
          selected = gene_val,
          options = mutation_selectize_options("Search genes..."),
          width = "100%"
        ),
        selectizeInput(
          paste0("codon_mut_", i, "_cohort1"),
          label = tags$div(style = "font-size: 12px;", "Codon (optional)"),
          choices = NULL,
          selected = if (!is.null(cached$codon)) cached$codon else "",
          options = mutation_selectize_options("Any codon"),
          width = "100%"
        ),
        selectizeInput(
          paste0("variant_mut_", i, "_cohort1"),
          label = tags$div(style = "font-size: 12px;", "Variant (optional)"),
          choices = NULL,
          selected = if (!is.null(cached$variant)) cached$variant else "",
          options = mutation_selectize_options("Any variant"),
          width = "100%"
        ),
        br()
      )
    })
  })

  output$mutation_rules_cohort2 <- renderUI({
    if (mut_row_counter$cohort2 == 0) return(NULL)

    lapply(seq_len(mut_row_counter$cohort2), function(i) {
      cached <- if (length(mut_rule_cache$cohort2) >= i) mut_rule_cache$cohort2[[i]] else NULL

      gene_val <- if (!is.null(cached) && !is.null(cached$gene)) cached$gene else ""
      state_val <- if (!is.null(cached)) cached$state else "Mutated"
      logic_val <- if (!is.null(cached)) cached$logic else "END"

      div(
        selectInput(
          paste0("logic_mut_", i, "_cohort2"),
          label = tags$div(style = "font-size: 12px;", "Logic"),
          choices = c("AND", "OR", "END"),
          selected = logic_val,
          width = "100%"
        ),
        selectInput(
          paste0("state_mut_", i, "_cohort2"),
          label = tags$div(style = "font-size: 12px;", "State"),
          choices = c("Mutated", "Not Mutated"),
          selected = state_val,
          width = "100%"
        ),
        selectizeInput(
          paste0("gene_mut_", i, "_cohort2"),
          label = tags$div(style = "font-size: 12px;", "Gene"),
          choices = NULL,
          selected = gene_val,
          options = mutation_selectize_options("Search genes..."),
          width = "100%"
        ),
        selectizeInput(
          paste0("codon_mut_", i, "_cohort2"),
          label = tags$div(style = "font-size: 12px;", "Codon (optional)"),
          choices = NULL,
          selected = if (!is.null(cached$codon)) cached$codon else "",
          options = mutation_selectize_options("Any codon"),
          width = "100%"
        ),
        selectizeInput(
          paste0("variant_mut_", i, "_cohort2"),
          label = tags$div(style = "font-size: 12px;", "Variant (optional)"),
          choices = NULL,
          selected = if (!is.null(cached$variant)) cached$variant else "",
          options = mutation_selectize_options("Any variant"),
          width = "100%"
        ),
        br()
      )
    })
  })

  .update_mutation_rule_choices <- function(cohort_id) {
    row_count <- if (identical(cohort_id, "cohort1")) {
      mut_row_counter$cohort1
    } else {
      mut_row_counter$cohort2
    }
    if (row_count < 1) return()

    cache <- if (identical(cohort_id, "cohort1")) {
      mut_rule_cache$cohort1
    } else {
      mut_rule_cache$cohort2
    }

    session$onFlushed(function() {
      for (i in seq_len(row_count)) {
        gene_id <- paste0("gene_mut_", i, "_", cohort_id)
        codon_id <- paste0("codon_mut_", i, "_", cohort_id)
        variant_id <- paste0("variant_mut_", i, "_", cohort_id)
        cached <- if (length(cache) >= i) cache[[i]] else NULL

        gene <- isolate(input[[gene_id]])
        if ((is.null(gene) || !length(gene) || !nzchar(gene[1])) && !is.null(cached$gene)) {
          gene <- cached$gene
        }
        gene_choices <- mutation_gene_choices()
        if (is.null(gene) || !length(gene) || !gene[1] %in% unname(gene_choices)) gene <- ""

        updateSelectizeInput(
          session,
          gene_id,
          choices = gene_choices,
          selected = gene,
          server = TRUE
        )

        codon <- isolate(input[[codon_id]])
        if ((is.null(codon) || !length(codon) || !nzchar(codon[1])) && !is.null(cached$codon)) {
          codon <- cached$codon
        }
        codon_choices <- mutation_codon_choices(gene)
        if (is.null(codon) || !length(codon) || !codon[1] %in% unname(codon_choices)) codon <- ""
        updateSelectizeInput(
          session, codon_id, choices = codon_choices, selected = codon, server = TRUE
        )

        variant <- isolate(input[[variant_id]])
        if ((is.null(variant) || !length(variant) || !nzchar(variant[1])) && !is.null(cached$variant)) {
          variant <- cached$variant
        }
        variant_choices <- mutation_variant_choices(gene, codon)
        if (is.null(variant) || !length(variant) || !variant[1] %in% unname(variant_choices)) variant <- ""
        updateSelectizeInput(
          session, variant_id, choices = variant_choices, selected = variant, server = TRUE
        )
      }
    }, once = TRUE)
  }

  observeEvent(mut_row_counter$cohort1, {
    .update_mutation_rule_choices("cohort1")
  }, ignoreInit = FALSE)

  observeEvent(mut_row_counter$cohort2, {
    .update_mutation_rule_choices("cohort2")
  }, ignoreInit = FALSE)

  observe({
    for (cohort_id in c("cohort1", "cohort2")) {
      row_count <- if (identical(cohort_id, "cohort1")) {
        mut_row_counter$cohort1
      } else {
        mut_row_counter$cohort2
      }
      if (row_count < 1) next

      for (i in seq_len(row_count)) {
        gene_id <- paste0("gene_mut_", i, "_", cohort_id)
        codon_id <- paste0("codon_mut_", i, "_", cohort_id)
        gene <- input[[gene_id]]
        codon <- isolate(input[[codon_id]])
        restoring <- isolate(mut_rule_restore[[cohort_id]])
        cache <- isolate(mut_rule_cache[[cohort_id]])
        cached <- if (restoring && length(cache) >= i) cache[[i]] else NULL

        codon_choices <- mutation_codon_choices(gene)
        if ((is.null(codon) || !length(codon) || !nzchar(codon[1])) &&
            !is.null(cached$gene) && identical(cached$gene, gene) &&
            !is.null(cached$codon)) {
          codon <- cached$codon
        }
        valid_codon <- if (!is.null(codon) && length(codon) && codon[1] %in% unname(codon_choices)) codon[1] else ""
        updateSelectizeInput(
          session, codon_id, choices = codon_choices, selected = valid_codon, server = TRUE
        )
      }
    }
  })

  observe({
    for (cohort_id in c("cohort1", "cohort2")) {
      row_count <- if (identical(cohort_id, "cohort1")) {
        mut_row_counter$cohort1
      } else {
        mut_row_counter$cohort2
      }
      if (row_count < 1) next

      for (i in seq_len(row_count)) {
        gene_id <- paste0("gene_mut_", i, "_", cohort_id)
        codon_id <- paste0("codon_mut_", i, "_", cohort_id)
        variant_id <- paste0("variant_mut_", i, "_", cohort_id)
        gene <- input[[gene_id]]
        codon <- input[[codon_id]]
        restoring <- isolate(mut_rule_restore[[cohort_id]])
        cache <- isolate(mut_rule_cache[[cohort_id]])
        cached <- if (restoring && length(cache) >= i) cache[[i]] else NULL

        variant_choices <- mutation_variant_choices(gene, codon)
        variant <- isolate(input[[variant_id]])
        if ((is.null(variant) || !length(variant) || !nzchar(variant[1])) &&
            !is.null(cached$gene) && identical(cached$gene, gene) &&
            !is.null(cached$codon) && identical(cached$codon, codon) &&
            !is.null(cached$variant)) {
          variant <- cached$variant
        }
        valid_variant <- if (!is.null(variant) && length(variant) && variant[1] %in% unname(variant_choices)) variant[1] else ""
        updateSelectizeInput(
          session, variant_id, choices = variant_choices, selected = valid_variant, server = TRUE
        )
      }
    }
  })

  observe({
    gene_choices <- c("None" = "", unique(rownames(bulkseq_tpm)))
    updateSelectizeInput(session, "gene_expr_search_cohort1",
                         choices   = gene_choices,
                         selected  = "",
                         server    = TRUE)

    updateSelectizeInput(session, "gene_expr_search_cohort2",
                         choices   = gene_choices,
                         selected  = "",
                         server    = TRUE)
  })

  observe({
    surv_var <- input$surv_variable_cohort1
    if (!is.null(surv_var) && surv_var %in% names(clinical_data)) {
      vals <- clinical_data[[surv_var]]
      min_val <- min(vals, na.rm = TRUE)
      max_val <- max(vals, na.rm = TRUE)

      updateNumericInput(session, "surv_threshold_min_value_cohort1", value = min_val)
      updateNumericInput(session, "surv_threshold_max_value_cohort1", value = max_val)
    }
  })
  observe({
    surv_var <- input$surv_variable_cohort2
    if (!is.null(surv_var) && surv_var %in% names(clinical_data)) {
      vals <- clinical_data[[surv_var]]
      min_val <- min(vals, na.rm = TRUE)
      max_val <- max(vals, na.rm = TRUE)

      updateNumericInput(session, "surv_threshold_min_value_cohort2", value = min_val)
      updateNumericInput(session, "surv_threshold_max_value_cohort2", value = max_val)
    }
  })

  # Filters interface
  # ---- Dynamic count labels inside pickers -------------------------------------

  # Map: inputId  -> (filter_key used by filter_cohort_data, clinical column name)
  .make_picker_map <- function(cohort_id) {
    # clinical
    clinical_map <- data.frame(
      inputId    = paste0(c("sex_filter_", "race_filter_", "genetic_ancestry_filter_", "stage_filter_", "risk_filter_", "cyto_risk_filter_", "cgs_risk_filter_",
                            "rna_subtype_filter_", "cna_subtype_filter_", "triplet_filter_", "asct_filter_", "regimen_filter_"),
                          cohort_id),
      filter_key = c("sex","race","genetic_ancestry","stage","risk","cyto_risk","cgs_risk",
                     "rna_subtype","cna_subtype","triplet","asct","regimen"),
      column     = c("Sex","Race","genetic_ancestry","ISS","IMWG_Risk_Class","Skerget_Cytogenetic_High_Risk","CGS_risk",
                     "Skerget_RNA_Subtype_Name","Skerget_CNA_Subtype_Name","Triplet_First","ASCT_First", "regimen"),
      stringsAsFactors = FALSE
    )

    # molecular
    molecular_map <- data.frame(
      inputId    = paste0(c("chr_1q21_gain_filter_","chr_1q21_amp_filter_","chr_13q14_del_filter_",
                            "chr_13q34_del_filter_","chr_17p13_del_filter_","diploidy_filter_",
                            "chromothripsis_filter_","t_11_14_filter_","t_4_14_filter_",
                            "maf_filter_","apobec_filter_","tp53_filter_","tp53_ns_filter_"),
                          cohort_id),
      filter_key = c("q21_gain","q21_amp","del13q14","del13q34","del17p13","diploidy",
                     "chromothripsis","t11_14","t4_14","maf","apobec","tp53","tp53_ns"),
      column     = c("chr_1q21_gain","chr_1q21_amp","chr_13q14_del","chr_13q34_del","chr_17p13_del",
                     "Hyperdiploidy","chromothripsis","t_11_14","t_4_14","MAF_MAFB","APOBEC",
                     "TP53_Funct_Copies","TP53_NS_Mut_Count"),
      stringsAsFactors = FALSE
    )

    rbind(clinical_map, molecular_map)
  }

  .update_picker_counts_for_cohort <- function(cohort_id) {
    # categorical filters
    filters <- c(
      get_cohort_filters(input, cohort_id, "clinical"),
      get_cohort_filters(input, cohort_id, "molecular")
    )

    mp <- .make_picker_map(cohort_id)

    for (i in seq_len(nrow(mp))) {
      input_id   <- mp$inputId[i]
      filter_key <- mp$filter_key[i]
      colname    <- mp$column[i]

      if (!(input_id %in% names(input))) next
      open_id <- isolate(input$picker_open_id)
      if (!is.null(open_id) && identical(open_id, input_id)) next

      filters_excl <- filters
      filters_excl[[filter_key]] <- NULL

      # Start from full clinical, apply remaining categorical filters
      preview <- filter_cohort_data(copy(clinical_data), filters_excl)

      # condition counts on non-categorical filters ----------
      # Mutation rules
      mut_rows <- if (cohort_id == "cohort1") mut_row_counter$cohort1 else mut_row_counter$cohort2
      ids_mut  <- get_mutation_filtered_ids(input, cohort_id, mut_rows)
      preview  <- preview[preview$Tumor_Sample_Barcode %in% ids_mut, ]

      # Gene expression
      gene <- input[[paste0("gene_expr_search_", cohort_id)]]
      if (!is.null(gene) && nzchar(gene)) {
        preview <- filter_by_gene_expression(
          clinical_data  = preview,
          gene           = gene,
          threshold_type = input[[paste0("expr_threshold_type_", cohort_id)]],
          min_value      = input[[paste0("gene_expr_min_", cohort_id)]],
          max_value      = input[[paste0("gene_expr_max_", cohort_id)]],
          min_percentile = input[[paste0("gene_expr_percentile_min_", cohort_id)]],
          max_percentile = input[[paste0("gene_expr_percentile_max_", cohort_id)]]
        )
      }

      # Survival
      preview <- .filter_survival_for_cohort(preview, cohort_id)

      # Respect loaded custom cohort IDs for live counts
      allowed <- if (cohort_id == "cohort1") user_cohorts$c1_tsb else user_cohorts$c2_tsb
      if (length(allowed)) {
        preview <- preview[preview$Tumor_Sample_Barcode %in% allowed, ]
      }
      # ------------------------------------------------------------------------

      # Use all possible values from the full dataset (so zero-count levels show)
      all_vals <- sort(unique(na.omit(as.character(clinical_data[[colname]]))))
      if (length(all_vals) == 0) next

      tbl <- table(as.character(preview[[colname]]), useNA = "no")
      labels <- vapply(all_vals, function(v) paste0(v, " (", .safe_count(tbl, v), ")"), FUN.VALUE = character(1))
      choices_named <- setNames(all_vals, labels)

      current_sel <- input[[input_id]]
      current_sel <- intersect(current_sel, all_vals)

      # change-detection cache to avoid dropdown-close jitter
      new_state <- list(labels = labels, selected = current_sel)
      old_state <- isolate(picker_state[[input_id]])
      if (!.inputs_identical(old_state, new_state)) {
        shinyWidgets::updatePickerInput(session, inputId = input_id, choices = choices_named, selected = current_sel)
        isolate(picker_state[[input_id]] <- new_state)
      }
    }
  }

  # Re-run the updater whenever relevant inputs change
  # (This listens broadly to that cohort's filter inputs.)
  observe({
    # Touch all cohort1 filter inputs to create reactivity
    dummy <- list(
      input$sex_filter_cohort1, input$race_filter_cohort1, input$genetic_ancestry_filter_cohort1, input$stage_filter_cohort1,
      input$risk_filter_cohort1, input$cyto_risk_filter_cohort1, input$cgs_risk_filter_cohort1,
      input$rna_subtype_filter_cohort1, input$cna_subtype_filter_cohort1,
      input$triplet_filter_cohort1, input$asct_filter_cohort1,
      input$regimen_filter_cohort1,
      input$chr_1q21_gain_filter_cohort1, input$chr_1q21_amp_filter_cohort1,
      input$chr_13q14_del_filter_cohort1, input$chr_13q34_del_filter_cohort1,
      input$chr_17p13_del_filter_cohort1, input$diploidy_filter_cohort1,
      input$chromothripsis_filter_cohort1, input$t_11_14_filter_cohort1,
      input$t_4_14_filter_cohort1, input$maf_filter_cohort1,
      input$apobec_filter_cohort1, input$tp53_filter_cohort1, input$tp53_ns_filter_cohort1,
      input$age_cohort1,  # include age so categorical counts reflect age subset, too
      input$gene_expr_search_cohort1, input$expr_threshold_type_cohort1,
      input$gene_expr_min_cohort1, input$gene_expr_max_cohort1,
      input$gene_expr_percentile_min_cohort1, input$gene_expr_percentile_max_cohort1,
      input$enable_survival_filter_cohort1, input$surv_threshold_type_cohort1,
      input$surv_variable_cohort1, input$require_surv_event_cohort1,
      input$surv_threshold_min_value_cohort1, input$surv_threshold_max_value_cohort1,
      input$surv_threshold_min_percentile_cohort1, input$surv_threshold_max_percentile_cohort1
    )
    .update_picker_counts_for_cohort("cohort1")
  })

  observe({
    # Same for cohort2
    dummy <- list(
      input$sex_filter_cohort2, input$race_filter_cohort2, input$genetic_ancestry_filter_cohort2, input$stage_filter_cohort2,
      input$risk_filter_cohort2, input$cyto_risk_filter_cohort2, input$cgs_risk_filter_cohort2,
      input$rna_subtype_filter_cohort2, input$cna_subtype_filter_cohort2,
      input$triplet_filter_cohort2, input$asct_filter_cohort2,
      input$regimen_filter_cohort2,
      input$chr_1q21_gain_filter_cohort2, input$chr_1q21_amp_filter_cohort2,
      input$chr_13q14_del_filter_cohort2, input$chr_13q34_del_filter_cohort2,
      input$chr_17p13_del_filter_cohort2, input$diploidy_filter_cohort2,
      input$chromothripsis_filter_cohort2, input$t_11_14_filter_cohort2,
      input$t_4_14_filter_cohort2, input$maf_filter_cohort2,
      input$apobec_filter_cohort2, input$tp53_filter_cohort2, input$tp53_ns_filter_cohort2,
      input$age_cohort2,
      input$gene_expr_search_cohort2, input$expr_threshold_type_cohort2,
      input$gene_expr_min_cohort2, input$gene_expr_max_cohort2,
      input$gene_expr_percentile_min_cohort2, input$gene_expr_percentile_max_cohort2,
      input$enable_survival_filter_cohort2, input$surv_threshold_type_cohort2,
      input$surv_variable_cohort2, input$require_surv_event_cohort2,
      input$surv_threshold_min_value_cohort2, input$surv_threshold_max_value_cohort2,
      input$surv_threshold_min_percentile_cohort2, input$surv_threshold_max_percentile_cohort2
    )
    .update_picker_counts_for_cohort("cohort2")
  })

  # Also refresh counts when users click "Clear All" for a cohort
  observeEvent(input$clear_cohort1, {
    shinyjs::reset("cohort1_filters")
    # Let the UI reset, then recompute counts
    shiny::invalidateLater(50, session)
    .update_picker_counts_for_cohort("cohort1")
  })

  observeEvent(input$clear_cohort2, {
    shinyjs::reset("cohort2_filters")
    shiny::invalidateLater(50, session)
    .update_picker_counts_for_cohort("cohort2")
  })

  # Cohort 1
  output$age_filter_cohort1 <- renderUI({
    min_age <- min(clinical_data$Age, na.rm = TRUE)
    max_age <- max(clinical_data$Age, na.rm = TRUE)
    sliderInput("age_cohort1", "Age", min = min_age, max = max_age, value = c(min_age, max_age))
  })

  # Clear Cohort 1 filters when the "Clear" button is clicked
  observeEvent(input$clear_cohort1, {
    shinyjs::reset("cohort1_filters") # Resets all inputs within the div for Cohort 1 filters
  })

  observeEvent(input$clear_cohort2, {
    shinyjs::reset("cohort2_filters")
  })

  .cohort_filters <- function(cohort_id) {
    c(
      get_cohort_filters(input, cohort_id, "clinical"),
      get_cohort_filters(input, cohort_id, "molecular"),
      get_cohort_filters(input, cohort_id, "gene")
    )
  }

  .has_active_mutation_filter <- function(cohort_id, row_count, read_input = function(id) input[[id]]) {
    if (row_count < 1) return(FALSE)

    for (i in seq_len(row_count)) {
      gene <- read_input(paste0("gene_mut_", i, "_", cohort_id))
      state <- read_input(paste0("state_mut_", i, "_", cohort_id))
      if (!is.null(gene) && nzchar(gene) &&
          !is.null(state) && state %in% c("Mutated", "Not Mutated")) {
        return(TRUE)
      }
    }
    FALSE
  }

  .cohort_eligible_universe <- function(cohort_id) {
    eligible <- copy(clinical_data)
    filters <- .cohort_filters(cohort_id)
    picker_map <- .make_picker_map(cohort_id)

    for (i in seq_len(nrow(picker_map))) {
      values <- filters[[picker_map$filter_key[i]]]
      if (is.null(values) || !length(values)) next

      column <- picker_map$column[i]
      available <- !is.na(eligible[[column]]) & nzchar(trimws(as.character(eligible[[column]])))
      eligible <- eligible[available, ]
    }

    age <- filters$age
    if (!is.null(age) && length(age) == 2) {
      eligible <- eligible[is.finite(eligible$Age), ]
    }

    mut_rows <- if (identical(cohort_id, "cohort1")) mut_row_counter$cohort1 else mut_row_counter$cohort2
    if (.has_active_mutation_filter(
      cohort_id,
      mut_rows,
      read_input = function(id) isolate(input[[id]])
    )) {
      eligible <- eligible[
        eligible$Tumor_Sample_Barcode %in% .mutation_data_available_ids(),
      ]
    }

    gene <- isolate(input[[paste0("gene_expr_search_", cohort_id)]])
    if (!is.null(gene) && nzchar(gene) && gene %in% rownames(bulkseq_tpm)) {
      expression <- bulkseq_tpm[gene, ]
      expression_ids <- names(expression)[is.finite(expression)]
      eligible <- eligible[eligible$Tumor_Sample_Barcode %in% expression_ids, ]
    }

    survival_enabled <- isTRUE(isolate(input[[paste0("enable_survival_filter_", cohort_id)]]))
    if (survival_enabled) {
      threshold_type <- isolate(input[[paste0("surv_threshold_type_", cohort_id)]])
      threshold_is_set <- if (identical(threshold_type, "percentile")) {
        values <- c(
          isolate(input[[paste0("surv_threshold_min_percentile_", cohort_id)]]),
          isolate(input[[paste0("surv_threshold_max_percentile_", cohort_id)]])
        )
        length(values) == 2 && all(is.finite(as.numeric(values)))
      } else {
        values <- c(
          isolate(input[[paste0("surv_threshold_min_value_", cohort_id)]]),
          isolate(input[[paste0("surv_threshold_max_value_", cohort_id)]])
        )
        length(values) == 2 && all(is.finite(as.numeric(values)))
      }

      if (threshold_is_set) {
        surv_var <- isolate(input[[paste0("surv_variable_", cohort_id)]])
        if (!is.null(surv_var) && surv_var %in% names(eligible)) {
          eligible <- eligible[is.finite(eligible[[surv_var]]), ]
        }

        if (isTRUE(isolate(input[[paste0("require_surv_event_", cohort_id)]]))) {
          event_var <- survival_event_column(surv_var)
          if (!is.null(event_var) && event_var %in% names(eligible)) {
            eligible <- eligible[!is.na(eligible[[event_var]]), ]
          }
        }
      }
    }

    eligible
  }

  .resolve_base_cohort <- function(cohort_id) {
    cohort_label <- if (identical(cohort_id, "cohort1")) "Cohort 1" else "Cohort 2"
    clinical_filtered <- filter_cohort_data(copy(clinical_data), .cohort_filters(cohort_id))
    allowed <- if (identical(cohort_id, "cohort1")) {
      user_cohorts$c1_tsb
    } else {
      user_cohorts$c2_tsb
    }
    if (length(allowed)) {
      clinical_filtered <- clinical_filtered[clinical_filtered$Tumor_Sample_Barcode %in% allowed, ]
    }

    clinical_error <- copy(clinical_filtered)

    # Gene mutation filter
    mut_rows <- if (identical(cohort_id, "cohort1")) mut_row_counter$cohort1 else mut_row_counter$cohort2
    mutated_ids <- get_mutation_filtered_ids(isolate(input), cohort_id, mut_rows)
    clinical_filtered <- clinical_filtered[clinical_filtered$Tumor_Sample_Barcode %in% mutated_ids, ]

    # Gene expression filter
    clinical_filtered <- filter_by_gene_expression(
      clinical_data = clinical_filtered,
      gene = isolate(input[[paste0("gene_expr_search_", cohort_id)]]),
      threshold_type = isolate(input[[paste0("expr_threshold_type_", cohort_id)]]),
      min_value = isolate(input[[paste0("gene_expr_min_", cohort_id)]]),
      max_value = isolate(input[[paste0("gene_expr_max_", cohort_id)]]),
      min_percentile = isolate(input[[paste0("gene_expr_percentile_min_", cohort_id)]]),
      max_percentile = isolate(input[[paste0("gene_expr_percentile_max_", cohort_id)]])
    )

    # Survival filter
    clinical_filtered <- .filter_survival_for_cohort(
      clinical_filtered,
      cohort_id,
      read_input = function(id) isolate(input[[id]])
    )

    if (nrow(clinical_filtered) == 0) {
      showNotification(sprintf("No patients in %s match the filters.", cohort_label), type = "error")
      return(clinical_error)
    }

    return(clinical_filtered)
  }

  .complement_of <- function(other_cohort_data, reference_cohort_id) {
    other_ids <- unique(other_cohort_data$Tumor_Sample_Barcode)
    eligible <- .cohort_eligible_universe(reference_cohort_id)
    copy(eligible[!(eligible$Tumor_Sample_Barcode %in% other_ids), ])
  }

  .empty_maf_subset <- function(patient_ids) {
    empty_maf <- maf_data
    empty_maf@data <- data.table::copy(maf_data@data[0, ])
    empty_maf@variants.per.sample <- data.table::copy(maf_data@variants.per.sample[0, ])
    empty_maf@variant.type.summary <- data.table::copy(maf_data@variant.type.summary[0, ])
    empty_maf@variant.classification.summary <- data.table::copy(maf_data@variant.classification.summary[0, ])
    empty_maf@gene.summary <- data.table::copy(maf_data@gene.summary[0, ])
    empty_maf@summary <- data.table::copy(maf_data@summary[0, ])
    empty_maf@maf.silent <- data.table::copy(maf_data@maf.silent[0, ])
    empty_maf@clinical.data <- data.table::copy(
      maf_data@clinical.data[Tumor_Sample_Barcode %in% patient_ids, ]
    )
    empty_maf
  }

  .safe_subset_maf <- function(patient_ids, cohort_label) {
    tryCatch(
      subsetMaf(maf = maf_data, tsb = patient_ids),
      error = function(e) {
        showNotification(
          sprintf("%s has no MAF variants after filtering; mutational plots may be empty.", cohort_label),
          type = "warning"
        )
        .empty_maf_subset(patient_ids)
      }
    )
  }

  .resolve_applied_cohorts <- function() {
    complement_mode <- isolate(input$complement_mode)
    if (is.null(complement_mode) || !complement_mode %in% c("none", "cohort1", "cohort2")) {
      complement_mode <- "none"
    }

    if (identical(complement_mode, "cohort1")) {
      base_cohort2 <- .resolve_base_cohort("cohort2")
      base_cohort1 <- .complement_of(base_cohort2, "cohort2")
    } else if (identical(complement_mode, "cohort2")) {
      base_cohort1 <- .resolve_base_cohort("cohort1")
      base_cohort2 <- .complement_of(base_cohort1, "cohort1")
    } else {
      base_cohort1 <- .resolve_base_cohort("cohort1")
      base_cohort2 <- .resolve_base_cohort("cohort2")
    }

    list(
      cohort1 = base_cohort1,
      cohort2 = base_cohort2,
      complement_mode = complement_mode
    )
  }


  # Cohort 2
  output$age_filter_cohort2 <- renderUI({
    min_age <- min(clinical_data$Age, na.rm = TRUE)
    max_age <- max(clinical_data$Age, na.rm = TRUE)
    sliderInput("age_cohort2", "Age", min = min_age, max = max_age, value = c(min_age, max_age))
  })

  # filtered_data that stores combined clinical data and cohort1, cohort2 data
  filtered_data <- reactiveValues(
    cohort1 = {
      data_cohort1 <- clinical_data
      data_cohort1$cohort <- "Cohort1"
      data_cohort1$cohort_label <- "Cohort 1"
      data_cohort1
    },
    cohort2 = {
      data_cohort2 <- clinical_data
      data_cohort2$cohort <- "Cohort2"
      data_cohort2$cohort_label <- "Cohort 2"
      data_cohort2
    },
    cohort1_maf = maf_data,
    cohort2_maf = maf_data,
    combined = {
      data_cohort1 <- clinical_data
      data_cohort1$cohort <- "Cohort1"
      data_cohort1$cohort_label <- "Cohort 1"
      data_cohort2 <- clinical_data
      data_cohort2$cohort <- "Cohort2"
      data_cohort2$cohort_label <- "Cohort 2"
      rbind(data_cohort1, data_cohort2)
    }
  )

  observeEvent(input$apply_filters, {
    applied_cohorts <- .resolve_applied_cohorts()
    data_cohort1 <- applied_cohorts$cohort1
    data_cohort2 <- applied_cohorts$cohort2
    complement_mode <- applied_cohorts$complement_mode

    if (nrow(data_cohort1) == 0 || nrow(data_cohort2) == 0) {
      empty_names <- c(
        if (nrow(data_cohort1) == 0) "Cohort 1",
        if (nrow(data_cohort2) == 0) "Cohort 2"
      )
      showNotification(
        sprintf("Cannot apply filters because %s would be empty.", paste(empty_names, collapse = " and ")),
        type = "error"
      )
      return()
    }

    labels <- c(Cohort1 = .cohort_name("cohort1"), Cohort2 = .cohort_name("cohort2"))
    if (identical(labels[["Cohort1"]], labels[["Cohort2"]])) {
      labels[["Cohort1"]] <- paste0(labels[["Cohort1"]], " (1)")
      labels[["Cohort2"]] <- paste0(labels[["Cohort2"]], " (2)")
    }
    cohort_metadata$labels <- labels
    cohort_metadata$desc1 <- if (identical(complement_mode, "cohort1")) {
      sprintf(
        "Complement of %s: samples with the data required by %s's filters that were not selected after its filters/uploads (%d samples). Own filters/uploads are ignored in this mode.",
        labels[["Cohort2"]],
        labels[["Cohort2"]],
        nrow(data_cohort1)
      )
    } else {
      .describe_cohort_selection("cohort1")
    }
    cohort_metadata$desc2 <- if (identical(complement_mode, "cohort2")) {
      sprintf(
        "Complement of %s: samples with the data required by %s's filters that were not selected after its filters/uploads (%d samples). Own filters/uploads are ignored in this mode.",
        labels[["Cohort1"]],
        labels[["Cohort1"]],
        nrow(data_cohort2)
      )
    } else {
      .describe_cohort_selection("cohort2")
    }

    data_cohort1$cohort <- "Cohort1"
    data_cohort2$cohort <- "Cohort2"
    data_cohort1$cohort_label <- labels[["Cohort1"]]
    data_cohort2$cohort_label <- labels[["Cohort2"]]

    n1 <- names(data_cohort1)
    n2 <- names(data_cohort2)
    common <- intersect(n1, n2)

    if (length(common) == 0L) {
      showNotification("No common columns between cohorts after filtering.", type = "error")
      return()
    }

    # keep only common columns, same order
    data_cohort1 <- as.data.frame(data_cohort1[, common, drop = FALSE])
    data_cohort2 <- as.data.frame(data_cohort2[, common, drop = FALSE])

    combined_clinical <- rbind(data_cohort1, data_cohort2)

    # Update reactiveValues
    filtered_data$cohort1 <- data_cohort1
    filtered_data$cohort2 <- data_cohort2
    filtered_data$combined <- combined_clinical

    patient_ids_cohort1 <- data_cohort1$Tumor_Sample_Barcode
    patient_ids_cohort2 <- data_cohort2$Tumor_Sample_Barcode
    filtered_data$cohort1_maf <- .safe_subset_maf(patient_ids_cohort1, labels[["Cohort1"]])
    filtered_data$cohort2_maf <- .safe_subset_maf(patient_ids_cohort2, labels[["Cohort2"]])
  })

  # Reactive expression for preprocessed bulk RNA-seq data
  preprocessed_bulkseq_data <- reactive({
    req(input$gene_search_bulk_distr)

    gene_interested <- input$gene_search_bulk_distr
    clinical_combined <- filtered_data$combined

    inter_samples <- intersect(clinical_combined$Tumor_Sample_Barcode, colnames(bulkseq_tpm))
    num_sample <- intersect(clinical_data$Tumor_Sample_Barcode, colnames(bulkseq_tpm))
    combined_bulkseq_tpm <- bulkseq_tpm[, inter_samples]
    clinical_combined <- clinical_combined[clinical_combined$Tumor_Sample_Barcode %in% inter_samples, ]

    list(
      num_sample = num_sample,
      gene_interested = gene_interested,
      clinical_combined = clinical_combined,
      combined_bulkseq_tpm = combined_bulkseq_tpm
    )
  })

  preprocessed_sc_meta <- reactive({
    clinical_combined <- filtered_data$combined[, c("Tumor_Sample_Barcode", "cohort")]

    # Derive public_id (MMRF_XXXX) from Tumor_Sample_Barcode like before
    patients_in_clinical <- sapply(
      strsplit(as.character(clinical_combined$Tumor_Sample_Barcode), "_"),
      function(x) paste(x[1], x[2], sep = "_")
    )
    clinical_combined$public_id <- patients_in_clinical

    # Intersect by public_id
    inter_public <- intersect(clinical_combined$public_id, sc_meta$public_id)

    # How many clinical patients have scRNA-seq data
    num_sample <- length(unique(sc_meta$d_visit_specimen_id))

    clinical_combined <- clinical_combined[clinical_combined$public_id %in% inter_public, ]

    cohort_info <- clinical_combined %>%
      dplyr::distinct(public_id, cohort)

    list(
      num_sample    = num_sample,
      inter_samples = inter_public,
      clinical_combined = clinical_combined,
      cohort_info   = cohort_info
    )
  })

  # Number of samples in each dataset after subsetting -------------------------
  # Clinical
  output$clinicalNum <- renderUI({
    num_total <- nrow(clinical_data)
    num_cohort1 <- nrow(filtered_data$cohort1)
    num_cohort2 <- nrow(filtered_data$cohort2)
    .count_card(num_total, num_cohort1, num_cohort2)
  })

  output$cohort_descriptions <- renderUI({
    labels <- .cohort_labels()
    HTML(sprintf(
      "<div class='cohort-definition-card cohort-definition-card-c1'><strong>%s</strong><p>%s</p></div>
       <div class='cohort-definition-card cohort-definition-card-c2'><strong>%s</strong><p>%s</p></div>",
      htmltools::htmlEscape(labels[["Cohort1"]]),
      htmltools::htmlEscape(cohort_metadata$desc1),
      htmltools::htmlEscape(labels[["Cohort2"]]),
      htmltools::htmlEscape(cohort_metadata$desc2)
    ))
  })

  # MAF
  output$mafNum <- renderUI({
    num_total <- uniqueN(maf_data@data$Tumor_Sample_Barcode)

    num_cohort1 <- nrow(filtered_data$cohort1_maf@clinical.data)
    num_cohort2 <- nrow(filtered_data$cohort2_maf@clinical.data)

    .count_card(num_total, num_cohort1, num_cohort2)
  })

  # Transcriptomics
  output$bulkNum <- renderUI({
    bulk_clinical <- preprocessed_bulkseq_data()$clinical_combined
    num_sample <- preprocessed_bulkseq_data()$num_sample
    num_total <- length(num_sample)
    num_cohort1 <- nrow(bulk_clinical[bulk_clinical$cohort == "Cohort1", ])
    num_cohort2 <- nrow(bulk_clinical[bulk_clinical$cohort == "Cohort2", ])
    .count_card(num_total, num_cohort1, num_cohort2)
  })

  # scRNA-seq
  output$scNum <- renderUI({
    sc_clinical <- preprocessed_sc_meta()$cohort_info
    inter_samples <- preprocessed_sc_meta()$inter_samples
    num_total <- preprocessed_sc_meta()$num_sample
    num_cohort1 <- nrow(sc_clinical[sc_clinical$cohort == "Cohort1", ])
    num_cohort2 <- nrow(sc_clinical[sc_clinical$cohort == "Cohort2", ])
    .count_card(num_total, num_cohort1, num_cohort2)
  })

  # Update gene choices for selectizeInput
  observe({
    updateSelectizeInput(session, "gene_search_lollipop_g1", choices = unique(maf_data@gene.summary$Hugo_Symbol), selected = "KRAS", server = TRUE)
    updateSelectizeInput(session, "gene_search_lollipop_g2", choices = unique(maf_data@gene.summary$Hugo_Symbol), selected = "KRAS", server = TRUE)
    updateSelectizeInput(session, "gene_search_maf", choices = unique(maf_data@gene.summary$Hugo_Symbol), server = TRUE)
    updateSelectizeInput(session, "gene_search_inter_g1", choices = unique(filtered_data[["cohort1_maf"]]@gene.summary$Hugo_Symbol), server = TRUE)
    updateSelectizeInput(session, "gene_search_inter_g2", choices = unique(filtered_data[["cohort2_maf"]]@gene.summary$Hugo_Symbol), server = TRUE)
    # updateSelectizeInput(session, "gene_search_bulk_heat", choices = unique(rownames(bulkseq_tpm)), server = TRUE)
    updateSelectizeInput(session, "gene_search_bulk_distr", choices = unique(rownames(bulkseq_tpm)), selected = "KRAS", server = TRUE)
  })

  # Summary ----------------------
  draw_summaryPlot_g1 <- function() {
    generate_summary_plot("cohort1", filtered_data) +
      ggplot2::labs(title = .cohort_labels()[["Cohort1"]])
  }

  draw_summaryPlot_g2 <- function() {
    generate_summary_plot("cohort2", filtered_data) +
      ggplot2::labs(title = .cohort_labels()[["Cohort2"]])
  }

  draw_survCompPlot_pfs <- function() {
    req(filtered_data$combined)
    plot_cohort_survival(
      data     = filtered_data$combined,
      time_col = "PFS",
      event_col = "PFS_event",
      xlab      = "Progression-Free Survival (Days)",
      break_by  = 500,
      legend_labs = unname(.cohort_labels())
    )
  }

  draw_survCompPlot_os <- function() {
    req(filtered_data$combined)
    plot_cohort_survival(
      data      = filtered_data$combined,
      time_col  = "OS",
      event_col = "OS_event",
      xlab      = "Overall Survival (Days)",
      break_by  = 500,
      legend_labs = unname(.cohort_labels())
    )
  }

  draw_survCompPlot_tt2Line_censored <- function() {
    req(filtered_data$combined)
    plot_cohort_survival(
      data      = filtered_data$combined,
      time_col  = "ttct2line",
      event_col = "censt2line",
      xlab      = "Time to Second Line (Days)",
      break_by  = 500,
      legend_labs = unname(.cohort_labels())
    )
  }

  draw_clin_distribution_ggplot <- function() {
    combined_data <- .with_display_cohort(filtered_data$combined)
    interested_feature <- input$clin_feature
    continuous_features <- get_continuous_features()

    create_distribution_ggplot(
      data = combined_data,
      feature = interested_feature,
      feature_label = input$clin_feature,
      continuous_features = continuous_features
    )
  }

  # Summary of clinical data
  output$summaryPlot_g1 <- renderPlot({
    draw_summaryPlot_g1()
  })

  output$summaryPlot_g2 <- renderPlot({
    draw_summaryPlot_g2()
  })

  # Draw survival curve comparison plot (PFS)
  output$survCompPlot_pfs <- renderPlot({
    draw_survCompPlot_pfs()
  })

  # Draw survival curve comparison plot (OS)
  output$survCompPlot_os <- renderPlot({
    draw_survCompPlot_os()
  })

  # Time to second line survival curve
  output$survCompPlot_tt2Line_censored <- renderPlot({
    draw_survCompPlot_tt2Line_censored()
  })

  # Distribution Plot
  output$clin_distribution <- renderPlotly({
    ggplotly(draw_clin_distribution_ggplot(), tooltip = "text") %>%
      layout(hovermode = "x")
  })

  register_plot_pdf_download(output, input, "summaryPlot_g1", draw_summaryPlot_g1)
  register_plot_pdf_download(output, input, "summaryPlot_g2", draw_summaryPlot_g2)
  register_plot_pdf_download(output, input, "survCompPlot_pfs", draw_survCompPlot_pfs)
  register_plot_pdf_download(output, input, "survCompPlot_os", draw_survCompPlot_os)
  register_plot_pdf_download(output, input, "survCompPlot_tt2Line_censored", draw_survCompPlot_tt2Line_censored)
  register_plot_pdf_download(output, input, "clin_distribution", draw_clin_distribution_ggplot)

  # Statistical Significance Table
  output$significance_table <- renderDT({
    combined_data <- filtered_data$combined
    clinical_features <- get_clinical_feature_choices(clinical_data)
    continuous_features <- get_continuous_features()

    significance_table <- create_significance_table(
      data = combined_data,
      clinical_features = clinical_features,
      continuous_features = continuous_features
    )

    return(significance_table)
  }, options = list(
    pageLength = 15,
    scrollX = TRUE,
    columnDefs = list(
      list(targets = 5, width = "80px", className = "dt-center") # Significance column
    )
  ))

  # --------------- Cox-PH ----------------
  # Build design frame for model based on inputs
  .build_cox_data <- function() {
    src <- input$cox_data_source
    df <- switch(src,
                 "Cohort 1" = filtered_data$cohort1,
                 "Cohort 2" = filtered_data$cohort2,
                 "Both"     = filtered_data$combined,
                 filtered_data$combined)

    df <- as.data.frame(df)  # avoid data.table class surprises
    if ("cohort" %in% names(df)) {
      df$cohort <- factor(.cohort_labels()[as.character(df$cohort)], levels = unname(.cohort_labels()))
    }

    ep <- .cox_endpoint_map(input$cox_endpoint)
    req(ep$time %in% names(df), ep$event %in% names(df))

    # TPM gene covariate
    gene_term <- NULL
    if (!is.null(input$cox_gene) && nzchar(input$cox_gene)) {
      g <- input$cox_gene
      inter <- intersect(df$Tumor_Sample_Barcode, colnames(bulkseq_tpm))
      df <- df[df$Tumor_Sample_Barcode %in% inter, , drop = FALSE]
      gvec <- as.numeric(bulkseq_tpm[g, df$Tumor_Sample_Barcode])

      if (identical(input$cox_gene_mode, "continuous")) {
        df$gene_expr <- log2(gvec + 1)
        gene_term <- "gene_expr"
      } else {
        med <- stats::median(gvec, na.rm = TRUE)
        df$gene_group <- factor(ifelse(gvec > med, "High", "Low"), levels = c("Low", "High"))
        gene_term <- "gene_group"
      }
    }

    # Covariates (keep only those that actually exist)
    covars <- input$cox_covars
    if (length(covars)) covars <- intersect(covars, names(df))
    allow_cohort <- identical(src, "Both")
    if (allow_cohort && isTRUE(input$cox_use_cohort)) covars <- unique(c("cohort", covars))
    if (!is.null(gene_term)) covars <- c(covars, gene_term)
    covars <- unique(covars)

    # rows complete for all needed columns
    needed <- unique(c(ep$time, ep$event, covars))
    needed <- needed[needed %in% names(df)]
    df <- df[stats::complete.cases(df[, needed, drop = FALSE]), , drop = FALSE]

    # coerce char/logical to factor
    df <- .coerce_for_cox(df, needed)

    list(df = df, time = ep$time, event = ep$event,
         covars = covars, gene_term = gene_term,
         ep_label = ep$label)
  }


  observe({
    covar_choices <- sort(get_clinical_feature_choices(clinical_data))
    updateSelectizeInput(session, "cox_covars", choices = covar_choices, server = TRUE)
    updateSelectizeInput(
      session, "cox_gene",
      choices = c("None" = "", sort(rownames(bulkseq_tpm))),
      server = TRUE
    )
    allow_cohort <- identical(input$cox_data_source, "Both")
    strata_choices <- c("None", covar_choices)
    if (allow_cohort) strata_choices <- c(strata_choices, "cohort")
    updateSelectInput(session, "cox_strata", choices = unique(strata_choices))
  })


  # When user isn't fitting on both cohorts, disable the cohort indicator
  observe({
    allow_cohort <- identical(input$cox_data_source, "Both")
    if (!allow_cohort && isTRUE(input$cox_use_cohort)) {
      updateCheckboxInput(session, "cox_use_cohort", value = FALSE)
    }
    shinyjs::toggleState("cox_use_cohort", condition = allow_cohort)
  })

  cox_state <- reactiveValues(fit = NULL, data = NULL, formula = NULL, ep_label = NULL)

  observeEvent(input$fit_cox, {
    dat <- .build_cox_data()

    # never pass "None" into strata()
    strata_var <- if (identical(input$cox_strata, "None")) NULL else input$cox_strata

    zero_var <- vapply(dat$df[, dat$covars, drop = FALSE], function(x) {
      ux <- unique(x); ux <- ux[!is.na(ux)]
      length(ux) <= 1
    }, logical(1))
    if (any(zero_var)) {
      dropped <- names(zero_var)[zero_var]
      dat$covars <- setdiff(dat$covars, dropped)
      showNotification(paste("Dropped constant covariates:", paste(dropped, collapse = ", ")),
                       type = "warning")
    }

    # --- prune covariates that cause separation/degeneracy
    pr <- .prune_separating_covars(dat$df, dat$covars, dat$event)
    if (length(pr$drop)) {
      showNotification(
        paste("Removed covariates with empty cells / zero variance:",
              paste(pr$drop, collapse = ", ")),
        type = "warning"
      )
    }
    dat$covars <- pr$keep
    validate(need(length(dat$covars) > 0 || !is.null(strata_var),
                  "No usable covariates left after pruning."))

    # re-compute required columns and keep complete rows only
    needed <- unique(c(dat$time, dat$event, dat$covars, if (!is.null(strata_var)) strata_var))
    dat$df  <- dat$df[stats::complete.cases(dat$df[, needed, drop = FALSE]), , drop = FALSE]

    dat$df <- droplevels(dat$df)

    validate(need(nrow(dat$df) >= 10, "Not enough rows after filtering to fit a Cox model."))
    validate(need(length(dat$covars) > 0 || !is.null(strata_var),
                  "No usable covariates left (all were constant/missing in this subset)."))

    has_event <- any(dat$df[[dat$event]] == 1, na.rm = TRUE)
    validate(need(has_event, "No events in the selected endpoint for this subset."))

    # build RHS
    rhs <- if (length(dat$covars)) paste(dat$covars, collapse = " + ") else "1"
    if (!is.null(strata_var)) rhs <- paste(rhs, paste0("+ strata(", strata_var, ")"))
    fml <- as.formula(paste0("survival::Surv(", dat$time, ", ", dat$event, ") ~ ", rhs))

    message("Fitting formula: ", deparse(fml))
    message("Rows: ", nrow(dat$df), " | Covars: ",
            if (length(dat$covars)) paste(dat$covars, collapse = ", ") else "(intercept only)")

    fit <- tryCatch(
      survival::coxph(fml, data = dat$df, ties = "efron", model = TRUE, x = TRUE),
      error = function(e) { showNotification(e$message, type = "error"); NULL }
    )
    validate(need(!is.null(fit), "Cox fit failed."))

    cox_state$fit      <- fit
    cox_state$data     <- dat$df
    cox_state$formula  <- fml
    cox_state$ep_label <- dat$ep_label
  })



  output$cox_formula <- renderText({
    req(cox_state$formula)
    paste("Model:", deparse(cox_state$formula))
  })

  output$cox_table <- renderDT({
    req(cox_state$fit)
    datatable(.cox_tidy_table(cox_state$fit),
              options = list(pageLength = 10, scrollX = TRUE))
  })

  output$download_cox_table <- downloadHandler(
    filename = function() paste0("cox_results_", Sys.Date(), ".csv"),
    content = function(file) {
      req(cox_state$fit)
      write.csv(.cox_tidy_table(cox_state$fit), file, row.names = FALSE)
    }
  )

  draw_cox_forest <- function() {
    req(cox_state$fit, cox_state$data)
    ok <- all(is.finite(coef(cox_state$fit)))

    if (ok) {
      p <- try(
        survminer::ggforest(
          cox_state$fit,
          data = cox_state$data,
          fontsize = 1.3   # increase table/label text inside ggforest
        ),
        silent = TRUE
      )

      if (!inherits(p, "try-error")) {
        p <- p +
          ggplot2::theme(
            text = ggplot2::element_text(size = 16),
            axis.text = ggplot2::element_text(size = 14),
            axis.title = ggplot2::element_text(size = 16),
            plot.title = ggplot2::element_text(size = 18, face = "bold")
          )
        return(p)
      }
    }

    # fallback if infinities or ggforest choked
    .safe_forest(cox_state$fit) +
      ggplot2::theme(
        text = ggplot2::element_text(size = 16),
        axis.text = ggplot2::element_text(size = 14),
        axis.title = ggplot2::element_text(size = 16),
        plot.title = ggplot2::element_text(size = 18, face = "bold")
      )
  }

  output$cox_forest <- renderPlot({
    draw_cox_forest()
  })

  register_plot_pdf_download(output, input, "cox_forest", draw_cox_forest)

  # WGS -------------------------
  draw_mafSummary_g1 <- function() {
    plotmafSummary(maf = filtered_data[["cohort1_maf"]], addStat = 'median', titvRaw = FALSE)
    invisible(NULL)
  }

  draw_mafSummary_g2 <- function() {
    plotmafSummary(maf = filtered_data[["cohort2_maf"]], addStat = 'median', titvRaw = FALSE)
    invisible(NULL)
  }

  draw_oncoplot_g1 <- function() {
    oncoplot(maf = filtered_data[["cohort1_maf"]], top = 10)
    invisible(NULL)
  }

  draw_oncoplot_g2 <- function() {
    oncoplot(maf = filtered_data[["cohort2_maf"]], top = 10)
    invisible(NULL)
  }

  draw_lollipopPlot_g1 <- function() {
    req(input$gene_search_lollipop_g1)
    req(nzchar(input$gene_search_lollipop_g1))
    lollipop_variant_plot(
      maf = filtered_data[["cohort1_maf"]],
      gene = input$gene_search_lollipop_g1,
      cohort_label = .cohort_labels()[["Cohort1"]]
    )
  }

  draw_lollipopPlot_g2 <- function() {
    req(input$gene_search_lollipop_g2)
    req(nzchar(input$gene_search_lollipop_g2))
    lollipop_variant_plot(
      maf = filtered_data[["cohort2_maf"]],
      gene = input$gene_search_lollipop_g2,
      cohort_label = .cohort_labels()[["Cohort2"]]
    )
  }

  draw_interactionPlot_g1 <- function() {
    cohort_selected <- "cohort1_maf"
    genes_g1 <- input$gene_search_inter_g1

    if (is.null(genes_g1)) {
      somaticInteractions(maf = filtered_data[[cohort_selected]], top = 25, pvalue = c(0.05, 0.1))
    }

    if (length(genes_g1) >= 1 && length(genes_g1) <= 4) {
      gene_summary <- getGeneSummary(filtered_data[[cohort_selected]])
      gene_summary <- gene_summary[order(-gene_summary$MutatedSamples), ]

      num_gene_need <- 5 - length(genes_g1)
      top_genes <- head(gene_summary$Hugo_Symbol, 5)
      tmp_genes <- head(setdiff(top_genes, genes_g1), num_gene_need)

      somaticInteractions(maf = filtered_data[[cohort_selected]], genes = c(genes_g1, tmp_genes), pvalue = c(0.05, 0.1))
    }

    if (length(genes_g1) >= 5) {
      somaticInteractions(maf = filtered_data[[cohort_selected]], genes = genes_g1, pvalue = c(0.05, 0.1))
    }
    invisible(NULL)
  }

  draw_interactionPlot_g2 <- function() {
    cohort_selected <- "cohort2_maf"
    genes_g2 <- input$gene_search_inter_g2

    if (is.null(genes_g2)) {
      somaticInteractions(maf = filtered_data[[cohort_selected]], top = 25, pvalue = c(0.05, 0.1))
    }

    if (length(genes_g2) >= 1 && length(genes_g2) <= 4) {
      gene_summary <- getGeneSummary(filtered_data[[cohort_selected]])
      gene_summary <- gene_summary[order(-gene_summary$MutatedSamples), ]

      num_gene_need <- 5 - length(genes_g2)
      top_genes <- head(gene_summary$Hugo_Symbol, 5)
      tmp_genes <- head(setdiff(top_genes, genes_g2), num_gene_need)

      somaticInteractions(maf = filtered_data[[cohort_selected]], genes = c(genes_g2, tmp_genes), pvalue = c(0.05, 0.1))
    }

    if (length(genes_g2) >= 5) {
      somaticInteractions(maf = filtered_data[[cohort_selected]], genes = genes_g2, pvalue = c(0.05, 0.1))
    }
    invisible(NULL)
  }

  draw_mafCompForestPlot <- function() {
    req(length(unique(filtered_data$combined$cohort)) == 2)
    labels <- .cohort_labels()

    g1.vs.g2 <- mafCompare(m1 = filtered_data$cohort1_maf, m2 = filtered_data$cohort2_maf,
                           m1Name = labels[["Cohort1"]], m2Name = labels[["Cohort2"]], minMut = 5)
    g1.vs.g2$results <- g1.vs.g2$results %>%
      arrange(pval)

    g1.vs.g2$results <- g1.vs.g2$results[1:10]
    forestPlot(mafCompareRes = g1.vs.g2,
               pVal = 0.05,
               color = c("blue", "red"))
    invisible(NULL)
  }

  draw_mafCompOncoPlot <- function() {
    req(input$gene_search_maf)
    labels <- .cohort_labels()
    genes <- input$gene_search_maf
    if (length(genes) > 0) {
      coOncoplot(m1 = filtered_data$cohort1_maf, m2 = filtered_data$cohort2_maf,
                 m1Name = labels[["Cohort1"]], m2Name = labels[["Cohort2"]], genes = genes, removeNonMutated = TRUE)
    }
    invisible(NULL)
  }

  draw_mafCompBarPlot <- function() {
    req(input$gene_search_maf)
    labels <- .cohort_labels()
    genes <- input$gene_search_maf
    if (length(genes) > 0) {
      coBarplot(m1 = filtered_data$cohort1_maf, m2 = filtered_data$cohort2_maf,
                m1Name = labels[["Cohort1"]], m2Name = labels[["Cohort2"]], genes = genes)
    }
    invisible(NULL)
  }

  # Draw MAF summary plot
  output$mafSummary_g1 <- renderPlot({
    draw_mafSummary_g1()
  })

  output$mafSummary_g2 <- renderPlot({
    draw_mafSummary_g2()
  })

  # Draw oncoplot
  output$oncoplot_g1 <- renderPlot({
    draw_oncoplot_g1()
  })

  output$oncoplot_g2 <- renderPlot({
    draw_oncoplot_g2()
  })

  # Draw lollipop plot based on gene search
  output$lollipopPlot_g1 <- renderPlotly({
    lollipop_variant_plotly(draw_lollipopPlot_g1())
  })
  output$lollipopPlot_g2 <- renderPlotly({
    lollipop_variant_plotly(draw_lollipopPlot_g2())
  })

  # Draw somatic interaction plot
  output$interactionPlot_g1 <- renderPlot({
    draw_interactionPlot_g1()
  })
  output$interactionPlot_g2 <- renderPlot({
    draw_interactionPlot_g2()
  })

  output$mafCompTable <- renderDataTable({
    req(length(unique(filtered_data$combined$cohort)) == 2)
    labels <- .cohort_labels()

    g1.vs.g2 <- mafCompare(m1 = filtered_data$cohort1_maf, m2 = filtered_data$cohort2_maf,
                           m1Name = labels[["Cohort1"]], m2Name = labels[["Cohort2"]], minMut = 5)
    g1.vs.g2$results <- g1.vs.g2$results %>%
      arrange(pval)

    g1.vs.g2$results
  })

  # MAF Comparison
  output$mafCompForestPlot <- renderPlot({
    draw_mafCompForestPlot()
  })

  output$mafCompOncoPlot <- renderPlot({
    draw_mafCompOncoPlot()
  })

  output$mafCompBarPlot <- renderPlot({
    draw_mafCompBarPlot()
  })

  register_plot_pdf_download(output, input, "mafSummary_g1", draw_mafSummary_g1)
  register_plot_pdf_download(output, input, "mafSummary_g2", draw_mafSummary_g2)
  register_plot_pdf_download(output, input, "oncoplot_g1", draw_oncoplot_g1)
  register_plot_pdf_download(output, input, "oncoplot_g2", draw_oncoplot_g2)
  register_plot_pdf_download(output, input, "lollipopPlot_g1", draw_lollipopPlot_g1)
  register_plot_pdf_download(output, input, "lollipopPlot_g2", draw_lollipopPlot_g2)
  register_plot_pdf_download(output, input, "interactionPlot_g1", draw_interactionPlot_g1)
  register_plot_pdf_download(output, input, "interactionPlot_g2", draw_interactionPlot_g2)
  register_plot_pdf_download(output, input, "mafCompForestPlot", draw_mafCompForestPlot)
  register_plot_pdf_download(output, input, "mafCompOncoPlot", draw_mafCompOncoPlot)
  register_plot_pdf_download(output, input, "mafCompBarPlot", draw_mafCompBarPlot)

  # BulkRNA-seq Distribution ----------
  draw_tpm_distr <- function() {
    data <- preprocessed_bulkseq_data()
    tpm_distr_dens(data$combined_bulkseq_tpm, .with_display_cohort(data$clinical_combined), data$gene_interested, "bulkRNAseq")
  }

  draw_tpm_distr_boxplot <- function() {
    data <- preprocessed_bulkseq_data()
    tpm_boxplot(data$combined_bulkseq_tpm, .with_display_cohort(data$clinical_combined), data$gene_interested, "bulkRNAseq")
  }

  draw_tpm_survCompPlot <- function(cohort_selected, cohorting_method) {
    req(input$gene_search_bulk_distr, cohorting_method)
    gene_interested <- input$gene_search_bulk_distr
    selected_clinical <- filtered_data[[cohort_selected]]
    selected_bulkseq_tpm <- bulkseq_tpm[, colnames(bulkseq_tpm) %in% selected_clinical$Tumor_Sample_Barcode, drop = FALSE]
    selected_clinical <- selected_clinical[selected_clinical$Tumor_Sample_Barcode %in% colnames(selected_bulkseq_tpm), ]
    gene_tpm <- selected_bulkseq_tpm[gene_interested, , drop = FALSE]

    tpm_distr_survival(gene_tpm, selected_clinical, cohorting_method)
  }

  # Distribution
  output$tpm_distr <- renderPlotly({
    ggplotly(draw_tpm_distr(), tooltip = "text") %>%
      layout(hovermode = "x")
  })

  # TPM Boxplot
  output$tpm_distr_boxplot <- renderPlot({
    draw_tpm_distr_boxplot()
  })

  # Quantile table
  output$quantile_table_cohort1 <- renderDT({
    req(input$gene_search_bulk_distr)
    cohort_selected <- "cohort1"
    clinical_selected <- filtered_data[[cohort_selected]]
    selected_bulkseq_tpm <- bulkseq_tpm[,colnames(bulkseq_tpm) %in% clinical_selected$Tumor_Sample_Barcode]
    gene_interested <- input$gene_search_bulk_distr

    # Generate table
    distr_table <- tpm_distr_table(selected_bulkseq_tpm, gene_interested)
    datatable(distr_table, caption = htmltools::tags$caption(.cohort_labels()[["Cohort1"]]), options = list(pageLength = 10, autoWidth = TRUE))
  })

  output$quantile_table_cohort2 <- renderDT({
    req(input$gene_search_bulk_distr)
    cohort_selected <- "cohort2"
    clinical_selected <- filtered_data[[cohort_selected]]
    selected_bulkseq_tpm <- bulkseq_tpm[,colnames(bulkseq_tpm) %in% clinical_selected$Tumor_Sample_Barcode]
    gene_interested <- input$gene_search_bulk_distr

    # Generate table
    distr_table <- tpm_distr_table(selected_bulkseq_tpm, gene_interested)
    datatable(distr_table, caption = htmltools::tags$caption(.cohort_labels()[["Cohort2"]]), options = list(pageLength = 10, autoWidth = TRUE))
  })

  # Survival curve for TPM quantile and mean
  output$tpm_survCompPlot_g1<- renderPlot({
    draw_tpm_survCompPlot("cohort1", input$cohorting_method_tpm_g1)
  })
  output$tpm_survCompPlot_g2<- renderPlot({
    draw_tpm_survCompPlot("cohort2", input$cohorting_method_tpm_g2)
  })

  register_plot_pdf_download(output, input, "tpm_distr", draw_tpm_distr)
  register_plot_pdf_download(output, input, "tpm_distr_boxplot", draw_tpm_distr_boxplot)
  register_plot_pdf_download(output, input, "tpm_survCompPlot_g1", function() draw_tpm_survCompPlot("cohort1", input$cohorting_method_tpm_g1))
  register_plot_pdf_download(output, input, "tpm_survCompPlot_g2", function() draw_tpm_survCompPlot("cohort2", input$cohorting_method_tpm_g2))

  # Reactive value to hold DESeq2 results
  deseq2_results <- reactiveValues(result = NULL)

  # Trigger DESeq2 analysis when the button is clicked
  observeEvent(input$start_deseq2, {
    withProgress(message = "Running DESeq2 analysis...", value = 0.5, {
      results <- process_deseq2(filtered_data, bulkseq, min_counts = 10, min_samples = 5)
      deseq2_results$result <- results
    })
  })

  # Reactive expression to apply user thresholds
  filtered_degs <- reactive({
    req(deseq2_results$result)
    degs <- deseq2_results$result
    degs$significant <- ifelse(degs$padj < input$p_threshold & degs$log2FoldChange > input$fc_threshold, "Up-regulated",
                               ifelse(degs$padj < input$p_threshold & degs$log2FoldChange < -input$fc_threshold, "Down-regulated", "Not Significant"))
    degs
  })

  draw_bulkVolcano <- function() {
    req(filtered_degs())
    deseq2_volcano_ggplot(filtered_degs(), input$p_threshold, input$fc_threshold)
  }

  # Render volcano plot
  # interactive volcano
  output$bulkVolcano <- plotly::renderPlotly({
    req(filtered_degs())
    deseq2_volcano_plotly(filtered_degs(), input$p_threshold, input$fc_threshold, src = "volcano")
  })

  register_plot_pdf_download(output, input, "bulkVolcano", draw_bulkVolcano)

  # Render DEGs table
  output$DEGs_table <- renderDT({
    req(filtered_degs())
    degs <- filtered_degs() %>%
      dplyr::select(baseMean, log2FoldChange, padj, significant)

    degs$baseMean <- round(degs$baseMean, 2)
    degs$log2FoldChange <- round(degs$log2FoldChange, 2)

    datatable(degs, options = list(pageLength = 10, autoWidth = TRUE))
  })

  output$download_DEGs <- downloadHandler(
    filename = function() {
      paste0("DEGs_table_", Sys.Date(), ".csv")
    },
    content = function(file) {
      degs <- filtered_degs() %>%
        dplyr::select(baseMean, log2FoldChange, padj, significant)

      write.csv(degs, file)
    }
  )

  # BulkRNA-seq Enrichment Analysis -----------
  ssgsea_data <- reactive({
    req(filtered_data)
    clinical_combined <- filtered_data$combined
    results <- compute_significant_gene_sets(ssgsea_result_ca, clinical_combined)
    return(results)
  })

  observe({
    results <- ssgsea_data()
    gene_sets <- unique(results$wilcox_results[order(results$wilcox_results$adjusted_p_value), ]$GeneSet)
    updateSelectizeInput(session, "selected_gene_sets",
                         choices = gene_sets,
                         selected = head(gene_sets, 3), # default selection
                         server = TRUE)
  })

  draw_ssgsea_violin <- function() {
    results <- ssgsea_data()
    req(input$selected_gene_sets)
    labels <- .cohort_labels()

    plot_data <- results$long_data %>%
      dplyr::filter(GeneSet %in% input$selected_gene_sets) %>%
      dplyr::mutate(Cohort = factor(labels[as.character(Cohort)], levels = unname(labels))) %>%
      droplevels()

    # Base plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Cohort, y = EnrichmentScore, fill = Cohort)) +
      ggplot2::geom_violin(trim = FALSE) +
      ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA) +
      ggplot2::facet_wrap(~GeneSet, scales = "free_y", ncol=4) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    # ---- Add asterisks + p-values from wilcoxon results ----
    if (!is.null(results$wilcox_results) && nrow(results$wilcox_results) > 0) {
      ann <- results$wilcox_results %>%
        dplyr::filter(GeneSet %in% input$selected_gene_sets) %>%
        dplyr::select(GeneSet, p = adjusted_p_value) %>%
        dplyr::left_join(
          plot_data %>%
            dplyr::group_by(GeneSet) %>%
            dplyr::summarise(
              y_top = max(EnrichmentScore, na.rm = TRUE),
              rng   = diff(range(EnrichmentScore, na.rm = TRUE)),
              .groups = "drop"
            ),
          by = "GeneSet"
        ) %>%
        dplyr::mutate(
          stars = cut(
            p,
            breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
            labels = c("***", "**", "*", "ns")
          ),
          label = paste0(stars, " (adj. p=", formatC(p, format = "e", digits = 2), ")"),
          x = 1.5,  # center between two cohorts on the x-axis
          y = y_top + 0.07 * ifelse(is.finite(rng) & rng > 0, rng, abs(y_top) + 1)
        ) %>%
        dplyr::filter(!is.na(p))

      p <- p + ggplot2::geom_text(
        data = ann,
        ggplot2::aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        size = 4.2,
        fontface = "bold"
      )
    }

    p
  }

  output$ssgsea_violin <- renderPlot({
    draw_ssgsea_violin()
  })

  register_plot_pdf_download(output, input, "ssgsea_violin", draw_ssgsea_violin)

  output$ssgsea_table <- renderDT({
    results <- ssgsea_data()
    table_data <- results$wilcox_results

    # Format numeric columns without losing tiny values
    numeric_cols <- sapply(table_data, is.numeric)
    table_data <- table_data[order(table_data$adjusted_p_value), ]

    table_data[numeric_cols] <- lapply(table_data[numeric_cols], function(x) {
      ifelse(x < 0.001, formatC(x, format = "e", digits = 2), formatC(x, format = "f", digits = 3))
    })
    datatable(table_data, options = list(pageLength = 10, autoWidth = TRUE))
  })

  # scRNA-seq ----------------------
  draw_sc_celltype_boxplot <- function() {
    preprocessed_sc_meta <- preprocessed_sc_meta()
    cohort_info <- .with_display_cohort(preprocessed_sc_meta$cohort_info)
    celltype_boxplot(cohort_info, sc_meta)
  }

  draw_sc_celltype_proportion <- function() {
    preprocessed_sc_meta <- preprocessed_sc_meta()
    cohort_info <- .with_display_cohort(preprocessed_sc_meta$cohort_info)
    celltype_proportion(cohort_info, sc_meta)
  }

  draw_sc_cellcycle_hist <- function() {
    preprocessed_sc_meta <- preprocessed_sc_meta()
    cohort_info <- .with_display_cohort(preprocessed_sc_meta$cohort_info)
    celltypes_interested <- input$celltypes_interested
    cell_cycle_hist(cohort_info, sc_meta, celltypes_interested)
  }

  output$sc_celltype_boxplot <- renderPlot({
    draw_sc_celltype_boxplot()
  })

  output$sc_celltype_proportion <- renderPlot({
    draw_sc_celltype_proportion()
  })

  output$sc_cellcycle_hist <- renderPlot({
    draw_sc_cellcycle_hist()
  })

  register_plot_pdf_download(output, input, "sc_celltype_boxplot", draw_sc_celltype_boxplot)
  register_plot_pdf_download(output, input, "sc_celltype_proportion", draw_sc_celltype_proportion)
  register_plot_pdf_download(output, input, "sc_cellcycle_hist", draw_sc_cellcycle_hist)

  # -------------------- PSEUDOBULK (by cell type) ----------------------------
  # Update gene list whenever cell type changes
  observe({
    req(input$pseudo_celltype)
    pb_norm <- .get_pb_matrix(pseudo_bulk_norm, input$pseudo_celltype)
    updateSelectizeInput(session, "gene_search_pseudo_distr",
                         choices = unique(rownames(pb_norm)), selected = "KRAS", server = TRUE)
  })

  preprocessed_pseudobulk_data <- reactive({
    req(input$pseudo_celltype)

    pb_norm   <- .get_pb_matrix(pseudo_bulk_norm,   input$pseudo_celltype)
    pb_counts <- .get_pb_matrix(pseudo_bulk_counts, input$pseudo_celltype)

    clinical_combined <- filtered_data$combined

    aln_norm   <- .align_pb_to_clinical(pb_norm,   clinical_combined)
    aln_counts <- .align_pb_to_clinical(pb_counts, clinical_combined)

    list(
      pb_norm   = aln_norm$pb_tpm,      # genes x TSB (normalized values for plots)
      pb_counts = aln_counts$pb_tpm,    # genes x TSB (integer counts for DE)
      clinical  = aln_norm$clinical_aligned,  # aligned clin rows matching columns
      celltype  = input$pseudo_celltype,
      total     = ncol(pb_norm)
    )
  })

  draw_pseudo_norm_distr <- function() {
    dat <- preprocessed_pseudobulk_data()
    req(input$gene_search_pseudo_distr)
    tpm_distr_dens(dat$pb_norm, .with_display_cohort(dat$clinical), input$gene_search_pseudo_distr, "scRNAseq")
  }

  draw_pseudo_norm_distr_boxplot <- function() {
    dat <- preprocessed_pseudobulk_data()
    req(input$gene_search_pseudo_distr)
    tpm_boxplot(dat$pb_norm, .with_display_cohort(dat$clinical), input$gene_search_pseudo_distr, "scRNAseq")
  }

  draw_pseudo_norm_survCompPlot <- function(cohort_selected, cohorting_method) {
    dat <- preprocessed_pseudobulk_data()
    req(input$gene_search_pseudo_distr, cohorting_method)
    clin <- filtered_data[[cohort_selected]]
    pb <- dat$pb_norm[, colnames(dat$pb_norm) %in% clin$Tumor_Sample_Barcode, drop = FALSE]
    clin <- clin[clin$Tumor_Sample_Barcode %in% colnames(pb), ]
    gene <- input$gene_search_pseudo_distr
    if (!gene %in% rownames(pb)) return(NULL)
    tpm_distr_survival(as.matrix(pb[gene, , drop = FALSE]), clin, cohorting_method)
  }

  # Counts summary card
  output$pseudoNum <- renderUI({
    dat <- preprocessed_pseudobulk_data()
    pb_clin <- dat$clinical
    num_total  <- dat$total
    num_cohort1 <- sum(pb_clin$cohort == "Cohort1")
    num_cohort2 <- sum(pb_clin$cohort == "Cohort2")
    .count_card(num_total, num_cohort1, num_cohort2, extra = paste0("<strong>Cell type:</strong> ", htmltools::htmlEscape(dat$celltype)))
  })

  # Distribution (density + histogram)
  output$pseudo_norm_distr <- renderPlotly({
    ggplotly(draw_pseudo_norm_distr(), tooltip = "text") %>% layout(hovermode = "x")
  })

  # Boxplot + Wilcoxon p‑value
  output$pseudo_norm_distr_boxplot <- renderPlot({
    draw_pseudo_norm_distr_boxplot()
  })


  # Quantile tables (split by cohort)
  output$pseudo_quantile_table_cohort1 <- renderDT({
    dat <- preprocessed_pseudobulk_data()
    req(input$gene_search_pseudo_distr)

    clin1 <- dat$clinical_combined[dat$clinical_combined$cohort == "Cohort1", ]
    pb1 <- dat$pb_tpm[, colnames(dat$pb_tpm) %in% clin1$Tumor_Sample_Barcode, drop = FALSE]

    tab <- tpm_distr_table(pb1, input$gene_search_pseudo_distr)
    datatable(tab, caption = htmltools::tags$caption(.cohort_labels()[["Cohort1"]]), options = list(pageLength = 10, autoWidth = TRUE))
  })

  output$pseudo_quantile_table_cohort1 <- renderDT({
    dat <- preprocessed_pseudobulk_data()
    req(input$gene_search_pseudo_distr)
    clin1 <- dat$clinical[dat$clinical$cohort == "Cohort1", ]
    pb1   <- dat$pb_norm[, colnames(dat$pb_norm) %in% clin1$Tumor_Sample_Barcode, drop = FALSE]
    tab <- tpm_distr_table(pb1, input$gene_search_pseudo_distr)
    datatable(tab, caption = htmltools::tags$caption(.cohort_labels()[["Cohort1"]]), options = list(pageLength = 10, autoWidth = TRUE))
  })

  output$pseudo_quantile_table_cohort2 <- renderDT({
    dat <- preprocessed_pseudobulk_data()
    req(input$gene_search_pseudo_distr)
    clin2 <- dat$clinical[dat$clinical$cohort == "Cohort2", ]
    pb2   <- dat$pb_norm[, colnames(dat$pb_norm) %in% clin2$Tumor_Sample_Barcode, drop = FALSE]
    tab <- tpm_distr_table(pb2, input$gene_search_pseudo_distr)
    datatable(tab, caption = htmltools::tags$caption(.cohort_labels()[["Cohort2"]]), options = list(pageLength = 10, autoWidth = TRUE))
  })

  # Survival curves by expression split
  output$pseudo_norm_survCompPlot_g1 <- renderPlot({
    draw_pseudo_norm_survCompPlot("cohort1", input$cohorting_method_pseudo_g1)
  })

  output$pseudo_norm_survCompPlot_g2 <- renderPlot({
    draw_pseudo_norm_survCompPlot("cohort2", input$cohorting_method_pseudo_g2)
  })

  register_plot_pdf_download(output, input, "pseudo_norm_distr", draw_pseudo_norm_distr)
  register_plot_pdf_download(output, input, "pseudo_norm_distr_boxplot", draw_pseudo_norm_distr_boxplot)
  register_plot_pdf_download(output, input, "pseudo_norm_survCompPlot_g1", function() draw_pseudo_norm_survCompPlot("cohort1", input$cohorting_method_pseudo_g1))
  register_plot_pdf_download(output, input, "pseudo_norm_survCompPlot_g2", function() draw_pseudo_norm_survCompPlot("cohort2", input$cohorting_method_pseudo_g2))

  # Differential analysis
  pseudo_diff_results <- reactiveValues(result = NULL)

  observeEvent(input$start_pseudo_diff, {
    dat <- preprocessed_pseudobulk_data()
    withProgress(message = sprintf("Running DESeq2 (pseudobulk: %s)...", dat$celltype), value = 0.5, {
      # DE on counts:
      .check_duplicate_samples(dat$pb_counts, dat$clinical, context = sprintf("pseudobulk (%s) DESeq2", dat$celltype))
      res <- pseudobulk_diff_counts(dat$pb_counts, dat$clinical)
      pseudo_diff_results$result <- res
    })
  })

  # Thresholded results
  filtered_pseudo_degs <- reactive({
    req(pseudo_diff_results$result)
    res <- pseudo_diff_results$result
    thr_p  <- input$p_threshold_pseudo
    thr_fc <- input$fc_threshold_pseudo
    res$significant <- ifelse(!is.na(res$padj) & res$padj < thr_p & res$log2FoldChange >  thr_fc, "Up-regulated",
                              ifelse(!is.na(res$padj) & res$padj < thr_p & res$log2FoldChange < -thr_fc, "Down-regulated", "Not Significant"))
    res
  })

  draw_pseudoVolcano <- function() {
    req(filtered_pseudo_degs())
    deseq2_volcano_ggplot(filtered_pseudo_degs(), input$p_threshold_pseudo, input$fc_threshold_pseudo)
  }

  # Volcano
  output$pseudoVolcano <- plotly::renderPlotly({
    req(filtered_pseudo_degs())
    deseq2_volcano_plotly(filtered_pseudo_degs(), input$p_threshold_pseudo, input$fc_threshold_pseudo, src = "pb_volcano")
  })

  register_plot_pdf_download(output, input, "pseudoVolcano", draw_pseudoVolcano)

  # Table + download
  output$pseudo_DEGs_table <- renderDT({
    req(filtered_pseudo_degs())
    d <- filtered_pseudo_degs()[, c("baseMean","log2FoldChange","padj","significant")]
    d$baseMean <- round(d$baseMean, 3)
    d$log2FoldChange <- round(d$log2FoldChange, 3)
    datatable(d, options = list(pageLength = 10, autoWidth = TRUE))
  })

  output$download_pseudo_DEGs <- downloadHandler(
    filename = function() paste0("PseudoBulk_DEGs_", input$pseudo_celltype, "_", Sys.Date(), ".csv"),
    content = function(file) {
      d <- filtered_pseudo_degs()[, c("baseMean", "log2FoldChange", "padj", "significant")]
      write.csv(d, file)
    }
  )

})
