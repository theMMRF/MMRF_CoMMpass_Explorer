# Server function ------
shinyServer(function(input, output, session) {
  shinyjs::useShinyjs()
  
  # --- user-defined cohorts  --------------------------------------------------
  user_cohorts <- reactiveValues(
    c1_public = character(), c2_public = character(),
    c1_tsb = character(), c2_tsb = character(),
    c1_unmatched = character(), c2_unmatched = character()
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
  
  # Status + previews
  output$cohort1_status <- renderText({
    sprintf("Matched: %d  |  Unmatched: %d", length(user_cohorts$c1_tsb), length(user_cohorts$c1_unmatched))
  })
  output$cohort2_status <- renderText({
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
  
  # Disable remove when no rules exist
  observe({
    toggleState("remove_mut_row_cohort1", condition = mut_row_counter$cohort1 > 0)
    toggleState("remove_mut_row_cohort2", condition = mut_row_counter$cohort2 > 0)
  })
  
  # Add new rows on click
  observeEvent(input$add_mut_row_cohort1, {
    isolate({
      mut_rule_cache$cohort1 <- lapply(1:mut_row_counter$cohort1, function(i) {
        list(
          gene = input[[paste0("gene_mut_", i, "_cohort1")]],
          state = input[[paste0("state_mut_", i, "_cohort1")]],
          logic = input[[paste0("logic_mut_", i, "_cohort1")]]
        )
      })
    })
    mut_row_counter$cohort1 <- mut_row_counter$cohort1 + 1
  })
  
  
  # Remove new rows
  observeEvent(input$remove_mut_row_cohort1, {
    if (mut_row_counter$cohort1 > 0) {
      # Save inputs from remaining rows before decrement
      isolate({
        mut_rule_cache$cohort1 <- lapply(1:(mut_row_counter$cohort1 - 1), function(i) {
          list(
            gene = input[[paste0("gene_mut_", i, "_cohort1")]],
            state = input[[paste0("state_mut_", i, "_cohort1")]],
            logic = input[[paste0("logic_mut_", i, "_cohort1")]]
          )
        })
      })
      mut_row_counter$cohort1 <- mut_row_counter$cohort1 - 1
    }
  })
  
  
  observeEvent(input$add_mut_row_cohort2, {
    isolate({
      mut_rule_cache$cohort2 <- lapply(1:mut_row_counter$cohort2, function(i) {
        list(
          gene = input[[paste0("gene_mut_", i, "_cohort2")]],
          state = input[[paste0("state_mut_", i, "_cohort2")]],
          logic = input[[paste0("logic_mut_", i, "_cohort2")]]
        )
      })
    })
    mut_row_counter$cohort2 <- mut_row_counter$cohort2 + 1
  })
  
  
  observeEvent(input$remove_mut_row_cohort2, {
    if (mut_row_counter$cohort2 > 0) {
      isolate({
        mut_rule_cache$cohort2 <- lapply(1:(mut_row_counter$cohort2 - 1), function(i) {
          list(
            gene = input[[paste0("gene_mut_", i, "_cohort2")]],
            state = input[[paste0("state_mut_", i, "_cohort2")]],
            logic = input[[paste0("logic_mut_", i, "_cohort2")]]
          )
        })
      })
      
      mut_row_counter$cohort2 <- mut_row_counter$cohort2 - 1
    }
  })
  
  # Render dynamic rule UI
  output$mutation_rules_cohort1 <- renderUI({
    if (mut_row_counter$cohort1 == 0) return(NULL)
    
    lapply(1:mut_row_counter$cohort1, function(i) {
      cached <- if (length(mut_rule_cache$cohort1) >= i) mut_rule_cache$cohort1[[i]] else NULL
      
      gene_val <- if (!is.null(cached)) cached$gene else NULL
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
          choices = unique(maf_data@data$Hugo_Symbol),
          selected = gene_val,
          options = list(maxOptions = 100),
          width = "100%"
        ),
        br()
      )
    })
  })
  
  output$mutation_rules_cohort2 <- renderUI({
    if (mut_row_counter$cohort2 == 0) return(NULL)
    
    lapply(1:mut_row_counter$cohort2, function(i) {
      cached <- if (length(mut_rule_cache$cohort2) >= i) mut_rule_cache$cohort2[[i]] else NULL
      
      gene_val <- if (!is.null(cached)) cached$gene else NULL
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
          choices = unique(maf_data@data$Hugo_Symbol),
          selected = gene_val,
          options = list(maxOptions = 100),
          width = "100%"
        ),
        br()
      )
    })
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
      inputId    = paste0(c("sex_filter_", "race_filter_", "stage_filter_", "risk_filter_", "cyto_risk_filter_",
                            "rna_subtype_filter_", "cna_subtype_filter_", "triplet_filter_", "asct_filter_"),
                          cohort_id),
      filter_key = c("sex","race","stage","risk","cyto_risk",
                     "rna_subtype","cna_subtype","triplet","asct"),
      column     = c("Sex","Race","ISS","IMWG_Risk_Class","Skerget_Cytogenetic_High_Risk",
                     "Skerget_RNA_Subtype_Name","Skerget_CNA_Subtype_Name","Triplet_First","ASCT_First"),
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
  
  # Core updater for one cohort
  .update_picker_counts_for_cohort <- function(cohort_id) {
    # current filters as entered (WITHOUT pressing Apply)
    filters <- c(
      get_cohort_filters(input, cohort_id, "clinical"),
      get_cohort_filters(input, cohort_id, "molecular")
    )
    
    # Build mapping
    mp <- .make_picker_map(cohort_id)
    
    # For every picker, recompute counts conditioned on OTHER filters
    for (i in seq_len(nrow(mp))) {
      input_id   <- mp$inputId[i]
      filter_key <- mp$filter_key[i]
      colname    <- mp$column[i]
      
      # Skip if this input doesn't exist yet (UI not mounted)
      if (!(input_id %in% names(input))) next
      
      # Exclude this one filter from the preview subset
      filters_excl <- filters
      filters_excl[[filter_key]] <- NULL
      
      preview <- filter_cohort_data(copy(clinical_data), filters_excl)
      
      # Use all possible values from the full dataset (so zero-count values still show)
      all_vals <- sort(unique(na.omit(as.character(clinical_data[[colname]]))))
      if (length(all_vals) == 0) next
      
      # Tally counts in preview subset
      tbl <- table(as.character(preview[[colname]]), useNA = "no")
      
      # Build labels "value (count)"
      labels <- vapply(all_vals,
                       function(v) paste0(v, " (", .safe_count(tbl, v), ")"),
                       FUN.VALUE = character(1))
      
      # Keep values the same, only change displayed labels
      choices_named <- setNames(all_vals, labels)
      
      # Preserve existing selections if still valid
      current_sel <- input[[input_id]]
      current_sel <- intersect(current_sel, all_vals)
      
      shinyWidgets::updatePickerInput(
        session,
        inputId = input_id,
        choices = choices_named,
        selected = current_sel
      )
    }
  }
  
  # Re-run the updater whenever relevant inputs change
  # (This listens broadly to that cohort's filter inputs.)
  observe({
    # Touch all cohort1 filter inputs to create reactivity
    dummy <- list(
      input$sex_filter_cohort1, input$race_filter_cohort1, input$stage_filter_cohort1,
      input$risk_filter_cohort1, input$cyto_risk_filter_cohort1,
      input$rna_subtype_filter_cohort1, input$cna_subtype_filter_cohort1,
      input$triplet_filter_cohort1, input$asct_filter_cohort1,
      input$chr_1q21_gain_filter_cohort1, input$chr_1q21_amp_filter_cohort1,
      input$chr_13q14_del_filter_cohort1, input$chr_13q34_del_filter_cohort1,
      input$chr_17p13_del_filter_cohort1, input$diploidy_filter_cohort1,
      input$chromothripsis_filter_cohort1, input$t_11_14_filter_cohort1,
      input$t_4_14_filter_cohort1, input$maf_filter_cohort1,
      input$apobec_filter_cohort1, input$tp53_filter_cohort1, input$tp53_ns_filter_cohort1,
      input$age_cohort1  # include age so categorical counts reflect age subset, too
    )
    .update_picker_counts_for_cohort("cohort1")
  })
  
  observe({
    # Same for cohort2
    dummy <- list(
      input$sex_filter_cohort2, input$race_filter_cohort2, input$stage_filter_cohort2,
      input$risk_filter_cohort2, input$cyto_risk_filter_cohort2,
      input$rna_subtype_filter_cohort2, input$cna_subtype_filter_cohort2,
      input$triplet_filter_cohort2, input$asct_filter_cohort2,
      input$chr_1q21_gain_filter_cohort2, input$chr_1q21_amp_filter_cohort2,
      input$chr_13q14_del_filter_cohort2, input$chr_13q34_del_filter_cohort2,
      input$chr_17p13_del_filter_cohort2, input$diploidy_filter_cohort2,
      input$chromothripsis_filter_cohort2, input$t_11_14_filter_cohort2,
      input$t_4_14_filter_cohort2, input$maf_filter_cohort2,
      input$apobec_filter_cohort2, input$tp53_filter_cohort2, input$tp53_ns_filter_cohort2,
      input$age_cohort2
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
  
  cohort1_filters <- reactive({
    c(
      get_cohort_filters(input, "cohort1", "clinical"),
      get_cohort_filters(input, "cohort1", "molecular"),
      get_cohort_filters(input, "cohort1", "gene")
    )
  })
  
  
  filtered_data_cohort1 <- eventReactive(input$apply_filters, {
    clinical_filtered <- filter_cohort_data(copy(clinical_data), cohort1_filters())
    allowed <- if (exists("user_cohorts") && length(user_cohorts$c1_tsb)) user_cohorts$c1_tsb else clinical_filtered$Tumor_Sample_Barcode
    clinical_filtered <- clinical_filtered[clinical_filtered$Tumor_Sample_Barcode %in% allowed, ]
    
    clinical_error <- copy(clinical_filtered)
    
    # Gene mutation filter
    mutated_ids_cohort1 <- get_mutation_filtered_ids(isolate(input), "cohort1", mut_row_counter$cohort1)
    clinical_filtered <- clinical_filtered[clinical_filtered$Tumor_Sample_Barcode %in% mutated_ids_cohort1, ]
    
    # Gene expression filter
    clinical_filtered <- filter_by_gene_expression(
      clinical_data = clinical_filtered,
      gene = isolate(input$gene_expr_search_cohort1),
      threshold_type = isolate(input$expr_threshold_type_cohort1),
      min_value = isolate(input$gene_expr_min_cohort1),
      max_value = isolate(input$gene_expr_max_cohort1),
      min_percentile = isolate(input$gene_expr_percentile_min_cohort1),
      max_percentile = isolate(input$gene_expr_percentile_max_cohort1)
    )
    
    # Survival filter
    if (isTRUE(isolate(input$enable_survival_filter_cohort1))) {
      threshold_type <- isolate(input$surv_threshold_type_cohort1)
      
      if (threshold_type == "percentile") {
        clinical_filtered <- filter_by_survival(
          clinical_filtered,
          surv_var = isolate(input$surv_variable_cohort1),
          threshold_type = "percentile",
          min_percentile = isolate(input$surv_threshold_min_percentile_cohort1),
          max_percentile = isolate(input$surv_threshold_max_percentile_cohort1)
        )
      } else if (threshold_type == "value") {
        clinical_filtered <- filter_by_survival(
          clinical_filtered,
          surv_var = isolate(input$surv_variable_cohort1),
          threshold_type = "value",
          min_value = isolate(input$surv_threshold_min_value_cohort1),
          max_value = isolate(input$surv_threshold_max_value_cohort1)
        )
      }
    }
    
    if (nrow(clinical_filtered) == 0) {
      showNotification("No patients in Cohort 1 match the filters.", type = "error")
      return(clinical_error)
    }
    
    return(clinical_filtered)
  })
  
  
  # Cohort 2
  output$age_filter_cohort2 <- renderUI({
    min_age <- min(clinical_data$Age, na.rm = TRUE)
    max_age <- max(clinical_data$Age, na.rm = TRUE)
    sliderInput("age_cohort2", "Age", min = min_age, max = max_age, value = c(min_age, max_age))
  })
  
  cohort2_filters <- reactive({
    c(
      get_cohort_filters(input, "cohort2", "clinical"),
      get_cohort_filters(input, "cohort2", "molecular"),
      get_cohort_filters(input, "cohort2", "gene")
    )
  })
  
  filtered_data_cohort2 <- eventReactive(input$apply_filters, {
    clinical_filtered <- filter_cohort_data(copy(clinical_data), cohort2_filters())
    allowed <- if (exists("user_cohorts") && length(user_cohorts$c2_tsb)) user_cohorts$c2_tsb else clinical_filtered$Tumor_Sample_Barcode
    clinical_filtered <- clinical_filtered[clinical_filtered$Tumor_Sample_Barcode %in% allowed, ]
    
    clinical_error <- copy(clinical_filtered)
    
    # Gene mutation filter
    mutated_ids_cohort2 <- get_mutation_filtered_ids(isolate(input), "cohort2", mut_row_counter$cohort2)
    clinical_filtered <- clinical_filtered[clinical_filtered$Tumor_Sample_Barcode %in% mutated_ids_cohort2, ]
    
    # Gene expression filter
    clinical_filtered <- filter_by_gene_expression(
      clinical_data = clinical_filtered,
      gene = isolate(input$gene_expr_search_cohort2),
      threshold_type = isolate(input$expr_threshold_type_cohort2),
      min_value = isolate(input$gene_expr_min_cohort2),
      max_value = isolate(input$gene_expr_max_cohort2),
      min_percentile = isolate(input$gene_expr_percentile_min_cohort2),
      max_percentile = isolate(input$gene_expr_percentile_max_cohort2)
    )
    
    # Survival filter
    if (isTRUE(isolate(input$enable_survival_filter_cohort2))) {
      threshold_type <- isolate(input$surv_threshold_type_cohort2)
      
      if (threshold_type == "percentile") {
        clinical_filtered <- filter_by_survival(
          clinical_filtered,
          surv_var = isolate(input$surv_variable_cohort2),
          threshold_type = "percentile",
          min_percentile = isolate(input$surv_threshold_min_percentile_cohort2),
          max_percentile = isolate(input$surv_threshold_max_percentile_cohort2)
        )
      } else if (threshold_type == "value") {
        clinical_filtered <- filter_by_survival(
          clinical_filtered,
          surv_var = isolate(input$surv_variable_cohort2),
          threshold_type = "value",
          min_value = isolate(input$surv_threshold_min_value_cohort2),
          max_value = isolate(input$surv_threshold_max_value_cohort2)
        )
      }
    }
    
    
    
    if (nrow(clinical_filtered) == 0) {
      showNotification("No patients in Cohort 2 match the filters.", type = "error")
      return(clinical_error)
    }
    
    return(clinical_filtered)
  })
  
  
  
  # filtered_data that stores combined clinical data and cohort1, cohort2 data
  filtered_data <- reactiveValues(
    cohort1 = {
      data_cohort1 <- clinical_data
      data_cohort1$cohort <- "Cohort1"
      data_cohort1
    },
    cohort2 = {
      data_cohort2 <- clinical_data
      data_cohort2$cohort <- "Cohort2"
      data_cohort2
    },
    cohort1_maf = maf_data,
    cohort2_maf = maf_data,
    combined = {
      data_cohort1 <- clinical_data
      data_cohort1$cohort <- "Cohort1"
      data_cohort2 <- clinical_data
      data_cohort2$cohort <- "Cohort2"
      rbind(data_cohort1, data_cohort2)
    }
  )
  
  observeEvent(input$apply_filters, {
    data_cohort1 <- filtered_data_cohort1()
    data_cohort2 <- filtered_data_cohort2()
    
    data_cohort1$cohort <- "Cohort1"
    data_cohort2$cohort <- "Cohort2"
    
    patient_ids_cohort1 <- data_cohort1$Tumor_Sample_Barcode
    patient_ids_cohort2 <- data_cohort2$Tumor_Sample_Barcode
    
    combined_clinical <- rbind(data_cohort1, data_cohort2)
    
    # Update reactiveValues
    filtered_data$cohort1 <- data_cohort1
    filtered_data$cohort2 <- data_cohort2
    filtered_data$combined <- combined_clinical
    filtered_data$cohort1_maf <- subsetMaf(maf = maf_data, tsb = patient_ids_cohort1)
    filtered_data$cohort2_maf <- subsetMaf(maf = maf_data, tsb = patient_ids_cohort2)
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
    
    # Get patients in clinical and create public_id column for clinical
    patients_in_clinical <- sapply(strsplit(as.character(clinical_combined$Tumor_Sample_Barcode), "_"),
                                   function(x) paste(x[1], x[2], sep="_"))
    clinical_combined$public_id <- patients_in_clinical
    
    # Intersect
    inter_samples <- intersect(clinical_combined$public_id, sc_meta$public_id)
    num_sample <- length(intersect(clinical_data$public_id, sc_meta$immune_analysis_id))
    clinical_combined <- clinical_combined[clinical_combined$public_id %in% inter_samples,]
    
    cohort_info <- clinical_combined %>% 
      distinct(public_id, cohort) # Select only unique combinations of public_id and cohort from the clinical_combined
    
    list(
      num_sample = num_sample,
      inter_samples = inter_samples,
      clinical_combined = clinical_combined,
      cohort_info = cohort_info
    )
  })
  
  # Number of samples in each dataset after subsetting -------------------------
  # Clinical
  output$clinicalNum <- renderUI({
    num_total <- nrow(clinical_data)
    num_cohort1 <- nrow(filtered_data$cohort1)
    num_cohort2 <- nrow(filtered_data$cohort2)
    text <- sprintf(
      "<div style='padding:10px; border: 1px solid #ccc; border-radius: 5px; background-color:#f5f5f5;'>
     <strong>Total Samples:</strong> %d<br>
     <span style='color: #E41A1C;'>Cohort 1:</span> %d (%.2f%%)<br>
     <span style='color: #4DBBD5;'>Cohort 2:</span> %d (%.2f%%)
   </div>",
      num_total, num_cohort1, num_cohort1 / num_total * 100, num_cohort2, num_cohort2 / num_total * 100
    )
    
    HTML(text)
  })
  
  # Download filtered clinical data
  output$download_clinical <- downloadHandler(
    filename = function() {
      paste0("Filtered_Clinical_Data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(filtered_data$combined, file, row.names = FALSE)
    }
  )
  
  # MAF
  output$mafNum <- renderUI({
    num_total <- uniqueN(maf_data@data$Tumor_Sample_Barcode)
    num_cohort1 <- nrow(filtered_data$cohort1_maf@clinical.data)
    num_cohort2 <- nrow(filtered_data$cohort2_maf@clinical.data)
    text <- sprintf(
      "<div style='padding:10px; border: 1px solid #ccc; border-radius: 5px; background-color:#f5f5f5;'>
     <strong>Total Samples:</strong> %d<br>
     <span style='color: #E41A1C;'>Cohort 1:</span> %d (%.2f%%)<br>
     <span style='color: #4DBBD5;'>Cohort 2:</span> %d (%.2f%%)
   </div>",
      num_total, num_cohort1, num_cohort1 / num_total * 100, num_cohort2, num_cohort2 / num_total * 100
    )
    HTML(text)
  })
  
  # Transcriptomics
  output$bulkNum <- renderUI({
    bulk_clinical <- preprocessed_bulkseq_data()$clinical_combined
    num_sample <- preprocessed_bulkseq_data()$num_sample
    num_total <- length(num_sample)
    num_cohort1 <- nrow(bulk_clinical[bulk_clinical$cohort == "Cohort1", ])
    num_cohort2 <- nrow(bulk_clinical[bulk_clinical$cohort == "Cohort2", ])
    text <- sprintf(
      "<div style='padding:10px; border: 1px solid #ccc; border-radius: 5px; background-color:#f5f5f5;'>
     <strong>Total Samples:</strong> %d<br>
     <span style='color: #E41A1C;'>Cohort 1:</span> %d (%.2f%%)<br>
     <span style='color: #4DBBD5;'>Cohort 2:</span> %d (%.2f%%)
   </div>",
      num_total, num_cohort1, num_cohort1 / num_total * 100, num_cohort2, num_cohort2 / num_total * 100
    )
    HTML(text)
  })
  
  # scRNA-seq
  output$scNum <- renderUI({
    sc_clinical <- preprocessed_sc_meta()$cohort_info
    inter_samples <- preprocessed_sc_meta()$inter_samples
    num_total <- preprocessed_sc_meta()$num_sample
    num_cohort1 <- nrow(sc_clinical[sc_clinical$cohort == "Cohort1", ])
    num_cohort2 <- nrow(sc_clinical[sc_clinical$cohort == "Cohort2", ])
    text <- sprintf(
      "<div style='padding:10px; border: 1px solid #ccc; border-radius: 5px; background-color:#f5f5f5;'>
     <strong>Total Samples:</strong> %d<br>
     <span style='color: #E41A1C;'>Cohort 1:</span> %d (%.2f%%)<br>
     <span style='color: #4DBBD5;'>Cohort 2:</span> %d (%.2f%%)
   </div>",
      num_total, num_cohort1, num_cohort1 / num_total * 100, num_cohort2, num_cohort2 / num_total * 100
    )
    HTML(text)
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
  # Summary of clinical data
  output$summaryPlot_g1 <- renderPlot({
    generate_summary_plot("cohort1", filtered_data)
  })
  
  output$summaryPlot_g2 <- renderPlot({
    generate_summary_plot("cohort2", filtered_data)
  })
  
  # Draw survival curve comparison plot (PFS_censored)
  output$survCompPlot_pfs_censored <- renderPlot({
    req(filtered_data$combined)
    combined_clinical <- filtered_data$combined
    
    combined_surv <- Surv(combined_clinical$PFS_censored, combined_clinical$PFS_event)
    
    fit_combined <- do.call(survfit, list(combined_surv ~ cohort, data = combined_clinical))
    
    ggsurvplot(fit_combined, data = combined_clinical,
               # Core aesthetics
               palette = c("#E41A1C", "#4DBBD5"),  #red for Cohort1, teal for Cohort2
               linetype = c("solid", "solid"),
               size = 1,           # Line thickness
               
               # Statistical elements
               conf.int = TRUE,
               pval = TRUE,
               pval.coord = c(500, 0.1),
               
               # Labels and titles
               title = "",  
               xlab = "Progression-Free Survival (Days)",
               ylab = "Survival Probability",
               legend.title = "",
               legend.labs = c("Cohort1", "Cohort2"),
               
               # Risk table configuration
               risk.table = TRUE,
               risk.table.height = 0.25,
               risk.table.title = "Number at risk",
               risk.table.fontsize = 3.5,
               tables.theme = theme_cleantable(),
               
               # Formatting and theme
               ggtheme = theme_bw() + theme(
                 panel.grid.minor = element_blank(),
                 axis.title = element_text(face = "bold", size = 12),
                 axis.text = element_text(size = 10),
                 legend.position = "top",
                 legend.text = element_text(size = 10),
                 plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
               ),
               
               # Customize the axes
               break.time.by = 500,  # X-axis tick marks every 500 days
               surv.scale = "percent" # Y-axis in percentage
    )
  })
  
  # Draw survival curve comparison plot (OS)
  output$survCompPlot_os_censored <- renderPlot({
    req(filtered_data$combined)
    combined_clinical <- filtered_data$combined
    
    combined_surv <- Surv(combined_clinical$OS_censored, combined_clinical$OS_event)
    
    fit_combined <- do.call(survfit, list(combined_surv ~ cohort, data = combined_clinical))
    
    ggsurvplot(fit_combined, data = combined_clinical,
               # Core aesthetics
               palette = c("#E41A1C", "#4DBBD5"),  # red for Cohort1, teal for Cohort2
               linetype = c("solid", "solid"),
               size = 1,           # Line thickness
               
               # Statistical elements
               conf.int = TRUE,    # Show 95% confidence intervals
               pval = TRUE,        # Show p-value
               pval.coord = c(500, 0.1), # Position of p-value
               
               # Labels and titles
               title = "",  
               xlab = "Overall Survival (Days)",
               ylab = "Survival Probability",
               legend.title = "",
               legend.labs = c("Cohort1", "Cohort2"),
               
               # Risk table configuration
               risk.table = TRUE,
               risk.table.height = 0.25,
               risk.table.title = "Number at risk",
               risk.table.fontsize = 3.5,
               tables.theme = theme_cleantable(),
               
               # Formatting and theme
               ggtheme = theme_bw() + theme(
                 panel.grid.minor = element_blank(),
                 axis.title = element_text(face = "bold", size = 12),
                 axis.text = element_text(size = 10),
                 legend.position = "top",
                 legend.text = element_text(size = 10),
                 plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
               ),
               
               # Customize the axes
               break.time.by = 500,  # X-axis tick marks every 500 days
               surv.scale = "percent" # Y-axis in percentage
    )
  })
  
  # Distribution Plot
  output$clin_distribution <- renderPlotly({
    combined_data <- filtered_data$combined
    interested_feature <- input$clin_feature
    continuous_features <- get_continuous_features()
    
    create_distribution_plot(
      data = combined_data,
      feature = interested_feature,
      feature_label = input$clin_feature,
      continuous_features = continuous_features
    )
  })
  
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
  
  
  
  # WGS -------------------------
  # Draw MAF summary plot
  output$mafSummary_g1 <- renderPlot({
    cohort_selected <- "cohort1"
    cohort_selected <- paste0(cohort_selected, "_maf")
    plotmafSummary(maf = filtered_data[[cohort_selected]], addStat = 'median', titvRaw = FALSE)
  })
  
  output$mafSummary_g2 <- renderPlot({
    cohort_selected <- "cohort2"
    cohort_selected <- paste0(cohort_selected, "_maf")
    plotmafSummary(maf = filtered_data[[cohort_selected]], addStat = 'median', titvRaw = FALSE)
  })
  
  # Draw oncoplot
  output$oncoplot_g1 <- renderPlot({
    cohort_selected <- "cohort1"
    cohort_selected <- paste0(cohort_selected, "_maf")
    oncoplot(maf = filtered_data[[cohort_selected]], top = 10)
  })
  
  output$oncoplot_g2 <- renderPlot({
    cohort_selected <- "cohort2"
    cohort_selected <- paste0(cohort_selected, "_maf")
    oncoplot(maf = filtered_data[[cohort_selected]], top = 10)
  })
  
  # Draw lollipop plot based on gene search
  output$lollipopPlot_g1 <- renderPlot({
    req(input$gene_search_lollipop_g1)
    cohort_selected <- "cohort1"
    cohort_selected <- paste0(cohort_selected, "_maf")
    
    if (input$gene_search_lollipop_g1 != "") {
      lollipopPlot(maf = filtered_data[[cohort_selected]], AACol = "HGVSp", gene = input$gene_search_lollipop_g1)
    }
  })
  output$lollipopPlot_g2 <- renderPlot({
    req(input$gene_search_lollipop_g2)
    cohort_selected <- "cohort2"
    cohort_selected <- paste0(cohort_selected, "_maf")
    
    if (input$gene_search_lollipop_g2 != "") {
      lollipopPlot(maf = filtered_data[[cohort_selected]], AACol = "HGVSp", gene = input$gene_search_lollipop_g2)
    }
  })
  
  # Draw somatic interaction plot
  output$interactionPlot_g1 <- renderPlot({
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
      
      somaticInteractions(maf = filtered_data[[cohort_selected]], genes = c(genes_g1, tmp_genes) , pvalue = c(0.05, 0.1))
    }
    
    if (length(genes_g1) >= 5) {
      somaticInteractions(maf = filtered_data[[cohort_selected]], genes = genes_g1, pvalue = c(0.05, 0.1))
    }
    
  })
  output$interactionPlot_g2 <- renderPlot({
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
      
      somaticInteractions(maf = filtered_data[[cohort_selected]], genes = c(genes_g2, tmp_genes) , pvalue = c(0.05, 0.1))
    }
    
    if (length(genes_g2) >= 5) {
      somaticInteractions(maf = filtered_data[[cohort_selected]], genes = genes_g2, pvalue = c(0.05, 0.1))
    }
  })
  
  output$mafCompTable <- renderDataTable({
    req(length(unique(filtered_data$combined$cohort)) == 2)
    
    g1.vs.g2 <- mafCompare(m1 = filtered_data$cohort1_maf, m2 = filtered_data$cohort2_maf,
                           m1Name = 'Cohort 1', m2Name = 'Cohort 2', minMut = 5)
    g1.vs.g2$results <- g1.vs.g2$results %>%
      arrange(pval)
    
    g1.vs.g2$results
  })
  
  # MAF Comparison
  output$mafCompForestPlot <- renderPlot({
    req(length(unique(filtered_data$combined$cohort)) == 2)
    
    g1.vs.g2 <- mafCompare(m1 = filtered_data$cohort1_maf, m2 = filtered_data$cohort2_maf,
                           m1Name = 'Cohort 1', m2Name = 'Cohort 2', minMut = 5)
    g1.vs.g2$results <- g1.vs.g2$results %>%
      arrange(pval)
    
    g1.vs.g2$results <- g1.vs.g2$results[1:10]
    forestPlot(mafCompareRes = g1.vs.g2, 
               pVal = 0.05, 
               color = c("blue", "red"))
    # forestPlot(mafCompareRes = g1.vs.g2, pVal = 0.05)
  })
  
  output$mafCompOncoPlot <- renderPlot({
    req(input$gene_search_maf)
    genes <- input$gene_search_maf
    if (length(genes) > 0) {
      coOncoplot(m1 = filtered_data$cohort1_maf, m2 = filtered_data$cohort2_maf,
                 m1Name = 'Cohort 1', m2Name = 'Cohort 2', genes = genes, removeNonMutated = TRUE)
    }
  })
  
  output$mafCompBarPlot <- renderPlot({
    req(input$gene_search_maf)
    genes <- input$gene_search_maf
    if (length(genes) > 0) {
      coBarplot(m1 = filtered_data$cohort1_maf, m2 = filtered_data$cohort2_maf,
                m1Name = 'Cohort 1', m2Name = 'Cohort 2', genes = genes)
    }
  })
  
  # BulkRNA-seq Distribution ----------
  # Distribution
  output$tpm_distr <- renderPlotly({
    data <- preprocessed_bulkseq_data()
    
    # Generate density plot
    distr_dens <- tpm_distr_dens(data$combined_bulkseq_tpm, data$clinical_combined, data$gene_interested, "bulkRNAseq")
    
    ggplotly(distr_dens, tooltip = "text") %>%
      layout(hovermode = "x")
  })
  
  # TPM Boxplot
  output$tpm_distr_boxplot <- renderPlot({
    data <- preprocessed_bulkseq_data()
    
    # Generate box plot
    tpm_box <- tpm_boxplot(data$combined_bulkseq_tpm, data$clinical_combined, data$gene_interested, "bulkRNAseq")
    print(tpm_box)
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
    datatable(distr_table, options = list(pageLength = 10, autoWidth = TRUE))
  })
  
  output$quantile_table_cohort2 <- renderDT({
    req(input$gene_search_bulk_distr)
    cohort_selected <- "cohort2"
    clinical_selected <- filtered_data[[cohort_selected]]
    selected_bulkseq_tpm <- bulkseq_tpm[,colnames(bulkseq_tpm) %in% clinical_selected$Tumor_Sample_Barcode]
    gene_interested <- input$gene_search_bulk_distr
    
    # Generate table
    distr_table <- tpm_distr_table(selected_bulkseq_tpm, gene_interested)
    datatable(distr_table, options = list(pageLength = 10, autoWidth = TRUE))
  })
  
  # Survival curve for TPM quantile and mean
  output$tpm_survCompPlot_g1<- renderPlot({
    req(input$gene_search_bulk_distr, input$cohorting_method_tpm_g1)
    cohort_selected <- "cohort1"
    cohorting_method <- input$cohorting_method_tpm_g1
    gene_interested <- input$gene_search_bulk_distr
    selected_clinical <- filtered_data[[cohort_selected]]
    selected_bulkseq_tpm <- bulkseq_tpm[,colnames(bulkseq_tpm) %in% selected_clinical$Tumor_Sample_Barcode]
    selected_clinical <- selected_clinical[selected_clinical$Tumor_Sample_Barcode %in% colnames(selected_bulkseq_tpm),]
    gene_tpm <- selected_bulkseq_tpm[gene_interested, ]
    
    tpm_distr_survival(gene_tpm, selected_clinical, cohorting_method)
  })
  output$tpm_survCompPlot_g2<- renderPlot({
    req(input$gene_search_bulk_distr, input$cohorting_method_tpm_g2)
    cohort_selected <- "cohort2"
    cohorting_method <- input$cohorting_method_tpm_g2
    gene_interested <- input$gene_search_bulk_distr
    selected_clinical <- filtered_data[[cohort_selected]]
    selected_bulkseq_tpm <- bulkseq_tpm[,colnames(bulkseq_tpm) %in% selected_clinical$Tumor_Sample_Barcode]
    selected_clinical <- selected_clinical[selected_clinical$Tumor_Sample_Barcode %in% colnames(selected_bulkseq_tpm),]
    gene_tpm <- selected_bulkseq_tpm[gene_interested, ]
    
    tpm_distr_survival(gene_tpm, selected_clinical, cohorting_method)
  })

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
  
  # Render volcano plot
  # interactive volcano
  output$bulkVolcano <- plotly::renderPlotly({
    req(filtered_degs())
    deseq2_volcano_plotly(filtered_degs(), input$p_threshold, input$fc_threshold, src = "volcano")
  })
  
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
  
  output$ssgsea_violin <- renderPlot({
    results <- ssgsea_data()
    req(input$selected_gene_sets)
    
    plot_data <- results$long_data %>%
      dplyr::filter(GeneSet %in% input$selected_gene_sets) %>%
      droplevels()
    
    # Base plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Cohort, y = EnrichmentScore, fill = Cohort)) +
      ggplot2::geom_violin(trim = FALSE) +
      ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA) +
      ggplot2::facet_wrap(~GeneSet, scales = "free_y", ncol=5) +
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
    
    print(p)
  })
  
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
  output$sc_celltype_boxplot <- renderPlot({
    preprocessed_sc_meta <- preprocessed_sc_meta()
    clinical_combined <- preprocessed_sc_meta$clinical_combined
    cohort_info <- preprocessed_sc_meta$cohort_info
    celltype_boxplot(cohort_info, sc_meta)
  })
  
  output$sc_celltype_proportion <- renderPlot({
    preprocessed_sc_meta <- preprocessed_sc_meta()
    clinical_combined <- preprocessed_sc_meta$clinical_combined
    cohort_info <- preprocessed_sc_meta$cohort_info
    
    celltype_proportion(cohort_info, sc_meta)
  })
  
  output$sc_cellcycle_hist <- renderPlot({
    preprocessed_sc_meta <- preprocessed_sc_meta()
    clinical_combined <- preprocessed_sc_meta$clinical_combined
    cohort_info <- preprocessed_sc_meta$cohort_info
    celltypes_interested <- input$celltypes_interested
    
    cell_cycle_hist(cohort_info, sc_meta, celltypes_interested)
  })
  
  # -------------------- PSEUDOBULK (by cell type) ----------------------------
  # Update gene list whenever cell type changes
  observe({
    req(input$pseudo_celltype)
    pb_mat <- .get_pb_matrix(pseudo_bulk, input$pseudo_celltype)
    updateSelectizeInput(session, "gene_search_pseudo_distr",
                         choices = unique(rownames(pb_mat)), selected = NULL, server = TRUE)
  })
  
  # Prepare pseudobulk data aligned to current filtered cohorts
  preprocessed_pseudobulk_data <- reactive({
    req(input$pseudo_celltype)
    pb_tpm <- .get_pb_matrix(pseudo_bulk, input$pseudo_celltype)
    
    # Use the same combined clinical
    clinical_combined <- filtered_data$combined
    aligned <- .align_pb_to_clinical(pb_tpm, clinical_combined)
    
    list(
      pb_tpm = aligned$pb_tpm,                       # genes x Tumor_Sample_Barcode
      clinical_combined = aligned$clinical_aligned,  # rows aligned to columns of pb_tpm
      celltype = input$pseudo_celltype,
      total = ncol(pb_tpm)
    )
  })
  
  # Counts summary card
  output$pseudoNum <- renderUI({
    dat <- preprocessed_pseudobulk_data()
    pb_clin <- dat$clinical_combined
    
    num_total <- dat$total
    num_cohort1 <- sum(pb_clin$cohort == "Cohort1")
    num_cohort2 <- sum(pb_clin$cohort == "Cohort2")
    
    text <- sprintf(
      "<div style='padding:10px; border: 1px solid #ccc; border-radius: 5px; background-color:#f5f5f5;'>
       <strong>Cell type:</strong> %s<br>
       <strong>Total Samples:</strong> %d<br>
       <span style='color: #E41A1C;'>Cohort 1:</span> %d (%.2f%%)<br>
       <span style='color: #4DBBD5;'>Cohort 2:</span> %d (%.2f%%)
       </div>",
      dat$celltype, num_total, num_cohort1, ifelse(num_total==0, 0, num_cohort1/num_total*100),
      num_cohort2, ifelse(num_total==0, 0, num_cohort2/num_total*100)
    )
    HTML(text)
  })
  
  # Distribution (density + histogram)
  output$pseudo_tpm_distr <- renderPlotly({
    dat <- preprocessed_pseudobulk_data()
    req(input$gene_search_pseudo_distr)
    
    p <- tpm_distr_dens(dat$pb_tpm, dat$clinical_combined, input$gene_search_pseudo_distr, "scRNAseq")
    ggplotly(p, tooltip = "text") %>% layout(hovermode = "x")
  })
  
  # Boxplot + Wilcoxon p‑value
  output$pseudo_tpm_distr_boxplot <- renderPlot({
    dat <- preprocessed_pseudobulk_data()
    req(input$gene_search_pseudo_distr)
    
    p <- tpm_boxplot(dat$pb_tpm, dat$clinical_combined, input$gene_search_pseudo_distr, "scRNAseq")
    print(p)
  })
  
  # Quantile tables (split by cohort)
  output$pseudo_quantile_table_cohort1 <- renderDT({
    dat <- preprocessed_pseudobulk_data()
    req(input$gene_search_pseudo_distr)
    
    clin1 <- dat$clinical_combined[dat$clinical_combined$cohort == "Cohort1", ]
    pb1 <- dat$pb_tpm[, colnames(dat$pb_tpm) %in% clin1$Tumor_Sample_Barcode, drop = FALSE]
    
    tab <- tpm_distr_table(pb1, input$gene_search_pseudo_distr)
    datatable(tab, options = list(pageLength = 10, autoWidth = TRUE))
  })
  
  output$pseudo_quantile_table_cohort2 <- renderDT({
    dat <- preprocessed_pseudobulk_data()
    req(input$gene_search_pseudo_distr)
    
    clin2 <- dat$clinical_combined[dat$clinical_combined$cohort == "Cohort2", ]
    pb2 <- dat$pb_tpm[, colnames(dat$pb_tpm) %in% clin2$Tumor_Sample_Barcode, drop = FALSE]
    
    tab <- tpm_distr_table(pb2, input$gene_search_pseudo_distr)
    datatable(tab, options = list(pageLength = 10, autoWidth = TRUE))
  })
  
  # Survival curves by expression split
  output$pseudo_tpm_survCompPlot_g1 <- renderPlot({
    dat <- preprocessed_pseudobulk_data()
    req(input$gene_search_pseudo_distr, input$cohorting_method_pseudo_g1)
    
    clin1 <- filtered_data$cohort1
    pb1 <- dat$pb_tpm[, colnames(dat$pb_tpm) %in% clin1$Tumor_Sample_Barcode, drop = FALSE]
    clin1 <- clin1[clin1$Tumor_Sample_Barcode %in% colnames(pb1), ]
    
    gene <- input$gene_search_pseudo_distr
    if (!gene %in% rownames(pb1)) return(NULL)
    
    tpm_distr_survival(as.matrix(pb1[gene, , drop = FALSE]), clin1, input$cohorting_method_pseudo_g1)
  })
  
  output$pseudo_tpm_survCompPlot_g2 <- renderPlot({
    dat <- preprocessed_pseudobulk_data()
    req(input$gene_search_pseudo_distr, input$cohorting_method_pseudo_g2)
    
    clin2 <- filtered_data$cohort2
    pb2 <- dat$pb_tpm[, colnames(dat$pb_tpm) %in% clin2$Tumor_Sample_Barcode, drop = FALSE]
    clin2 <- clin2[clin2$Tumor_Sample_Barcode %in% colnames(pb2), ]
    
    gene <- input$gene_search_pseudo_distr
    if (!gene %in% rownames(pb2)) return(NULL)
    
    tpm_distr_survival(as.matrix(pb2[gene, , drop = FALSE]), clin2, input$cohorting_method_pseudo_g2)
  })
  
  # Differential analysis (Wilcoxon + log2FC)
  pseudo_diff_results <- reactiveValues(result = NULL)
  
  observeEvent(input$start_pseudo_diff, {
    dat <- preprocessed_pseudobulk_data()
    withProgress(message = paste0("Running differential analysis (", dat$celltype, ")..."), value = 0.5, {
      res <- pseudobulk_diff(dat$pb_tpm, dat$clinical_combined)
      pseudo_diff_results$result <- res
    })
  })
  
  # Thresholded results
  filtered_pseudo_degs <- reactive({
    req(pseudo_diff_results$result)
    res <- pseudo_diff_results$result
    thr_p <- input$p_threshold_pseudo
    thr_fc <- input$fc_threshold_pseudo
    
    res$significant <- ifelse(!is.na(res$padj) & res$padj < thr_p & res$log2FoldChange > thr_fc, "Up-regulated",
                              ifelse(!is.na(res$padj) & res$padj < thr_p & res$log2FoldChange < -thr_fc, "Down-regulated", "Not Significant"))
    res
  })
  
  # Volcano
  output$pseudoVolcano <- renderPlot({
    req(filtered_pseudo_degs())
    print(deseq2_volcano(filtered_pseudo_degs(), input$p_threshold_pseudo, input$fc_threshold_pseudo))
  })
  
  # Table + download
  output$pseudo_DEGs_table <- renderDT({
    req(filtered_pseudo_degs())
    d <- filtered_pseudo_degs()[, c("baseMean", "log2FoldChange", "padj", "significant")]
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

