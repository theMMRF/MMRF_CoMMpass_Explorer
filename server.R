# Server function ------
shinyServer(function(input, output, session) {
  shinyjs::useShinyjs()
  
  # Load available genes for selection
  observe({
    gene_choices <- unique(maf_data@data$Hugo_Symbol)
    updateSelectizeInput(session, "gene_include_filter", choices = gene_choices, server = TRUE)
  })
  
  observe({
    gene_choices <- unique(rownames(bulkseq_tpm))
    updateSelectizeInput(session, "gene_expr_search", choices = gene_choices, server = TRUE)
  })
  
  observeEvent(input$mut_number, {
    if (!is.na(input$mut_number) && input$mut_number %% 1 != 0) {
      updateNumericInput(session, "mut_number", value = round(input$mut_number))
    }
  })
  
  # Filters interface
  # Group 1
  output$age_filter_group1 <- renderUI({
    min_age <- min(clinical_data$Age, na.rm = TRUE)
    max_age <- max(clinical_data$Age, na.rm = TRUE)
    sliderInput("age_group1", "Age", min = min_age, max = max_age, value = c(min_age, max_age))
  })
  
  # Clear Group 1 filters when the "Clear" button is clicked
  observeEvent(input$clear_group1, {
    shinyjs::reset("group1_filters") # Resets all inputs within the div for Group 1 filters
  })
  
  observeEvent(input$clear_group2, {
    shinyjs::reset("group2_filters")
  })
  
  group1_filters <- reactive({
    c(
      get_group_filters(input, "group1", "clinical"),
      get_group_filters(input, "group1", "molecular"),
      get_group_filters(input, "group1", "gene")
    )
  })
  
  filtered_data_group1 <- reactive({
    clinical_filtered <- filter_group_data(clinical_data, group1_filters())
    clinical_error <- copy(clinical_filtered)
    
    include_genes <- input$gene_include_filter
    clinical_filtered <- subset_by_gene_mutations(clinical_filtered, include_genes, "group1", input$mut_logic)

    if (dim(clinical_filtered)[1] == 0) {
      showNotification("No patients in Group 1 match the selected genes. Please reset the filter and select different genes.", type = "error")
      clinical_filtered <- copy(clinical_error)
    }
    
    gene_selected <- input$gene_expr_search
    threshold <- input$gene_expr_threshold
    clinical_filtered <- filter_by_gene_expression(clinical_filtered, gene_selected, threshold, "group1")
    
    if (isTRUE(input$enable_survival_filter)) {
      clinical_filtered <- filter_by_survival(
        clinical_filtered,
        input$surv_variable,
        input$surv_threshold_type,
        if (input$surv_threshold_type == "value") input$surv_threshold_value else input$surv_threshold_percentile,
        "group1"
      )
    }
    
    return(clinical_filtered)
  })
  
  # Group 2
  output$age_filter_group2 <- renderUI({
    min_age <- min(clinical_data$Age, na.rm = TRUE)
    max_age <- max(clinical_data$Age, na.rm = TRUE)
    sliderInput("age_group2", "Age", min = min_age, max = max_age, value = c(min_age, max_age))
  })
  
  group2_filters <- reactive({
    c(
      get_group_filters(input, "group2", "clinical"),
      get_group_filters(input, "group2", "molecular"),
      get_group_filters(input, "group2", "gene")
    )
  })
  
  filtered_data_group2 <- reactive({
    clinical_filtered <- filter_group_data(clinical_data, group2_filters())
    clinical_error <- copy(clinical_filtered)
    
    include_genes <- input$gene_include_filter
    clinical_filtered <- subset_by_gene_mutations(clinical_filtered, include_genes, "group2", input$mut_logic)
    
    if (dim(clinical_filtered)[1] == 0) {
      showNotification("No patients in Group 2 match the selected genes. Please reset the filter and select different genes.", type = "error")
      clinical_filtered <- copy(clinical_error)
    }
    
    gene_selected <- input$gene_expr_search
    threshold <- input$gene_expr_threshold
    clinical_filtered <- filter_by_gene_expression(clinical_filtered, gene_selected, threshold, "group2")
    
    if (isTRUE(input$enable_survival_filter)) {
      clinical_filtered <- filter_by_survival(
        clinical_filtered,
        input$surv_variable,
        input$surv_threshold_type,
        if (input$surv_threshold_type == "value") input$surv_threshold_value else input$surv_threshold_percentile,
        "group2"
      )
    }
    
    return(clinical_filtered)
  })
  
  # filtered_data that stores combined clinical data and group1, group2 data
  filtered_data <- reactiveValues(
    group1 = NULL,
    group2 = NULL,
    group1_maf = NULL,
    group2_maf = NULL,
    combined = NULL,
  )
  
  observeEvent({
    list(filtered_data_group1(), filtered_data_group2())
  }, {
    # Group 1
    data_group1 <- filtered_data_group1()
    patient_ids_group1 <- data_group1$Tumor_Sample_Barcode
    data_group1$group <- "Group1"
    
    # Group 2
    data_group2 <- filtered_data_group2()
    patient_ids_group2 <- data_group2$Tumor_Sample_Barcode
    data_group2$group <- "Group2"
    
    # Combine the two groups (# of samples may be greater than # of samples in total)
    combined_clinical <- rbind(data_group1, data_group2)
    
    # Update reactive values
    filtered_data$combined <- combined_clinical
    filtered_data$group1 <- data_group1
    filtered_data$group2 <- data_group2
    
    # MAF data processing for Group 1 and Group 2
    filtered_data$group1_maf <- subsetMaf(maf = maf_data, tsb = patient_ids_group1)
    filtered_data$group2_maf <- subsetMaf(maf = maf_data, tsb = patient_ids_group2)
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
    clinical_combined <- filtered_data$combined[, c("Tumor_Sample_Barcode", "group")]
    
    # Get patients in clinical and create public_id column for clinical
    patients_in_clinical <- sapply(strsplit(as.character(clinical_combined$Tumor_Sample_Barcode), "_"),
                                   function(x) paste(x[1], x[2], sep="_"))
    clinical_combined$public_id <- patients_in_clinical
    
    # Intersect
    inter_samples <- intersect(clinical_combined$public_id, sc_meta$public_id)
    num_sample <- intersect(clinical_data$Tumor_Sample_Barcode, colnames(sc_meta))
    clinical_combined <- clinical_combined[clinical_combined$public_id %in% inter_samples,]
    
    group_info <- clinical_combined %>% 
      distinct(public_id, group) # Select only unique combinations of public_id and group from the clinical_combined
    
    list(
      inter_samples = inter_samples,
      clinical_combined = clinical_combined,
      group_info = group_info
    )
  })
  
  # Number of samples in each dataset after subsetting -------------------------
  # Clinical
  output$clinicalNum <- renderUI({
    num_total <- nrow(clinical_data)
    num_group1 <- nrow(filtered_data$group1)
    num_group2 <- nrow(filtered_data$group2)
    text <- sprintf(
      "Total number of samples: %d.<br>Group 1 (n = %d (%0.2f%%)) vs. Group 2 (n = %d (%0.2f%%))", 
      num_total,
      num_group1, 
      num_group1 / num_total * 100, 
      num_group2, 
      num_group2 / num_total * 100
    )
    HTML(text)
  })
  
  output$clinicalNum <- renderUI({
    num_total <- nrow(clinical_data)
    num_group1 <- nrow(filtered_data$group1)
    num_group2 <- nrow(filtered_data$group2)
    text <- sprintf(
      "Total number of samples: %d.<br>Group 1 (n = %d (%0.2f%%)) vs. Group 2 (n = %d (%0.2f%%))", 
      num_total,
      num_group1, 
      num_group1 / num_total * 100, 
      num_group2, 
      num_group2 / num_total * 100
    )
    HTML(text)
  })
  
  # WGS
  output$mafNum <- renderUI({
    num_total <- nrow(clinical_data)
    num_group1 <- nrow(filtered_data$group1)
    num_group2 <- nrow(filtered_data$group2)
    text <- sprintf(
      "Total number of samples: %d.<br>Group 1 (n = %d (%0.2f%%)) vs. Group 2 (n = %d (%0.2f%%))", 
      num_total,
      num_group1, 
      num_group1 / num_total * 100, 
      num_group2, 
      num_group2 / num_total * 100
    )
    HTML(text)
  })
  
  # Transcriptomics
  output$bulkNum <- renderUI({
    bulk_clinical <- preprocessed_bulkseq_data()$clinical_combined
    num_sample <- preprocessed_bulkseq_data()$num_sample
    num_total <- length(num_sample)
    num_group1 <- nrow(bulk_clinical[bulk_clinical$group == "Group1", ])
    num_group2 <- nrow(bulk_clinical[bulk_clinical$group == "Group2", ])
    text <- sprintf(
      "Total number of samples: %d.<br>Group 1 (n = %d (%0.2f%%)) vs. Group 2 (n = %d (%0.2f%%))", 
      num_total,
      num_group1, 
      num_group1 / num_total * 100, 
      num_group2, 
      num_group2 / num_total * 100
    )
    HTML(text)
  })
  
  # scRNA-seq
  output$scNum <- renderUI({
    sc_clinical <- preprocessed_sc_meta()$group_info
    inter_samples <- preprocessed_sc_meta()$inter_samples
    num_total <- length(inter_samples)
    num_group1 <- nrow(sc_clinical[sc_clinical$group == "Group1", ])
    num_group2 <- nrow(sc_clinical[sc_clinical$group == "Group2", ])
    text <- sprintf(
      "Total number of samples: %d.<br>Group 1 (n = %d (%0.2f%%)) vs. Group 2 (n = %d (%0.2f%%))", 
      num_total,
      num_group1, 
      num_group1 / num_total * 100, 
      num_group2, 
      num_group2 / num_total * 100
    )
    HTML(text)
  })
  
  # Update gene choices for selectizeInput
  observe({
    updateSelectizeInput(session, "gene_search_lollipop_g1", choices = unique(maf_data@gene.summary$Hugo_Symbol), selected = "KRAS", server = TRUE)
    updateSelectizeInput(session, "gene_search_lollipop_g2", choices = unique(maf_data@gene.summary$Hugo_Symbol), selected = "KRAS", server = TRUE)
    updateSelectizeInput(session, "gene_search_maf", choices = unique(maf_data@gene.summary$Hugo_Symbol), server = TRUE)
    updateSelectizeInput(session, "gene_search_inter_g1", choices = unique(filtered_data[["group1_maf"]]@gene.summary$Hugo_Symbol), server = TRUE)
    updateSelectizeInput(session, "gene_search_inter_g2", choices = unique(filtered_data[["group2_maf"]]@gene.summary$Hugo_Symbol), server = TRUE)
    # updateSelectizeInput(session, "gene_search_bulk_heat", choices = unique(rownames(bulkseq_tpm)), server = TRUE)
    updateSelectizeInput(session, "gene_search_bulk_distr", choices = unique(rownames(bulkseq_tpm)), selected = "KRAS", server = TRUE)
  })
  
  # Summary ----------------------
  # Summary of clinical data
  output$summaryPlot_g1 <- renderPlot({
    generate_summary_plot("group1", filtered_data)
  })
  
  output$summaryPlot_g2 <- renderPlot({
    generate_summary_plot("group2", filtered_data)
  })
  
  # Draw survival curve comparison plot (PFS_censored)
  output$survCompPlot_pfs_censored <- renderPlot({
    req(filtered_data$combined)
    combined_clinical <- filtered_data$combined
    
    combined_surv <- Surv(combined_clinical$PFS_censored, combined_clinical$PFS_event)
    
    fit_combined <- do.call(survfit, list(combined_surv ~ group, data = combined_clinical))
    
    ggsurvplot(fit_combined, data = combined_clinical,
               # Core aesthetics
               palette = c("#E41A1C", "#4DBBD5"),  # Professional color scheme (red for Group1, teal for Group2)
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
               legend.labs = c("Group1", "Group2"),
               
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
    
    fit_combined <- do.call(survfit, list(combined_surv ~ group, data = combined_clinical))
    
    ggsurvplot(fit_combined, data = combined_clinical,
               # Core aesthetics
               palette = c("#E41A1C", "#4DBBD5"),  # Professional color scheme (red for Group1, teal for Group2)
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
               legend.labs = c("Group1", "Group2"),
               
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
  
  # Distribution
  output$clin_distribution <- renderPlotly({
    combined_data <- filtered_data$combined
    interested_feature <- input$clin_feature
    
    # Handle column names with numbers or special characters
    interested_feature <- paste0("`", interested_feature, "`")
    
    p <- ggplot(combined_data, aes_string(x = interested_feature, alpha = "group", fill = "group")) +
      geom_bar(position = "dodge") +
      labs(title = "", x = input$clin_feature, y = "Count") +
      scale_alpha_manual(values = c(Group1 = 1, Group2 = 0.5)) +
      scale_fill_manual(values = c(Group1 = "#E87D72", Group2 = "#5BAEB0")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggplotly(p) %>%
      layout(hovermode = "x")
  })
  
  # WGS -------------------------
  # Draw MAF summary plot
  output$mafSummary_g1 <- renderPlot({
    group_selected <- "group1"
    group_selected <- paste0(group_selected, "_maf")
    plotmafSummary(maf = filtered_data[[group_selected]], addStat = 'median', titvRaw = FALSE)
  })
  
  output$mafSummary_g2 <- renderPlot({
    group_selected <- "group2"
    group_selected <- paste0(group_selected, "_maf")
    plotmafSummary(maf = filtered_data[[group_selected]], addStat = 'median', titvRaw = FALSE)
  })
  
  # Draw oncoplot
  output$oncoplot_g1 <- renderPlot({
    group_selected <- "group1"
    group_selected <- paste0(group_selected, "_maf")
    oncoplot(maf = filtered_data[[group_selected]], top = 10)
  })
  
  output$oncoplot_g2 <- renderPlot({
    group_selected <- "group2"
    group_selected <- paste0(group_selected, "_maf")
    oncoplot(maf = filtered_data[[group_selected]], top = 10)
  })
  
  # Draw lollipop plot based on gene search
  output$lollipopPlot_g1 <- renderPlot({
    req(input$gene_search_lollipop_g1)
    group_selected <- "group1"
    group_selected <- paste0(group_selected, "_maf")
    
    if (input$gene_search_lollipop_g1 != "") {
      lollipopPlot(maf = filtered_data[[group_selected]], AACol = "HGVSp", gene = input$gene_search_lollipop_g1)
    }
  })
  output$lollipopPlot_g2 <- renderPlot({
    req(input$gene_search_lollipop_g2)
    group_selected <- "group2"
    group_selected <- paste0(group_selected, "_maf")
    
    if (input$gene_search_lollipop_g2 != "") {
      lollipopPlot(maf = filtered_data[[group_selected]], AACol = "HGVSp", gene = input$gene_search_lollipop_g2)
    }
  })
  
  
  # # Draw Ti/Tv ratio plot
  # output$titvPlot <- renderPlot({
  #   maf_titv <- titv(maf = filtered_data$selected_maf, plot = FALSE, useSyn = TRUE)
  #   plotTiTv(res = maf_titv)
  # })
  
  # # Draw Variant allele frequency plot
  # output$vafPlot <- renderPlot({
  #   group_selected <- filtered_data$group_selected
  #   group_selected <- paste0(group_selected, "_maf")
  #   plotVaf(maf = filtered_data[[group_selected]], vafCol = 'i_TumorVAF_WU')
  # })
  
  # Draw somatic interaction plot
  output$interactionPlot_g1 <- renderPlot({
    group_selected <- "group1_maf"
    genes_g1 <- input$gene_search_inter_g1
    
    if (is.null(genes_g1)) {
      somaticInteractions(maf = filtered_data[[group_selected]], top = 25, pvalue = c(0.05, 0.1))
    }
    
    if (length(genes_g1) >= 1 && length(genes_g1) <= 4) {
      gene_summary <- getGeneSummary(filtered_data[[group_selected]])
      gene_summary <- gene_summary[order(-gene_summary$MutatedSamples), ]
      
      num_gene_need <- 5 - length(genes_g1)
      top_genes <- head(gene_summary$Hugo_Symbol, 5)
      tmp_genes <- head(setdiff(top_genes, genes_g1), num_gene_need)
      
      somaticInteractions(maf = filtered_data[[group_selected]], genes = c(genes_g1, tmp_genes) , pvalue = c(0.05, 0.1))
    }
    
    if (length(genes_g1) >= 5) {
      somaticInteractions(maf = filtered_data[[group_selected]], genes = genes_g1, pvalue = c(0.05, 0.1))
    }
    
  })
  output$interactionPlot_g2 <- renderPlot({
    group_selected <- "group2_maf"
    genes_g2 <- input$gene_search_inter_g2
    
    if (is.null(genes_g2)) {
      somaticInteractions(maf = filtered_data[[group_selected]], top = 25, pvalue = c(0.05, 0.1))
    }
    
    if (length(genes_g2) >= 1 && length(genes_g2) <= 4) {
      gene_summary <- getGeneSummary(filtered_data[[group_selected]])
      gene_summary <- gene_summary[order(-gene_summary$MutatedSamples), ]
      
      num_gene_need <- 5 - length(genes_g2)
      top_genes <- head(gene_summary$Hugo_Symbol, 5)
      tmp_genes <- head(setdiff(top_genes, genes_g2), num_gene_need)
      
      somaticInteractions(maf = filtered_data[[group_selected]], genes = c(genes_g2, tmp_genes) , pvalue = c(0.05, 0.1))
    }
    
    if (length(genes_g2) >= 5) {
      somaticInteractions(maf = filtered_data[[group_selected]], genes = genes_g2, pvalue = c(0.05, 0.1))
    }
  })
  
  # MAF Comparison
  output$mafCompForestPlot <- renderPlot({
    req(length(unique(filtered_data$combined$group)) == 2)
    
    g1.vs.g2 <- mafCompare(m1 = filtered_data$group1_maf, m2 = filtered_data$group2_maf,
                           m1Name = 'Group 1', m2Name = 'Group 2', minMut = 5)
    g1.vs.g2$results <- g1.vs.g2$results[1:10]
    forestPlot(mafCompareRes = g1.vs.g2, pVal = 0.05)
  })
  
  output$mafCompOncoPlot <- renderPlot({
    req(input$gene_search_maf)
    genes <- input$gene_search_maf
    if (length(genes) > 0) {
      coOncoplot(m1 = filtered_data$group1_maf, m2 = filtered_data$group2_maf,
                 m1Name = 'Group 1', m2Name = 'Group 2', genes = genes, removeNonMutated = TRUE)
    }
  })
  
  output$mafCompBarPlot <- renderPlot({
    req(input$gene_search_maf)
    genes <- input$gene_search_maf
    if (length(genes) > 0) {
      coBarplot(m1 = filtered_data$group1_maf, m2 = filtered_data$group2_maf,
                m1Name = 'Group 1', m2Name = 'Group 2', genes = genes)
    }
  })
  
  # BulkRNA-seq Distribution ----------
  # Distribution
  output$tpm_distr <- renderPlotly({
    data <- preprocessed_bulkseq_data()
    
    # Generate density plot
    distr_dens <- tpm_distr_dens(data$combined_bulkseq_tpm, data$clinical_combined, data$gene_interested)
    
    ggplotly(distr_dens, tooltip = "text") %>%
      layout(hovermode = "x")
  })
  
  # TPM Boxplot
  output$tpm_distr_boxplot <- renderPlot({
    data <- preprocessed_bulkseq_data()
    
    # Generate box plot
    tpm_box <- tpm_boxplot(data$combined_bulkseq_tpm, data$clinical_combined, data$gene_interested)
    print(tpm_box)
  })
  
  # Quantile table
  output$quantile_table_group1 <- renderDT({
    req(input$gene_search_bulk_distr)
    group_selected <- "group1"
    clinical_selected <- filtered_data[[group_selected]]
    selected_bulkseq_tpm <- bulkseq_tpm[,colnames(bulkseq_tpm) %in% clinical_selected$Tumor_Sample_Barcode]
    gene_interested <- input$gene_search_bulk_distr
    
    # Generate table
    distr_table <- tpm_distr_table(selected_bulkseq_tpm, gene_interested)
    datatable(distr_table, options = list(pageLength = 10, autoWidth = TRUE))
  })
  
  output$quantile_table_group2 <- renderDT({
    req(input$gene_search_bulk_distr)
    group_selected <- "group2"
    clinical_selected <- filtered_data[[group_selected]]
    selected_bulkseq_tpm <- bulkseq_tpm[,colnames(bulkseq_tpm) %in% clinical_selected$Tumor_Sample_Barcode]
    gene_interested <- input$gene_search_bulk_distr
    
    # Generate table
    distr_table <- tpm_distr_table(selected_bulkseq_tpm, gene_interested)
    datatable(distr_table, options = list(pageLength = 10, autoWidth = TRUE))
  })
  
  # Survival curve for TPM quantile and mean
  output$tpm_survCompPlot_g1<- renderPlot({
    req(input$gene_search_bulk_distr, input$grouping_method_tpm_g1)
    group_selected <- "group1"
    grouping_method <- input$grouping_method_tpm_g1
    gene_interested <- input$gene_search_bulk_distr
    selected_clinical <- filtered_data[[group_selected]]
    selected_bulkseq_tpm <- bulkseq_tpm[,colnames(bulkseq_tpm) %in% selected_clinical$Tumor_Sample_Barcode]
    selected_clinical <- selected_clinical[selected_clinical$Tumor_Sample_Barcode %in% colnames(selected_bulkseq_tpm),]
    gene_tpm <- selected_bulkseq_tpm[gene_interested, ]
    
    tpm_distr_survival(gene_tpm, selected_clinical, grouping_method)
  })
  output$tpm_survCompPlot_g2<- renderPlot({
    req(input$gene_search_bulk_distr, input$grouping_method_tpm_g2)
    group_selected <- "group2"
    grouping_method <- input$grouping_method_tpm_g2
    gene_interested <- input$gene_search_bulk_distr
    selected_clinical <- filtered_data[[group_selected]]
    selected_bulkseq_tpm <- bulkseq_tpm[,colnames(bulkseq_tpm) %in% selected_clinical$Tumor_Sample_Barcode]
    selected_clinical <- selected_clinical[selected_clinical$Tumor_Sample_Barcode %in% colnames(selected_bulkseq_tpm),]
    gene_tpm <- selected_bulkseq_tpm[gene_interested, ]
    
    tpm_distr_survival(gene_tpm, selected_clinical, grouping_method)
  })
  
  # BulkRNA-seq DESeq2 and Heatmap ----------
  # # Heatmap for selected genes in bulkRNAseq
  # output$bulkHeat <- renderPlot({
  #   req(input$gene_search_bulk_heat)
  #   # filtered_data$combined has already been sorted by "group" (selected, unselected)
  #   clinical_sorted <- filtered_data$combined[filtered_data$combined$Tumor_Sample_Barcode %in% colnames(bulkseq_tpm),]
  #   
  #   genes <- input$gene_search_bulk_heat
  #   
  #   # Get expression data of checkpoints
  #   checkpoint_data <- bulkseq_tpm[genes,]
  #   # log2 transfromation for better visualization
  #   checkpoint_data_log <- log2(checkpoint_data + 1)
  #   checkpoint_data_log <- checkpoint_data_log[, clinical_sorted$Tumor_Sample_Barcode]
  #   
  #   # Calculate the percentage of samples expressing each gene
  #   gene_expression_percentage <- rowMeans(checkpoint_data > 0) * 100
  #   rownames(checkpoint_data_log) <- paste0(genes, " (", round(gene_expression_percentage, 2), "%)")
  #   
  #   annotation_col <- data.frame(Group = clinical_sorted$group)
  #   row.names(annotation_col) <- clinical_sorted$Tumor_Sample_Barcode
  # 
  #   # Plot the heatmap
  #   pheatmap(checkpoint_data_log, 
  #            scale = "row",  # Normalize
  #            color = colorRampPalette(c("blue", "white", "red"))(50),
  #            main = "Gene Expression Heatmap",
  #            show_colnames = FALSE,
  #            # cluster_cols = FALSE, # set to TRUE If we want to order by "group"
  #            annotation_col = annotation_col)
  # })
  
  # Perform DESeq2 analysis
  deseq2_results <- reactiveValues(result = NULL)
  observeEvent(input$start_deseq2, {
    results <- process_deseq2(filtered_data, bulkseq, min_counts=10, min_samples=5) # Default p=0.05, logFC=1.5
    deseq2_results$result <- results
  })
  
  # Observe changes in threshold inputs and update volcano plot
  observeEvent({
    input$p_threshold
    input$fc_threshold
  }, {
    # Volcano plot
    if (!is.null(deseq2_results$result)) {
      output$bulkVolcano <- renderPlot({
        degs <- deseq2_results$result
        volcano_plot <- deseq2_volcano(degs, input$p_threshold, input$fc_threshold)
        print(volcano_plot)
      })
    }
    
    # DEGs table
    if (!is.null(deseq2_results$result)) {
      degs <- deseq2_results$result %>%
        dplyr::select(c(1,2,6,7))
      
      degs$baseMean <- round(degs$baseMean, digits=2)
      degs$log2FoldChange <- round(degs$log2FoldChange, digits=2)
      
      degs$significant <- ifelse(degs$padj < input$p_threshold & degs$log2FoldChange > input$fc_threshold, "Up-regulated",
                                 ifelse(degs$padj < input$p_threshold & degs$log2FoldChange < -input$fc_threshold, "Down-regulated", "Not Significant"))
      
      output$DEGs_table <- renderDT({
        datatable(degs, options = list(pageLength = 10, autoWidth = TRUE))
      })
    }
  })
  
  # BulkRNA-seq Enrichment Analysis -----------
  ssgsea_data <- reactive({
    req(input$genesets, filtered_data)
    clinical_combined <- filtered_data$combined

    if (input$genesets == "Acharya") {
      # Use genesets defined by Chaitanya Acharya
      results <- compute_significant_gene_sets(ssgsea_result_ca, clinical_combined)
      
    } else if (input$genesets == "MSigDB-C2") {
      # Use MSigDB-C2 genesets
      results <- compute_significant_gene_sets(ssgsea_result_c2, clinical_combined)
    }
    return(results)
  })
  
  output$ssgsea_violin <- renderPlot({
    results <- ssgsea_data()
    p <- create_violin_plot(results$wilcox_results, results$long_data)
    print(p)
  })
  
  output$ssgsea_table <- renderDT({
    results <- ssgsea_data()
    datatable(results$wilcox_results, options = list(pageLength = 10, autoWidth = TRUE))
  })
  
  # scRNA-seq ----------------------
  output$sc_celltype_boxplot <- renderPlot({
    preprocessed_sc_meta <- preprocessed_sc_meta()
    clinical_combined <- preprocessed_sc_meta$clinical_combined
    group_info <- preprocessed_sc_meta$group_info
    celltype_boxplot(group_info, sc_meta)
  })
  
  output$sc_celltype_proportion <- renderPlot({
    preprocessed_sc_meta <- preprocessed_sc_meta()
    clinical_combined <- preprocessed_sc_meta$clinical_combined
    group_info <- preprocessed_sc_meta$group_info
    
    celltype_proportion(group_info, sc_meta)
  })
  
  output$sc_cellcycle_hist <- renderPlot({
    preprocessed_sc_meta <- preprocessed_sc_meta()
    clinical_combined <- preprocessed_sc_meta$clinical_combined
    group_info <- preprocessed_sc_meta$group_info
    celltypes_interested <- input$celltypes_interested
    
    cell_cycle_hist(group_info, sc_meta, celltypes_interested)
  })
  
  # Feature Distribution Analysis
  output$feature_distribution_plot <- renderPlot({
    req(input$feature_x, input$feature_y)
    
    create_distribution_stacked_barplot(
      data = clinical_data,
      x_feature = input$feature_x,
      y_features = input$feature_y
    )
  })
  
})