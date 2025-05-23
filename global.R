# Load dependencies -------
packages <- c("shiny", "shinydashboard", "shinyWidgets", "shinyjs",
              "tidyr", "dplyr", "DT", "data.table", "tibble", "gridExtra", "readxl",
              "ggplot2", "pheatmap", "ggrepel", "plotly",
              "maftools", "survival", "survminer", "limma", "DESeq2",
              "clusterProfiler", "org.Hs.eg.db", "biomaRt")

lapply(packages, library, character.only = TRUE)

# Load data ------
bulkseq <- readRDS("data/bulkseq_baseline.rds")
bulkseq_tpm <- readRDS("data/bulkseq_tpm_baseline.rds")
# clinical_data <- readRDS("data/clinical_data_cleaned.rds")
# clinical_data <- readRDS("data/clinical_data_cleaned_new_surv.rds")
clinical_data <- readRDS("data/clinical_data_n1411.rds")
maf_data <- readRDS("data/maf_data.rds")
# maf_data <- readRDS("data/maf_data_n937.rds")


sc_meta <- readRDS("data/scRNAseq_metadata.rds")
ssgsea_result_ca <- readRDS("data/ssgsea_result_ca.rds")
# ssgsea_result_c2 <- readRDS("data/ssgsea_result_c2.rds")
# bulkseq <- readRDS("../../data/commpass_explorer/bulkseq_baseline.rds")
# bulkseq_tpm <- readRDS("../../data/commpass_explorer/bulkseq_tpm_baseline.rds")
# clinical_data <- readRDS("../../data/commpass_explorer/clinical_data_cleaned.rds")
# maf_data <- readRDS("../../data/commpass_explorer/maf_data.rds")
# sc_meta <- readRDS("../../data/commpass_explorer/scRNAseq_metadata.rds")
# ssgsea_result_ca <- readRDS("../../data/commpass_explorer/ssgsea_result_ca.rds")
# ssgsea_result_c2 <- readRDS("../../data/commpass_explorer/ssgsea_result_c2.rds")

clinical_data$PFS_censored <- as.numeric(as.character(clinical_data$PFS_censored))
clinical_data$PFS_event <- as.numeric(as.character(clinical_data$PFS_event))
clinical_data$OS_censored <- as.numeric(as.character(clinical_data$OS_censored))
clinical_data$OS_event <- as.numeric(as.character(clinical_data$OS_event))

# Helper functions -----------
subset_by_gene_mutations <- function(clinical_data, include_genes = NULL, group, logic) {
  if (is.null(include_genes)) return(clinical_data)
  maf_data_table <- maf_data@data
  
  if (logic == "Or") {
    matched <- maf_data_table[Hugo_Symbol %in% include_genes, unique(Tumor_Sample_Barcode)]
    clinical_data <- if (group == "group1") {
      clinical_data[clinical_data$Tumor_Sample_Barcode %in% matched,]
    } else {
      clinical_data[!clinical_data$Tumor_Sample_Barcode %in% matched,]
    }
  }
  
  if (logic == "And") {
    sample_gene_count <- maf_data_table %>%
      filter(Hugo_Symbol %in% include_genes) %>%
      distinct(Tumor_Sample_Barcode, Hugo_Symbol) %>%
      group_by(Tumor_Sample_Barcode) %>%
      summarise(gene_count = n(), .groups = "drop")
    
    matched <- sample_gene_count %>%
      filter(gene_count == length(include_genes)) %>%
      pull(Tumor_Sample_Barcode)
    
    clinical_data <- if (group == "group1") {
      clinical_data[clinical_data$Tumor_Sample_Barcode %in% matched,]
    } else {
      clinical_data[!clinical_data$Tumor_Sample_Barcode %in% matched,]
    }
  }
  
  return(clinical_data)
}

filter_by_gene_expression <- function(clinical_data, gene=NULL, threshold=NULL, group) {
  if (is.null(gene) || is.na(threshold)) return(clinical_data)
  
  if (gene %in% unique(rownames(bulkseq_tpm)) && (group == "group1")) {
    patient_ids <- colnames(bulkseq_tpm[, bulkseq_tpm[gene, ] >= threshold, drop = FALSE])
    clinical_data <- clinical_data[clinical_data$Tumor_Sample_Barcode %in% patient_ids, ]
  }
  
  if (gene %in% unique(rownames(bulkseq_tpm)) && (group == "group2")) {
    patient_ids <- colnames(bulkseq_tpm[, bulkseq_tpm[gene, ] < threshold, drop = FALSE])
    clinical_data <- clinical_data[clinical_data$Tumor_Sample_Barcode %in% patient_ids, ]
  }
  
  return(clinical_data)
}


filter_by_survival <- function(clinical_data, surv_var, threshold_type, threshold_value, group) {
  if (!(surv_var %in% colnames(clinical_data))) return(clinical_data)
  
  surv_data <- clinical_data[[surv_var]]
  
  if (threshold_type == "percentile") {
    threshold <- quantile(surv_data, probs = threshold_value / 100, na.rm = TRUE)
  } else {
    threshold <- threshold_value
  }
  
  if (group == "group1") {
    filtered <- clinical_data[surv_data >= threshold, ]
  } else if (group == "group2") {
    filtered <- clinical_data[surv_data < threshold, ]
  } else {
    filtered <- clinical_data
  }
  
  return(filtered)
}


create_picker_input <- function(inputId, label, choices) {
  choices <- choices[!is.na(choices)]
  choices <- choices[choices!=""]
  pickerInput(
    inputId = inputId,
    label = label,
    choices = choices,
    options = list(
      `actions-box` = TRUE,
      `deselect-all-text` = "None",
      `select-all-text` = "Select All",
      `none-selected-text` = "None Selected"
    ),
    multiple = TRUE,
    selected = "All"
  )
}

create_group_filters_ui <- function(group_id, category) {
  tagList(
    div(id = paste0(group_id, "_", category, "_filters"),
        if (category == "clinical") {
          return(list(
            # Demographics
            create_picker_input(paste0("sex_filter_", group_id), "Sex", sort(unique(clinical_data$Sex))),
            create_picker_input(paste0("race_filter_", group_id), "Race", sort(unique(clinical_data$Race))),
            uiOutput(paste0("age_filter_", group_id)),
            
            # Clinical classification
            create_picker_input(paste0("stage_filter_", group_id), "ISS Stage", sort(unique(clinical_data$ISS))),
            create_picker_input(paste0("risk_filter_", group_id), "IMWG Risk Classification", sort(unique(clinical_data$IMWG_Risk_Class))),
            create_picker_input(paste0("cyto_risk_filter_", group_id), "Cytogenetic High Risk (Skerget)", sort(unique(clinical_data$Skerget_Cytogenetic_High_Risk))),
            
            # Subtypes
            create_picker_input(paste0("rna_subtype_filter_", group_id), "RNA Subtype (Skerget)", sort(unique(clinical_data$Skerget_RNA_Subtype_Name))),
            create_picker_input(paste0("cna_subtype_filter_", group_id), "CNA Subtype (Skerget)", sort(unique(clinical_data$Skerget_CNA_Subtype_Name))),
            
            # Treatment
            create_picker_input(paste0("triplet_filter_", group_id), "Triplet Firstline", sort(unique(clinical_data$Triplet_First))),
            create_picker_input(paste0("asct_filter_", group_id), "ASCT Firstline", sort(unique(clinical_data$ASCT_First)))
          ))
        } else if (category == "molecular") {
          return(list(
            # Chromosomal abnormalities
            create_picker_input(paste0("chr_1q21_gain_filter_", group_id), "1q21 Gain", sort(unique(clinical_data$chr_1q21_gain))),
            create_picker_input(paste0("chr_1q21_amp_filter_", group_id), "1q21 Amplification", sort(unique(clinical_data$chr_1q21_amp))),
            create_picker_input(paste0("chr_13q14_del_filter_", group_id), "13q14 Deletion", sort(unique(clinical_data$chr_13q14_del))),
            create_picker_input(paste0("chr_13q34_del_filter_", group_id), "13q34 Deletion", sort(unique(clinical_data$chr_13q34_del))),
            create_picker_input(paste0("chr_17p13_del_filter_", group_id), "17p13 Deletion", sort(unique(clinical_data$chr_17p13_del))),
            create_picker_input(paste0("diploidy_filter_", group_id), "Hyperdiploidy", sort(unique(clinical_data$Hyperdiploidy))),
            create_picker_input(paste0("chromothripsis_filter_", group_id), "Chromothripsis", sort(unique(clinical_data$chromothripsis))),
            
            # Translocations
            create_picker_input(paste0("t_11_14_filter_", group_id), "t(11;14)", sort(unique(clinical_data$t_11_14))),
            create_picker_input(paste0("t_4_14_filter_", group_id), "t(4;14)", sort(unique(clinical_data$t_4_14))),
            
            # Mutational markers
            create_picker_input(paste0("maf_filter_", group_id), "MAF/MAFB", sort(unique(clinical_data$MAF_MAFB))),
            create_picker_input(paste0("apobec_filter_", group_id), "APOBEC", sort(unique(clinical_data$APOBEC))),
            create_picker_input(paste0("tp53_filter_", group_id), "TP53 Functional Copies", sort(unique(clinical_data$TP53_Funct_Copies))),
            create_picker_input(paste0("tp53_ns_filter_", group_id), "TP53 Non-Synonymous Mutation Count", sort(unique(clinical_data$TP53_NS_Mut_Count)))
          ))
        }
    )
  )
}


get_group_filters <- function(input, group_id, category) {
  if (category == "clinical") {
    return(list(
      sex = input[[paste0("sex_filter_", group_id)]],
      race = input[[paste0("race_filter_", group_id)]],
      age = input[[paste0("age_", group_id)]],
      
      stage = input[[paste0("stage_filter_", group_id)]],
      risk = input[[paste0("risk_filter_", group_id)]],
      cyto_risk = input[[paste0("cyto_risk_filter_", group_id)]],
      
      rna_subtype = input[[paste0("rna_subtype_filter_", group_id)]],
      cna_subtype = input[[paste0("cna_subtype_filter_", group_id)]],
      
      triplet = input[[paste0("triplet_filter_", group_id)]],
      asct = input[[paste0("asct_filter_", group_id)]]
    ))
  } else if (category == "molecular") {
    return(list(
      q21_gain = input[[paste0("chr_1q21_gain_filter_", group_id)]],
      q21_amp = input[[paste0("chr_1q21_amp_filter_", group_id)]],
      del13q14 = input[[paste0("chr_13q14_del_filter_", group_id)]],
      del13q34 = input[[paste0("chr_13q34_del_filter_", group_id)]],
      del17p13 = input[[paste0("chr_17p13_del_filter_", group_id)]],
      
      diploidy = input[[paste0("diploidy_filter_", group_id)]],
      chromothripsis = input[[paste0("chromothripsis_filter_", group_id)]],
      
      t11_14 = input[[paste0("t_11_14_filter_", group_id)]],
      t4_14 = input[[paste0("t_4_14_filter_", group_id)]],
      
      maf = input[[paste0("maf_filter_", group_id)]],
      apobec = input[[paste0("apobec_filter_", group_id)]],
      tp53 = input[[paste0("tp53_filter_", group_id)]],
      tp53_ns = input[[paste0("tp53_ns_filter_", group_id)]]
    ))
  }
}

filter_group_data <- function(data, filters) {
  data <- as_tibble(data) # Ensures compatibility with dplyr functions
  
  data <- data %>%
    # Demographic filters
    {
      if (!is.null(filters$sex) && length(filters$sex) > 0)
        filter(., Sex %in% filters$sex) else .
    } %>%
    {
      if (!is.null(filters$race) && length(filters$race) > 0)
        filter(., Race %in% filters$race) else .
    } %>%
    {
      if (!is.null(filters$age) && length(filters$age) > 0) {
        min_age <- min(.$Age, na.rm = TRUE)
        max_age <- max(.$Age, na.rm = TRUE)
        if (filters$age[1] != min_age | filters$age[2] != max_age)
          filter(., Age >= filters$age[1], Age <= filters$age[2]) else .
      } else .
    } %>%
    
    # Clinical classification
    {
      if (!is.null(filters$stage) && length(filters$stage) > 0)
        filter(., ISS %in% filters$stage) else .
    } %>%
    {
      if (!is.null(filters$risk) && length(filters$risk) > 0)
        filter(., IMWG_Risk_Class %in% filters$risk) else .
    } %>%
    {
      if (!is.null(filters$cyto_risk) && length(filters$cyto_risk) > 0)
        filter(., Skerget_Cytogenetic_High_Risk %in% filters$cyto_risk) else .
    } %>%
    
    # Subtypes
    {
      if (!is.null(filters$rna_subtype) && length(filters$rna_subtype) > 0)
        filter(., Skerget_RNA_Subtype_Name %in% filters$rna_subtype) else .
    } %>%
    {
      if (!is.null(filters$cna_subtype) && length(filters$cna_subtype) > 0)
        filter(., Skerget_CNA_Subtype_Name %in% filters$cna_subtype) else .
    } %>%
    
    # Treatment
    {
      if (!is.null(filters$triplet) && length(filters$triplet) > 0)
        filter(., Triplet_First %in% filters$triplet) else .
    } %>%
    {
      if (!is.null(filters$asct) && length(filters$asct) > 0)
        filter(., ASCT_First %in% filters$asct) else .
    } %>%
    
    # Chromosomal abnormalities
    {
      if (!is.null(filters$q21_gain) && length(filters$q21_gain) > 0)
        filter(., chr_1q21_gain %in% filters$q21_gain) else .
    } %>%
    {
      if (!is.null(filters$q21_amp) && length(filters$q21_amp) > 0)
        filter(., chr_1q21_amp %in% filters$q21_amp) else .
    } %>%
    {
      if (!is.null(filters$del13q14) && length(filters$del13q14) > 0)
        filter(., chr_13q14_del %in% filters$del13q14) else .
    } %>%
    {
      if (!is.null(filters$del13q34) && length(filters$del13q34) > 0)
        filter(., chr_13q34_del %in% filters$del13q34) else .
    } %>%
    {
      if (!is.null(filters$del17p13) && length(filters$del17p13) > 0)
        filter(., chr_17p13_del %in% filters$del17p13) else .
    } %>%
    {
      if (!is.null(filters$diploidy) && length(filters$diploidy) > 0)
        filter(., Hyperdiploidy %in% filters$diploidy) else .
    } %>%
    {
      if (!is.null(filters$chromothripsis) && length(filters$chromothripsis) > 0)
        filter(., chromothripsis %in% filters$chromothripsis) else .
    } %>%
    
    # Translocations
    {
      if (!is.null(filters$t11_14) && length(filters$t11_14) > 0)
        filter(., t_11_14 %in% filters$t11_14) else .
    } %>%
    {
      if (!is.null(filters$t4_14) && length(filters$t4_14) > 0)
        filter(., t_4_14 %in% filters$t4_14) else .
    } %>%
    
    # Mutational markers
    {
      if (!is.null(filters$maf) && length(filters$maf) > 0)
        filter(., MAF_MAFB %in% filters$maf) else .
    } %>%
    {
      if (!is.null(filters$apobec) && length(filters$apobec) > 0)
        filter(., APOBEC %in% filters$apobec) else .
    } %>%
    {
      if (!is.null(filters$tp53) && length(filters$tp53) > 0)
        filter(., TP53_Funct_Copies %in% filters$tp53) else .
    } %>%
    {
      if (!is.null(filters$tp53_ns) && length(filters$tp53_ns) > 0)
        filter(., TP53_NS_Mut_Count %in% filters$tp53_ns) else .
    }
  
  return(data)
}

# Determine selected group
get_group_selected <- function(filtered_data, group_selected) {
  if (group_selected == "group1") {
    return(filtered_data$group1)
  } else {
    return(filtered_data$group2)
  }
}

# Summary of clinical data
generate_summary_plot <- function(group_selected, filtered_data) {
  selected_data <- filtered_data[[group_selected]]
  features <- c("Race", "Sex", "Age_range", "ISS", "ASCT_First")
  
  unique_values_counts_list <- lapply(features, function(feature) {
    data.frame(Value = names(table(selected_data[[feature]])),
               Count = as.integer(table(selected_data[[feature]])),
               Feature = feature)
  })
  
  unique_values_counts_df <- do.call(rbind, unique_values_counts_list)
  
  ggplot(unique_values_counts_df, aes(x = Value, y = Count, fill = Feature)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_wrap(~ Feature, scales = "free") +
    theme_minimal() +
    labs(title = "", x = "", y = "Count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(),
          legend.position = "none") +
    geom_text(aes(label = Count), vjust = 0.5, size = 3)
}

create_categorical_plot <- function(data, feature, feature_label) {
  plot_data <- data %>%
    group_by(group, !!sym(feature)) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(group) %>%
    mutate(
      percentage = n / sum(n) * 100,
      label = paste0("Group: ", group,
                     "<br>", feature_label, ": ", !!sym(feature),
                     "<br>Percentage: ", sprintf("%.2f", percentage), "%")
    )
  
  ggplot(plot_data, aes_string(x = feature, y = "percentage", fill = "group", alpha = "group", text = "label")) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = feature_label, y = "Percentage") +
    scale_alpha_manual(values = c(Group1 = 0.9, Group2 = 0.9)) +
    scale_fill_manual(values = c(Group1 = "#E87D72", Group2 = "#5BAEB0")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Create distribution plot data for continuous variables
#' @param data The combined dataset
#' @param feature The feature name to analyze
#' @param feature_label The display label for the feature
#' @return ggplot object
create_continuous_plot <- function(data, feature, feature_label) {
  data <- data %>%
    mutate(label = paste0("Group: ", group,
                          "<br>", feature_label, ": ", sprintf("%.2f", !!sym(feature))))
  
  ggplot(data, aes_string(x = "group", y = feature, fill = "group", text = "label")) +
    geom_boxplot(outlier.alpha = 0.4) +
    labs(x = "Group", y = feature_label) +
    scale_fill_manual(values = c(Group1 = "#E87D72", Group2 = "#5BAEB0")) +
    theme_minimal()
}

#' Main function to create distribution plot
#' @param data The combined dataset
#' @param feature The feature name to analyze
#' @param feature_label The display label for the feature
#' @param continuous_features Vector of continuous feature names
#' @return plotly object
create_distribution_plot <- function(data, feature, feature_label, continuous_features) {
  is_continuous <- feature %in% continuous_features
  
  if (!is_continuous) {
    p <- create_categorical_plot(data, feature, feature_label)
  } else {
    p <- create_continuous_plot(data, feature, feature_label)
  }
  
  ggplotly(p, tooltip = "text") %>%
    layout(hovermode = "x")
}

# =============================================================================
# STATISTICAL TESTING FUNCTIONS
# =============================================================================

#' Perform statistical test for continuous variables
#' @param group1_vals Vector of values for group 1
#' @param group2_vals Vector of values for group 2
#' @return List with test results
test_continuous_variable <- function(group1_vals, group2_vals) {
  # Remove NA values
  group1_vals <- group1_vals[!is.na(group1_vals)]
  group2_vals <- group2_vals[!is.na(group2_vals)]
  
  if (length(group1_vals) <= 2 || length(group2_vals) <= 2) {
    return(list(
      test_name = "Insufficient data",
      statistic = "NA",
      p_value = NA,
      group1_summary = "Insufficient data",
      group2_summary = "Insufficient data"
    ))
  }
  
  # Check normality
  normal1 <- tryCatch(shapiro.test(group1_vals)$p.value > 0.05, error = function(e) FALSE)
  normal2 <- tryCatch(shapiro.test(group2_vals)$p.value > 0.05, error = function(e) FALSE)
  
  if (normal1 && normal2 && length(group1_vals) >= 3 && length(group2_vals) >= 3) {
    # Use t-test
    test_result <- tryCatch(t.test(group1_vals, group2_vals), error = function(e) NULL)
    test_name <- "t-test"
    if (!is.null(test_result)) {
      statistic <- paste0("t = ", round(test_result$statistic, 3))
      p_value <- test_result$p.value
    } else {
      statistic <- "NA"
      p_value <- NA
    }
  } else {
    # Use Mann-Whitney U test
    test_result <- tryCatch(wilcox.test(group1_vals, group2_vals), error = function(e) NULL)
    test_name <- "Mann-Whitney U"
    if (!is.null(test_result)) {
      statistic <- paste0("W = ", round(test_result$statistic, 3))
      p_value <- test_result$p.value
    } else {
      statistic <- "NA"
      p_value <- NA
    }
  }
  
  # Summary statistics
  group1_summary <- paste0("Mean: ", round(mean(group1_vals, na.rm = TRUE), 2),
                           " (SD: ", round(sd(group1_vals, na.rm = TRUE), 2), ")")
  group2_summary <- paste0("Mean: ", round(mean(group2_vals, na.rm = TRUE), 2),
                           " (SD: ", round(sd(group2_vals, na.rm = TRUE), 2), ")")
  
  return(list(
    test_name = test_name,
    statistic = statistic,
    p_value = p_value,
    group1_summary = group1_summary,
    group2_summary = group2_summary
  ))
}

#' Perform statistical test for categorical variables
#' @param data The combined dataset
#' @param feature The feature name to analyze
#' @return List with test results
test_categorical_variable <- function(data, feature) {
  contingency_table <- tryCatch(table(data$group, data[[feature]], useNA = "no"),
                                error = function(e) NULL)
  
  if (is.null(contingency_table) || sum(contingency_table) == 0) {
    return(list(
      test_name = "Insufficient data",
      statistic = "NA",
      p_value = NA,
      group1_summary = "Insufficient data",
      group2_summary = "Insufficient data"
    ))
  }
  
  # Check for very low frequencies that could cause issues
  min_cell_count <- min(contingency_table)
  total_cells <- length(contingency_table)
  cells_less_than_5 <- sum(contingency_table < 5)
  
  # Skip analysis if too many cells have low counts or if any cell has 0 counts
  if (min_cell_count == 0 || cells_less_than_5 > (total_cells/2)) {
    return(list(
      test_name = "Low frequency categories",
      statistic = "NA",
      p_value = NA,
      group1_summary = "Low frequency",
      group2_summary = "Low frequency"
    ))
  }
  
  # Check if any expected frequencies are < 5
  expected_freq <- tryCatch(chisq.test(contingency_table)$expected, error = function(e) NULL)
  use_fisher <- !is.null(expected_freq) && any(expected_freq < 5)
  
  if (use_fisher) {
    # For tables larger than 2x2, Fisher's exact test can be computationally intensive
    if (nrow(contingency_table) == 2 && ncol(contingency_table) == 2) {
      test_result <- tryCatch(fisher.test(contingency_table), error = function(e) NULL)
      test_name <- "Fisher's exact"
      if (!is.null(test_result)) {
        statistic <- paste0("OR = ", round(test_result$estimate, 3))
        p_value <- test_result$p.value
      } else {
        statistic <- "NA"
        p_value <- NA
      }
    } else {
      # Use simulation for larger tables
      test_result <- tryCatch(chisq.test(contingency_table, simulate.p.value = TRUE, B = 2000), 
                              error = function(e) NULL)
      test_name <- "Chi-square (simulated)"
      if (!is.null(test_result)) {
        statistic <- paste0("χ² = ", round(test_result$statistic, 3))
        p_value <- test_result$p.value
      } else {
        statistic <- "NA"
        p_value <- NA
      }
    }
  } else {
    test_result <- tryCatch(chisq.test(contingency_table), error = function(e) NULL)
    test_name <- "Chi-square"
    if (!is.null(test_result)) {
      statistic <- paste0("χ² = ", round(test_result$statistic, 3))
      p_value <- test_result$p.value
    } else {
      statistic <- "NA"
      p_value <- NA
    }
  }
  
  # Summary - proportions
  if ("Group1" %in% rownames(contingency_table) && "Group2" %in% rownames(contingency_table)) {
    group1_counts <- contingency_table["Group1", ]
    group2_counts <- contingency_table["Group2", ]
    group1_total <- sum(group1_counts)
    group2_total <- sum(group2_counts)
    
    group1_summary <- paste(names(group1_counts), 
                            paste0(group1_counts, " (", round(group1_counts/group1_total*100, 1), "%)"),
                            sep = ": ", collapse = "; ")
    group2_summary <- paste(names(group2_counts),
                            paste0(group2_counts, " (", round(group2_counts/group2_total*100, 1), "%)"),
                            sep = ": ", collapse = "; ")
  } else {
    group1_summary <- "Data unavailable"
    group2_summary <- "Data unavailable"
  }
  
  return(list(
    test_name = test_name,
    statistic = statistic,
    p_value = p_value,
    group1_summary = group1_summary,
    group2_summary = group2_summary
  ))
}

#' Determine significance level from p-value
#' @param p_value The p-value
#' @return Significance indicator string
get_significance_level <- function(p_value) {
  if (is.na(p_value)) {
    return("NA")
  } else if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("NS")
  }
}

#' Format p-value for display
#' @param p_value The p-value
#' @return Formatted p-value string
format_p_value <- function(p_value) {
  if (is.na(p_value)) {
    return("NA")
  } else if (p_value < 0.001) {
    return("<0.001")
  } else {
    return(sprintf("%.4f", p_value))
  }
}

#' Create complete significance table
#' @param data The combined dataset
#' @param clinical_features Vector of clinical feature names to test
#' @param continuous_features Vector of continuous feature names
#' @return Data frame with significance results
create_significance_table <- function(data, clinical_features, continuous_features) {
  results <- data.frame(
    Feature = character(),
    Type = character(),
    Test_Used = character(),
    P_Value = numeric(),
    Statistic = character(),
    Group1_Summary = character(),
    Group2_Summary = character(),
    Significance = character(),
    stringsAsFactors = FALSE
  )
  
  for (feature in clinical_features) {
    # Skip if feature has too many missing values
    if (sum(!is.na(data[[feature]])) < 10) next
    
    # Determine if continuous or categorical
    is_continuous <- feature %in% continuous_features
    
    if (is_continuous) {
      group1_vals <- data[data$group == "Group1", feature]
      group2_vals <- data[data$group == "Group2", feature]
      test_results <- test_continuous_variable(group1_vals, group2_vals)
    } else {
      test_results <- test_categorical_variable(data, feature)
    }
    
    # Determine significance level
    significance <- get_significance_level(test_results$p_value)
    
    # Add to results
    results <- rbind(results, data.frame(
      Feature = feature,
      Type = ifelse(is_continuous, "Continuous", "Categorical"),
      Test_Used = test_results$test_name,
      P_Value = test_results$p_value,
      Statistic = test_results$statistic,
      Group1_Summary = test_results$group1_summary,
      Group2_Summary = test_results$group2_summary,
      Significance = significance,
      stringsAsFactors = FALSE
    ))
  }
  
  # Format p-values for display
  results$P_Value_Display <- sapply(results$P_Value, format_p_value)
  
  # Create final display table (without group summaries)
  display_table <- results[, c("Feature", "Type", "Test_Used", "P_Value_Display", 
                               "Statistic", "Significance")]
  colnames(display_table) <- c("Clinical Feature", "Type", "Statistical Test", "P-Value", 
                               "Test Statistic", "Significance")
  
  # Sort by p-value (significant first)
  display_table <- display_table[order(results$P_Value, na.last = TRUE), ]
  
  return(display_table)
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

#' Get clinical feature choices for selectInput
#' @param clinical_data The clinical dataset
#' @param exclude_cols Vector of column names to exclude
#' @return Vector of feature choices
get_clinical_feature_choices <- function(clinical_data, exclude_cols = NULL) {
  default_exclude <- c("public_id", "Tumor_Sample_Barcode", "Age", "Tx",
                       "PFS", "PFS_event", "PFS_censored", "OS", "OS_censored", 
                       "OS_event", "PFS_1", "PFS_1_censored", "PFS_1_event")
  
  if (!is.null(exclude_cols)) {
    exclude_cols <- c(default_exclude, exclude_cols)
  } else {
    exclude_cols <- default_exclude
  }
  
  setdiff(colnames(clinical_data), exclude_cols)
}

#' Define continuous features (customize this for your dataset)
#' @return Vector of continuous feature names
get_continuous_features <- function() {
  c("BMI", "Serum_B2M", "Serum_LDH", "Creatinine")
}


# Function to process DESeq2 -----
process_deseq2 <- function(filtered_data, bulkseq, min_counts=5, min_samples=5) {
  clinical_combined <- filtered_data$combined
  clinical_combined <- clinical_combined[clinical_combined$Tumor_Sample_Barcode %in% colnames(bulkseq),]
  count_data <- bulkseq[, clinical_combined$Tumor_Sample_Barcode]
  count_data <- round(count_data)
  
  # Create metadata
  metadata <- data.frame(
    row.names = clinical_combined$Tumor_Sample_Barcode,
    condition = as.factor(clinical_combined$group)
  )
  
  # Filter genes
  keep_genes <- rowSums(count_data >= min_counts) >= min_samples
  count_data <- count_data[keep_genes, ]
  
  # DESeq2 analysis
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~ condition)
  dds <- dds[rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)
  results <- results(dds, contrast = c("condition", "Group1", "Group2"))
  
  results <- as.data.frame(results)
  results$significant <- ifelse(results$padj < 0.05 & results$log2FoldChange > 1.5, "Up-regulated",
                                ifelse(results$padj < 0.05 & results$log2FoldChange < -1.5, "Down-regulated", "Not Significant"))
  results <- results[!is.na(results$significant), ]
  
  return(results)
}

# Function for volcano plot
deseq2_volcano <- function(deseq2_result, p_thr, logfc_thr) {
  deseq2_result$significant <- ifelse(deseq2_result$padj < p_thr & deseq2_result$log2FoldChange > logfc_thr, "Up-regulated",
                                      ifelse(deseq2_result$padj < p_thr & deseq2_result$log2FoldChange < -logfc_thr, "Down-regulated", "Not Significant"))
  
  top_genes <- deseq2_result %>%
    dplyr::arrange(padj) %>%
    dplyr::filter(significant != "Not Significant") %>%
    dplyr::slice(1:10)
  
  volcano_plot <- ggplot(deseq2_result, aes(x = log2FoldChange, y = -log10(padj), color = significant, shape = significant, fill = significant)) +
    geom_point(alpha = 0.8, size = 2, stroke = 0.5) +
    scale_color_manual(values = c("Not Significant" = "grey", "Up-regulated" = "red", "Down-regulated" = "blue")) +
    scale_fill_manual(values = c("Not Significant" = "grey", "Up-regulated" = "red", "Down-regulated" = "blue")) +
    scale_shape_manual(values = c("Not Significant" = 16, "Up-regulated" = 17, "Down-regulated" = 25)) +
    xlab("Log2 Fold Change") +
    ylab("-Log10 Adjusted P-value") +
    ggtitle("") +
    theme_minimal() +
    geom_text_repel(data = top_genes, aes(label = rownames(top_genes)), 
                    box.padding = 0.3, point.padding = 0.5, segment.color = 'grey50', 
                    show.legend = FALSE)
  
  return(volcano_plot)
}

# Function to create distribution density plot for gene of interest
tpm_distr_dens <- function(count_data_tpm, clinical_combined, gene_interested) {
  # Merge the datasets
  merged_data <- count_data_tpm %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "Tumor_Sample_Barcode") %>%
    inner_join(clinical_combined, by = "Tumor_Sample_Barcode")
  
  # Ensure 'group' is a factor
  merged_data$group <- as.factor(merged_data$group)
  
  # Filter for the gene of interest
  gene_data <- merged_data %>%
    dplyr::select(Tumor_Sample_Barcode, group, all_of(gene_interested)) %>%
    rename(TPM = all_of(gene_interested))
  
  # Calculate median TPM for each group
  median_tpm <- gene_data %>%
    group_by(group) %>%
    summarise(median_TPM = median(TPM, na.rm = TRUE))
  
  # Calculate the number of bins using the Freedman-Diaconis rule
  bin_width <- 2 * IQR(gene_data$TPM, na.rm = TRUE) / (length(gene_data$TPM)^(1/3))
  num_bins <- max(1, round((max(gene_data$TPM, na.rm = TRUE) - min(gene_data$TPM, na.rm = TRUE)) / bin_width))
  if (bin_width == 0) {
    num_bins <- 30
    bin_width <- max(1, round((max(gene_data$TPM, na.rm = TRUE) - min(gene_data$TPM, na.rm = TRUE)) / num_bins))
  }
  
  # Create a temp p for stats
  p_temp <- ggplot(gene_data, aes(x = TPM)) +
    geom_histogram(bins = num_bins) +
    geom_density(adjust = 1.5)
  
  # Extract data
  p_build <- ggplot_build(p_temp)
  histogram_data <- p_build$data[[1]]
  density_data <- p_build$data[[2]]
  
  scale_factor <- max(histogram_data$count) / max(density_data$density)
  p <- ggplot(gene_data, aes(x = TPM, fill = group, color = group)) +
    geom_histogram(aes(y = after_stat(count),
                       text = paste("TPM range:", round(after_stat(x) - after_stat(width)/2, 2), "-", round(after_stat(x) + after_stat(width)/2, 2), "<br>Count:", after_stat(count))),
                   position = "identity", bins = num_bins, fill = NA, alpha = 0) +
    geom_density(aes(y = after_stat(density) * scale_factor), alpha = 0.3, adjust = 1.5) +
    geom_vline(data = median_tpm, aes(xintercept = median_TPM, color = group, text = paste("Median TPM:", round(median_TPM, 2))), linetype = "dashed", linewidth = 1) +
    labs(x = "TPM",
         y = "Count",
         fill = "Group",
         color = "Group") +
    theme_minimal()
  
  return(p)
}

# Function to draw boxplot of tpm values and use Wilcoxon test for p values
tpm_boxplot <- function(count_data_tpm, clinical_combined, gene_interested) {
  # Merge the datasets
  merged_data <- count_data_tpm %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "Tumor_Sample_Barcode") %>%
    inner_join(clinical_combined, by = "Tumor_Sample_Barcode")
  
  # Ensure 'group' is a factor
  merged_data$group <- as.factor(merged_data$group)
  
  # Filter for the gene of interest
  gene_data <- merged_data %>%
    dplyr::select(Tumor_Sample_Barcode, group, all_of(gene_interested)) %>%
    rename(TPM = all_of(gene_interested))
  
  # Create the boxplot
  p <- ggboxplot(gene_data, x = "group", y = "TPM",
                 color = "group", add = "jitter") +
    stat_compare_means(method = "wilcox.test", label = "p.format", label.y = max(gene_data$TPM) * 1.1) +
    labs(x = "Group",
         y = paste("TPM of", gene_interested))
  
  return(p)
}

# Function to create distribution table for gene of interest
tpm_distr_table <- function(count_data_tpm, gene_interested) {
  gene_tpm <- as.numeric(count_data_tpm[gene_interested, ])
  
  # TPM quantile distribution for gene of interest
  quartiles <- quantile(gene_tpm, probs = c(0, 0.25, 0.5, 0.75, 1))
  quartile_table <- data.frame(
    Quantile = c("0-25%", "25-50%", "50-75%", "75-100%"),
    Expression_Range = c(
      paste0(round(quartiles[1], 2), "-", round(quartiles[2], 2)),
      paste0(round(quartiles[2], 2), "-", round(quartiles[3], 2)),
      paste0(round(quartiles[3], 2), "-", round(quartiles[4], 2)),
      paste0(round(quartiles[4], 2), "-", round(quartiles[5], 2))
    ),
    Number_of_Samples = c(
      sum(gene_tpm >= quartiles[1] & gene_tpm < quartiles[2]),
      sum(gene_tpm >= quartiles[2] & gene_tpm < quartiles[3]),
      sum(gene_tpm >= quartiles[3] & gene_tpm < quartiles[4]),
      sum(gene_tpm >= quartiles[4] & gene_tpm <= quartiles[5])
    )
  )
  return(quartile_table)
}

tpm_distr_survival <- function(gene_tpm, selected_clinical, grouping_method) {
  gene_interested <- rownames(gene_tpm)
  all_samples <- colnames(gene_tpm)
  gene_tpm <- as.numeric(gene_tpm)
  
  if (grouping_method == "quartiles") {
    # Group by quartiles
    quartiles <- quantile(gene_tpm, probs = c(0, 0.25, 0.5, 0.75, 1))
    
    # Prevent 'breaks' not unique error
    epsilon <- .Machine$double.eps
    quartiles <- quartiles + cumsum(duplicated(quartiles)) * epsilon
    
    # Create group
    group <- cut(gene_tpm,
                 breaks = quartiles,
                 labels = c("Q1", "Q2", "Q3", "Q4"),
                 include.lowest = TRUE)
    
  } else if (grouping_method == "median") {
    # Group by median
    median_value <- median(gene_tpm, na.rm = TRUE)
    group <- ifelse(gene_tpm <= median_value, "Below Median", "Above Median")
  }
  
  # Map group to clinical data
  group_df <- data.frame(Tumor_Sample_Barcode = all_samples,
                         group = group)
  
  selected_clinical <- merge(selected_clinical, group_df, by = "Tumor_Sample_Barcode")
  
  # Some samples does not have PFS and PFS_event, so the num of samples used in survival curve and selected_clinical is not consistent
  # Num of NAs: PFS, PFS_event; 308, 12. 
  # Survival curves
  surv_object <- Surv(time = selected_clinical$PFS_censored, event = selected_clinical$PFS_event)
  fit <- do.call(survfit, list(surv_object ~ group, data = selected_clinical))
  ggsurvplot(fit, data = selected_clinical, pval = TRUE,
             risk.table = TRUE, risk.table.col = "strata",
             ggtheme = theme_minimal(),
             tables.theme = theme_void(),
             xlab = "Days",
             ylab = "Progression-Free Survival",
             title = paste("Survival Curve by", gene_interested, "Expression by", if (grouping_method == "quartiles") "Quartiles" else "median"))
}

compute_significant_gene_sets <- function(ssgsea_result, clinical_combined) {
  # Intersect
  inter_samples <- intersect(colnames(ssgsea_result), clinical_combined$Tumor_Sample_Barcode)
  clinical_combined <- clinical_combined[clinical_combined$Tumor_Sample_Barcode %in% inter_samples, c("Tumor_Sample_Barcode", "group")]
  
  req(length(unique(clinical_combined$group)) == 2)
  
  ssgsea_result <- ssgsea_result[,inter_samples]
  ssgsea_result <- t(scale(t(ssgsea_result)))
  
  # Transpose ssGSEA result
  ssgsea_result_t <- as.data.frame(t(ssgsea_result))
  ssgsea_result_t$Tumor_Sample_Barcode <- rownames(ssgsea_result_t)
  
  # Merge
  merged_data <- merge(clinical_combined, ssgsea_result_t, by = "Tumor_Sample_Barcode")
  
  # Transform to long data
  long_data <- pivot_longer(
    merged_data,
    cols = -c(Tumor_Sample_Barcode, group),
    names_to = "GeneSet",
    values_to = "EnrichmentScore"
  )
  colnames(long_data)[1:2] <- c("Sample", "Group")
  
  
  unique_gene_sets <- unique(long_data$GeneSet)
  wilcox_results <- data.frame(GeneSet = character(), p_value = numeric(), stringsAsFactors = FALSE)
  
  for (gene_set in unique_gene_sets) {
    subset_data <- filter(long_data, GeneSet == gene_set)
    wilcox_result <- wilcox.test(EnrichmentScore ~ Group, data = subset_data)
    wilcox_results <- rbind(wilcox_results, data.frame(GeneSet = gene_set, p_value = wilcox_result$p.value))
  }
  
  # Adjust p-values using Benjamini-Hochberg method
  wilcox_results$adjusted_p_value <- p.adjust(wilcox_results$p_value, method = "BH")
  
  # Add significance stars based on p-value
  wilcox_results <- wilcox_results %>%
    mutate(Significance = case_when(
      adjusted_p_value < 0.001 ~ "***",
      adjusted_p_value < 0.01 ~ "**",
      adjusted_p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ))
  
  # Add direction
  wilcox_results$Direction <- ''
  for (gene_set in unique(wilcox_results$GeneSet)) {
    
    group1_mean <- mean(as.numeric(unlist(long_data[long_data$GeneSet == gene_set & long_data$Group == 'Group1', 'EnrichmentScore'])), na.rm = TRUE)
    group2_mean <- mean(as.numeric(unlist(long_data[long_data$GeneSet == gene_set & long_data$Group == 'Group2', 'EnrichmentScore'])), na.rm = TRUE)
    
    if (group1_mean > group2_mean) {
      direction <- 'Up'
    } else {
      direction <- 'Down'
    }
    
    wilcox_results[wilcox_results$GeneSet == gene_set, 'Direction'] <- direction
  }
  
  return(list(wilcox_results = wilcox_results, long_data = long_data))
}

create_violin_plot <- function(wilcox_results, long_data) {
  # Sort significant gene sets by p_value and select top 10
  top_significant_gene_sets <- wilcox_results %>% 
    filter(Significance != "ns") %>%
    arrange(adjusted_p_value) %>%
    head(10)
  
  # Filter data for top significant gene sets
  top_significant_data <- filter(long_data, GeneSet %in% top_significant_gene_sets$GeneSet)
  
  # Merge significance stars with top_significant_data
  top_significant_data <- merge(top_significant_data, top_significant_gene_sets, by = "GeneSet")
  
  # Violin plot
  dodge <- position_dodge(width=0.6)
  p <- ggplot(top_significant_data, aes(x = GeneSet, y = EnrichmentScore, fill = Group)) +
    geom_violin(trim = FALSE, width = 0.7, position = dodge) +
    geom_boxplot(width = 0.1, position = dodge) +
    stat_summary(fun = median, geom = "point", size = 2, color = "black", position = dodge) +
    geom_text(aes(x = GeneSet, y = max(EnrichmentScore) + 0.1, label = Significance), size = 5, vjust = 0) +
    theme_minimal() +
    labs(title = "Violin Plots for Top Significant Gene Sets",
         x = "Gene Set",
         y = "Enrichment Score") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # scale_x_discrete(labels = function(x) strwrap(x, width = 10)) + # What's causing the `$<-.data.frame`(`*tmp*`,... error
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  
  return(p)
}

celltype_boxplot <- function(group_info, sc_meta) {
  # Calculate cell type abundance for each patient
  sc_meta <- sc_meta[, c("public_id", "celltypes", "barcode")]
  patient_abundance <- sc_meta %>%
    group_by(public_id, celltypes) %>%
    summarise(count = n()) %>%
    mutate(proportion = count / sum(count))
  
  # Merge patient abundance with clinical data
  merged_data <- merge(patient_abundance, group_info, by = "public_id")
  
  req(length(unique(merged_data$group)) == 2)
  # Compare cell type abundance between groups
  comparison_results <- merged_data %>%
    group_by(celltypes) %>%
    summarise(
      p_value = t.test(proportion ~ group)$p.value
    ) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"))
  
  comparison_results <- comparison_results %>%
    mutate(Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ))
  
  merged_data <- merge(merged_data, comparison_results, by = "celltypes")
  
  # Create the boxplot
  ggplot(merged_data, aes(x = celltypes, y = proportion, fill = group)) +
    geom_boxplot() +
    labs(title = "Cell Type Abundance by Group",
         x = "Cell Type",
         y = "Proportion",
         fill = "Group") +
    theme_minimal() +
    geom_text(aes(x = celltypes, y = max(proportion) + 0.1, label = Significance), size = 5, vjust = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

celltype_proportion <- function(group_info, sc_meta) {
  # Calculate cell type abundance by group
  sc_meta <- sc_meta[, c("public_id", "celltypes", "barcode")]
  merged_meta <- merge(sc_meta, group_info, by = "public_id")
  cell_proportions <- merged_meta %>%
    group_by(group, celltypes) %>%
    summarise(count = n()) %>%
    mutate(proportion = count / sum(count) * 100)
  
  colors <- c(
    "CD4+" = "#7b3294", 
    "CD8+" = "#c2a5cf", 
    "B_Cells" = "#a6dba0",
    "Monocytes" = "#008837",  
    "NK" = "#fdae61", 
    "PlasmaCells" = "#e66101", 
    "UNK" = "#b2182b"
  )
  # Set the order of the group factor levels
  cell_proportions$group <- factor(cell_proportions$group, levels = c("Group2", "Group1"))
  # Reverse the order of factors for correct stacking
  cell_proportions$celltypes <- factor(cell_proportions$celltypes, levels = rev(levels(factor(cell_proportions$celltypes))))
  ggplot(cell_proportions, aes(y = group, x = proportion, fill = celltypes)) +
    geom_bar(stat = "identity", color = "black") + # Add border to the bars
    geom_label(aes(label = sprintf("%s\n%.2f%%", celltypes, proportion)),
               position = position_stack(vjust = 0.5, reverse = TRUE), size = 3.5,
               fill = "white", # Background for the label
               label.padding = unit(0.15, "lines"),
               color = "black", # Set text color to black for clarity
               fontface = "bold") + 
    scale_fill_manual(values = colors) +
    labs(title = "Cell Proportions by Group", y = "Group", x = "Percentage (%)") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )
}

cell_cycle_hist <- function(group_info, sc_meta, celltypes) {
  if (!('All' %in% celltypes)) {
    sc_meta <- sc_meta[sc_meta$celltypes %in% celltypes, ]
  }
  
  if (length(celltypes) == 0){
    showNotification("No cell types selected. Please select at least one.", type = "error")
    return(NULL)
  }
  
  sc_meta <- merge(sc_meta, group_info, by = "public_id")
  
  phase_group_counts <- table(sc_meta$Phase, sc_meta$group)
  phase_group_df <- as.data.frame(phase_group_counts)
  colnames(phase_group_df) <- c("Phase", "Group", "Count")
  
  # Calculate the percentage of cells in each phase within each group
  phase_group_summary <- phase_group_df %>%
    group_by(Group) %>%
    mutate(Total = sum(Count),
           Percentage = (Count / Total) * 100)
  
  # Create the histogram plot by group
  ggplot(phase_group_summary, aes(x = Phase, y = Percentage, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(title = "Cell Cycle Phase Distribution by Group (Percentage)",
         x = "Cell Cycle Phase",
         y = "Percentage of Cells") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


create_distribution_stacked_barplot <- function(data, x_feature, y_features) {
  # Function to create a single stacked barplot for one y feature
  create_single_plot <- function(data, x_feature, y_feature) {
    # Ensure factors are properly ordered
    data[[x_feature]] <- factor(data[[x_feature]])
    
    # Handle NA values in a cleaner way
    data[[y_feature]] <- ifelse(is.na(data[[y_feature]]), "NA", as.character(data[[y_feature]]))
    data[[y_feature]] <- factor(data[[y_feature]])
    
    # Count occurrences of each combination
    prop_data <- data %>%
      group_by(!!sym(x_feature), !!sym(y_feature)) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(!!sym(x_feature)) %>%
      mutate(proportion = count / sum(count))
    
    # Create plot
    p <- ggplot(prop_data, aes(x = !!sym(x_feature), y = proportion, fill = !!sym(y_feature))) +
      geom_bar(stat = "identity", position = "stack", color = "black") +
      labs(title = paste(y_feature, "by", x_feature),
           x = x_feature,
           y = "Proportion",
           fill = y_feature) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    return(p)
  }
  
  # Create a list of plots
  plot_list <- lapply(y_features, function(y) {
    create_single_plot(data, x_feature, y)
  })
  
  # Determine grid layout based on number of plots
  n_plots <- length(plot_list)
  n_cols <- min(3, n_plots)  # Maximum 3 columns
  n_rows <- ceiling(n_plots / n_cols)
  
  # Arrange plots in a grid
  grid_plot <- do.call(grid.arrange, c(plot_list, ncol = n_cols))
  
  return(grid_plot)
}
