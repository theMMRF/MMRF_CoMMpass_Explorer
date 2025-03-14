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
clinical_data <- readRDS("data/clinical_data_cleaned.rds")
maf_data <- readRDS("data/maf_data.rds")
sc_meta <- readRDS("data/scRNAseq_metadata.rds")
ssgsea_result_ca <- readRDS("data/ssgsea_result_ca.rds")
ssgsea_result_c2 <- readRDS("data/ssgsea_result_c2.rds")
# bulkseq <- readRDS("../../data/commpass_explorer/bulkseq_baseline.rds")
# bulkseq_tpm <- readRDS("../../data/commpass_explorer/bulkseq_tpm_baseline.rds")
# clinical_data <- readRDS("../../data/commpass_explorer/clinical_data_cleaned.rds")
# maf_data <- readRDS("../../data/commpass_explorer/maf_data.rds")
# sc_meta <- readRDS("../../data/commpass_explorer/scRNAseq_metadata.rds")
# ssgsea_result_ca <- readRDS("../../data/commpass_explorer/ssgsea_result_ca.rds")
# ssgsea_result_c2 <- readRDS("../../data/commpass_explorer/ssgsea_result_c2.rds")

# Helper functions -----------
# NA processing
fill_na_with_same_patient <- function(data) {
  unique_patients <- unique(sapply(strsplit(as.character(data$Tumor_Sample_Barcode), "_"),
                                   function(x) paste(x[1], x[2], sep="_")))
  
  for (patient in unique_patients) {
    # print(paste("Processing patient:", patient))
    patient_rows <- grep(patient, data$Tumor_Sample_Barcode)
    patient_data <- data[patient_rows, ]
    
    if (nrow(patient_data) > 1) {
      for (col in names(data)) {
        if (col != "Tumor_Sample_Barcode") {
          
          # Use first non-NA value of a column to fill in this column of other samples
          all_non_na <- unlist(patient_data[!is.na(patient_data[[col]]), ..col])
          
          # If > 1, then there are multiple values for this clinical parameter of the same patient
          if (length(unique(all_non_na)) > 1) {
            stop("Error: conflicting clinical parameters for samples of the same patient.")
          }
          
          non_na_value <- all_non_na[1]
          
          if (!is.na(non_na_value)) {
            data[patient_rows, (col) := fifelse(is.na(get(..col)), non_na_value, get(..col))]
          }
        }
      }
    }
  }
  
  return(data)
}

remove_all_na_rows <- function(data) {
  cleaned_data <- data[!(rowSums(is.na(data)) == ncol(data)-1), ]
  return(cleaned_data)
}

subset_by_gene_mutations <- function(clinical_data, include_genes = NULL, group, logic) {
  if (is.null(include_genes)) return(clinical_data)
  maf_data_table <- maf_data@data
  
  ## "Or" logic
  if (logic == "Or") {
    # Filter for included genes if specified
    if (group == "group1") {
      included_patients <- unique(maf_data_table[Hugo_Symbol %in% include_genes, .(Tumor_Sample_Barcode)])
      clinical_data <- clinical_data[clinical_data$Tumor_Sample_Barcode %in% included_patients$Tumor_Sample_Barcode]
    }
    
    # Filter out excluded genes if specified
    if (group == "group2") {
      excluded_patients <- unique(maf_data_table[Hugo_Symbol %in% include_genes, .(Tumor_Sample_Barcode)])
      clinical_data <- clinical_data[!clinical_data$Tumor_Sample_Barcode %in% excluded_patients$Tumor_Sample_Barcode]
    }
  }
  
  ## "And" logic
  # Filter for included genes if specified
  if (logic == "And") {
    if (group == "group1") {
      sample_gene_count <- maf_data_table %>%
        filter(Hugo_Symbol %in% include_genes) %>%
        distinct(Tumor_Sample_Barcode, Hugo_Symbol) %>%
        group_by(Tumor_Sample_Barcode) %>%
        summarise(gene_count = n(), .groups = "drop")
      
      samples_with_all_genes <- sample_gene_count %>%
        filter(gene_count == length(include_genes)) %>%
        pull(Tumor_Sample_Barcode)
      
      clinical_data <- clinical_data[clinical_data$Tumor_Sample_Barcode %in% samples_with_all_genes]
    }
    
    # Filter out excluded genes if specified
    if (group == "group2") {
      sample_gene_count <- maf_data_table %>%
        filter(Hugo_Symbol %in% include_genes) %>%
        distinct(Tumor_Sample_Barcode, Hugo_Symbol) %>%
        group_by(Tumor_Sample_Barcode) %>%
        summarise(gene_count = n(), .groups = "drop")
      
      samples_with_all_genes <- sample_gene_count %>%
        filter(gene_count == length(include_genes)) %>%
        pull(Tumor_Sample_Barcode)
      
      clinical_data <- clinical_data[!clinical_data$Tumor_Sample_Barcode %in% samples_with_all_genes]
    }
  }
  return(clinical_data)
}

filter_by_gene_expression <- function(clinical_data, gene=NULL, threshold=NULL, group) {
  if (is.null(gene) || is.na(threshold)) return(clinical_data)
  print(gene)
  print(threshold)
  if (gene %in% unique(rownames(bulkseq_tpm)) && (group == "group1")) {
    patient_ids <- colnames(bulkseq_tpm[, bulkseq_tpm[gene, ] >= threshold, drop = FALSE])
    print("patient_ids")
    print(length(patient_ids))
    print(patient_ids)
    clinical_data <- clinical_data[clinical_data$Tumor_Sample_Barcode %in% patient_ids, ]
  }
  
  if (gene %in% unique(rownames(bulkseq_tpm)) && (group == "group2")) {
    patient_ids <- colnames(bulkseq_tpm[, bulkseq_tpm[gene, ] < threshold, drop = FALSE])
    print("patient_ids")
    print(length(patient_ids))
    clinical_data <- clinical_data[clinical_data$Tumor_Sample_Barcode %in% patient_ids, ]
  }
  
  # if clinical_data
  
  return(clinical_data)
}


create_picker_input <- function(inputId, label, choices) {
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
            create_picker_input(paste0("risk_filter_", group_id), "Risk Group (Davies-based)", sort(unique(clinical_data$Risk))),
            uiOutput(paste0("age_filter_", group_id)),
            create_picker_input(paste0("race_filter_", group_id), "Race", sort(unique(clinical_data$Race))),
            create_picker_input(paste0("ethnicity_filter_", group_id), "Ethnicity", sort(unique(clinical_data$Ethnicity))),
            create_picker_input(paste0("sex_filter_", group_id), "Sex", sort(unique(clinical_data$Sex))),
            create_picker_input(paste0("stage_filter_", group_id), "Stage (ISS)", sort(unique(clinical_data$ISS))),
            create_picker_input(paste0("asct_filter_", group_id), "ASCT Firstline", sort(unique(clinical_data$ASCT))),
            create_picker_input(paste0("triplet_filter_", group_id), "Triplet Firstline", sort(unique(clinical_data$TRIP_FirstLine)))
          ))
        } else if (category == "molecular") {
          return(list(
            create_picker_input(paste0("diploidy_filter_", group_id), "Hyperdiploidy", sort(unique(clinical_data$Hyperdiploidy))),
            create_picker_input(paste0("chromothripsis_filter_", group_id), "Chromothripsis", sort(unique(clinical_data$chromothripsis))),
            create_picker_input(paste0("t(11;14)_filter_", group_id), "t(11;14)", sort(unique(clinical_data$`t(11;14)`))),
            create_picker_input(paste0("t(4;14)_filter_", group_id), "t(4;14)", sort(unique(clinical_data$`t(4;14)`))),
            create_picker_input(paste0("1q21_amp_filter_", group_id), "1q21_amp", sort(unique(clinical_data$`1q21_amp`))),
            create_picker_input(paste0("1q21_gain_filter_", group_id), "1q21_gain", sort(unique(clinical_data$`1q21_gain`))),
            create_picker_input(paste0("13q14_del_filter_", group_id), "13q14_del", sort(unique(clinical_data$`13q14_del`))),
            create_picker_input(paste0("13q34_del_filter_", group_id), "13q34_del", sort(unique(clinical_data$`13q34_del`))),
            create_picker_input(paste0("17p13_del_filter_", group_id), "17p13_del", sort(unique(clinical_data$`17p13_del`))),
            create_picker_input(paste0("apobec_filter_", group_id), "Apobec", sort(unique(clinical_data$APOBEC))),
            create_picker_input(paste0("maf_filter_", group_id), "MAF/MAFB", sort(unique(clinical_data$`MAF/MAFB`))),
            create_picker_input(paste0("tp53_filter_", group_id), "TP53 inactivation", sort(unique(clinical_data$`TP53 inactivation`)))
          ))
        } else if (category == "gene") {
          
        }
    )
  )
  
}

get_group_filters <- function(input, group_id, category) {
  if (category == "clinical") {
    return(list(
      risk = input[[paste0("risk_filter_", group_id)]],
      age = input[[paste0("age_", group_id)]],
      race = input[[paste0("race_filter_", group_id)]],
      ethnicity = input[[paste0("ethnicity_filter_", group_id)]],
      sex = input[[paste0("sex_filter_", group_id)]],
      stage = input[[paste0("stage_filter_", group_id)]],
      asct = input[[paste0("asct_filter_", group_id)]],
      triplet = input[[paste0("triplet_filter_", group_id)]]
    ))
  } else if (category == "molecular") {
    return(list(
      diploidy = input[[paste0("diploidy_filter_", group_id)]],
      chromothripsis = input[[paste0("chromothripsis_filter_", group_id)]],
      t11_14 = input[[paste0("t(11;14)_filter_", group_id)]],
      t4_14 = input[[paste0("t(4;14)_filter_", group_id)]],
      q21_amp = input[[paste0("1q21_amp_filter_", group_id)]],
      q21_gain = input[[paste0("1q21_gain_filter_", group_id)]],
      del13q14 = input[[paste0("13q14_del_filter_", group_id)]],
      del13q34 = input[[paste0("13q34_del_filter_", group_id)]],
      del17p13 = input[[paste0("17p13_del_filter_", group_id)]]
    ))
  } else if (category == "gene") {
    return(list(
      apobec = input[[paste0("apobec_filter_", group_id)]],
      maf = input[[paste0("maf_filter_", group_id)]],
      tp53 = input[[paste0("tp53_filter_", group_id)]]
    ))
  }
}

filter_group_data <- function(data, filters) {
  # Apply clinical filters
  if (!is.null(filters$risk) && length(filters$risk) > 0) {
    data <- data[Risk %in% filters$risk]
  }
  
  if (!is.null(filters$age) && length(filters$age) > 0) {
    min_age <- min(clinical_data$Age, na.rm = TRUE)
    max_age <- max(clinical_data$Age, na.rm = TRUE)
    if (filters$age[1] != min_age | filters$age[2] != max_age) {
      data <- data[Age >= filters$age[1] & Age <= filters$age[2]]
    }
  }
  
  if (!is.null(filters$race) && length(filters$race) > 0) {
    data <- data[Race %in% filters$race]
  }
  if (!is.null(filters$ethnicity) && length(filters$ethnicity) > 0) {
    data <- data[Ethnicity %in% filters$ethnicity]
  }
  if (!is.null(filters$sex) && length(filters$sex) > 0) {
    data <- data[Sex %in% filters$sex]
  }
  if (!is.null(filters$stage) && length(filters$stage) > 0) {
    data <- data[ISS %in% filters$stage]
  }
  if (!is.null(filters$asct) && length(filters$asct) > 0) {
    data <- data[ASCT %in% filters$asct]
  }
  if (!is.null(filters$triplet) && length(filters$triplet) > 0) {
    data <- data[TRIP_FirstLine %in% filters$triplet]
  }
  
  # Apply molecular filters
  if (!is.null(filters$diploidy) && length(filters$diploidy) > 0) {
    data <- data[Hyperdiploidy %in% filters$diploidy]
  }
  if (!is.null(filters$chromothripsis) && length(filters$chromothripsis) > 0) {
    data <- data[chromothripsis %in% filters$chromothripsis]
  }
  if (!is.null(filters$t11_14) && length(filters$t11_14) > 0) {
    data <- data[`t(11;14)` %in% filters$t11_14]
  }
  if (!is.null(filters$t4_14) && length(filters$t4_14) > 0) {
    data <- data[`t(4;14)` %in% filters$t4_14]
  }
  if (!is.null(filters$q21_amp) && length(filters$q21_amp) > 0) {
    data <- data[`1q21_amp` %in% filters$q21_amp]
  }
  if (!is.null(filters$q21_gain) && length(filters$q21_gain) > 0) {
    data <- data[`1q21_gain` %in% filters$q21_gain]
  }
  if (!is.null(filters$del13q14) && length(filters$del13q14) > 0) {
    data <- data[`13q14_del` %in% filters$del13q14]
  }
  if (!is.null(filters$del13q34) && length(filters$del13q34) > 0) {
    data <- data[`13q34_del` %in% filters$del13q34]
  }
  if (!is.null(filters$del17p13) && length(filters$del17p13) > 0) {
    data <- data[`17p13_del` %in% filters$del17p13]
  }
  
  # Apply gene mutation filters
  if (!is.null(filters$apobec) && length(filters$apobec) > 0) {
    data <- data[APOBEC %in% filters$apobec]
  }
  if (!is.null(filters$maf) && length(filters$maf) > 0) {
    data <- data[`MAF/MAFB` %in% filters$maf]
  }
  if (!is.null(filters$tp53) && length(filters$tp53) > 0) {
    data <- data[`TP53 inactivation` %in% filters$tp53]
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
  features <- c("Race", "Ethnicity", "Sex", "Age_range", "ISS", "ASCT")
  
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
  surv_object <- Surv(time = selected_clinical$PFS, event = selected_clinical$PFS_event)
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
  long_data <- melt(merged_data, id.vars = c("Tumor_Sample_Barcode", "group"))
  colnames(long_data) <- c("Sample", "Group", "GeneSet", "EnrichmentScore")
  
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
