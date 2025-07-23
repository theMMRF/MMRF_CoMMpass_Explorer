source("global.R")
# User interface ------
dashboardPage(
  # skin = "purple",
  # Header
  dashboardHeader(title = "CoMMpass Explorer"),
  
  # Sidebar
  dashboardSidebar(
    useShinyjs(),
    sidebarMenu(
      id = "sidebarMenu",
      
      # Filters
      menuItem("Cohort 1 Filters", icon = icon("filter"),
               div(id = "cohort1_filters",
                   menuItem("Clinical Features", icon = icon("heartbeat"),
                            create_cohort_filters_ui("cohort1", "clinical")),
                   menuItem("Molecular Features", icon = icon("dna"),
                            create_cohort_filters_ui("cohort1", "molecular")),
                   menuItem("Gene Mutations", icon = icon("vial"),
                            actionButton("add_mut_row_cohort1", "Add Mutation Rule"),
                            actionButton("remove_mut_row_cohort1", "Remove Last Rule"),
                            uiOutput("mutation_rules_cohort1")
                   ),
                   
                   menuItem("Gene Expression", icon = icon("chart-line"),
                            selectizeInput("gene_expr_search_cohort1", "Gene", choices = NULL),
                            radioButtons("expr_threshold_type_cohort1", "Threshold Type", choices = c("value", "percentile"), selected = "percentile"),
                            conditionalPanel(
                              condition = "input.expr_threshold_type_cohort1 == 'percentile'",
                              fluidRow(
                                column(6, numericInput("gene_expr_percentile_min_cohort1", "Min Percentile", value = 0, min = 0, max = 100)),
                                column(6, numericInput("gene_expr_percentile_max_cohort1", "Max Percentile", value = 100, min = 0, max = 100))
                              )
                            ),
                            conditionalPanel(
                              condition = "input.expr_threshold_type_cohort1 == 'value'",
                              fluidRow(
                                column(6, numericInput("gene_expr_min_cohort1", "Min Expression", value = NULL)),
                                column(6, numericInput("gene_expr_max_cohort1", "Max Expression", value = NULL))
                              )
                            )
                   ),
                   
                   menuItem("Survival Filter", icon = icon("hourglass-half"),
                            checkboxInput("enable_survival_filter_cohort1", "Enable", FALSE),
                            conditionalPanel(
                              condition = "input.enable_survival_filter_cohort1 == true",
                              radioButtons("surv_threshold_type_cohort1", "Type", choices = c("value", "percentile"), selected = "value"),
                              selectInput("surv_variable_cohort1", "Variable", choices = c("PFS_censored", "OS_censored"), selected = "OS_censored"),
                              
                              conditionalPanel(
                                condition = "input.surv_threshold_type_cohort1 == 'value'",
                                fluidRow(
                                  column(6, div(style = "width: 100%;", numericInput("surv_threshold_min_value_cohort1", "Min Days", value = NULL))),
                                  column(6, div(style = "width: 100%;", numericInput("surv_threshold_max_value_cohort1", "Max Days", value = NULL)))
                                )
                              ),
                              conditionalPanel(
                                condition = "input.surv_threshold_type_cohort1 == 'percentile'",
                                fluidRow(
                                  column(6, div(style = "width: 100%;", numericInput("surv_threshold_min_percentile_cohort1", "Min Percentile", value = 0, min = 0, max = 100))),
                                  column(6, div(style = "width: 100%;", numericInput("surv_threshold_max_percentile_cohort1", "Max Percentile", value = 100, min = 0, max = 100)))
                                )
                              )
                            )
                   )
               ),
               actionButton("clear_cohort1", "Clear All Cohort 1 Filters", icon = icon("times")) # Clear button for Cohort 1
      ),
      
      menuItem("Cohort 2 Filters", icon = icon("filter"),
               div(id = "cohort2_filters",
                   menuItem("Clinical Features", icon = icon("heartbeat"),
                            create_cohort_filters_ui("cohort2", "clinical")),
                   menuItem("Molecular Features", icon = icon("dna"),
                            create_cohort_filters_ui("cohort2", "molecular")),
                   menuItem("Gene Mutations", icon = icon("vial"),
                            actionButton("add_mut_row_cohort2", "Add Mutation Rule"),
                            actionButton("remove_mut_row_cohort2", "Remove Last Rule"),
                            uiOutput("mutation_rules_cohort2")
                   ),
                   
                   menuItem("Gene Expression", icon = icon("chart-line"),
                            selectizeInput("gene_expr_search_cohort2", "Gene", choices = NULL),
                            radioButtons("expr_threshold_type_cohort2", "Threshold Type", choices = c("value", "percentile"), selected = "percentile"),
                            conditionalPanel(
                              condition = "input.expr_threshold_type_cohort2 == 'percentile'",
                              fluidRow(
                                column(6, numericInput("gene_expr_percentile_min_cohort2", "Min Percentile", value = 0, min = 0, max = 100)),
                                column(6, numericInput("gene_expr_percentile_max_cohort2", "Max Percentile", value = 100, min = 0, max = 100))
                              )
                            ),
                            conditionalPanel(
                              condition = "input.expr_threshold_type_cohort2 == 'value'",
                              fluidRow(
                                column(6, numericInput("gene_expr_min_cohort2", "Min Expression", value = NULL)),
                                column(6, numericInput("gene_expr_max_cohort2", "Max Expression", value = NULL))
                              )
                            )
                   ),
                   menuItem("Survival Filter", icon = icon("hourglass-half"),
                            checkboxInput("enable_survival_filter_cohort2", "Enable", FALSE),
                            conditionalPanel(
                              condition = "input.enable_survival_filter_cohort2 == true",
                              radioButtons("surv_threshold_type_cohort2", "Type", choices = c("value", "percentile"), selected = "value"),
                              selectInput("surv_variable_cohort2", "Variable", choices = c("PFS_censored", "OS_censored"), selected = "OS_censored"),
                              
                              conditionalPanel(
                                condition = "input.surv_threshold_type_cohort2 == 'value'",
                                fluidRow(
                                  column(6, div(style = "width: 100%;", numericInput("surv_threshold_min_value_cohort2", "Min Days", value = NULL))),
                                  column(6, div(style = "width: 100%;", numericInput("surv_threshold_max_value_cohort2", "Max Days", value = NULL)))
                                )
                              ),
                              conditionalPanel(
                                condition = "input.surv_threshold_type_cohort2 == 'percentile'",
                                fluidRow(
                                  column(6, div(style = "width: 100%;", numericInput("surv_threshold_min_percentile_cohort2", "Min Percentile", value = 0, min = 0, max = 100))),
                                  column(6, div(style = "width: 100%;", numericInput("surv_threshold_max_percentile_cohort2", "Max Percentile", value = 100, min = 0, max = 100)))
                                )
                              )
                            )
                   )
               ),
               actionButton("clear_cohort2", "Clear All Cohort 2 Filters", icon = icon("times"))
      ),
      
      actionButton("apply_filters", "Apply Filters", icon = icon("play"), class = "btn-primary")
    )
  ),
  
  # Dashboard
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
    ),
    tabsetPanel(
      tabPanel("Overall Summary",
               fluidRow(
                 column(6, htmlOutput("clinicalNum")),
                 column(6, align = "right", downloadButton("download_clinical", "Download Filtered Clinical Data"))
               ),
               
               fluidRow(
                 box(title = "Cohort 1", width = 6,
                     plotOutput("summaryPlot_g1")
                 ),
                 
                 box(title = "Cohort 2", width = 6,
                     plotOutput("summaryPlot_g2")
                 )
               ),
               
               fluidRow(
                 box(title = "Cohort1 vs Cohort2 Survival Curve (PFS)", width = 6,
                     plotOutput("survCompPlot_pfs_censored")
                 ),
                 
                 box(title = "Cohort1 vs Cohort2 Survival Curve (OS)", width = 6,
                     plotOutput("survCompPlot_os_censored")
                 )
               ),
               
               fluidRow(
                 box(title = "Distribution by Cohort", width = 12,
                     selectInput("clin_feature", 
                                 "Enter a clinical feature", 
                                 choices = get_clinical_feature_choices(clinical_data), 
                                 width = "200px"),
                     plotlyOutput("clin_distribution")
                 )
               ),
               fluidRow(
                 box(title = "Clinical Feature Significance Table", width = 12,
                     p("Statistical comparison of clinical features between Cohort 1 and Cohort 2."),
                     p(HTML("<strong>Significance levels:</strong> *** p<0.001, ** p<0.01, * p<0.05, NS = Not Significant")),
                     DTOutput("significance_table")
                 )
               )
               
      ),
      
      tabPanel("Mutational Profile",
               fluidRow(
                 column(12, htmlOutput("mafNum"))
               ),
               
               # MAF summary
               fluidRow(
                 box(title = "Cohort 1", width = 6,
                     plotOutput("mafSummary_g1")
                 ),
                 box(title = "Cohort 2", width = 6,
                     plotOutput("mafSummary_g2")
                 )
               ),
               
               # Oncoplot
               fluidRow(
                 box(title = "", width = 6,
                     plotOutput("oncoplot_g1")
                 ),
                 box(title = "", width = 6,
                     plotOutput("oncoplot_g2")
                 )
               ),
               
               # Lollipop plot
               fluidRow(
                 box(title = "", width = 6,
                     selectizeInput("gene_search_lollipop_g1", "Enter Gene for Lollipop Plot", choices = NULL, selected = "KRAS", options = list(create = TRUE, placeholder = 'Search for genes'),  width = "200px"),
                     plotOutput("lollipopPlot_g1")
                 ),
                 box(title = "", width = 6,
                     selectizeInput("gene_search_lollipop_g2", "Enter Gene for Lollipop Plot", choices = NULL, selected = "KRAS", options = list(create = TRUE, placeholder = 'Search for genes'),  width = "200px"),
                     plotOutput("lollipopPlot_g2")
                 ),
                 # box(title = "Variant Allele Frequency", width = 6,
                 #     plotOutput("vafPlot")
                 # )
               ),
               
               # Interaction
               fluidRow(
                 box(title = "", width = 6,
                     selectizeInput("gene_search_inter_g1", "Enter 5 or more Genes for Interaction Plot", choices = NULL,
                                    multiple = TRUE, options = list(create = TRUE, placeholder = 'Search for genes')),
                     plotOutput("interactionPlot_g1")
                 ),
                 box(title = "", width = 6,
                     selectizeInput("gene_search_inter_g2", "Enter 5 or more Genes for Interaction Plot", choices = NULL,
                                    multiple = TRUE, options = list(create = TRUE, placeholder = 'Search for genes')),
                     plotOutput("interactionPlot_g2")
                 )
               ),
               
               # CoOncoplot and Barplot
               fluidRow(
                 box(title = "", width = 6,
                     selectizeInput("gene_search_maf", "Enter Genes for MAF Comparison", choices = NULL,
                                    multiple = TRUE, options = list(create = TRUE, placeholder = 'Search for genes')),
                     plotOutput("mafCompOncoPlot")
                 ),
                 box(title = "", width = 6,
                     plotOutput("mafCompBarPlot")
                 )
               ),
               
               fluidRow(
                 box(title = "", width = 6,
                     plotOutput("mafCompForestPlot")
                 )
               ),
      ),
      
      tabPanel("Tumor Profile",
               fluidRow(
                 column(12, htmlOutput("bulkNum"))
               ),
               fluidRow(
                 box(title = "TPM Distribution", width = 12,
                     selectizeInput("gene_search_bulk_distr", "Enter Gene", choices = NULL, selected = "KRAS", options = list(create = TRUE, placeholder = 'Search for genes'), width = "300px"),
                     plotlyOutput("tpm_distr")
                 )
               ),
               
               fluidRow(
                 box(title = "", width = 12,
                     plotOutput("tpm_distr_boxplot")
                 )
               ),
               
               fluidRow(
                 box(title = "Cohort 1", width = 6,
                     # p("This plot visualizes the survival probability of patient samples of the selected cohort based on the expression levels of a gene selected in TPM Distribution box."),
                     selectInput("cohorting_method_tpm_g1", "Cohorting Method", choices = c("quartiles", "median"), width = "150px"),
                     plotOutput("tpm_survCompPlot_g1")
                 ),
                 box(title = "Cohort 2", width = 6,
                     selectInput("cohorting_method_tpm_g2", "Cohorting Method", choices = c("quartiles", "median"), width = "150px"),
                     plotOutput("tpm_survCompPlot_g2")
                 )
               ),
               
               fluidRow(
                 box(title = "", width = 6,
                     DTOutput("quantile_table_cohort1"),
                 ),
                 
                 box(title = "", width = 6,
                     DTOutput("quantile_table_cohort2"),
                 )
               ),
               
               fluidRow(
                 box(title = "DESeq2", width = 12,
                     p("The DESeq2 analysis will perform differential expression analysis on Cohort 1 and Cohort 2,"),
                     p("comparing their gene expression profiles to identify significant differences. Click 'Start DESeq2 Analysis' to proceed."),
                     p("Each run would take several minutes depending on the sample size."),
                     fluidRow(
                       column(6, numericInput("p_threshold", "P-value Threshold", value = 0.05, min = 0, max = 1, step = 0.01)),
                       column(6, numericInput("fc_threshold", "Log2 Fold Change Threshold", value = 1.5, min = 0, max = 10, step = 0.1))
                     ),
                     verbatimTextOutput("deseq_status"),
                     actionButton("start_deseq2", "Start DESeq2 Analysis"),
                     plotOutput("bulkVolcano")
                 ),
                 # box(title = "Heatmap", width = 6,
                 #     selectizeInput("gene_search_bulk_heat", "Enter Genes",choices = NULL,
                 #                    multiple = TRUE, options = list(create = TRUE, placeholder = 'Search for genes')),
                 #     plotOutput("bulkHeat"),
                 #     textOutput("portionExpressed")
                 # ),
               ),
               
               fluidRow(
                 box(title = "DEGs", width = 12,
                     DTOutput("DEGs_table"),
                 ),
                 downloadButton("download_DEGs", "Download Full Table")
               ),
               
               # fluidRow(
               #   box(title = "ssGSEA", width = 12,
               #       plotOutput("ssgsea_violin"),
               #   ),
               #   box(title = "Gene set table", width = 12,
               #       DTOutput("ssgsea_table")
               #   )
               # )
               
               fluidRow(
                 box(title = "ssGSEA", width = 12,
                     selectizeInput("selected_gene_sets", "Select Gene Sets to Display:",
                                    choices = NULL, multiple = TRUE),
                     plotOutput("ssgsea_violin")
                 ),
                 box(title = "Gene set table", width = 12,
                     DTOutput("ssgsea_table")
                 )
               )
               
               
               
      ),
      
      tabPanel("Immune Microenvironment",
               fluidRow(
                 column(12, htmlOutput("scNum"))
               ),
               fluidRow(
                 box(title = "Cell Type Abundance", width = 12,
                     plotOutput("sc_celltype_boxplot")
                 )
               ),
               
               fluidRow(
                 box(title = "Cell Type Proportion", width = 12,
                     plotOutput("sc_celltype_proportion")
                 )
               ),
               
               fluidRow(
                 box(title = "Cell Cycle Distribution", width = 12,
                     selectizeInput(
                       "celltypes_interested", 
                       "Enter Cell Types", 
                       choices = c("All", unique(sc_meta$celltypes)), 
                       selected = "All",
                       multiple = TRUE, 
                       options = list(create = TRUE, placeholder = 'Search for cell types')
                     ),
                     plotOutput("sc_cellcycle_hist")
                 )
               )
      )
    ),
    tags$div(
      style = "text-align: center; padding: 30px 10px 20px 10px; border-top: 1px solid #ccc; background-color: #f9f9f9;",
      
      tags$img(src = "MMRF_sign.png", height = "40px", style = "margin-bottom: 10px;"),
      
      tags$p(
        strong("Reference: "), 
        "Zhang W, Acharya C.R. (Chuck), ", 
        em("Multiple Myeloma Research Foundation (MMRF)"),
        strong(", "),
        " (2025). doi: ", 
        tags$a("[Link]", href = "", target = "_blank"),
        style = "font-size: 90%; color: #771544;"
      ),
      
      tags$p(
        HTML("&copy; 2025 Multiple Myeloma Research Foundation"),
        style = "font-size: 85%; color: #555;"
      )
    )
  )
)