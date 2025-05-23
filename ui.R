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
      menuItem("Group 1 Filters", icon = icon("filter"),
               div(id = "group1_filters",
                   menuItem("Clinical Features", icon = icon("heartbeat"),
                            create_group_filters_ui("group1", "clinical")),
                   menuItem("Molecular Features", icon = icon("dna"),
                            create_group_filters_ui("group1", "molecular"))
                   # menuItem("Gene Mutations", icon = icon("vial"),
                   #          create_group_filters_ui("group1", "gene"))
               ),
               actionButton("clear_group1", "Clear All Group 1 Filters", icon = icon("times")) # Clear button for Group 1
      ),
      
      menuItem("Group 2 Filters", icon = icon("filter"),
               div(id = "group2_filters",
                   menuItem("Clinical Features", icon = icon("heartbeat"),
                            create_group_filters_ui("group2", "clinical")),
                   menuItem("Molecular Features", icon = icon("dna"),
                            create_group_filters_ui("group2", "molecular"))
                   # menuItem("Gene Mutations", icon = icon("vial"),
                   #          create_group_filters_ui("group2", "gene"))
               ),
               actionButton("clear_group2", "Clear All Group 2 Filters", icon = icon("times"))
      ),
      
      menuItem("Gene Mutations", icon = icon("filter"),
               selectizeInput("mut_logic", "Filter Logic", choices = c("And", "Or"), selected = "And"),
               selectizeInput("gene_include_filter", "Input Genes", choices = NULL, multiple = TRUE, options = list(create = TRUE, placeholder = 'Search for genes'))
      ),
      
      menuItem("Gene Expression", icon = icon("filter"),
               selectizeInput("gene_expr_search", "Input a Gene", choices = NULL, multiple = FALSE, options = list(create = TRUE, placeholder = 'Search for a gene')),
               numericInput("gene_expr_threshold", "Expression Threshold", value = NULL, min = 0, step = 0.1)
      ),
      
      menuItem("Survival Filter", icon = icon("hourglass-half"),
               checkboxInput("enable_survival_filter", "Enable Survival Filtering", value = FALSE),
               
               conditionalPanel(
                 condition = "input.enable_survival_filter == true",
                 
                 selectInput("surv_variable", "Select Survival Variable",
                             choices = c("PFS_censored", "OS_censored"), 
                             selected = "OS_censored"),
                 
                 radioButtons("surv_threshold_type", "Threshold Type",
                              choices = c("Days" = "value", "Percentile" = "percentile"),
                              selected = "value", inline = TRUE),
                 
                 conditionalPanel(
                   condition = "input.surv_threshold_type == 'value'",
                   numericInput("surv_threshold_value", "Survival Threshold (Days)", value = 365, min = 0)
                 ),
                 
                 conditionalPanel(
                   condition = "input.surv_threshold_type == 'percentile'",
                   sliderInput("surv_threshold_percentile", "Percentile Threshold", min = 0, max = 100, value = 50)
                 )
               )
      )
      
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
                 column(12, htmlOutput("clinicalNum"))
               ),
               fluidRow(
                 box(title = "Group 1", width = 6,
                     plotOutput("summaryPlot_g1")
                 ),
                 
                 box(title = "Group 2", width = 6,
                     plotOutput("summaryPlot_g2")
                 )
               ),
               
               fluidRow(
                 box(title = "Group1 vs Group2 Survival Curve (PFS)", width = 6,
                     plotOutput("survCompPlot_pfs_censored")
                 ),
                 
                 box(title = "Group1 vs Group2 Survival Curve (OS)", width = 6,
                     plotOutput("survCompPlot_os_censored")
                 )
               ),
               
               fluidRow(
                 box(title = "Distribution by Group", width = 12,
                     selectInput("clin_feature", 
                                 "Enter a clinical feature", 
                                 choices = get_clinical_feature_choices(clinical_data), 
                                 width = "200px"),
                     plotlyOutput("clin_distribution")
                 )
               ),
               fluidRow(
                 box(title = "Clinical Feature Significance Table", width = 12,
                     p("Statistical comparison of clinical features between Group 1 and Group 2."),
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
                 box(title = "Group 1", width = 6,
                     plotOutput("mafSummary_g1")
                 ),
                 box(title = "Group 2", width = 6,
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
                 box(title = "TPM Distribution", width = 12,
                     htmlOutput("bulkNum"),
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
                 box(title = "Group 1", width = 6,
                     # p("This plot visualizes the survival probability of patient samples of the selected group based on the expression levels of a gene selected in TPM Distribution box."),
                     selectInput("grouping_method_tpm_g1", "Grouping Method", choices = c("quartiles", "median"), width = "150px"),
                     plotOutput("tpm_survCompPlot_g1")
                 ),
                 box(title = "Group 2", width = 6,
                     selectInput("grouping_method_tpm_g2", "Grouping Method", choices = c("quartiles", "median"), width = "150px"),
                     plotOutput("tpm_survCompPlot_g2")
                 )
               ),
               
               fluidRow(
                 box(title = "", width = 6,
                     DTOutput("quantile_table_group1"),
                 ),
                 
                 box(title = "", width = 6,
                     DTOutput("quantile_table_group2"),
                 )
               ),
               
               fluidRow(
                 box(title = "DESeq2", width = 12,
                     p("The DESeq2 analysis will perform differential expression analysis on Group 1 and Group 2,"),
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
                 box(title = "Cell Type Abundance", width = 12,
                     htmlOutput("scNum"),
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
      ),
      
      tabPanel("Feature Distribution Analysis",
               fluidRow(
                 box(title = "Select Features", width = 12,
                     column(6, 
                            selectInput("feature_x", "Primary Feature (X-axis)", 
                                        choices = setdiff(colnames(clinical_data), 
                                                          c("Tumor_Sample_Barcode", "Age", "PFS", "PFS_event", "PFS_censored", 
                                                            "OS", "OS_censored", "OS_event")),
                                        selected = "chr_1q21_gain")
                     ),
                     column(6, 
                            selectInput("feature_y", "Secondary Features to Display", 
                                        choices = setdiff(colnames(clinical_data), 
                                                          c("Tumor_Sample_Barcode", "Age", "PFS", "PFS_event", "PFS_censored", 
                                                            "OS", "OS_censored", "OS_event")),
                                        multiple = TRUE,
                                        selected = c("Sex", "Ethnicity", "ISS"))
                     )
                 )
               ),
               fluidRow(
                 box(title = "Distribution Plots", width = 12,
                     plotOutput("feature_distribution_plot", height = "800px")
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
          " (2024). doi: ", 
          tags$a("[Link]", href = "", target = "_blank"),
          style = "font-size: 90%; color: #771544;"
        ),
        
        tags$p(
          HTML("&copy; 2025 Multiple Myeloma Research Foundation"),
          style = "font-size: 85%; color: #555;"
        )
      )
    ),
  )
)