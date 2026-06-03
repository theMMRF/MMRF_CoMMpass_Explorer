# User interface ------
dashboardPage(
  # skin = "purple",
  # Header
  dashboardHeader(
    title = "CoMMpass Explorer",
    tags$li(class = "dropdown",
            actionButton("start_tour", "Interactive Tour", icon = icon("question-circle"))
    )
  ),

  # Sidebar
  dashboardSidebar(
    useShinyjs(),
    div(id = "filters_sidebar",
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
                                  selectInput("surv_variable_cohort1", "Variable",
                                              choices = c(
                                                "PFS (days)"               = "PFS",
                                                "OS (days)"                = "OS",
                                                "Time to Second Line (days)"  = "ttct2line"
                                              ),
                                              selected = "PFS"),
                                  checkboxInput("require_surv_event_cohort1", "Require event = 1", FALSE),

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
                                  selectInput("surv_variable_cohort2", "Variable",
                                              choices = c(
                                                "PFS (days)"               = "PFS",
                                                "OS (days)"                = "OS",
                                                "Time to Second Line (days)"  = "ttct2line"
                                              ),
                                              selected = "PFS"),
                                  checkboxInput("require_surv_event_cohort2", "Require event = 1", FALSE),

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
          div(id = "apply_filters_wrap",
              actionButton("apply_filters", "Apply Filters", icon = icon("play"), class = "btn-primary")
          )
        )
    )

  ),

  # Dashboard
  dashboardBody(
    rintrojs::introjsUI(),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
    ),
    tabsetPanel(
      id = "main_tabs",
      tabPanel("Upload Cohorts", value = "Custom Cohorts",
               fluidRow(
                 box(title = "How it works", width = 12, id = "how_it_works",
                     HTML("<ul style='margin-bottom:0'>
             <li>Upload a CSV/TXT <em>or</em> paste <code>public_id</code>s (e.g., <code>MMRF_1013_1</code>), one per line or comma-separated.</li>
             <li>Click <strong>Load</strong> for each cohort.</li>
             <li>Click <strong>Apply Filters</strong> on the left menu to enable your loaded cohorts.</li>
             </ul>")
                 )
               ),
               fluidRow(
                 column(6,
                        div(id = "cohort1_upload",
                            box(title = "Cohort 1", width = 12,
                                fileInput("upload_cohort1", "Upload Cohort 1 public_ids (csv/txt)", accept = c(".csv", ".txt", ".tsv")),
                                textAreaInput("paste_cohort1", "Or paste public_ids", rows = 5, placeholder = "MMRF_1007_1\nMMRF_1013_1\n..."),
                                actionButton("load_cohort1", "Load Cohort 1 IDs"),
                                downloadButton("download_unmatched_c1", "Download unmatched IDs"),
                                tags$hr(),
                                verbatimTextOutput("cohort1_status"),
                                DTOutput("cohort1_preview")
                            )
                        )
                 ),
                 column(6,
                        div(id = "cohort2_upload",
                            box(title = "Cohort 2", width = 12,
                                fileInput("upload_cohort2", "Upload Cohort 2 public_ids (csv/txt)", accept = c(".csv", ".txt", ".tsv")),
                                textAreaInput("paste_cohort2", "Or paste public_ids", rows = 5, placeholder = "MMRF_1020_1\nMMRF_1055_1\n..."),
                                actionButton("load_cohort2", "Load Cohort 2 IDs"),
                                downloadButton("download_unmatched_c2", "Download unmatched IDs"),
                                tags$hr(),
                                verbatimTextOutput("cohort2_status"),
                                DTOutput("cohort2_preview")
                            )
                        )
                 )
               )
      ),
      tabPanel("Overall Summary", value = "Overall Summary",
               fluidRow(
                 column(6, div(id = "counts_card", htmlOutput("clinicalNum"))),
                 column(6, align = "right", downloadButton("download_clinical", "Download Filtered Clinical Data"))
               ),

               fluidRow(
                 box(title = "Cohort 1", width = 6, id = "summary_g1",
                     plotOutput("summaryPlot_g1"),
                     plot_export_controls_ui("summaryPlot_g1")),
                 box(title = "Cohort 2", width = 6, id = "summary_g2",
                     plotOutput("summaryPlot_g2"),
                     plot_export_controls_ui("summaryPlot_g2"))

               ),

               fluidRow(
                 box(title = "Survival Curve (PFS)", width = 4,
                     plotOutput("survCompPlot_pfs"),
                     plot_export_controls_ui("survCompPlot_pfs")
                 ),

                 box(title = "Survival Curve (Time to Second Line)", width = 4,
                     plotOutput("survCompPlot_tt2Line_censored"),
                     plot_export_controls_ui("survCompPlot_tt2Line_censored")
                 ),

                 box(title = "Survival Curve (OS)", width = 4,
                     plotOutput("survCompPlot_os"),
                     plot_export_controls_ui("survCompPlot_os")
                 )
               ),

               fluidRow(
                 box(title = "Distribution by Cohort", width = 12,
                     selectInput("clin_feature",
                                 "Enter a clinical feature",
                                 choices = get_clinical_feature_choices(clinical_data),
                                 width = "200px"),
                     plotlyOutput("clin_distribution"),
                     plot_export_controls_ui("clin_distribution")
                 )
               ),
               fluidRow(
                 box(title = "Clinical Feature Significance Table", width = 12,
                     p("Statistical comparison of clinical features between Cohort 1 and Cohort 2."),
                     p(HTML("<strong>Significance levels:</strong> *** p<0.001, ** p<0.01, * p<0.05, NS = Not Significant")),
                     DTOutput("significance_table")
                 )
               ),

               ###################### COX-PH ######################
               fluidRow(
                 box(title = "Model setup", width = 12,
                     fluidRow(
                       column(3,
                              radioButtons(
                                "cox_data_source", "Data for model",
                                choices = c("Cohort 1", "Cohort 2", "Both"),
                                selected = "Both"
                              ),

                              selectInput(
                                "cox_endpoint", "Endpoint",
                                choices = c(
                                  "PFS"                    = "PFS",
                                  "OS"                     = "OS",
                                  "Time to Second Line"    = "TT2L"
                                ),
                                selected = "OS"
                              ),

                              checkboxInput("cox_use_cohort", "Include cohort indicator", TRUE),
                              selectInput("cox_strata", "Stratify by",
                                          choices = c("None"), selected = "None")
                       ),
                       column(5,
                              selectizeInput(
                                "cox_covars",
                                "Covariates (clinical / molecular)",
                                choices = NULL, multiple = TRUE,
                                options = list(placeholder = 'Pick variables to adjust for…',
                                               plugins = list("remove_button")),
                                width = "100%"
                              )
                       ),
                       column(4,
                              selectizeInput(
                                "cox_gene",
                                "Add bulk TPM gene covariate",
                                choices = NULL,
                                multiple = FALSE,
                                selected = "",
                                options = list(placeholder = 'None')
                              ),

                              radioButtons("cox_gene_mode", "Gene encoding",
                                           choices = c("continuous (log2(TPM+1))" = "continuous",
                                                       "median split (Low/High)" = "median"),
                                           selected = "continuous")
                       )
                     ),
                     actionButton("fit_cox", "Fit Cox model", icon = icon("play"), class = "btn-primary"),
                     tags$br(), tags$br(),
                     verbatimTextOutput("cox_formula")
                 )
               ),
               fluidRow(
                 box(title = "Cox-PH Results", width = 12,
                     DTOutput("cox_table"),
                     downloadButton("download_cox_table", "Download HR table")
                 )
               ),
               fluidRow(
                 box(title = "Forest plot", width = 12,
                     plotOutput("cox_forest", height = "600px"),
                     plot_export_controls_ui("cox_forest")
                 )
               )
      ),

      tabPanel("Mutational Profile", value = "Mutational Profile",
               fluidRow(
                 column(12, htmlOutput("mafNum"))
               ),

               # MAF summary
               fluidRow(
                 box(title = "Cohort 1", width = 6,
                     plotOutput("mafSummary_g1"),
                     plot_export_controls_ui("mafSummary_g1")
                 ),
                 box(title = "Cohort 2", width = 6,
                     plotOutput("mafSummary_g2"),
                     plot_export_controls_ui("mafSummary_g2")
                 )
               ),

               # Oncoplot
               fluidRow(
                 box(title = "", width = 6,
                     plotOutput("oncoplot_g1"),
                     plot_export_controls_ui("oncoplot_g1")
                 ),
                 box(title = "", width = 6,
                     plotOutput("oncoplot_g2"),
                     plot_export_controls_ui("oncoplot_g2")
                 )
               ),

               # Lollipop plot
               fluidRow(
                 box(title = "", width = 6,
                     selectizeInput("gene_search_lollipop_g1", "Enter Gene for Lollipop Plot", choices = NULL, selected = "KRAS", options = list(create = TRUE, placeholder = 'Search for genes'),  width = "200px"),
                     plotOutput("lollipopPlot_g1"),
                     plot_export_controls_ui("lollipopPlot_g1")
                 ),
                 box(title = "", width = 6,
                     selectizeInput("gene_search_lollipop_g2", "Enter Gene for Lollipop Plot", choices = NULL, selected = "KRAS", options = list(create = TRUE, placeholder = 'Search for genes'),  width = "200px"),
                     plotOutput("lollipopPlot_g2"),
                     plot_export_controls_ui("lollipopPlot_g2")
                 )
               ),

               # Interaction
               fluidRow(
                 box(title = "", width = 6,
                     selectizeInput("gene_search_inter_g1", "Enter 5 or more Genes for Interaction Plot", choices = NULL,
                                    multiple = TRUE, options = list(create = TRUE, placeholder = 'Search for genes')),
                     plotOutput("interactionPlot_g1"),
                     plot_export_controls_ui("interactionPlot_g1")
                 ),
                 box(title = "", width = 6,
                     selectizeInput("gene_search_inter_g2", "Enter 5 or more Genes for Interaction Plot", choices = NULL,
                                    multiple = TRUE, options = list(create = TRUE, placeholder = 'Search for genes')),
                     plotOutput("interactionPlot_g2"),
                     plot_export_controls_ui("interactionPlot_g2")
                 )
               ),

               # CoOncoplot and Barplot
               fluidRow(
                 box(title = "", width = 6,
                     selectizeInput("gene_search_maf", "Enter Genes for MAF Comparison", choices = NULL,
                                    multiple = TRUE, options = list(create = TRUE, placeholder = 'Search for genes')),
                     plotOutput("mafCompOncoPlot"),
                     plot_export_controls_ui("mafCompOncoPlot")
                 ),
                 box(title = "", width = 6,
                     plotOutput("mafCompBarPlot"),
                     plot_export_controls_ui("mafCompBarPlot")
                 )
               ),

               fluidRow(
                 box(title = "", width = 12,
                     dataTableOutput("mafCompTable")
                 )
               ),

               fluidRow(
                 box(title = "", width = 6,
                     plotOutput("mafCompForestPlot"),
                     plot_export_controls_ui("mafCompForestPlot")
                 )
               ),
      ),

      tabPanel("Tumor Profile", value = "Tumor Profile",
               fluidRow(
                 column(12, htmlOutput("bulkNum"))
               ),
               fluidRow(
                 box(title = "Distribution", width = 12,
                     selectizeInput("gene_search_bulk_distr", "Enter Gene", choices = NULL, selected = "KRAS", options = list(create = TRUE, placeholder = 'Search for genes'), width = "300px"),
                     plotlyOutput("tpm_distr"),
                     plot_export_controls_ui("tpm_distr")
                 )
               ),

               fluidRow(
                 box(title = "", width = 12,
                     plotOutput("tpm_distr_boxplot"),
                     plot_export_controls_ui("tpm_distr_boxplot")
                 )
               ),

               fluidRow(
                 box(title = "Cohort 1", width = 6,
                     # p("This plot visualizes the survival probability of patient samples of the selected cohort based on the expression levels of a gene selected in TPM Distribution box."),
                     selectInput("cohorting_method_tpm_g1", "Cohorting Method", choices = c("quartiles", "median"), width = "150px"),
                     plotOutput("tpm_survCompPlot_g1"),
                     plot_export_controls_ui("tpm_survCompPlot_g1")
                 ),
                 box(title = "Cohort 2", width = 6,
                     selectInput("cohorting_method_tpm_g2", "Cohorting Method", choices = c("quartiles", "median"), width = "150px"),
                     plotOutput("tpm_survCompPlot_g2"),
                     plot_export_controls_ui("tpm_survCompPlot_g2")
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

               box(title = "DESeq2", width = 12,
                   p("The DESeq2 analysis will perform differential expression ..."),
                   fluidRow(
                     column(6, numericInput("p_threshold", "P-value Threshold", value = 0.05, min = 0, max = 1, step = 0.01)),
                     column(6, numericInput("fc_threshold", "Log2 Fold Change Threshold", value = 1.5, min = 0, max = 10, step = 0.1))
                   ),
                   verbatimTextOutput("deseq_status"),
                   actionButton("start_deseq2", "Start DESeq2 Analysis"),
                   plotlyOutput("bulkVolcano", height = "520px"),
                   plot_export_controls_ui("bulkVolcano")
               ),

               fluidRow(
                 box(title = "DEGs", width = 12,
                     DTOutput("DEGs_table"),
                     downloadButton("download_DEGs", "Download Full Table")
                 ),
               ),

               fluidRow(
                 box(title = "ssGSEA", width = 12,
                     selectizeInput("selected_gene_sets", "Select Gene Sets to Display:",
                                    choices = NULL, multiple = TRUE),
                     plotOutput("ssgsea_violin"),
                     plot_export_controls_ui("ssgsea_violin")
                 ),
                 box(title = "Gene set table", width = 12,
                     DTOutput("ssgsea_table")
                 )
               )
      ),

      tabPanel("Immune Microenvironment", value = "Immune Microenvironment",
               tabsetPanel(
                 id = "immune_tabs",

                 # ---- scRNAseq TAB ----
                 tabPanel("scRNAseq", value = "scRNAseq",
                          fluidRow(
                            column(12, htmlOutput("scNum"))
                          ),
                          fluidRow(
                            box(title = "Cell Type Abundance", width = 12,
                                plotOutput("sc_celltype_boxplot"),
                                plot_export_controls_ui("sc_celltype_boxplot")
                            )
                          ),
                          fluidRow(
                            box(title = "Cell Type Proportion", width = 12,
                                plotOutput("sc_celltype_proportion"),
                                plot_export_controls_ui("sc_celltype_proportion")
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
                                plotOutput("sc_cellcycle_hist"),
                                plot_export_controls_ui("sc_cellcycle_hist")
                            )
                          )
                 ),

                 # ---- PSEUDO-BULK TAB ----
                 tabPanel("Pseudo-Bulk", value = "Pseudo-Bulk",
                          fluidRow(
                            column(12, htmlOutput("pseudoNum"))
                          ),

                          fluidRow(
                            box(title = "Distribution", width = 12,
                                selectInput("pseudo_celltype", "Cell type", choices = names(pseudo_bulk_norm), width = "240px"),
                                selectizeInput("gene_search_pseudo_distr", "Enter Gene", choices = NULL, selected = "KRAS",
                                               options = list(create = TRUE, placeholder = 'Search for genes'), width = "300px"),
                                plotlyOutput("pseudo_norm_distr"),
                                plot_export_controls_ui("pseudo_norm_distr")
                            )
                          ),

                          fluidRow(
                            box(title = NULL, width = 12,
                                plotOutput("pseudo_norm_distr_boxplot"),
                                plot_export_controls_ui("pseudo_norm_distr_boxplot")
                            )
                          ),

                          fluidRow(
                            box(title = "Cohort 1", width = 6,
                                selectInput("cohorting_method_pseudo_g1", "Cohorting Method", choices = c("quartiles", "median"), width = "150px"),
                                plotOutput("pseudo_norm_survCompPlot_g1"),
                                plot_export_controls_ui("pseudo_norm_survCompPlot_g1")
                            ),
                            box(title = "Cohort 2", width = 6,
                                selectInput("cohorting_method_pseudo_g2", "Cohorting Method", choices = c("quartiles", "median"), width = "150px"),
                                plotOutput("pseudo_norm_survCompPlot_g2"),
                                plot_export_controls_ui("pseudo_norm_survCompPlot_g2")
                            )
                          ),

                          fluidRow(
                            box(title = NULL, width = 6, DTOutput("pseudo_quantile_table_cohort1")),
                            box(title = NULL, width = 6, DTOutput("pseudo_quantile_table_cohort2"))
                          ),

                          fluidRow(
                            box(title = "Differential analysis (pseudobulk)", width = 12,
                                p("Using DESeq2."),
                                fluidRow(
                                  column(6, numericInput("p_threshold_pseudo", "P-value Threshold", value = 0.05, min = 0, max = 1, step = 0.01)),
                                  column(6, numericInput("fc_threshold_pseudo", "Log2 Fold Change Threshold", value = 1.5, min = 0, max = 10, step = 0.1))
                                ),
                                verbatimTextOutput("pseudo_deseq_status"),
                                actionButton("start_pseudo_diff", "Run differential analysis"),
                                plotlyOutput("pseudoVolcano", height = "520px"),
                                plot_export_controls_ui("pseudoVolcano")
                            )
                          ),

                          fluidRow(
                            box(title = "DEGs", width = 12, DTOutput("pseudo_DEGs_table")),
                            downloadButton("download_pseudo_DEGs", "Download Full Table")
                          )
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
