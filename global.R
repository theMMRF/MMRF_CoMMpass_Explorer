# Load dependencies -------
packages <- c("shiny", "shinydashboard", "shinyWidgets", "shinyjs", "rintrojs",
              "tidyr", "dplyr", "DT", "data.table", "tibble", "gridExtra", "readxl",
              "ggplot2", "pheatmap", "ggrepel", "plotly",
              "maftools", "survival", "survminer", "limma", "DESeq2",
              "clusterProfiler", "org.Hs.eg.db", "biomaRt")

lapply(packages, library, character.only = TRUE)

# Load data ------
bulkseq <- readRDS("data/bulkseq_baseline_cleaned.rds")
bulkseq_tpm <- readRDS("data/bulkseq_tpm_baseline_cleaned.rds")
clinical_data <- readRDS("data/clinical_data_n1141.rds")
maf_data <- readRDS("data/maf_data.rds")
sc_meta <- readRDS("data/sc_meta_manuscript_baseline.rds")
ssgsea_result_ca <- readRDS("data/ssgsea_result_ca.rds")
pseudo_bulk_counts <- readRDS("data/pseudobulk_data_manuscript_Counts.rds")
pseudo_bulk_norm <- readRDS("data/pseudobulk_data_manuscript_LogNormalize.rds")

clinical_data$PFS <- as.numeric(as.character(clinical_data$PFS))
clinical_data$PFS_event <- as.numeric(as.character(clinical_data$PFS_event))
clinical_data$OS <- as.numeric(as.character(clinical_data$OS))
clinical_data$OS_event <- as.numeric(as.character(clinical_data$OS_event))
clinical_data$ttct2line <- as.numeric(as.character(clinical_data$ttct2line))
clinical_data$censt2line <- as.numeric(as.character(clinical_data$censt2line))

clinical_data$CGS_risk <- factor(clinical_data$CGS_risk,
                                 levels = c("Standard", "High"))

# Safe lookup of counts from a table
.safe_count <- function(tbl, v) {
  n <- unname(tbl[as.character(v)])
  if (is.na(n)) 0L else as.integer(n)
}

.sanitize_percentile_range <- function(min_percentile, max_percentile) {
  min_percentile <- suppressWarnings(as.numeric(min_percentile))
  max_percentile <- suppressWarnings(as.numeric(max_percentile))
  if (!is.finite(min_percentile) || !is.finite(max_percentile)) return(NULL)

  min_percentile <- max(0, min(100, min_percentile))
  max_percentile <- max(0, min(100, max_percentile))
  if (min_percentile > max_percentile) {
    tmp <- min_percentile
    min_percentile <- max_percentile
    max_percentile <- tmp
  }

  list(min = min_percentile, max = max_percentile)
}

# ---- user's cohort helpers ----------------------------------------------------
.parse_public_ids_from_file <- function(file_input) {
  if (is.null(file_input)) return(character())
  path <- file_input$datapath
  if (!nzchar(path)) return(character())
  # Try fread for flexible csv/tsv; fallback to readLines for raw txt
  ids <- tryCatch({
    dt <- data.table::fread(path, header = TRUE, sep = ",", fill = TRUE, data.table = FALSE)
    # if single column without header, fread will name it V1
    cols <- names(dt)
    cand <- intersect(cols, c("public_id", "Public_ID", "PublicId", "id", "ID", "V1"))
    if (length(cand) == 0) cand <- cols[1]
    as.character(dt[[cand[1]]])
  }, error = function(e) {
    # raw text: one per line / comma separated
    txt <- readLines(path, warn = FALSE)
    unlist(strsplit(paste(txt, collapse = ","), "[\n,;\t ]+"))
  })
  ids <- trimws(ids)
  ids[nzchar(ids)]
}

.parse_public_ids_from_text <- function(txt) {
  if (is.null(txt) || !nzchar(txt)) return(character())
  ids <- unlist(strsplit(txt, "[\n,;\t ]+"))
  ids <- trimws(ids)
  ids[nzchar(ids)]
}

# Map public_ids to Tumor_Sample_Barcode; allow fallback by base id (drop trailing _#)
.map_public_ids_to_tsb <- function(public_ids, clinical) {
  if (length(public_ids) == 0) return(list(tsb = character(), matched_public = character(), unmatched = character()))
  public_ids <- unique(public_ids)
  baseify <- function(x) sub("(_\\d+)?$", "", x)

  exact_match <- clinical$Tumor_Sample_Barcode[match(public_ids, clinical$public_id)]
  matched_exact <- !is.na(exact_match)

  # Fallback: match by base id (MMRF_1234)
  need_fallback <- which(!matched_exact)
  tsb_fb <- character(length(public_ids)); tsb_fb[] <- NA_character_
  if (length(need_fallback)) {
    clin_base <- baseify(as.character(clinical$public_id))
    pid_base  <- baseify(public_ids[need_fallback])
    idx <- match(pid_base, clin_base)
    tsb_fb[need_fallback] <- clinical$Tumor_Sample_Barcode[idx]
  }

  tsb_all <- ifelse(matched_exact, exact_match, tsb_fb)
  matched <- !is.na(tsb_all)

  list(
    tsb = unique(tsb_all[matched]),
    matched_public = unique(public_ids[matched]),
    unmatched = unique(public_ids[!matched])
  )
}


# Helper functions -----------
.mutation_data_available_ids <- function() {
  ids <- character()
  if ("clinical.data" %in% slotNames(maf_data) &&
      "Tumor_Sample_Barcode" %in% names(maf_data@clinical.data)) {
    ids <- as.character(maf_data@clinical.data$Tumor_Sample_Barcode)
  }
  if (!length(ids) && "Tumor_Sample_Barcode" %in% names(maf_data@data)) {
    ids <- as.character(maf_data@data$Tumor_Sample_Barcode)
  }

  clinical_ids <- unique(as.character(clinical_data$Tumor_Sample_Barcode))
  base::intersect(clinical_ids, unique(ids[!is.na(ids) & nzchar(ids)]))
}

get_mutation_filtered_ids <- function(input, cohort_id, row_count) {
  if (row_count == 0) {
    # No rules, don't filter, return all IDs
    return(unique(clinical_data$Tumor_Sample_Barcode))
  }

  rows <- row_count
  assayed_samples <- .mutation_data_available_ids()
  get_ids_for_rule <- function(selection, state) {
    selector <- .parse_mutation_selector(selection)
    if (is.null(selector)) return(character())

    if (identical(selector$type, "gene")) {
      mutated_ids <- maf_data@data[
        as.character(Hugo_Symbol) == selector$gene,
        unique(as.character(Tumor_Sample_Barcode))
      ]
    } else if (identical(selector$type, "codon")) {
      mutated_ids <- .mutation_variant_index[
        gene == selector$gene & ref_aa == selector$ref_aa &
          aa_position == selector$aa_position,
        unique(sample_id)
      ]
    } else {
      mutated_ids <- .mutation_variant_index[
        gene == selector$gene & short_change == selector$short_change,
        unique(sample_id)
      ]
    }

    if (state == "Mutated") return(base::intersect(assayed_samples, mutated_ids))
    base::setdiff(assayed_samples, mutated_ids)
  }

  result_ids <- NULL
  for (i in 1:rows) {
    selection <- .mutation_rule_selection(input, cohort_id, i)
    state <- input[[paste0("state_mut_", i, "_", cohort_id)]]
    logic <- input[[paste0("logic_mut_", i, "_", cohort_id)]]

    if (is.null(selection) || !nzchar(selection) ||
        is.null(state) || !state %in% c("Mutated", "Not Mutated")) next

    ids <- get_ids_for_rule(selection, state)

    if (is.null(result_ids)) {
      result_ids <- ids
    } else {
      if (logic == "AND") {
        result_ids <- intersect(result_ids, ids)
      } else if (logic == "OR") {
        result_ids <- union(result_ids, ids)
      }
    }

    if (logic == "END") break
  }

  if (is.null(result_ids)) {
    return(unique(clinical_data$Tumor_Sample_Barcode))
  }
  result_ids
}


filter_by_gene_expression <- function(clinical_data, gene = NULL,
                                      threshold_type = "value",
                                      min_value = NULL, max_value = NULL,
                                      min_percentile = 0, max_percentile = 100) {
  if (is.null(gene) || gene == "" || !gene %in% rownames(bulkseq_tpm)) return(clinical_data)

  gene_expr <- bulkseq_tpm[gene, ]
  names(gene_expr) <- colnames(bulkseq_tpm)
  gene_expr <- gene_expr[is.finite(gene_expr)]
  if (!length(gene_expr)) return(clinical_data[0, , drop = FALSE])

  if (threshold_type == "value") {
    if (is.null(min_value) || is.null(max_value) || is.na(min_value) || is.na(max_value)) return(clinical_data)
    keep_ids <- names(gene_expr)[gene_expr >= min_value & gene_expr <= max_value]
  } else {
    # Percentile range logic
    pct <- .sanitize_percentile_range(min_percentile, max_percentile)
    if (is.null(pct)) return(clinical_data)
    lower_cutoff <- quantile(gene_expr, probs = pct$min / 100)
    upper_cutoff <- quantile(gene_expr, probs = pct$max / 100)
    keep_ids <- names(gene_expr)[gene_expr >= lower_cutoff & gene_expr <= upper_cutoff]
  }

  clinical_data <- clinical_data[clinical_data$Tumor_Sample_Barcode %in% keep_ids, ]
  return(clinical_data)
}

filter_by_survival <- function(clinical_data, surv_var, threshold_type,
                               min_value = NULL, max_value = NULL,
                               min_percentile = NULL, max_percentile = NULL,
                               require_event = FALSE, event_var = NULL) {
  if (is.null(surv_var) || !(surv_var %in% colnames(clinical_data))) return(clinical_data)

  surv_data <- clinical_data[[surv_var]]

  if (threshold_type == "percentile") {
    if (is.null(min_percentile) || is.null(max_percentile) ||
        is.na(min_percentile) || is.na(max_percentile)) return(clinical_data)
    pct <- .sanitize_percentile_range(min_percentile, max_percentile)
    if (is.null(pct)) return(clinical_data)
    min_thresh <- quantile(surv_data, probs = pct$min / 100, na.rm = TRUE)
    max_thresh <- quantile(surv_data, probs = pct$max / 100, na.rm = TRUE)
  } else {
    if (is.null(min_value) || is.null(max_value) || is.na(min_value) || is.na(max_value)) return(clinical_data)
    min_thresh <- min_value
    max_thresh <- max_value
  }

  filtered <- clinical_data[surv_data >= min_thresh & surv_data <= max_thresh, ]
  if (isTRUE(require_event) && !is.null(event_var) && event_var %in% colnames(filtered)) {
    filtered <- filtered[filtered[[event_var]] == 1 & !is.na(filtered[[event_var]]), ]
  }
  return(filtered)
}

survival_event_column <- function(surv_var) {
  switch(
    surv_var,
    PFS = "PFS_event",
    OS = "OS_event",
    ttct2line = "censt2line",
    NULL
  )
}

.short_protein_change <- function(hgvsp) {
  aa_codes <- c(
    Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C",
    Gln = "Q", Glu = "E", Gly = "G", His = "H", Ile = "I",
    Leu = "L", Lys = "K", Met = "M", Phe = "F", Pro = "P",
    Ser = "S", Thr = "T", Trp = "W", Tyr = "Y", Val = "V",
    Ter = "*", Sec = "U"
  )

  result <- sub("^.*p\\.", "", as.character(hgvsp))
  for (code in names(aa_codes)) {
    result <- gsub(code, aa_codes[[code]], result, fixed = TRUE)
  }
  result
}

.protein_position <- function(hgvsp) {
  change <- sub("^.*p\\.", "", as.character(hgvsp))
  has_position <- grepl("[0-9]+", change)
  position <- rep(NA_real_, length(change))
  position[has_position] <- suppressWarnings(as.numeric(
    sub("^.*?([0-9]+).*$", "\\1", change[has_position])
  ))
  position
}

.parse_mutation_selector <- function(value) {
  if (is.null(value) || !length(value) || !nzchar(value[1])) return(NULL)
  value <- as.character(value[1])
  parts <- strsplit(value, "::", fixed = TRUE)[[1]]

  if (length(parts) == 2 && identical(parts[1], "gene")) {
    return(list(type = "gene", gene = parts[2]))
  }
  if (length(parts) == 4 && identical(parts[1], "codon")) {
    position <- suppressWarnings(as.numeric(parts[4]))
    if (!is.finite(position)) return(NULL)
    return(list(
      type = "codon",
      gene = parts[2],
      ref_aa = parts[3],
      aa_position = position
    ))
  }
  if (length(parts) == 3 && identical(parts[1], "variant")) {
    return(list(type = "variant", gene = parts[2], short_change = parts[3]))
  }

  # Preserve compatibility with gene-only values from an existing session.
  list(type = "gene", gene = value)
}

.mutation_selector_display <- function(value) {
  selector <- .parse_mutation_selector(value)
  if (is.null(selector)) return("")
  if (identical(selector$type, "gene")) return(selector$gene)
  if (identical(selector$type, "codon")) {
    return(paste0(selector$gene, " ", selector$ref_aa, as.integer(selector$aa_position)))
  }
  paste0(selector$gene, " ", selector$short_change)
}

.mutation_rule_selection <- function(input, cohort_id, row_index) {
  suffix <- paste0("_", row_index, "_", cohort_id)
  gene <- input[[paste0("gene_mut", suffix)]]
  codon_value <- input[[paste0("codon_mut", suffix)]]
  variant_value <- input[[paste0("variant_mut", suffix)]]

  if (is.null(gene) || !length(gene) || !nzchar(gene[1])) return(NULL)
  gene <- as.character(gene[1])

  codon <- .parse_mutation_selector(codon_value)
  codon_is_valid <- !is.null(codon) && identical(codon$type, "codon") &&
    identical(codon$gene, gene)

  variant <- .parse_mutation_selector(variant_value)
  variant_is_valid <- !is.null(variant) && identical(variant$type, "variant") &&
    identical(variant$gene, gene)
  if (variant_is_valid && codon_is_valid) {
    variant_position <- .protein_position(variant$short_change)
    variant_is_valid <- length(variant_position) == 1 && is.finite(variant_position) &&
      identical(substr(variant$short_change, 1, 1), codon$ref_aa) &&
      identical(as.numeric(variant_position), as.numeric(codon$aa_position))
  }

  if (variant_is_valid) return(as.character(variant_value[1]))
  if (codon_is_valid) return(as.character(codon_value[1]))
  paste("gene", gene, sep = "::")
}

.build_mutation_filter_catalog <- function() {
  assayed_samples <- .mutation_data_available_ids()
  gene_counts <- maf_data@data[
    as.character(Tumor_Sample_Barcode) %in% assayed_samples,
    .(
    patients = data.table::uniqueN(as.character(Tumor_Sample_Barcode))
    ),
    by = .(gene = as.character(Hugo_Symbol))
  ]
  gene_counts <- gene_counts[!is.na(gene) & nzchar(gene)]
  gene_counts[, `:=`(
    type_rank = 1L,
    type = "gene",
    ref_aa = NA_character_,
    aa_position = NA_real_,
    short_change = NA_character_,
    value = paste("gene", gene, sep = "::"),
    label = sprintf("%s - any mutation (%d patients)", gene, patients)
  )]

  codon_counts <- .mutation_variant_index[, .(
    patients = data.table::uniqueN(sample_id)
  ), by = .(gene, ref_aa, aa_position)]
  codon_counts[, `:=`(
    type_rank = 2L,
    type = "codon",
    short_change = NA_character_,
    value = paste("codon", gene, ref_aa, as.integer(aa_position), sep = "::"),
    label = sprintf(
      "%s %s%d - any change at codon (%d patients)",
      gene, ref_aa, as.integer(aa_position), patients
    )
  )]

  variant_counts <- .mutation_variant_index[, .(
    patients = data.table::uniqueN(sample_id)
  ), by = .(gene, short_change)]
  variant_counts[, `:=`(
    type_rank = 3L,
    type = "variant",
    ref_aa = substr(short_change, 1, 1),
    aa_position = .protein_position(short_change),
    value = paste("variant", gene, short_change, sep = "::"),
    label = sprintf("%s %s (%d patients)", gene, short_change, patients)
  )]

  catalog <- data.table::rbindlist(list(
    gene_counts[, .(type, gene, ref_aa, aa_position, short_change, value, label, patients, type_rank)],
    codon_counts[, .(type, gene, ref_aa, aa_position, short_change, value, label, patients, type_rank)],
    variant_counts[, .(type, gene, ref_aa, aa_position, short_change, value, label, patients, type_rank)]
  ))
  catalog <- unique(catalog, by = "value")
  catalog <- catalog[order(-patients, type_rank, label)]
  catalog
}

.mutation_variant_index <- data.table::copy(maf_data@data[, .(
  gene = as.character(Hugo_Symbol),
  sample_id = as.character(Tumor_Sample_Barcode),
  hgvsp = as.character(HGVSp)
)])
.mutation_variant_index[, `:=`(
  short_change = .short_protein_change(hgvsp),
  aa_position = .protein_position(hgvsp)
)]
.mutation_variant_index[, ref_aa := substr(short_change, 1, 1)]
.mutation_variant_index <- .mutation_variant_index[
  !is.na(gene) & nzchar(gene) &
    !is.na(sample_id) & nzchar(sample_id) &
    sample_id %in% .mutation_data_available_ids() &
    is.finite(aa_position) & grepl("^[A-Z*][0-9]+", short_change)
]

mutation_filter_catalog <- .build_mutation_filter_catalog()

mutation_gene_choices <- function() {
  choices <- mutation_filter_catalog[type == "gene"]
  labels <- sprintf("%s (%d patients)", choices$gene, choices$patients)
  c("Select a gene" = "", stats::setNames(choices$gene, labels))
}

mutation_codon_choices <- function(gene) {
  if (is.null(gene) || !length(gene) || !nzchar(gene[1])) {
    return(c("Any codon" = ""))
  }
  target_gene <- as.character(gene[1])
  choices <- mutation_filter_catalog[type == "codon" & gene == target_gene]
  labels <- sprintf(
    "%s%d (%d patients)", choices$ref_aa, as.integer(choices$aa_position), choices$patients
  )
  c("Any codon" = "", stats::setNames(choices$value, labels))
}

mutation_variant_choices <- function(gene, codon_value = NULL) {
  if (is.null(gene) || !length(gene) || !nzchar(gene[1])) {
    return(c("Any variant" = ""))
  }

  target_gene <- as.character(gene[1])
  choices <- mutation_filter_catalog[type == "variant" & gene == target_gene]
  codon <- .parse_mutation_selector(codon_value)
  if (!is.null(codon) && identical(codon$type, "codon") && identical(codon$gene, target_gene)) {
    choices <- choices[
      ref_aa == codon$ref_aa & aa_position == codon$aa_position
    ]
  }
  labels <- sprintf("%s (%d patients)", choices$short_change, choices$patients)
  c("Any variant" = "", stats::setNames(choices$value, labels))
}

summarize_lollipop_variants <- function(maf, gene) {
  required <- c(
    "Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification",
    "HGVSp", "HGVSc"
  )
  if (!all(required %in% names(maf@data))) return(data.table())

  variants <- data.table::copy(maf@data[as.character(Hugo_Symbol) == gene])
  variants <- variants[
    !is.na(HGVSp) & nzchar(trimws(as.character(HGVSp))) &
      !is.na(Tumor_Sample_Barcode)
  ]
  if (!nrow(variants)) return(data.table())

  variants[, aa_position := .protein_position(HGVSp)]
  variants <- variants[is.finite(aa_position)]
  if (!nrow(variants)) return(data.table())

  cohort_ids <- character()
  if ("Tumor_Sample_Barcode" %in% names(maf@clinical.data)) {
    cohort_ids <- unique(as.character(maf@clinical.data$Tumor_Sample_Barcode))
  }
  if (!length(cohort_ids)) {
    cohort_ids <- unique(as.character(maf@data$Tumor_Sample_Barcode))
  }
  cohort_size <- length(cohort_ids[!is.na(cohort_ids) & nzchar(cohort_ids)])

  summarized <- variants[, .(
    sample_count = data.table::uniqueN(as.character(Tumor_Sample_Barcode)),
    classification = paste(
      sort(unique(as.character(Variant_Classification))),
      collapse = ", "
    ),
    coding_change = paste(
      sort(unique(as.character(HGVSc[!is.na(HGVSc) & nzchar(HGVSc)]))),
      collapse = ", "
    )
  ), by = .(
    gene = as.character(Hugo_Symbol),
    protein_change = as.character(HGVSp),
    aa_position
  )]

  summarized[, `:=`(
    short_change = .short_protein_change(protein_change),
    cohort_size = cohort_size,
    cohort_percent = if (cohort_size > 0) sample_count / cohort_size * 100 else NA_real_
  )]
  summarized[, hover_text := sprintf(
    paste0(
      "<b>%s %s</b>",
      "<br>Samples: %d of %d",
      "<br>Cohort frequency: %.1f%%",
      "<br>Protein position: %d",
      "<br>Classification: %s",
      "%s"
    ),
    gene,
    short_change,
    sample_count,
    cohort_size,
    cohort_percent,
    as.integer(aa_position),
    gsub("_", " ", classification, fixed = TRUE),
    ifelse(nzchar(coding_change), paste0("<br>Coding change: ", coding_change), "")
  )]

  summarized[order(aa_position, -sample_count)]
}

.lollipop_domains <- function(gene, fallback_length) {
  domains <- tryCatch(
    suppressMessages(
      getFromNamespace(".getdomains", "maftools")(geneID = gene)
    ),
    error = function(e) NULL
  )

  if (is.null(domains) || !nrow(domains)) {
    return(list(data = data.table(), protein_length = fallback_length, transcript = NULL))
  }

  domains <- data.table::as.data.table(domains)
  domains <- domains[is.finite(Start) & is.finite(End)]
  protein_length <- suppressWarnings(max(as.numeric(domains$aa.length), na.rm = TRUE))
  if (!is.finite(protein_length)) protein_length <- fallback_length
  transcript <- unique(as.character(domains$refseq.ID))
  transcript <- transcript[!is.na(transcript) & nzchar(transcript)]

  list(
    data = unique(domains[, .(Start, End, Label)]),
    protein_length = max(protein_length, fallback_length),
    transcript = if (length(transcript)) transcript[1] else NULL
  )
}

lollipop_variant_plot <- function(maf, gene, cohort_label) {
  variants <- summarize_lollipop_variants(maf, gene)
  if (!nrow(variants)) {
    return(.empty_plot_message(paste("No protein variants found for", gene)))
  }

  max_count <- max(variants$sample_count)
  max_position <- max(variants$aa_position)
  domain_info <- .lollipop_domains(gene, fallback_length = max_position)
  domains <- domain_info$data
  protein_length <- max(domain_info$protein_length, max_position)
  domain_height <- max(0.18, max_count * 0.035)
  domain_gap <- domain_height * 0.18
  domain_bottom <- -(domain_height + domain_gap)

  variant_classes <- unique(variants$classification)
  class_palette <- c(
    Missense_Mutation = "#238B8D",
    Nonsense_Mutation = "#D1495B",
    Frame_Shift_Del = "#7A5195",
    Frame_Shift_Ins = "#955196",
    In_Frame_Del = "#E59F3A",
    In_Frame_Ins = "#F2C14E",
    Splice_Site = "#4D648D",
    Nonstop_Mutation = "#8F2D56"
  )
  missing_classes <- setdiff(variant_classes, names(class_palette))
  if (length(missing_classes)) {
    fallback_colors <- grDevices::hcl.colors(length(missing_classes), "Dark 3")
    class_palette <- c(class_palette, stats::setNames(fallback_colors, missing_classes))
  }

  plot <- ggplot() +
    geom_hline(yintercept = 0, color = "#68717B", linewidth = 0.5) +
    geom_segment(
      data = variants,
      aes(x = aa_position, xend = aa_position, y = 0, yend = sample_count),
      color = "#AAB1B8",
      linewidth = 0.7
    )

  if (nrow(domains)) {
    domain_labels <- unique(as.character(domains$Label))
    domain_palette <- stats::setNames(
      grDevices::hcl.colors(length(domain_labels), "Set 2"),
      domain_labels
    )
    domains[, label_text := gsub("_", " ", as.character(Label), fixed = TRUE)]
    domains[, label_x := (Start + End) / 2]
    domains[, domain_track := seq_len(.N)]
    domains[, domain_ymax := -(domain_track - 1) * (domain_height + domain_gap) - domain_gap]
    domains[, domain_ymin := domain_ymax - domain_height]
    domain_bottom <- min(domains$domain_ymin) - domain_gap

    plot <- plot +
      geom_rect(
        data = domains,
        aes(xmin = Start, xmax = End, ymin = domain_ymin, ymax = domain_ymax, fill = Label),
        color = "#4E5964",
        linewidth = 0.35,
        alpha = 0.9
      ) +
      geom_text(
        data = domains,
        aes(x = label_x, y = (domain_ymin + domain_ymax) / 2, label = label_text),
        size = 3,
        check_overlap = TRUE
      ) +
      scale_fill_manual(values = domain_palette, guide = "none")
  }

  subtitle <- sprintf(
    "%s | %d samples with mutation data%s",
    cohort_label,
    unique(variants$cohort_size)[1],
    if (!is.null(domain_info$transcript)) paste0(" | ", domain_info$transcript) else ""
  )

  point_layer <- suppressWarnings(
    geom_point(
      data = variants,
      aes(
        x = aa_position,
        y = sample_count,
        color = classification,
        text = hover_text
      ),
      size = 2,
      alpha = 0.6
    )
  )

  plot +
    point_layer +
    scale_color_manual(
      values = class_palette,
      breaks = variant_classes,
      labels = gsub("_", " ", variant_classes, fixed = TRUE)
    ) +
    scale_x_continuous(
      limits = c(0, protein_length),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_y_continuous(
      breaks = pretty(c(0, max_count)),
      limits = c(domain_bottom * 1.12, max_count * 1.12),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = paste(gene, "protein variants"),
      subtitle = subtitle,
      x = "Amino acid position",
      y = "Samples",
      color = "Variant classification"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(color = "#56616C"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )
}

lollipop_variant_plotly <- function(plot) {
  widget <- suppressWarnings(ggplotly(plot, tooltip = "text"))

  for (i in seq_along(widget$x$data)) {
    trace <- widget$x$data[[i]]
    if (identical(trace$fill, "toself")) {
      widget$x$data[[i]]$showlegend <- FALSE
      widget$x$data[[i]]$hoverinfo <- "skip"
    } else if (identical(trace$mode, "markers") && !is.null(trace$name)) {
      clean_name <- sub("^\\((.*),1\\)$", "\\1", trace$name)
      clean_name <- gsub("_", " ", clean_name, fixed = TRUE)
      clean_name <- switch(
        clean_name,
        "Missense Mutation" = "Missense",
        "Nonsense Mutation" = "Nonsense",
        "In Frame Del" = "In-frame deletion",
        "In Frame Ins" = "In-frame insertion",
        "Frame Shift Del" = "Frameshift deletion",
        "Frame Shift Ins" = "Frameshift insertion",
        clean_name
      )
      widget$x$data[[i]]$name <- clean_name
      widget$x$data[[i]]$legendgroup <- clean_name
    } else if (!identical(trace$mode, "markers")) {
      widget$x$data[[i]]$hoverinfo <- "skip"
    }
  }

  widget %>%
    layout(
      hovermode = "closest",
      legend = list(
        orientation = "h",
        x = 0,
        xanchor = "left",
        y = -0.22,
        yanchor = "top",
        title = list(text = ""),
        font = list(size = 10)
      ),
      margin = list(b = 95)
    ) %>%
    config(displaylogo = FALSE)
}

.panel_help_text <- function(id) {
  help_text <- c(
    # how_it_works = "Load optional public ID lists, name cohorts, then apply sidebar filters to define analysis groups.",
    cohort1_upload = "Name this cohort, load optional public IDs, or clear loaded IDs to remove that constraint.",
    cohort2_upload = "Name this cohort, load optional public IDs, or clear loaded IDs to remove that constraint.",
    cohort_complement = "Optionally define one cohort as all eligible samples not selected by the other cohort.",
    cohort_definitions = "Shows cohort names and filter rules captured the last time filters were applied.",
    summaryPlot_g1 = "Summarizes demographics and key clinical features for this cohort. Use counts to check cohort composition.",
    summaryPlot_g2 = "Summarizes demographics and key clinical features for this cohort. Use counts to check cohort composition.",
    survCompPlot_pfs = "Progression-free survival comparison. Wider curve separation suggests different PFS between cohorts.",
    survCompPlot_tt2Line_censored = "Time-to-second-line comparison. Curves show time until second line.",
    survCompPlot_os = "Overall survival comparison. Wider curve separation suggests different OS between cohorts.",
    clin_distribution = "Shows the selected clinical feature by cohort. Compare distributions before interpreting p-values.",
    significance_table = "Tests clinical feature differences between cohorts.",
    cox_model = "Select endpoint, cohorts, covariates, and stratification before fitting the Cox model.",
    cox_results = "Cox model hazard ratios. HR greater than 1 means higher hazard after selected adjustments.",
    cox_forest = "Forest plot of Cox hazard ratios.",
    mafSummary_g1 = "Mutation summary for this cohort, including variant classes and sample-level mutation burden.",
    mafSummary_g2 = "Mutation summary for this cohort, including variant classes and sample-level mutation burden.",
    oncoplot_g1 = "Top mutated genes in this cohort. Columns are samples and colors mark mutation types.",
    oncoplot_g2 = "Top mutated genes in this cohort. Columns are samples and colors mark mutation types.",
    lollipopPlot_g1 = "Protein-position mutation plot for this cohort. Hover over a point for the variant and cohort frequency.",
    lollipopPlot_g2 = "Protein-position mutation plot for this cohort. Hover over a point for the variant and cohort frequency.",
    interactionPlot_g1 = "Gene mutation co-occurrence and exclusivity patterns. Significant pairs are non-random.",
    interactionPlot_g2 = "Gene mutation co-occurrence and exclusivity patterns. Significant pairs are non-random.",
    mafCompOncoPlot = "Side-by-side mutation view for selected genes across cohorts.",
    mafCompBarPlot = "Compares selected gene mutation frequencies between cohorts.",
    mafCompTable = "Genes with different mutation frequencies between cohorts from maftools comparison.",
    mafCompForestPlot = "Mutation enrichment forest plot. Direction reflects relative mutation frequency between cohorts.",
    tpm_distr = "Bulk RNA expression distribution for the selected gene by cohort.",
    tpm_distr_boxplot = "Bulk RNA expression boxplot with cohort comparison p-value.",
    tpm_survCompPlot_g1 = "Survival within this cohort after splitting patients by selected gene expression.",
    tpm_survCompPlot_g2 = "Survival within this cohort after splitting patients by selected gene expression.",
    quantile_table_cohort1 = "Expression cutoffs and sample counts for the selected gene in this cohort.",
    quantile_table_cohort2 = "Expression cutoffs and sample counts for the selected gene in this cohort.",
    bulkVolcano = "Bulk DESeq2 volcano plot. Far left or right points with high y-values are genes down-regulated or up-regulated in the first cohot (Cohort 1).",
    DEGs_table = "Thresholded bulk DESeq2 results. Adjust p-value and fold-change cutoffs above.",
    ssgsea_violin = "Gene-set enrichment scores by cohort. Higher values indicate stronger pathway activity.",
    ssgsea_table = "Statistical table for gene-set enrichment differences between cohorts. Up means this pathway is up-regulated in the first cohort (Cohort 1), vice versa.",
    sc_celltype_boxplot = "Patient-level cell type proportions by cohort. Stars mark differential abundance tests.",
    sc_celltype_proportion = "Stacked cell type composition by cohort for broad immune shifts.",
    sc_cellcycle_hist = "Cell-cycle phase distribution for selected cell types by cohort.",
    pseudo_norm_distr = "Pseudobulk expression distribution for the selected cell type and gene.",
    pseudo_norm_distr_boxplot = "Pseudobulk expression boxplot with cohort comparison p-value.",
    pseudo_norm_survCompPlot_g1 = "Survival within this cohort after splitting by pseudobulk expression.",
    pseudo_norm_survCompPlot_g2 = "Survival within this cohort after splitting by pseudobulk expression.",
    pseudo_quantile_table_cohort1 = "Pseudobulk expression cutoffs and sample counts for this cohort.",
    pseudo_quantile_table_cohort2 = "Pseudobulk expression cutoffs and sample counts for this cohort.",
    pseudoVolcano = "Cell-type-specific pseudobulk DESeq2 volcano plot.",
    pseudo_DEGs_table = "Thresholded pseudobulk differential expression results for the selected cell type."
  )
  text <- unname(help_text[id])
  if (is.na(text)) {
    "Explains the inputs or results in this panel. Use it as context for interpretation."
  } else {
    text
  }
}

panel_help_ui <- function(id, text = NULL) {
  if (is.null(text)) text <- .panel_help_text(id)
  div(
    class = "panel-help-control",
    `data-help` = text,
    tags$button(
      type = "button",
      class = "btn btn-default btn-xs panel-help-toggle",
      `aria-label` = text,
      icon("question-circle")
    )
  )
}

plot_export_controls_ui <- function(id, width = 7, height = 5, units = "in") {
  div(
    class = "plot-export-controls",
    div(
      class = "plot-export-dropdown dropdown",
      tags$button(
        type = "button",
        class = "btn btn-default btn-xs dropdown-toggle plot-export-toggle",
        `data-toggle` = "dropdown",
        `aria-haspopup` = "true",
        `aria-expanded` = "false",
        title = "Download PDF",
        icon("download")
      ),
      tags$ul(
        class = "dropdown-menu dropdown-menu-right plot-export-menu",
        onclick = "event.stopPropagation();",
        tags$li(
          tags$div(
            class = "plot-export-menu-body",
            numericInput(paste0(id, "_pdf_width"), "Width", value = width, min = 0.1, step = 0.5),
            numericInput(paste0(id, "_pdf_height"), "Height", value = height, min = 0.1, step = 0.5),
            selectInput(paste0(id, "_pdf_units"), "Units", choices = c("in", "cm", "mm"), selected = units),
            downloadButton(paste0("download_", id, "_pdf"), "Download PDF", class = "btn-primary btn-block")
          )
        )
      )
    ),
    panel_help_ui(id)
  )
}

.pdf_dim_to_inches <- function(value, units) {
  value <- as.numeric(value)
  if (!is.finite(value) || value <= 0) value <- 7
  switch(
    units,
    cm = value / 2.54,
    mm = value / 25.4,
    value
  )
}

render_plot_for_pdf <- function(file, draw_plot, width = 7, height = 5, units = "in") {
  width_in <- .pdf_dim_to_inches(width, units)
  height_in <- .pdf_dim_to_inches(height, units)
  grDevices::pdf(file, width = width_in, height = height_in, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  p <- draw_plot()
  if (is.null(p)) return(invisible(NULL))
  print(p)
  invisible(p)
}

register_plot_pdf_download <- function(output, input, id, draw_plot, filename = NULL) {
  output[[paste0("download_", id, "_pdf")]] <- downloadHandler(
    filename = function() {
      if (is.null(filename)) {
        paste0(id, "_", Sys.Date(), ".pdf")
      } else if (is.function(filename)) {
        filename()
      } else {
        filename
      }
    },
    content = function(file) {
      render_plot_for_pdf(
        file = file,
        draw_plot = draw_plot,
        width = input[[paste0(id, "_pdf_width")]],
        height = input[[paste0(id, "_pdf_height")]],
        units = input[[paste0(id, "_pdf_units")]]
      )
    }
  )
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
      `none-selected-text` = "All"
    ),
    multiple = TRUE,
    selected = NULL
  )
}

apply_filter <- function(data, column, values) {
  if (!is.null(values) && length(values) > 0) {
    data %>% filter(!!sym(column) %in% values)
  } else {
    data  # If nothing selected, keep all (default "All")
  }
}

create_cohort_filters_ui <- function(cohort_id, category) {
  tagList(
    div(id = paste0(cohort_id, "_", category, "_filters"),
        if (category == "clinical") {
          return(list(
            # Demographics
            create_picker_input(paste0("sex_filter_", cohort_id), "Sex", sort(unique(clinical_data$Sex))),
            create_picker_input(paste0("race_filter_", cohort_id), "Race (self-reported)", sort(unique(clinical_data$Race))),
            create_picker_input(paste0("genetic_ancestry_filter_", cohort_id), "Genetic Ancestry", sort(unique(clinical_data$genetic_ancestry))),
            uiOutput(paste0("age_filter_", cohort_id)),

            # Clinical classification
            create_picker_input(paste0("stage_filter_", cohort_id), "ISS Stage", sort(unique(clinical_data$ISS))),
            create_picker_input(paste0("risk_filter_", cohort_id), "IMWG Risk Classification", sort(unique(clinical_data$IMWG_Risk_Class))),
            create_picker_input(paste0("cyto_risk_filter_", cohort_id), "Cytogenetic High Risk (Skerget)", sort(unique(clinical_data$Skerget_Cytogenetic_High_Risk))),
            create_picker_input(paste0("cgs_risk_filter_", cohort_id), "CGS Risk (Consensus Genomic Staging)", sort(unique(clinical_data$CGS_risk))),

            # Subtypes
            create_picker_input(paste0("rna_subtype_filter_", cohort_id), "RNA Subtype (Skerget)", sort(unique(clinical_data$Skerget_RNA_Subtype_Name))),
            create_picker_input(paste0("cna_subtype_filter_", cohort_id), "CNA Subtype (Skerget)", sort(unique(clinical_data$Skerget_CNA_Subtype_Name))),

            # Treatment
            create_picker_input(paste0("triplet_filter_", cohort_id), "Triplet Firstline", sort(unique(clinical_data$Triplet_First))),
            create_picker_input(paste0("asct_filter_", cohort_id), "ASCT Firstline", sort(unique(clinical_data$ASCT_First))),

            # Regimen
            create_picker_input(paste0("regimen_filter_", cohort_id), "Regimen Firstline", sort(unique(clinical_data$regimen)))
          ))
        } else if (category == "molecular") {
          return(list(
            # Chromosomal abnormalities
            create_picker_input(paste0("chr_1q21_gain_filter_", cohort_id), "1q21 Gain", sort(unique(clinical_data$chr_1q21_gain))),
            create_picker_input(paste0("chr_1q21_amp_filter_", cohort_id), "1q21 Amplification", sort(unique(clinical_data$chr_1q21_amp))),
            create_picker_input(paste0("chr_13q14_del_filter_", cohort_id), "13q14 Deletion", sort(unique(clinical_data$chr_13q14_del))),
            create_picker_input(paste0("chr_13q34_del_filter_", cohort_id), "13q34 Deletion", sort(unique(clinical_data$chr_13q34_del))),
            create_picker_input(paste0("chr_17p13_del_filter_", cohort_id), "17p13 Deletion", sort(unique(clinical_data$chr_17p13_del))),
            create_picker_input(paste0("diploidy_filter_", cohort_id), "Hyperdiploidy", sort(unique(clinical_data$Hyperdiploidy))),
            create_picker_input(paste0("chromothripsis_filter_", cohort_id), "Chromothripsis", sort(unique(clinical_data$chromothripsis))),

            # Translocations
            create_picker_input(paste0("t_11_14_filter_", cohort_id), "t(11;14)", sort(unique(clinical_data$t_11_14))),
            create_picker_input(paste0("t_4_14_filter_", cohort_id), "t(4;14)", sort(unique(clinical_data$t_4_14))),

            # Mutational markers
            create_picker_input(paste0("maf_filter_", cohort_id), "MAF/MAFB", sort(unique(clinical_data$MAF_MAFB))),
            create_picker_input(paste0("apobec_filter_", cohort_id), "APOBEC", sort(unique(clinical_data$APOBEC))),
            create_picker_input(paste0("tp53_filter_", cohort_id), "TP53 Functional Copies", sort(unique(clinical_data$TP53_Funct_Copies))),
            create_picker_input(paste0("tp53_ns_filter_", cohort_id), "TP53 Non-Synonymous Mutation Count", sort(unique(clinical_data$TP53_NS_Mut_Count)))
          ))
        }
    )
  )
}

get_cohort_filters <- function(input, cohort_id, category) {
  if (category == "clinical") {
    return(list(
      sex = input[[paste0("sex_filter_", cohort_id)]],
      race = input[[paste0("race_filter_", cohort_id)]],
      genetic_ancestry = input[[paste0("genetic_ancestry_filter_", cohort_id)]],
      age = input[[paste0("age_", cohort_id)]],

      stage = input[[paste0("stage_filter_", cohort_id)]],
      risk = input[[paste0("risk_filter_", cohort_id)]],
      cyto_risk = input[[paste0("cyto_risk_filter_", cohort_id)]],
      cgs_risk  = input[[paste0("cgs_risk_filter_", cohort_id)]],

      rna_subtype = input[[paste0("rna_subtype_filter_", cohort_id)]],
      cna_subtype = input[[paste0("cna_subtype_filter_", cohort_id)]],

      triplet = input[[paste0("triplet_filter_", cohort_id)]],
      asct = input[[paste0("asct_filter_", cohort_id)]],

      regimen = input[[paste0("regimen_filter_", cohort_id)]]
    ))
  } else if (category == "molecular") {
    return(list(
      q21_gain = input[[paste0("chr_1q21_gain_filter_", cohort_id)]],
      q21_amp = input[[paste0("chr_1q21_amp_filter_", cohort_id)]],
      del13q14 = input[[paste0("chr_13q14_del_filter_", cohort_id)]],
      del13q34 = input[[paste0("chr_13q34_del_filter_", cohort_id)]],
      del17p13 = input[[paste0("chr_17p13_del_filter_", cohort_id)]],

      diploidy = input[[paste0("diploidy_filter_", cohort_id)]],
      chromothripsis = input[[paste0("chromothripsis_filter_", cohort_id)]],

      t11_14 = input[[paste0("t_11_14_filter_", cohort_id)]],
      t4_14 = input[[paste0("t_4_14_filter_", cohort_id)]],

      maf = input[[paste0("maf_filter_", cohort_id)]],
      apobec = input[[paste0("apobec_filter_", cohort_id)]],
      tp53 = input[[paste0("tp53_filter_", cohort_id)]],
      tp53_ns = input[[paste0("tp53_ns_filter_", cohort_id)]]
    ))
  }
}

filter_cohort_data <- function(data, filters) {
  data <- as_tibble(data)  # For compatibility with dplyr pipes

  data <- data %>%
    # Demographic filters
    apply_filter("Sex", filters$sex) %>%
    apply_filter("Race", filters$race) %>%
    apply_filter("genetic_ancestry", filters$genetic_ancestry) %>%
    {
      if (!is.null(filters$age) && length(filters$age) == 2) {
        filter(., Age >= filters$age[1], Age <= filters$age[2])
      } else {
        .
      }
    } %>%

    # Clinical classification
    apply_filter("ISS", filters$stage) %>%
    apply_filter("IMWG_Risk_Class", filters$risk) %>%
    apply_filter("Skerget_Cytogenetic_High_Risk", filters$cyto_risk) %>%
    apply_filter("CGS_risk", filters$cgs_risk) %>%

    # Subtypes
    apply_filter("Skerget_RNA_Subtype_Name", filters$rna_subtype) %>%
    apply_filter("Skerget_CNA_Subtype_Name", filters$cna_subtype) %>%

    # Treatment
    apply_filter("Triplet_First", filters$triplet) %>%
    apply_filter("ASCT_First", filters$asct) %>%
    apply_filter("regimen", filters$regimen) %>%

    # Chromosomal abnormalities
    apply_filter("chr_1q21_gain", filters$q21_gain) %>%
    apply_filter("chr_1q21_amp", filters$q21_amp) %>%
    apply_filter("chr_13q14_del", filters$del13q14) %>%
    apply_filter("chr_13q34_del", filters$del13q34) %>%
    apply_filter("chr_17p13_del", filters$del17p13) %>%
    apply_filter("Hyperdiploidy", filters$diploidy) %>%
    apply_filter("chromothripsis", filters$chromothripsis) %>%

    # Translocations
    apply_filter("t_11_14", filters$t11_14) %>%
    apply_filter("t_4_14", filters$t4_14) %>%

    # Mutational markers
    apply_filter("MAF_MAFB", filters$maf) %>%
    apply_filter("APOBEC", filters$apobec) %>%
    apply_filter("TP53_Funct_Copies", filters$tp53) %>%
    apply_filter("TP53_NS_Mut_Count", filters$tp53_ns)

  return(data)
}

# Determine selected cohort
get_cohort_selected <- function(filtered_data, cohort_selected) {
  if (cohort_selected == "cohort1") {
    return(filtered_data$cohort1)
  } else {
    return(filtered_data$cohort2)
  }
}

# Summary of clinical data
.empty_plot_message <- function(message) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = message, size = 4) +
    theme_void()
}

generate_summary_plot <- function(cohort_selected, filtered_data) {
  selected_data <- filtered_data[[cohort_selected]]
  features <- c("Race", "Sex", "Age_range", "IMWG_Risk_Class",
                "CGS_risk", "ASCT_First", "Triplet_First", "Hyperdiploidy", "chromothripsis")

  if (is.null(selected_data) || nrow(selected_data) == 0) {
    return(.empty_plot_message("No samples available for this cohort."))
  }

  unique_values_counts_list <- lapply(features, function(feature) {
    if (!feature %in% names(selected_data)) {
      return(data.frame(Value = "No data", Count = 0L, Feature = feature))
    }
    tbl <- table(selected_data[[feature]], useNA = "no")
    if (length(tbl) == 0) {
      return(data.frame(Value = "No data", Count = 0L, Feature = feature))
    }
    data.frame(Value = names(tbl),
               Count = as.integer(tbl),
               Feature = feature)
  })

  unique_values_counts_df <- do.call(rbind, unique_values_counts_list)
  # Ensure Feature is a factor with specified order
  unique_values_counts_df$Feature <- factor(unique_values_counts_df$Feature, levels = features)

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
  cohort_levels <- unique(as.character(data$cohort))
  cohort_colors <- setNames(c("#E87D72", "#5BAEB0")[seq_along(cohort_levels)], cohort_levels)
  plot_data <- data %>%
    group_by(cohort, !!sym(feature)) %>%
    summarise(n = n(), .cohorts = "drop") %>%
    group_by(cohort) %>%
    mutate(
      percentage = n / sum(n) * 100,
      label = paste0("Cohort: ", cohort,
                     "<br>", feature_label, ": ", !!sym(feature),
                     "<br>Percentage: ", sprintf("%.2f", percentage), "%")
    )

  ggplot(plot_data, aes_string(x = feature, y = "percentage", fill = "cohort", alpha = "cohort", text = "label")) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = feature_label, y = "Percentage") +
    scale_alpha_manual(values = setNames(rep(0.9, length(cohort_levels)), cohort_levels)) +
    scale_fill_manual(values = cohort_colors) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Create distribution plot data for continuous variables
create_continuous_plot <- function(data, feature, feature_label) {
  cohort_levels <- unique(as.character(data$cohort))
  cohort_colors <- setNames(c("#E87D72", "#5BAEB0")[seq_along(cohort_levels)], cohort_levels)
  data <- data %>%
    dplyr::mutate(
      label = paste0(
        "Cohort: ", cohort,
        "<br>", feature_label, ": ", sprintf("%.2f", !!rlang::sym(feature))
      )
    )

  ggplot(data, aes_string(x = "cohort", y = feature, fill = "cohort")) +
    # Cleaner boxplot
    geom_boxplot(width = 0.5, outlier.shape = NA, color = "#2b2b2b", alpha = 0.9) +
    # Light jitter to reveal distribution
    geom_jitter(aes_string(text = "label", color = "cohort"),
                width = 0.12, size = 1.6, alpha = 0.45, show.legend = FALSE) +
    # Emphasize the median
    stat_summary(fun = median, geom = "point",
                 shape = 21, size = 3.2, fill = "white", color = "black") +
    scale_fill_manual(values = cohort_colors) +
    scale_color_manual(values = cohort_colors) +
    labs(x = "Cohort", y = feature_label) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(margin = margin(t = 4)),
      axis.title.x = element_text(margin = margin(t = 6)),
      axis.title.y = element_text(margin = margin(r = 6))
    )
}

# Main function to create distribution plot
create_distribution_ggplot <- function(data, feature, feature_label, continuous_features) {
  is_continuous <- feature %in% continuous_features

  if (!is_continuous) {
    p <- create_categorical_plot(data, feature, feature_label)
  } else {
    p <- create_continuous_plot(data, feature, feature_label)
  }

  p
}

create_distribution_plot <- function(data, feature, feature_label, continuous_features) {
  p <- create_distribution_ggplot(data, feature, feature_label, continuous_features)
  ggplotly(p, tooltip = "text") %>%
    layout(hovermode = "x")
}

# -------------------------------------
# Statistical testing
safe_dir <- function(delta, up_lbl = "Cohort1 ↑", down_lbl = "Cohort2 ↑") {
  if (is.na(delta)) return("NA")
  if (delta > 0) up_lbl else if (delta < 0) down_lbl else "Tie"
}

cohens_d <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  n1 <- length(x); n2 <- length(y)
  if (n1 < 2 || n2 < 2) return(NA_real_)
  s1 <- stats::sd(x); s2 <- stats::sd(y)
  sp <- sqrt(((n1 - 1)*s1^2 + (n2 - 1)*s2^2) / (n1 + n2 - 2))
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (mean(x) - mean(y)) / sp
}

# Perform statistical test for continuous variables
test_continuous_variable <- function(cohort1_vals, cohort2_vals) {
  cohort1_vals <- cohort1_vals[!is.na(cohort1_vals)]
  cohort2_vals <- cohort2_vals[!is.na(cohort2_vals)]

  if (length(cohort1_vals) <= 2 || length(cohort2_vals) <= 2) {
    return(list(
      test_name = "Insufficient data",
      statistic = "NA",
      p_value = NA,
      cohort1_summary = "Insufficient data",
      cohort2_summary = "Insufficient data",
      direction = "NA",
      effect = "NA"
    ))
  }

  # normality
  normal1 <- tryCatch(shapiro.test(cohort1_vals)$p.value > 0.05, error = function(e) FALSE)
  normal2 <- tryCatch(shapiro.test(cohort2_vals)$p.value > 0.05, error = function(e) FALSE)

  # defaults
  effect_num <- NA_real_
  effect_label <- "NA"
  direction <- "NA"

  if (normal1 && normal2 && length(cohort1_vals) >= 3 && length(cohort2_vals) >= 3) {
    # t-test
    test_result <- tryCatch(t.test(cohort1_vals, cohort2_vals), error = function(e) NULL)
    test_name <- "t-test"
    if (!is.null(test_result)) {
      statistic <- paste0("t = ", round(test_result$statistic, 3))
      p_value <- test_result$p.value
      # Direction & effect: mean difference (+ Cohen's d)
      m1 <- mean(cohort1_vals); m2 <- mean(cohort2_vals)
      effect_num <- m1 - m2
      d <- cohens_d(cohort1_vals, cohort2_vals)
      direction <- safe_dir(effect_num)
      effect_label <- paste0("Δmean = ", round(effect_num, 2),
                             if (!is.na(d)) paste0(" (Cohen's d = ", round(d, 2), ")") else "")
    } else {
      statistic <- "NA"; p_value <- NA
    }
  } else {
    # Mann-Whitney
    test_result <- tryCatch(wilcox.test(cohort1_vals, cohort2_vals), error = function(e) NULL)
    test_name <- "Mann-Whitney U"
    if (!is.null(test_result)) {
      statistic <- paste0("W = ", round(test_result$statistic, 3))
      p_value <- test_result$p.value
      # Direction & effect: median difference
      med1 <- stats::median(cohort1_vals); med2 <- stats::median(cohort2_vals)
      effect_num <- med1 - med2
      direction <- safe_dir(effect_num)
      effect_label <- paste0("Δmedian = ", round(effect_num, 2))
    } else {
      statistic <- "NA"; p_value <- NA
    }
  }

  cohort1_summary <- paste0("Mean: ", round(mean(cohort1_vals), 2),
                            " (SD: ", round(sd(cohort1_vals), 2), ")")
  cohort2_summary <- paste0("Mean: ", round(mean(cohort2_vals), 2),
                            " (SD: ", round(sd(cohort2_vals), 2), ")")

  list(
    test_name = test_name,
    statistic = statistic,
    p_value = p_value,
    cohort1_summary = cohort1_summary,
    cohort2_summary = cohort2_summary,
    direction = direction,
    effect = effect_label
  )
}

# --------------------------------
# Modify: test_categorical_variable()
test_categorical_variable <- function(data, feature) {
  contingency_table <- tryCatch(table(data$cohort, data[[feature]], useNA = "no"),
                                error = function(e) NULL)

  if (is.null(contingency_table) || sum(contingency_table) == 0) {
    return(list(
      test_name = "Insufficient data",
      statistic = "NA",
      p_value = NA,
      cohort1_summary = "Insufficient data",
      cohort2_summary = "Insufficient data",
      direction = "NA",
      effect = "NA"
    ))
  }

  min_cell_count <- min(contingency_table)
  total_cells <- length(contingency_table)
  cells_less_than_5 <- sum(contingency_table < 5)

  if (min_cell_count == 0 || cells_less_than_5 > (total_cells/2)) {
    return(list(
      test_name = "Low frequency categories",
      statistic = "NA",
      p_value = NA,
      cohort1_summary = "Low frequency",
      cohort2_summary = "Low frequency",
      direction = "NA",
      effect = "NA"
    ))
  }

  expected_freq <- tryCatch(chisq.test(contingency_table)$expected, error = function(e) NULL)
  use_fisher <- !is.null(expected_freq) && any(expected_freq < 5)

  if (use_fisher) {
    if (nrow(contingency_table) == 2 && ncol(contingency_table) == 2) {
      test_result <- tryCatch(fisher.test(contingency_table), error = function(e) NULL)
      test_name <- "Fisher's exact"
      if (!is.null(test_result)) {
        statistic <- paste0("OR = ", round(test_result$estimate, 3))
        p_value <- test_result$p.value
      } else {
        statistic <- "NA"; p_value <- NA
      }
    } else {
      test_result <- tryCatch(chisq.test(contingency_table, simulate.p.value = TRUE, B = 2000),
                              error = function(e) NULL)
      test_name <- "Chi-square (simulated)"
      if (!is.null(test_result)) {
        statistic <- paste0("\u03C7\u00B2 = ", round(test_result$statistic, 3))
        p_value <- test_result$p.value
      } else {
        statistic <- "NA"; p_value <- NA
      }
    }
  } else {
    test_result <- tryCatch(chisq.test(contingency_table), error = function(e) NULL)
    test_name <- "Chi-square"
    if (!is.null(test_result)) {
      statistic <- paste0("\u03C7\u00B2 = ", round(test_result$statistic, 3))
      p_value <- test_result$p.value
    } else {
      statistic <- "NA"; p_value <- NA
    }
  }

  # Summaries
  if ("Cohort1" %in% rownames(contingency_table) && "Cohort2" %in% rownames(contingency_table)) {
    cohort1_counts <- contingency_table["Cohort1", ]
    cohort2_counts <- contingency_table["Cohort2", ]
    cohort1_total <- sum(cohort1_counts)
    cohort2_total <- sum(cohort2_counts)

    cohort1_summary <- paste(names(cohort1_counts),
                             paste0(cohort1_counts, " (", round(cohort1_counts/cohort1_total*100, 1), "%)"),
                             sep = ": ", collapse = "; ")
    cohort2_summary <- paste(names(cohort2_counts),
                             paste0(cohort2_counts, " (", round(cohort2_counts/cohort2_total*100, 1), "%)"),
                             sep = ": ", collapse = "; ")

    # Direction: level with largest absolute proportion difference
    p1 <- cohort1_counts / cohort1_total
    p2 <- cohort2_counts / cohort2_total
    diffs <- p1 - p2
    idx <- which.max(abs(diffs))
    top_level <- names(diffs)[idx]
    delta_p <- diffs[idx] * 100
    direction <- if (delta_p > 0) {
      paste0("Cohort1 \u2191 for ", top_level)
    } else if (delta_p < 0) {
      paste0("Cohort2 \u2191 for ", top_level)
    } else {
      "Tie"
    }
    effect_label <- paste0(top_level, ": \u0394p = ", round(delta_p, 1), "%")
  } else {
    cohort1_summary <- "Data unavailable"
    cohort2_summary <- "Data unavailable"
    direction <- "NA"
    effect_label <- "NA"
  }

  list(
    test_name = test_name,
    statistic = statistic,
    p_value = p_value,
    cohort1_summary = cohort1_summary,
    cohort2_summary = cohort2_summary,
    direction = direction,
    effect = effect_label
  )
}

# Determine significance level from p-value
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

# Format p-value for display
format_p_value <- function(p_value) {
  if (is.na(p_value)) {
    return("NA")
  } else if (p_value < 0.001) {
    return("<0.001")
  } else {
    return(sprintf("%.4f", p_value))
  }
}

# Create complete significance table
create_significance_table <- function(data, clinical_features, continuous_features) {
  results <- data.frame(
    Feature = character(),
    Type = character(),
    Test_Used = character(),
    P_Value = numeric(),
    Statistic = character(),
    Cohort1_Summary = character(),
    Cohort2_Summary = character(),
    Direction = character(),
    Effect = character(),
    Significance = character(),
    stringsAsFactors = FALSE
  )

  for (feature in clinical_features) {
    if (sum(!is.na(data[[feature]])) < 10) next
    is_continuous <- feature %in% continuous_features

    if (is_continuous) {
      cohort1_vals <- data[data$cohort == "Cohort1", feature]
      cohort2_vals <- data[data$cohort == "Cohort2", feature]
      test_results <- test_continuous_variable(cohort1_vals, cohort2_vals)
    } else {
      test_results <- test_categorical_variable(data, feature)
    }

    significance <- get_significance_level(test_results$p_value)

    results <- rbind(results, data.frame(
      Feature = feature,
      Type = ifelse(is_continuous, "Continuous", "Categorical"),
      Test_Used = test_results$test_name,
      P_Value = test_results$p_value,
      Statistic = test_results$statistic,
      Cohort1_Summary = test_results$cohort1_summary,
      Cohort2_Summary = test_results$cohort2_summary,
      Direction = test_results$direction,
      Effect = test_results$effect,
      Significance = significance,
      stringsAsFactors = FALSE
    ))
  }

  results$P_Value_Display <- sapply(results$P_Value, format_p_value)

  display_table <- results[, c("Feature", "Type", "Test_Used", "P_Value_Display",
                               "Statistic", "Direction", "Effect", "Significance")]
  colnames(display_table) <- c("Clinical Feature", "Type", "Statistical Test", "P-Value",
                               "Test Statistic", "Direction", "Effect", "Significance")

  display_table <- display_table[order(results$P_Value, na.last = TRUE), ]
  return(display_table)
}

# ------------------
# Utils
# Get clinical feature choices for selectInput
get_clinical_feature_choices <- function(clinical_data, exclude_cols = NULL) {
  default_exclude <- c("public_id", "Tumor_Sample_Barcode", "Tx",
                       "PFS_event", "PFS", "OS",
                       "OS_event", "PFS_1", "PFS_1_censored", "PFS_1_event", "ttct2line", "censt2line")

  if (!is.null(exclude_cols)) {
    exclude_cols <- c(default_exclude, exclude_cols)
  } else {
    exclude_cols <- default_exclude
  }

  setdiff(colnames(clinical_data), exclude_cols)
}

# Define continuous features
get_continuous_features <- function() {
  c("Age", "BMI", "Serum_B2M", "Serum_LDH", "Creatinine")
}


# Generic cohort survival curve helper -------------------------------------
plot_cohort_survival <- function(data,
                                 time_col,   # string: column name for time
                                 event_col,  # string: column name for event (0/1)
                                 xlab,
                                 break_by = 500,
                                 legend_labs = c("Cohort1", "Cohort2")) {
  stopifnot(all(c("cohort", time_col, event_col) %in% names(data)))
  data <- data[!is.na(data[[time_col]]) & !is.na(data[[event_col]]), ]
  if (nrow(data) == 0 || length(unique(data$cohort)) < 2) {
    return(.empty_plot_message("Not enough survival data for this comparison."))
  }

  group_counts <- table(data$cohort)
  show_interval <- all(group_counts >= 2)
  show_pval <- all(group_counts >= 2)

  surv_obj <- survival::Surv(time = data[[time_col]],
                             event = data[[event_col]])

  fit <- do.call(survival::survfit,
                 list(surv_obj ~ cohort, data = data))

  xmax <- max(data[[time_col]], na.rm = TRUE)
  pval_x <- if (is.finite(xmax)) xmax * 0.6 else break_by

  survminer::ggsurvplot(
    fit, data = data,
    # Core aesthetics
    palette = c("#E41A1C", "#4DBBD5"),  # red = Cohort1, teal = Cohort2
    linetype = c("solid", "solid"),
    size = 1,

    conf.int = show_interval,
    pval = show_pval,
    pval.coord = c(pval_x, 0.1),

    title = "",
    xlab = xlab,
    ylab = "Survival Probability",
    legend.title = "",
    legend.labs = legend_labs,

    # Risk table
    risk.table = TRUE,
    risk.table.height = 0.25,
    risk.table.title = "Number at risk",
    risk.table.fontsize = 3.5,
    tables.theme = survminer::theme_cleantable(),

    ggtheme = ggplot2::theme_bw() + ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(face = "bold", size = 12),
      axis.text = ggplot2::element_text(size = 10),
      legend.position = "top",
      legend.text = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5)
    ),

    break.time.by = break_by,
    surv.scale = "percent"
  )
}


############### Cox-PH ###############
# Pretty HR table from a fitted coxph
.cox_tidy_table <- function(fit) {
  sm <- summary(fit)
  co <- as.data.frame(sm$coefficients)
  ci <- as.data.frame(sm$conf.int)
  # Match rows
  out <- data.frame(
    Term = rownames(co),
    HR = ci$`exp(coef)`,
    CI_low = ci$`lower .95`,
    CI_high = ci$`upper .95`,
    z = co$`z`,
    p = co$`Pr(>|z|)`,
    row.names = NULL,
    check.names = FALSE
  )
  out$HR <- round(out$HR, 3)
  out$CI_low <- round(out$CI_low, 3)
  out$CI_high <- round(out$CI_high, 3)
  out$z <- round(out$z, 3)
  out$p <- signif(out$p, 3)
  out$`HR (95% CI)` <- paste0(out$HR, " (", out$CI_low, "-", out$CI_high, ")")
  out[, c("Term", "HR (95% CI)", "z", "p")]
}

# Maps endpoint to time/event columns
.cox_endpoint_map <- function(ep) {
  if (identical(ep, "OS")) {
    list(
      time  = "OS",
      event = "OS_event",
      label = "Overall Survival"
    )
  } else if (identical(ep, "TT2L")) {
    list(
      time  = "ttct2line",
      event = "censt2line",
      label = "Time to Second Line"
    )
  } else {  # default to PFS
    list(
      time  = "PFS",
      event = "PFS_event",
      label = "Progression-Free Survival"
    )
  }
}

# Converts character columns to factor; leaves numeric as numeric
.coerce_for_cox <- function(df, cols) {
  for (nm in cols) {
    if (!nm %in% names(df)) next
    if (is.character(df[[nm]]) || is.logical(df[[nm]])) df[[nm]] <- factor(df[[nm]])
    # leave numeric/integer as is
  }
  df
}

# Does a factor level produce 0 events or 0 non-events? (separation)
.has_zero_cell <- function(x, event) {
  if (!is.factor(x)) return(FALSE)
  tab <- table(x, event)
  any(tab == 0)  # any level has 0 in a cell
}

# For numeric: zero variance overall or within event/non-event strata
.bad_numeric <- function(x, event) {
  if (!is.numeric(x)) return(FALSE)
  ux <- unique(x[!is.na(x)])
  if (length(ux) <= 1) return(TRUE)
  v0 <- suppressWarnings(stats::var(x[event == 1], na.rm = TRUE))
  v1 <- suppressWarnings(stats::var(x[event == 0], na.rm = TRUE))
  (is.na(v0) || v0 == 0) || (is.na(v1) || v1 == 0)
}

# Drop covariates that would cause infinite estimates
.prune_separating_covars <- function(df, covars, event_col) {
  bad <- vapply(covars, function(v) {
    if (!v %in% names(df)) return(TRUE)
    x <- df[[v]]
    ev <- df[[event_col]]
    if (is.factor(x)) .has_zero_cell(x, ev) else .bad_numeric(x, ev)
  }, logical(1))
  list(keep = covars[!bad], drop = covars[bad])
}

.safe_forest <- function(fit) {
  sm <- summary(fit)
  co <- as.data.frame(sm$coefficients)
  ci <- as.data.frame(sm$conf.int)

  # Build tidy frame (drop non-finite)
  df <- data.frame(
    term     = rownames(co),
    hr       = ci$`exp(coef)`,
    lo       = ci$`lower .95`,
    hi       = ci$`upper .95`,
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$hr) & is.finite(df$lo) & is.finite(df$hi), , drop = FALSE]
  if (!nrow(df)) {
    # nothing plottable
    return(ggplot2::ggplot() + ggplot2::geom_blank() +
             ggplot2::labs(title = "No finite HRs to plot"))
  }

  # clip extreme CIs to keep axis sane; adjust limits from remaining finite values
  lo_min <- max(min(df$lo[df$lo > 0 & is.finite(df$lo)], na.rm = TRUE), 1e-3)
  hi_max <- min(max(df$hi[is.finite(df$hi)], na.rm = TRUE), 1e3)

  df$term <- factor(df$term, levels = rev(df$term))

  ggplot2::ggplot(df, ggplot2::aes(x = term, y = hr)) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = pmax(lo, lo_min), ymax = pmin(hi, hi_max)), width = 0.2) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_log10(limits = c(lo_min, hi_max)) +
    ggplot2::labs(x = NULL, y = "Hazard Ratio (log scale)")
}

####################################################

.check_duplicate_samples <- function(count_mat, clin, sample_col = "Tumor_Sample_Barcode", context = "bulkRNA-seq") {
  # what DESeq2 will use as rownames
  ids_from_counts <- colnames(count_mat)
  ids_from_clin   <- as.character(clin[[sample_col]])

  dup_counts <- unique(ids_from_counts[duplicated(ids_from_counts)])
  dup_clin   <- unique(ids_from_clin[duplicated(ids_from_clin)])

  if (length(dup_counts) || length(dup_clin)) {
    # Build a short, user-friendly message
    list_preview <- function(x) {
      if (length(x) == 0) return("None")
      if (length(x) <= 10) paste(x, collapse = ", ")
      else paste(paste(x[1:10], collapse = ", "), sprintf("… (+%d more)", length(x) - 10))
    }
    msg <- paste0(
      "Differential gene expression analysis failed because at least one same sample detected in both cohorts.<br>",
      "Please resolve and try again."
    )
    showNotification(
      htmltools::HTML(msg),
      type = "error",
      duration = 10
    )
    validate(need(FALSE, msg))
  }
}

# Function to process DESeq2 -----
process_deseq2 <- function(filtered_data, bulkseq, min_counts=5, min_samples=5) {
  clinical_combined <- filtered_data$combined
  clinical_combined <- clinical_combined[clinical_combined$Tumor_Sample_Barcode %in% colnames(bulkseq),]
  count_data <- bulkseq[, clinical_combined$Tumor_Sample_Barcode]
  count_data <- round(count_data)

  .check_duplicate_samples(count_data, clinical_combined, context = "bulk RNA-seq (DESeq2)")

  # Create metadata
  metadata <- data.frame(
    row.names = clinical_combined$Tumor_Sample_Barcode,
    condition = as.factor(clinical_combined$cohort)
  )

  # Filter genes
  keep_genes <- rowSums(count_data >= min_counts) >= min_samples
  count_data <- count_data[keep_genes, ]

  # DESeq2 analysis
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~ condition)
  dds <- dds[rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)
  results <- results(dds, contrast = c("condition", "Cohort1", "Cohort2"))

  results <- as.data.frame(results)
  results$significant <- ifelse(results$padj < 0.05 & results$log2FoldChange > 1.5, "Up-regulated",
                                ifelse(results$padj < 0.05 & results$log2FoldChange < -1.5, "Down-regulated", "Not Significant"))
  results <- results[!is.na(results$significant), ]

  return(results)
}

.deseq2_volcano_df <- function(deseq2_result, p_thr, logfc_thr) {
  stopifnot(all(c("log2FoldChange","padj") %in% colnames(deseq2_result)))

  deseq2_result %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(
      padj = as.numeric(padj),
      padj_safe = dplyr::if_else(!is.finite(padj) | is.na(padj) | padj == 0, 1, padj),
      neglog10padj = -log10(padj_safe),
      significant = dplyr::case_when(
        padj < p_thr & log2FoldChange >  logfc_thr ~ "Up-regulated",
        padj < p_thr & log2FoldChange < -logfc_thr ~ "Down-regulated",
        TRUE ~ "Not significant"
      )
    )
}

deseq2_volcano_ggplot <- function(deseq2_result, p_thr, logfc_thr, n_labels = 10,
                                  palette = c(
                                    "Not significant" = "#B3B3B3",
                                    "Up-regulated"    = "#D55E00",
                                    "Down-regulated"  = "#0072B2"
                                  )) {
  df <- .deseq2_volcano_df(deseq2_result, p_thr, logfc_thr)
  ythr <- -log10(p_thr)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FoldChange, y = neglog10padj, color = significant)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.6) +
    ggplot2::geom_vline(xintercept = c(-logfc_thr, logfc_thr), linetype = "dotted", color = "grey40") +
    ggplot2::geom_hline(yintercept = ythr, linetype = "dotted", color = "grey40") +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::labs(x = "log2 fold change", y = "-log10 adjusted p-value", color = "Status") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top")

  lab_df <- df %>%
    dplyr::filter(significant != "Not significant" & is.finite(padj)) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice(seq_len(min(n_labels, dplyr::n())))

  if (nrow(lab_df) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = lab_df,
      ggplot2::aes(label = gene),
      size = 3,
      max.overlaps = Inf,
      show.legend = FALSE
    )
  }

  p
}

# Function for volcano plot
deseq2_volcano_plotly <- function(deseq2_result, p_thr, logfc_thr, src = "volcano",
                                  n_labels = 10, point_alpha = 0.7, sig_size = 7,
                                  nonsig_size = 5, palette = c(
                                    "Not significant" = "#B3B3B3",  # grey70
                                    "Up-regulated"    = "#D55E00",  # vermilion
                                    "Down-regulated"  = "#0072B2"   # blue
                                  )) {
  df <- .deseq2_volcano_df(deseq2_result, p_thr, logfc_thr) %>%
    dplyr::mutate(
      point_size = dplyr::if_else(significant == "Not significant", nonsig_size, sig_size),
      tooltip = paste0(
        "<b>", gene, "</b>",
        "<br>log2FC: ", signif(log2FoldChange, 3),
        "<br>padj: ", ifelse(is.na(padj), "NA", signif(padj, 3)),
        "<br>status: ", significant
      )
    )

  # axis limits for threshold lines (symmetric x-range feels cleaner)
  xmax <- max(abs(df$log2FoldChange), na.rm = TRUE)
  ymax <- max(df$neglog10padj, na.rm = TRUE)
  ythr <- -log10(p_thr)

  p <- plotly::plot_ly(
    data = df, source = src,
    x = ~log2FoldChange, y = ~neglog10padj,
    type = "scattergl", mode = "markers",
    color = ~significant, colors = palette,
    text = ~tooltip, hoverinfo = "text",
    key = ~gene,
    marker = list(
      size = ~point_size,
      opacity = point_alpha,
      line = list(width = 0.4, color = "rgba(0,0,0,0.25)")
    )
  ) %>%
    plotly::layout(
      template = "plotly_white",
      hovermode = "closest",
      margin = list(l = 60, r = 20, t = 35, b = 55),
      font = list(family = "Helvetica, Arial, sans-serif", size = 13),
      xaxis = list(
        title = "log<sub>2</sub> fold change",
        zeroline = FALSE,
        range = c(-xmax, xmax),
        gridcolor = "rgba(0,0,0,0.07)"
      ),
      yaxis = list(
        title = "-log<sub>10</sub> adjusted p-value",
        rangemode = "tozero",
        gridcolor = "rgba(0,0,0,0.07)"
      ),
      legend = list(
        title = list(text = "<b>Status</b>"),
        orientation = "h", x = 0, y = 1.12, xanchor = "left"
      ),
      shapes = list(
        # vertical FC thresholds
        list(type = "line", x0 =  logfc_thr, x1 =  logfc_thr, xref = "x",
             y0 = 0, y1 = ymax, yref = "y",
             line = list(color = "rgba(0,0,0,0.35)", dash = "dot", width = 1.5)),
        list(type = "line", x0 = -logfc_thr, x1 = -logfc_thr, xref = "x",
             y0 = 0, y1 = ymax, yref = "y",
             line = list(color = "rgba(0,0,0,0.35)", dash = "dot", width = 1.5)),
        # horizontal p-value threshold
        list(type = "line", x0 = -xmax, x1 = xmax, xref = "x",
             y0 = ythr, y1 = ythr, yref = "y",
             line = list(color = "rgba(0,0,0,0.35)", dash = "dot", width = 1.5))
      ),
      annotations = list(
        # tiny labels on threshold lines for clarity
        list(x =  logfc_thr, y = ymax, xref = "x", yref = "y",
             text = paste0("|log2FC| ≥ ", signif(logfc_thr, 3)),
             ax = 10, ay = -20, showarrow = TRUE, arrowwidth = 1, arrowsize = 0.6,
             font = list(size = 11), bgcolor = "rgba(255,255,255,0.7)", bordercolor = "rgba(0,0,0,0.2)"),
        list(x =  0, y = ythr, xref = "x", yref = "y",
             text = paste0("padj ≤ ", format(p_thr, digits = 2, scientific = TRUE)),
             ax = 0, ay = -25, showarrow = TRUE, arrowwidth = 1, arrowsize = 0.6,
             font = list(size = 11), bgcolor = "rgba(255,255,255,0.7)", bordercolor = "rgba(0,0,0,0.2)")
      )
    ) %>%
    # cleaner hover text
    plotly::style(hoverlabel = list(bgcolor = "white", bordercolor = "rgba(0,0,0,0.2)", font = list(size = 12))) %>%
    # nicer export button
    plotly::config(
      displaylogo = FALSE,
      modeBarButtonsToRemove = c("autoScale2d","lasso2d","select2d","toggleSpikelines"),
      toImageButtonOptions = list(format = "png", filename = "volcano_plot")
    )

  # labels for top significant genes (by padj)
  lab_df <- df %>%
    dplyr::filter(significant != "Not significant" & is.finite(padj)) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice(seq_len(min(n_labels, dplyr::n()))) %>%
    dplyr::mutate(label = htmltools::htmlEscape(as.character(gene)))

  if (nrow(lab_df) > 0) {
    p <- p %>%
      plotly::add_trace(
        data = lab_df,
        x = ~log2FoldChange, y = ~neglog10padj,
        type = "scatter", mode = "text",
        text = ~label, textposition = "top center",
        textfont = list(size = 11, color = "black"),
        hoverinfo = "skip", showlegend = FALSE
      )
  }
  p
}

# Function to create distribution density plot for gene of interest
tpm_distr_dens <- function(count_data_tpm, clinical_combined, gene_interested, data_type) {
  # Merge the datasets
  value_type <- "TPM"
  if (data_type=="scRNAseq") {
    value_type <- "Log-normalized Expression"
  } else {
    if (data_type=="bulkRNAseq") {
      value_type <- "TPM"
    }
  }

  merged_data <- count_data_tpm %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "Tumor_Sample_Barcode") %>%
    inner_join(clinical_combined, by = "Tumor_Sample_Barcode")

  # Ensure 'cohort' is a factor
  merged_data$cohort <- as.factor(merged_data$cohort)

  # Filter for the gene of interest
  gene_data <- merged_data %>%
    dplyr::select(Tumor_Sample_Barcode, cohort, all_of(gene_interested)) %>%
    rename(TPM = all_of(gene_interested))

  # Calculate median TPM for each cohort
  median_tpm <- gene_data %>%
    group_by(cohort) %>%
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
  p <- ggplot(gene_data, aes(x = TPM, fill = cohort, color = cohort)) +
    geom_histogram(aes(y = after_stat(count),
                       text = paste(paste0(value_type, " range:"), round(after_stat(x) - after_stat(width)/2, 2), "-", round(after_stat(x) + after_stat(width)/2, 2), "<br>Count:", after_stat(count))),
                   position = "identity", bins = num_bins, fill = NA, alpha = 0) +
    geom_density(aes(y = after_stat(density) * scale_factor), alpha = 0.3, adjust = 1.5) +
    geom_vline(data = median_tpm, aes(xintercept = median_TPM, color = cohort, text = paste(paste0("Median ", value_type), round(median_TPM, 2))), linetype = "dashed", linewidth = 1) +
    labs(x = value_type,
         y = "Count",
         fill = "Cohort",
         color = "Cohort") +
    theme_minimal()

  return(p)
}

# Function to draw boxplot of tpm values and use Wilcoxon test for p values
tpm_boxplot <- function(count_data_tpm, clinical_combined, gene_interested, data_type) {

  value_label <- if (data_type == "scRNAseq") "Log-normalized Expression" else "TPM"
  value_col   <- value_label

  y_mapping <- if (make.names(value_col) != value_col) paste0("`", value_col, "`") else value_col

  merged_data <- count_data_tpm %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "Tumor_Sample_Barcode") %>%
    inner_join(clinical_combined, by = "Tumor_Sample_Barcode")

  merged_data$cohort <- as.factor(merged_data$cohort)

  gene_data <- merged_data %>%
    dplyr::select(Tumor_Sample_Barcode, cohort, all_of(gene_interested)) %>%
    dplyr::rename(!!value_col := all_of(gene_interested))

  max_y <- suppressWarnings(max(gene_data[[value_col]], na.rm = TRUE))

  p <- ggpubr::ggboxplot(gene_data, x = "cohort", y = y_mapping,
                         color = "cohort", add = "jitter") +
    labs(x = "Cohort",
         y = paste0(value_label, " of ", gene_interested))

  if (is.finite(max_y)) {
    p <- p + ggpubr::stat_compare_means(method = "wilcox.test",
                                        label = "p.format",
                                        label.y = max_y * 1.1)
  }
  p
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

tpm_distr_survival <- function(gene_tpm, selected_clinical, cohorting_method) {
  gene_interested <- rownames(gene_tpm)
  all_samples <- colnames(gene_tpm)
  gene_tpm <- as.numeric(gene_tpm)

  if (cohorting_method == "quartiles") {
    # Cohort by quartiles
    quartiles <- quantile(gene_tpm, probs = c(0, 0.25, 0.5, 0.75, 1))

    # Prevent 'breaks' not unique error
    epsilon <- .Machine$double.eps
    quartiles <- quartiles + cumsum(duplicated(quartiles)) * epsilon

    # Create cohort
    cohort <- cut(gene_tpm,
                  breaks = quartiles,
                  labels = c("Q1", "Q2", "Q3", "Q4"),
                  include.lowest = TRUE)

  } else if (cohorting_method == "median") {
    # Cohort by median
    median_value <- median(gene_tpm, na.rm = TRUE)
    cohort <- ifelse(gene_tpm <= median_value, "Below Median", "Above Median")
  }

  # Map cohort to clinical data
  cohort_df <- data.frame(Tumor_Sample_Barcode = all_samples,
                          cohort = cohort)

  selected_clinical <- merge(selected_clinical, cohort_df, by = "Tumor_Sample_Barcode")

  # Some samples does not have PFS and PFS_event, so the num of samples used in survival curve and selected_clinical is not consistent
  # Num of NAs: PFS, PFS_event; 308, 12.
  # Survival curves
  surv_object <- Surv(time = selected_clinical$PFS, event = selected_clinical$PFS_event)
  fit <- do.call(survfit, list(surv_object ~ cohort, data = selected_clinical))
  ggsurvplot(fit, data = selected_clinical, pval = TRUE,
             risk.table = TRUE, risk.table.col = "strata",
             ggtheme = theme_minimal(),
             tables.theme = theme_void(),
             xlab = "Days",
             ylab = "Progression-Free Survival",
             title = paste("Survival Curve by", gene_interested, "Expression by", if (cohorting_method == "quartiles") "Quartiles" else "median"))
}

compute_significant_gene_sets <- function(ssgsea_result, clinical_combined) {
  # Intersect
  inter_samples <- intersect(colnames(ssgsea_result), clinical_combined$Tumor_Sample_Barcode)
  clinical_combined <- clinical_combined[clinical_combined$Tumor_Sample_Barcode %in% inter_samples, c("Tumor_Sample_Barcode", "cohort")]

  req(length(unique(clinical_combined$cohort)) == 2)

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
    cols = -c(Tumor_Sample_Barcode, cohort),
    names_to = "GeneSet",
    values_to = "EnrichmentScore"
  )
  colnames(long_data)[1:2] <- c("Sample", "Cohort")


  unique_gene_sets <- unique(long_data$GeneSet)
  wilcox_results <- data.frame(GeneSet = character(), p_value = numeric(), stringsAsFactors = FALSE)

  for (gene_set in unique_gene_sets) {
    subset_data <- filter(long_data, GeneSet == gene_set)
    wilcox_result <- wilcox.test(EnrichmentScore ~ Cohort, data = subset_data)
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

    cohort1_mean <- mean(as.numeric(unlist(long_data[long_data$GeneSet == gene_set & long_data$Cohort == 'Cohort1', 'EnrichmentScore'])), na.rm = TRUE)
    cohort2_mean <- mean(as.numeric(unlist(long_data[long_data$GeneSet == gene_set & long_data$Cohort == 'Cohort2', 'EnrichmentScore'])), na.rm = TRUE)

    if (cohort1_mean > cohort2_mean) {
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
  p <- ggplot(top_significant_data, aes(x = GeneSet, y = EnrichmentScore, fill = Cohort)) +
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

celltype_boxplot <- function(cohort_info, sc_meta) {
  # Calculate cell type abundance for each patient
  sc_meta <- sc_meta[, c("public_id", "celltypes", "barcode")]
  patient_abundance <- sc_meta %>%
    group_by(public_id, celltypes) %>%
    summarise(count = n()) %>%
    mutate(proportion = count / sum(count))

  # Merge patient abundance with clinical data
  merged_data <- merge(patient_abundance, cohort_info, by = "public_id")

  req(length(unique(merged_data$cohort)) == 2)
  # Compare cell type abundance between cohorts
  comparison_results <- merged_data %>%
    group_by(celltypes) %>%
    summarise(
      p_value = wilcox.test(proportion ~ cohort)$p.value
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
  ggplot(merged_data, aes(x = celltypes, y = proportion, fill = cohort)) +
    geom_boxplot() +
    labs(title = "Cell Type Abundance by Cohort",
         x = "Cell Type",
         y = "Proportion",
         fill = "Cohort") +
    theme_minimal() +
    geom_text(aes(x = celltypes, y = max(proportion) + 0.1, label = Significance), size = 5, vjust = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

celltype_proportion <- function(cohort_info, sc_meta) {
  # Calculate cell type abundance by cohort
  sc_meta <- sc_meta[, c("public_id", "celltypes", "barcode")]
  merged_meta <- merge(sc_meta, cohort_info, by = "public_id")
  cell_proportions <- merged_meta %>%
    group_by(cohort, celltypes) %>%
    summarise(count = n()) %>%
    mutate(proportion = count / sum(count) * 100)

  colors <- c(
    "CD4_T_cells" = "#7b3294",
    "CD8_T_cells" = "#c2a5cf",
    "B_cells" = "#a6dba0",
    "Monocytes" = "#008837",
    "NK_cells" = "#fdae61",
    "Plasma_cells" = "#e66101"
  )
  cohort_levels <- rev(unique(as.character(cell_proportions$cohort)))
  cell_proportions$cohort <- factor(cell_proportions$cohort, levels = cohort_levels)
  # Reverse the order of factors for correct stacking
  cell_proportions$celltypes <- factor(cell_proportions$celltypes, levels = rev(levels(factor(cell_proportions$celltypes))))
  ggplot(cell_proportions, aes(y = cohort, x = proportion, fill = celltypes)) +
    geom_bar(stat = "identity", color = "black") + # Add border to the bars
    geom_label(aes(label = sprintf("%s\n%.2f%%", celltypes, proportion)),
               position = position_stack(vjust = 0.5, reverse = TRUE), size = 3.5,
               fill = "white", # Background for the label
               label.padding = unit(0.15, "lines"),
               color = "black", # Set text color to black for clarity
               fontface = "bold") +
    scale_fill_manual(values = colors) +
    labs(title = "Cell Proportions by Cohort", y = "Cohort", x = "Percentage (%)") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )
}

cell_cycle_hist <- function(cohort_info, sc_meta, celltypes) {
  if (!('All' %in% celltypes)) {
    sc_meta <- sc_meta[sc_meta$celltypes %in% celltypes, ]
  }

  if (length(celltypes) == 0){
    showNotification("No cell types selected. Please select at least one.", type = "error")
    return(NULL)
  }

  sc_meta <- merge(sc_meta, cohort_info, by = "public_id")

  phase_cohort_counts <- table(sc_meta$Phase, sc_meta$cohort)
  phase_cohort_df <- as.data.frame(phase_cohort_counts)
  colnames(phase_cohort_df) <- c("Phase", "Cohort", "Count")

  # Calculate the percentage of cells in each phase within each cohort
  phase_cohort_summary <- phase_cohort_df %>%
    group_by(Cohort) %>%
    mutate(Total = sum(Count),
           Percentage = (Count / Total) * 100)

  # Create the histogram plot by cohort
  ggplot(phase_cohort_summary, aes(x = Phase, y = Percentage, fill = Cohort)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(title = "Cell Cycle Phase Distribution by Cohort (Percentage)",
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
      summarise(count = n(), .cohorts = 'drop') %>%
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


# --- PSEUDOBULK: load + helpers ---------------------------------------------
# Convert IDs like "MMRF_1013_1" → "MMRF_1013" (matches pseudobulk columns)
.to_base_id <- function(x) sub("(_[0-9]+)?$", "", as.character(x))

# Make sure a pseudobulk matrix is numeric and well‑formed
.get_pb_matrix <- function(pb_list, celltype) {
  if (is.null(pb_list) || !is.list(pb_list)) stop("Invalid pseudobulk object")
  if (!(celltype %in% names(pb_list))) stop(sprintf("Cell type '%s' not found.", celltype))
  mat <- pb_list[[celltype]]
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  mat
}

# Align a pseudobulk TPM matrix to clinical samples by patient (base_id)
# Returns: list(pb_tpm, clinical_aligned, inter_base_ids)
.align_pb_to_clinical <- function(pb_tpm, clinical_combined) {
  clinical_combined$base_id <- .to_base_id(clinical_combined$public_id)
  inter <- intersect(colnames(pb_tpm), clinical_combined$base_id)
  if (length(inter) < 2) {
    warning("Very small overlap between pseudobulk and clinical")
  }

  # Keep one clinical row per base_id
  clin_one <- clinical_combined[clinical_combined$base_id %in% inter, , drop = FALSE]

  # Reorder columns to the clinical order and relabel columns to Tumor_Sample_Barcode
  pb_tpm2 <- pb_tpm[, inter, drop = FALSE]
  pb_tpm2 <- pb_tpm2[, clin_one$base_id, drop = FALSE]
  colnames(pb_tpm2) <- clin_one$Tumor_Sample_Barcode

  list(pb_tpm = pb_tpm2, clinical_aligned = clin_one, inter_base_ids = inter)
}

pseudobulk_diff_counts <- function(pb_counts, clin) {
  stopifnot(ncol(pb_counts) == nrow(clin))
  grp <- factor(clin$cohort)
  if (length(levels(grp)) != 2) stop("Need exactly two cohorts for differential analysis")

  # Round to integers just in case
  cnt <- round(pb_counts)
  keep <- rowSums(cnt) > 1
  cnt <- cnt[keep, , drop = FALSE]

  coldata <- data.frame(row.names = colnames(cnt), condition = grp)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = cnt, colData = coldata, design = ~ condition)
  dds <- DESeq2::DESeq(dds)
  res <- as.data.frame(DESeq2::results(dds, contrast = c("condition", levels(grp)[1], levels(grp)[2])))
  res$baseMean <- if (!"baseMean" %in% names(res)) rowMeans(cnt) else res$baseMean
  res <- res[, intersect(c("baseMean","log2FoldChange","pvalue","padj"), names(res))]
  res
}
