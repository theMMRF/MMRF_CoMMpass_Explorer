setwd("~/CMU/MMRF/RShiny/MMRF_CoMMpass_Explorer")
source("global.R")

# Load CID meta data
cid_meta <- read_xlsx("data/CID_metadata_all_may2024.xlsx")

# scTCR-seq
cid_meta_tcr <- cid_meta[cid_meta$data_type=="scTCR-seq", ]
tcr_samples <- unique(cid_meta_tcr$sample_id)

# scRNA-seq
cid_meta_scrna <- cid_meta[cid_meta$data_type=="scRNA-seq", ]
cid_meta_scrna <- cid_meta_scrna[cid_meta_scrna$dataset=="CID-005", ]
scrna_samples <- unique(cid_meta_scrna$sample_id)

# Sample intersect
length(intersect(tcr_samples, scrna_samples))

# cid_meta_bl <- cid_meta[grep("_1$", cid_meta$sample_id), ]
# cid_meta_bl <- cid_meta_bl[cid_meta_bl["data_type"]=="scRNA-seq",] # Baseline, scRNA-seq

# Table for unique values in each column
create_unique_count_table <- function(cid_meta_data) {
  unique_counts <- lapply(cid_meta_data[,-c(1:5)], function(column) {
    table(column)
  })
  return(unique_counts)
}

unique_count <- create_unique_count_table(cid_meta)
unique_count

dim(unique(cid_meta_bl['sample_id'])) # num of unique sample ids

################
#MAF subsetting#
################
clinical_data <- readRDS("data/clinical_data_cleaned.rds")
maf_data <- readRDS("data/maf_data.rds")

# Subset based on the genes interested
genes_interested <- c("KRAS")

maf_data_table <- maf_data@data
mutated_patients <- maf_data_table[Hugo_Symbol %in% genes_interested, .(Tumor_Sample_Barcode)]
mutated_patients <- unique(mutated_patients)


library(maftools)

# Get full sample IDs
sample_ids <- maf_data@clinical.data$Tumor_Sample_Barcode

# Extract patient IDs (first 10 characters)
patient_ids <- substr(sample_ids, 1, 9)

# Create a mapping of full sample ID to patient ID
barcode_df <- data.frame(
  full_id = sample_ids,
  patient_id = patient_ids,
  stringsAsFactors = FALSE
)
# Filter for excellent outcome barcodes
excellent_barcodes <- barcode_df$full_id[barcode_df$patient_id %in% excellent_outcome_ids]

# Same for poor outcome
poor_barcodes <- barcode_df$full_id[barcode_df$patient_id %in% poor_outcome_ids]

excellent_maf <- subsetMaf(maf = maf_data, tsb = excellent_barcodes)
poor_maf <- subsetMaf(maf = maf_data, tsb = poor_barcodes)

# Summarize cohort sizes
cat("Excellent outcome cohort size:", length(excellent_barcodes), "\n")
cat("Poor outcome cohort size:", length(poor_barcodes), "\n")

# Basic mutation summaries for each group
excellent_summary <- getSampleSummary(excellent_maf)
poor_summary <- getSampleSummary(poor_maf)

# Compare average mutations per sample
mean_muts_excellent <- mean(excellent_summary$total)
mean_muts_poor <- mean(poor_summary$total)
cat("Mean mutations per sample (Excellent):", mean_muts_excellent, "\n")
cat("Mean mutations per sample (Poor):", mean_muts_poor, "\n")

# Mutation types distribution
excellent_types <- getGeneSummary(excellent_maf)
poor_types <- getGeneSummary(poor_maf)

# Compare mutation frequencies between groups
maf_comparison <- mafCompare(m1 = excellent_maf, m2 = poor_maf, 
                             m1Name = "Excellent Outcome", 
                             m2Name = "Poor Outcome")

# Get the significant genes
sig_genes <- maf_comparison$results[maf_comparison$results$adjPval < 0.05,]
print(sig_genes)

# Create forest plot of significantly different mutations
forestPlot(mafCompareRes = maf_comparison, pVal = 0.05)

# Oncoplots for top mutated genes in each group
oncoplot(maf = excellent_maf, top = 20, titleText = "Excellent Outcome")
oncoplot(maf = poor_maf, top = 20, titleText = "Poor Outcome")

# Oncoplot with genes sorted by statistical significance
sig_genes_list <- sig_genes$Hugo_Symbol
oncoplot(maf = maf_data, genes = sig_genes_list, 
         clinicalFeatures = "Outcome_Group", 
         sortByAnnotation = TRUE,
         titleText = "Significant Genes by Outcome")


#############
# scRNA-seq #
#############
# scRNA-seq Meta data
library("ggplot2")
# sc_meta <- read.table("data/IA_baseline_meta_data.txt")
# sc_meta$barcode <- rownames(sc_meta)
# saveRDS(sc_meta, "scRNAseq_metadata.rds")
readRDS("data/scRNAseq_metadata.rds")


# get the counts
phase_counts <- table(sc_meta$Phase)
phase_df <- as.data.frame(phase_counts)
colnames(phase_df) <- c("Phase", "Count")

# calculate percentage
total_cells <- sum(phase_df$Count)
phase_df$Percentage <- (phase_df$Count / total_cells) * 100

ggplot(phase_df, aes(x = Phase, y = Percentage, fill = Phase)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Cell Cycle Phase Distribution (Percentage)",
       x = "Cell Cycle Phase",
       y = "Percentage of Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Clean and standardize Tx column in clinical data
cleaned_Tx <- sapply(clinical_data$Tx, function(x) {
  if (is.na(x)) {
    return(NA)
  } else {
    paste(sort(unlist(strsplit(x, "/"))), collapse = "/")
  }
})

unique(cleaned_Tx)

# Focus on IMID, Proteasome, and Steroid
simplified_Tx <- sapply(cleaned_Tx, function(x) {
  if (is.na(x)) return(NA)
  regimes <- unlist(strsplit(x, "/"))
  relevant <- regimes[regimes %in% c("IMID", "Proteasome", "Steroid")]
  if (length(relevant) > 0) paste(sort(unique(relevant)), collapse = "/") else NA
})

simplified_unique <- unique(simplified_Tx)

simplified_unique


barplot(table(clinical_data$Tx),
        main = "Distribution of Tx",
        xlab = "Treatment Regimes (Tx)",
        ylab = "Frequency",
        las = 2,                # Rotate x-axis labels for better readability
        cex.names = 0.7)

#########
## Orig
#########
library(maftools)
library(sigminer)
library(NMF)

clinical <- read.csv("orig_data/IA_clinical_data.txt", sep = "\t")[,-1]

maf <- read.table("orig_data/commpass_multi_omic_summary_flat_files_MMRF_CoMMpass_IA22_exome_vcfmerger2_All_Canonical_Variants.tsv", 
                  header = TRUE,
                  sep = "\t",
                  stringsAsFactors = FALSE,
                  check.names = FALSE)

dim(maf)

#########
## Clinical
#########
patients_in_clinical <- data.frame(sapply(strsplit(as.character(clinical_data$Tumor_Sample_Barcode), "_"),
                                   function(x) paste(x[1], x[2], sep="_")))

patients_in_clinical['visit'] <- sapply(strsplit(as.character(clinical_data$Tumor_Sample_Barcode), "_"),
                                        function(x) x[3])

patients_in_clinical['type'] <- sapply(strsplit(as.character(clinical_data$Tumor_Sample_Barcode), "_"),
                                       function(x) x[4])

unique(patients_in_clinical['type'])

clin_baseline <- read.csv("data/clinical_data_cleaned_new_surv.csv")
saveRDS(clin_baseline, "data/clinical_data_cleaned_new_surv.rds")

# For survival


#########
## Original clinical data for WGS and WES
#########
IA_clinical_data <- read.csv("orig_clinical/IA_clinical_data.txt", sep="\t")[-1]
NDMM_clinical <- read.csv("orig_clinical/CoMMpass_NDMM_clinical.txt", sep="\t")[-1]
clinical_BM_NDMM <- readRDS("orig_clinical/CoMMpass_clinical_BM_NDMM.rds")
NDMM_BM_clinical <- readRDS("orig_clinical/NDMM_BM_clinical.rds")


##########
### Mutation subsetting
##########
maf <- maf_data@data
include_genes <- c('KRAS', 'NLRP10', 'NRAS')

# or
uniqueN(maf[Hugo_Symbol %in% include_genes, Tumor_Sample_Barcode])

# And
# length(intersect(unique(maf[Hugo_Symbol == 'KRAS', Tumor_Sample_Barcode]), unique(maf[Hugo_Symbol == 'NRAS', Tumor_Sample_Barcode])))
sample_gene_count <- maf %>%
  filter(Hugo_Symbol %in% include_genes) %>%
  distinct(Tumor_Sample_Barcode, Hugo_Symbol) %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(gene_count = n(), .groups = "drop")

samples_with_all_genes <- sample_gene_count %>%
  filter(gene_count == length(include_genes)) %>%
  pull(Tumor_Sample_Barcode)

length(samples_with_all_genes)
tmp <- clinical_data[clinical_data$Tumor_Sample_Barcode %in% samples_with_all_genes]
dim(tmp)[1]


# Interaction Plot
somaticInteractions(maf = maf_data, top = 25, genes = c("ACKR1", "KRAS", "NRAS", "BRAF", "PRKD2"), pvalue = c(0.05, 0.1))
gene_summary <- getGeneSummary(maf_data)
gene_summary <- gene_summary[order(-gene_summary$MutatedSamples), ]
head(gene_summary$Hugo_Symbol, 25)

somaticInteractions(maf = clinical_data, top = 25, genes = c("ACKR1", "KRAS", "NRAS", "BRAF", "PRKD2"), pvalue = c(0.05, 0.1))


#### Convert them into csv
bulkseq <- readRDS("data/bulkseq_baseline.rds")
bulkseq_tpm <- readRDS("data/bulkseq_tpm_baseline.rds")
clinical_data <- readRDS("data/clinical_data_cleaned.rds")
maf_data <- readRDS("data/maf_data.rds")
sc_meta <- readRDS("data/scRNAseq_metadata.rds")
ssgsea_result_ca <- readRDS("data/ssgsea_result_ca.rds")
ssgsea_result_c2 <- readRDS("data/ssgsea_result_c2.rds")

write.csv(clinical_data, "clinical_data_cleaned.csv", quote = F, row.names = F)
write.csv(bulkseq_tpm, "bulkseq_tpm_baseline.csv", quote = F)
write.csv(bulkseq, "bulkseq_baseline.csv", quote = F)
write.csv(sc_meta, "scRNAseq_metadata.csv", quote = F)
write.csv(maf_data@data, "maf_data.csv", quote = F, row.names = F)


#####
# Clinical processing
#####
colnames(clinical_data)[20] <- "TP53_inactivation"
colnames(clinical_data)[22] <- "MAF_MAFB"
colnames(clinical_data)[23] <- "t_11_14"
colnames(clinical_data)[24] <- "t_4_14"

colnames(clinical_data)[15] <- "chr_1q21_amp"
colnames(clinical_data)[16] <- "chr_13q14_del"
colnames(clinical_data)[17] <- "chr_13q34_del"
colnames(clinical_data)[18] <- "chr_17p13_del"
colnames(clinical_data)[21] <- "chr_1q21_gain"


saveRDS(clinical_data, "data/clinical_data_cleaned_new_surv.rds")

#####
# For new clinical data (n=1141), firstline, BM
#####
clin <- read.csv("data/CoMMpass_baseline_clinical_n1141.txt", sep="\t")
# Recode from numeric to categorical
clin$Gender <- factor(clin$Gender, 
                      levels = c(1, 2),
                      labels = c("Male", "Female"))

clin$Ethnicity <- factor(clin$Ethnicity, 
                         levels = c(1, 2, 3, 4, 6),
                         labels = c("Caucasian", "Black", "American_Indian", "Asian", "Other_Unknown"))

# Reorganize clinical data columns by category
clin <- clin[, c(
  # Identifiers
  "public_id", "Visit_ID", "Specimen_ID", 
  
  # Demographics
  "Gender", "Ethnicity", "Age", "BMI", 
  
  # Clinical Params
  "ISS_Stage", "IMWG_Risk_Class", "Clinical_IgH", "Clinical_IgL", "Clinical_Ig_Status", 
  
  # Lab Values
  "Serum_Albumin", "Serum_B2M", "Serum_LDH", "Creatinine", 
  
  # Molecular Characteristics
  "Cytogenetic_High_Risk", "RNA_Subtype_Name", "CNA_Subtype_Name", 
  "TP53_Funct_Copies", "TP53_NS_Mut_Count", "Hyperdiploid_Call", 
  
  # Treatment Related
  "ASCT_First", "Trip_First", "Duration_First_Resp", "Tx", 
  
  # Outcomes
  "PD", "PFS", "PFS_censored", "PFS_event", "PFS_1", "PFS_1_censored", 
  "PFS_1_event", "OS", "OS_censored", "OS_event"
)]

# Check the structure of the recoded data
str(clin[, c("Gender", "Ethnicity")])

table(clin$Gender, useNA = "always")
table(clin$Ethnicity, useNA = "always")

# Old clinical data (n=885)
clin_old <- readRDS("data/clinical_data_cleaned.rds")
table(clin_old$Sex, useNA = "always")
table(clin_old$Race, useNA = "always")

# Generate Tumor_Sample_Barcode, Skerget risk, Sex->Gender, Generate Age Range,
# ISS->ISS_Stage, ASCT->ASCT_Fistline, Hyperdiploidy->Hyperdiploid_Call,
# Trip_Firstline -> Trip_First,
# "1q21_amp", "13q14_del", "13q34_del", "17p13_del", "1q21_gain", "t(11;14)", "t(4;14)"
# Chromothripsis, "APOBEC", "MAF/MAFB"

# Add columns to new clinical data






