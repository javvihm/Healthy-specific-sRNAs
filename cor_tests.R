# ===================
# INSTALL AND LOAD PACKAGES
# ===================
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

# Libraries
library(tidyverse)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)
library(stringr) 

# Define some variables
project <- 'LUAD'
cohort <- 'CPTAC'

# ===================
# GET GE DATA FROM CPTAC
# FILTER ONLY FOR LUAD-HEALTHY SAMPLES
# ===================

# Query with all CPTAC-3 data
query <- GDCquery(
  project = "CPTAC-3",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

saveRDS(query, file= 'rds_cptac3_query.rds')

# Get the sample sheet
sample_sheet <- getResults(query)


getManifest(query,save = T)
system(paste0("mv gdc_manifest.txt manifest_file_",cohort,"_",project,".txt"))

clinical <- GDCquery_clinic(project = "CPTAC-3", type = "clinical")
clinical <- subset(clinical, select=which(!duplicated(colnames(clinical)))) 

#LUAD
if (project == "LUAD") {
  tissue <- clinical |> 
    filter(primary_site == "Bronchus and lung") |> 
    filter(disease_type == "Adenomas and Adenocarcinomas")
}

# Filter the sample sheet
sample_sheet_filt <- sample_sheet %>%
  filter(cases.submitter_id %in% tissue$submitter_id) %>%
  filter(sample_type == 'Solid Tissue Normal')

# Filter the query
query_filt <- query
query_filt$results[[1]] <- query$results[[1]] %>%
  filter(sample.submitter_id %in% sample_sheet_filt$sample.submitter_id)

# Prepare the data
GDCdownload(query_filt)
data <- GDCprepare(query_filt)

# =================
# GET MIRNA AND LINC DATA FOR COR TEST
# =================

# Get the GE data and extract the LINC counts
count_matrix <- assay(data)
linc <- count_matrix[ 'ENSG00000259471.1', , drop = FALSE]
linc <- data.frame(submitter_id = colnames(linc),
                     linc_counts = linc[,])

# So far we got the x-data, now we have to import the y-data
mirna <- read.delim("./Analysis/AGGGCTGAACCTAGAATT_greps_cptac.txt",
                    header = FALSE)
colnames(mirna) <- c('Counts', 'Sequence', 'id')
# Some of the counts came from neigbouring kmers, so we add them
mirna <- mirna %>%
  group_by(id) %>%
  summarise(Counts = sum(Counts))

# We need to import the metadata and play around to get submitter_id from filenames
long_metadata <- read.delim("./metadata_CPTAC_LUAD/clean_metadata_long_CPTAC_LUAD.txt",
                            header = TRUE, sep = " ")
# Get the correspondent id for each file name
tmp <- long_metadata %>%
  select(file_name, sample.submitter_id)
colnames(tmp) <- c('id', 'submitter_id')

tmp <- tmp %>%
  mutate(id = sub("_.*", "", id))

# Add the submitter_id so they match with x-data
mirna <- mirna %>%
  left_join(tmp, by = 'id') %>%
  mutate(submitter_id = str_sub(submitter_id, 1, 12))
rm(tmp)
rm(long_metadata)

# Join the data
df <- linc %>%
  left_join(mirna, by= 'submitter_id')
colnames(df)[colnames(df) == "Counts"] <- "mirna_counts"
# Replace NA by 0
df$mirna_counts[is.na(df$mirna_counts)] <- 0

# =================
# COR TEST LINC VS MIRNA
# EXTEND TO ALL GENES
# =================
# Spearman cor test
res <- cor.test(df$linc_counts, df$mirna_counts, method = 'spearman')
print(res)

# Extension: all genes

# Sanity check: order of samples
if (all(colnames(count_matrix) == df$submitter_id)){
  print('You may continue')
}

# Correlation tests!
cor_results <- apply(count_matrix, 1, function(x) {
  test <- cor.test(x,df$mirna_counts, method = "spearman")
  c(rho = test$estimate, p = test$p.value)
})
cor_results <- as.data.frame(t(cor_results))

cor_results$adj.p <- p.adjust(cor_results$p, method = 'BH')
cor_results$gene_id <- colnames(count_matrix)

# Save the results
write_delim(x = cor_results,
            file = 'correlation_test_cptac.tsv',
            delim = "\t")
# =================
# PRECURSOR ANALYSIS
# =================

scatter_plot <- ggplot(df, aes(x = linc_counts, y = mirna_counts)) +
  geom_point(alpha = 0.6, color = "#1f77b4") +
  geom_smooth(method = "lm", color = "darkred", se = FALSE, linetype = "dashed") + # línea recta guía
  labs(
    title = "Spearman correlation between lincRNA and miRNA counts",
    subtitle = paste0("ρ = ", rho, ", p = ", pval),
    x = "lincRNA counts",
    y = "miRNA counts"
  ) +
  theme_minimal(base_size = 13)

ggsave(plot = scatter_plot, path = "./Plots", 
       filename = "spearman_rho_novel_miRNA.png",
       dpi = 600, device = "png", width = 8, height = 6)

# =================
# TARGET ANALYSIS
# =================
target <- count_matrix[ 'ENSG00000182158.15', , drop = FALSE]
target <- data.frame(submitter_id = colnames(target),
                   target = target[,])

df <- df %>%
  left_join(target, by= 'submitter_id')

res <- cor.test(df$mirna_counts, df$target, method = 'spearman')
rho <- round(res$estimate, 3)
pval <- signif(res$p.value, 3)

scatter_plot <- ggplot(df, aes(x = mirna_counts, y = target)) +
  geom_point(alpha = 0.6, color = "#1f77b4") +
  geom_smooth(method = "lm", color = "darkred", se = FALSE, linetype = "dashed") + 
  labs(
    title = "Spearman correlation between miRNA and CREB3L2 counts",
    subtitle = paste0("ρ = ", rho, ", p = ", pval),
    x = "miRNA counts",
    y = "CREB3L2 counts"
  ) +
  theme_minimal(base_size = 13)

ggsave(plot = scatter_plot, path = "./Plots", 
       filename = "spearman_rho_novel_miRNA_target.png",
       dpi = 600, device = "png", width = 8, height = 6)


