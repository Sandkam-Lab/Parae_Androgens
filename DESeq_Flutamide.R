# Set working directory
setwd("")

# Load libraries
{
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(heatmap3)
library(RColorBrewer)
library(ggalt)
library(VennDiagram)
library(ComplexUpset)
library(eulerr)
library(reshape2)
library(dplyr)
library(broom)
}

# Clear and reload data
{
rm(list = ls())
graphics.off()

# Set global parameters
coveragethrd <- 0.0  # Min average number of reads per sample (set to 0 to disable filtering)
Wald_FDR <- 0.05     # Max FDR value for significance for Wald tests

# Load the counts table
df <- read.table("Combined_Counts.tsv", header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE, strip.white = TRUE)
head(df)

# Load the sample information
dkey <- read.table("Sample_Info.tsv", header = TRUE, sep = "\t", stringsAsFactors = TRUE)
dkey
}

#######################################################################################
# Start Analyses 
#######################################################################################

################################
###   Brain - Male  - parae  ###
################################

# Subset dkey for Brain Male_parae
dkey_Brain_Male_parae <- subset(dkey, Tissue == "Brain" & Sex == "Male" & Morph == "parae")
print("Sample information for Brain Male parae:")
dkey_Brain_Male_parae

# Get sample IDs for Brain Male_parae
samples_Brain_Male_parae <- as.character(dkey_Brain_Male_parae$Sample)

# Subset the counts table for Brain Male_parae
df_Brain_Male_parae <- df[, samples_Brain_Male_parae]
print("Head of counts table for Brain Male_parae:")
head(df_Brain_Male_parae)

# Check if sample names match without modifying any names
samples_match_Brain_Male_parae <- all(sort(dkey_Brain_Male_parae$Sample) == sort(colnames(df_Brain_Male_parae)))
print("Do sample names match for Brain Male_parae?")
print(samples_match_Brain_Male_parae)

# Filter out low expression genes
columns_Brain_Male_parae <- ncol(df_Brain_Male_parae)
nrow_before_filter_Brain_Male_parae <- nrow(df_Brain_Male_parae)
df_Brain_Male_parae <- df_Brain_Male_parae[rowSums(df_Brain_Male_parae) >= (columns_Brain_Male_parae * coveragethrd), ]
nrow_after_filter_Brain_Male_parae <- nrow(df_Brain_Male_parae)
cat("Brain Male_parae - Number of genes before filtering:", nrow_before_filter_Brain_Male_parae, "\n")
cat("Brain Male_parae - Number of genes after filtering:", nrow_after_filter_Brain_Male_parae, "\n")

# Create DESeqDataSet
dds_Brain_Male_parae <- DESeqDataSetFromMatrix(countData = df_Brain_Male_parae,
                                               colData = dkey_Brain_Male_parae,
                                               design = ~ Treatment)

# Relevel Treatment to have 'veh' as reference
dds_Brain_Male_parae$Treatment <- relevel(dds_Brain_Male_parae$Treatment, ref = "veh")

# Run DESeq
dds_Brain_Male_parae <- DESeq(dds_Brain_Male_parae)

# Get results (unshrunken log2 fold changes)
res_Brain_Male_parae <- results(dds_Brain_Male_parae, alpha = Wald_FDR)
res_Brain_Male_parae <- res_Brain_Male_parae[order(res_Brain_Male_parae$padj), ]
n_significant_Brain_Male_parae <- sum(res_Brain_Male_parae$padj <= Wald_FDR, na.rm = TRUE)
print("Summary of DESeq2 results for Brain Male_parae:")
cat("Brain Male_parae - Number of significant DEGs (padj <= ", Wald_FDR, "): ", n_significant_Brain_Male_parae, "\n")
summary(res_Brain_Male_parae)

################################
###   Brain - Male  - immac  ###
################################

# Subset dkey for Brain Male_immac
dkey_Brain_Male_immac <- subset(dkey, Tissue == "Brain" & Sex == "Male" & Morph == "immac")
print("Sample information for Brain Male immac:")
dkey_Brain_Male_immac

# Get sample IDs for Brain Male_immac
samples_Brain_Male_immac <- as.character(dkey_Brain_Male_immac$Sample)

# Subset the counts table for Brain Male_immac
df_Brain_Male_immac <- df[, samples_Brain_Male_immac]
print("Head of counts table for Brain Male_immac:")
head(df_Brain_Male_immac)

# Check if sample names match without modifying any names
samples_match_Brain_Male_immac <- all(sort(dkey_Brain_Male_immac$Sample) == sort(colnames(df_Brain_Male_immac)))
print("Do sample names match for Brain Male_immac?")
print(samples_match_Brain_Male_immac)

# Filter out low expression genes
columns_Brain_Male_immac <- ncol(df_Brain_Male_immac)
nrow_before_filter_Brain_Male_immac <- nrow(df_Brain_Male_immac)
df_Brain_Male_immac <- df_Brain_Male_immac[rowSums(df_Brain_Male_immac) >= (columns_Brain_Male_immac * coveragethrd), ]
nrow_after_filter_Brain_Male_immac <- nrow(df_Brain_Male_immac)
cat("Brain Male_immac - Number of genes before filtering:", nrow_before_filter_Brain_Male_immac, "\n")
cat("Brain Male_immac - Number of genes after filtering:", nrow_after_filter_Brain_Male_immac, "\n")

# Create DESeqDataSet
dds_Brain_Male_immac <- DESeqDataSetFromMatrix(countData = df_Brain_Male_immac,
                                               colData = dkey_Brain_Male_immac,
                                               design = ~ Treatment)

# Relevel Treatment to have 'veh' as reference
dds_Brain_Male_immac$Treatment <- relevel(dds_Brain_Male_immac$Treatment, ref = "veh")

# Run DESeq
dds_Brain_Male_immac <- DESeq(dds_Brain_Male_immac)

# Get results (unshrunken log2 fold changes)
res_Brain_Male_immac <- results(dds_Brain_Male_immac, alpha = Wald_FDR)
res_Brain_Male_immac <- res_Brain_Male_immac[order(res_Brain_Male_immac$padj), ]
n_significant_Brain_Male_immac <- sum(res_Brain_Male_immac$padj <= Wald_FDR, na.rm = TRUE)
print("Summary of DESeq2 results for Brain Male_immac:")
cat("Brain Male_immac - Number of significant DEGs (padj <= ", Wald_FDR, "): ", n_significant_Brain_Male_immac, "\n")
summary(res_Brain_Male_immac)

##############################
###   Brain - Male  - mel  ###
##############################

# Subset dkey for Brain Male_mel
dkey_Brain_Male_mel <- subset(dkey, Tissue == "Brain" & Sex == "Male" & Morph == "mel")
print("Sample information for Brain Male mel:")
dkey_Brain_Male_mel

# Get sample IDs for Brain Male_mel
samples_Brain_Male_mel <- as.character(dkey_Brain_Male_mel$Sample)

# Subset the counts table for Brain Male_mel
df_Brain_Male_mel <- df[, samples_Brain_Male_mel]
print("Head of counts table for Brain Male_mel:")
head(df_Brain_Male_mel)

# Check if sample names match without modifying any names
samples_match_Brain_Male_mel <- all(sort(dkey_Brain_Male_mel$Sample) == sort(colnames(df_Brain_Male_mel)))
print("Do sample names match for Brain Male_mel?")
print(samples_match_Brain_Male_mel)

# Filter out low expression genes
columns_Brain_Male_mel <- ncol(df_Brain_Male_mel)
nrow_before_filter_Brain_Male_mel <- nrow(df_Brain_Male_mel)
df_Brain_Male_mel <- df_Brain_Male_mel[rowSums(df_Brain_Male_mel) >= (columns_Brain_Male_mel * coveragethrd), ]
nrow_after_filter_Brain_Male_mel <- nrow(df_Brain_Male_mel)
cat("Brain Male_mel - Number of genes before filtering:", nrow_before_filter_Brain_Male_mel, "\n")
cat("Brain Male_mel - Number of genes after filtering:", nrow_after_filter_Brain_Male_mel, "\n")

# Create DESeqDataSet
dds_Brain_Male_mel <- DESeqDataSetFromMatrix(countData = df_Brain_Male_mel,
                                             colData = dkey_Brain_Male_mel,
                                             design = ~ Treatment)

# Relevel Treatment to have 'veh' as reference
dds_Brain_Male_mel$Treatment <- relevel(dds_Brain_Male_mel$Treatment, ref = "veh")

# Run DESeq
dds_Brain_Male_mel <- DESeq(dds_Brain_Male_mel)

# Get results (unshrunken log2 fold changes)
res_Brain_Male_mel <- results(dds_Brain_Male_mel, alpha = Wald_FDR)
res_Brain_Male_mel <- res_Brain_Male_mel[order(res_Brain_Male_mel$padj), ]
n_significant_Brain_Male_mel <- sum(res_Brain_Male_mel$padj <= Wald_FDR, na.rm = TRUE)
print("Summary of DESeq2 results for Brain Male_mel:")
cat("Brain Male_mel - Number of significant DEGs (padj <= ", Wald_FDR, "): ", n_significant_Brain_Male_mel, "\n")
summary(res_Brain_Male_mel)

#################################
###   Muscle - Male  - parae  ###
#################################

# Subset dkey for Muscle Male_parae
dkey_Muscle_Male_parae <- subset(dkey, Tissue == "Muscle" & Sex == "Male" & Morph == "parae")
print("Sample information for Muscle Male parae:")
dkey_Muscle_Male_parae

# Get sample IDs for Muscle Male_parae
samples_Muscle_Male_parae <- as.character(dkey_Muscle_Male_parae$Sample)

# Subset the counts table for Muscle Male_parae
df_Muscle_Male_parae <- df[, samples_Muscle_Male_parae]
print("Head of counts table for Muscle Male_parae:")
head(df_Muscle_Male_parae)

# Check if sample names match without modifying any names
samples_match_Muscle_Male_parae <- all(sort(dkey_Muscle_Male_parae$Sample) == sort(colnames(df_Muscle_Male_parae)))
print("Do sample names match for Muscle Male_parae?")
print(samples_match_Muscle_Male_parae)

# Filter out low expression genes
columns_Muscle_Male_parae <- ncol(df_Muscle_Male_parae)
nrow_before_filter_Muscle_Male_parae <- nrow(df_Muscle_Male_parae)
df_Muscle_Male_parae <- df_Muscle_Male_parae[rowSums(df_Muscle_Male_parae) >= (columns_Muscle_Male_parae * coveragethrd), ]
nrow_after_filter_Muscle_Male_parae <- nrow(df_Muscle_Male_parae)
cat("Muscle Male_parae - Number of genes before filtering:", nrow_before_filter_Muscle_Male_parae, "\n")
cat("Muscle Male_parae - Number of genes after filtering:", nrow_after_filter_Muscle_Male_parae, "\n")

# Create DESeqDataSet
dds_Muscle_Male_parae <- DESeqDataSetFromMatrix(countData = df_Muscle_Male_parae,
                                                colData = dkey_Muscle_Male_parae,
                                                design = ~ Treatment)

# Relevel Treatment to have 'veh' as reference
dds_Muscle_Male_parae$Treatment <- relevel(dds_Muscle_Male_parae$Treatment, ref = "veh")

# Run DESeq
dds_Muscle_Male_parae <- DESeq(dds_Muscle_Male_parae)

# Get results (unshrunken log2 fold changes)
res_Muscle_Male_parae <- results(dds_Muscle_Male_parae, alpha = Wald_FDR)
res_Muscle_Male_parae <- res_Muscle_Male_parae[order(res_Muscle_Male_parae$padj), ]
n_significant_Muscle_Male_parae <- sum(res_Muscle_Male_parae$padj <= Wald_FDR, na.rm = TRUE)
print("Summary of DESeq2 results for Muscle Male_parae:")
cat("Muscle Male_parae - Number of significant DEGs (padj <= ", Wald_FDR, "): ", n_significant_Muscle_Male_parae, "\n")
summary(res_Muscle_Male_parae)

################################
###   Muscle - Male  - immac ###
################################

# Subset dkey for Muscle Male_immac
dkey_Muscle_Male_immac <- subset(dkey, Tissue == "Muscle" & Sex == "Male" & Morph == "immac")
print("Sample information for Muscle Male immac:")
dkey_Muscle_Male_immac

# Get sample IDs for Muscle Male_immac
samples_Muscle_Male_immac <- as.character(dkey_Muscle_Male_immac$Sample)

# Subset the counts table for Muscle Male_immac
df_Muscle_Male_immac <- df[, samples_Muscle_Male_immac]
print("Head of counts table for Muscle Male_immac:")
head(df_Muscle_Male_immac)

# Check if sample names match without modifying any names
samples_match_Muscle_Male_immac <- all(sort(dkey_Muscle_Male_immac$Sample) == sort(colnames(df_Muscle_Male_immac)))
print("Do sample names match for Muscle Male_immac?")
print(samples_match_Muscle_Male_immac)

# Filter out low expression genes
columns_Muscle_Male_immac <- ncol(df_Muscle_Male_immac)
nrow_before_filter_Muscle_Male_immac <- nrow(df_Muscle_Male_immac)
df_Muscle_Male_immac <- df_Muscle_Male_immac[rowSums(df_Muscle_Male_immac) >= (columns_Muscle_Male_immac * coveragethrd), ]
nrow_after_filter_Muscle_Male_immac <- nrow(df_Muscle_Male_immac)
cat("Muscle Male_immac - Number of genes before filtering:", nrow_before_filter_Muscle_Male_immac, "\n")
cat("Muscle Male_immac - Number of genes after filtering:", nrow_after_filter_Muscle_Male_immac, "\n")

# Create DESeqDataSet
dds_Muscle_Male_immac <- DESeqDataSetFromMatrix(countData = df_Muscle_Male_immac,
                                                colData = dkey_Muscle_Male_immac,
                                                design = ~ Treatment)

# Relevel Treatment to have 'veh' as reference
dds_Muscle_Male_immac$Treatment <- relevel(dds_Muscle_Male_immac$Treatment, ref = "veh")

# Run DESeq
dds_Muscle_Male_immac <- DESeq(dds_Muscle_Male_immac)

# Get results (unshrunken log2 fold changes)
res_Muscle_Male_immac <- results(dds_Muscle_Male_immac, alpha = Wald_FDR)
res_Muscle_Male_immac <- res_Muscle_Male_immac[order(res_Muscle_Male_immac$padj), ]
n_significant_Muscle_Male_immac <- sum(res_Muscle_Male_immac$padj <= Wald_FDR, na.rm = TRUE)
print("Summary of DESeq2 results for Muscle Male_immac:")
cat("Muscle Male_immac - Number of significant DEGs (padj <= ", Wald_FDR, "): ", n_significant_Muscle_Male_immac, "\n")
summary(res_Muscle_Male_immac)

################################
###   Muscle - Male  - mel   ###
################################

# Subset dkey for Muscle Male_mel
dkey_Muscle_Male_mel <- subset(dkey, Tissue == "Muscle" & Sex == "Male" & Morph == "mel")
print("Sample information for Muscle Male mel:")
dkey_Muscle_Male_mel

# Get sample IDs for Muscle Male_mel
samples_Muscle_Male_mel <- as.character(dkey_Muscle_Male_mel$Sample)

# Subset the counts table for Muscle Male_mel
df_Muscle_Male_mel <- df[, samples_Muscle_Male_mel]
print("Head of counts table for Muscle Male_mel:")
head(df_Muscle_Male_mel)

# Check if sample names match without modifying any names
samples_match_Muscle_Male_mel <- all(sort(dkey_Muscle_Male_mel$Sample) == sort(colnames(df_Muscle_Male_mel)))
print("Do sample names match for Muscle Male_mel?")
print(samples_match_Muscle_Male_mel)

# Filter out low expression genes
columns_Muscle_Male_mel <- ncol(df_Muscle_Male_mel)
nrow_before_filter_Muscle_Male_mel <- nrow(df_Muscle_Male_mel)
df_Muscle_Male_mel <- df_Muscle_Male_mel[rowSums(df_Muscle_Male_mel) >= (columns_Muscle_Male_mel * coveragethrd), ]
nrow_after_filter_Muscle_Male_mel <- nrow(df_Muscle_Male_mel)
cat("Muscle Male_mel - Number of genes before filtering:", nrow_before_filter_Muscle_Male_mel, "\n")
cat("Muscle Male_mel - Number of genes after filtering:", nrow_after_filter_Muscle_Male_mel, "\n")

# Create DESeqDataSet
dds_Muscle_Male_mel <- DESeqDataSetFromMatrix(countData = df_Muscle_Male_mel,
                                              colData = dkey_Muscle_Male_mel,
                                              design = ~ Treatment)

# Relevel Treatment to have 'veh' as reference
dds_Muscle_Male_mel$Treatment <- relevel(dds_Muscle_Male_mel$Treatment, ref = "veh")

# Run DESeq
dds_Muscle_Male_mel <- DESeq(dds_Muscle_Male_mel)

# Get results (unshrunken log2 fold changes)
res_Muscle_Male_mel <- results(dds_Muscle_Male_mel, alpha = Wald_FDR)
res_Muscle_Male_mel <- res_Muscle_Male_mel[order(res_Muscle_Male_mel$padj), ]
n_significant_Muscle_Male_mel <- sum(res_Muscle_Male_mel$padj <= Wald_FDR, na.rm = TRUE)
print("Summary of DESeq2 results for Muscle Male_mel:")
cat("Muscle Male_mel - Number of significant DEGs (padj <= ", Wald_FDR, "): ", n_significant_Muscle_Male_mel, "\n")
summary(res_Muscle_Male_mel)

#########################
### Results Summaries ###
#########################

print("Summary of DESeq2 results for Brain Male_parae:")
cat("Brain Male_parae - Number of significant DEGs (padj <= ", Wald_FDR, "): ", n_significant_Brain_Male_parae, "\n")
summary(res_Brain_Male_parae)

print("Summary of DESeq2 results for Brain Male_immac:")
cat("Brain Male_immac - Number of significant DEGs (padj <= ", Wald_FDR, "): ", n_significant_Brain_Male_immac, "\n")
summary(res_Brain_Male_immac)

print("Summary of DESeq2 results for Brain Male_mel:")
cat("Brain Male_mel - Number of significant DEGs (padj <= ", Wald_FDR, "): ", n_significant_Brain_Male_mel, "\n")
summary(res_Brain_Male_mel)

print("Summary of DESeq2 results for Muscle Male_parae:")
cat("Muscle Male_parae - Number of significant DEGs (padj <= ", Wald_FDR, "): ", n_significant_Muscle_Male_parae, "\n")
summary(res_Muscle_Male_parae)

print("Summary of DESeq2 results for Muscle Male_immac:")
cat("Muscle Male_immac - Number of significant DEGs (padj <= ", Wald_FDR, "): ", n_significant_Muscle_Male_immac, "\n")
summary(res_Muscle_Male_immac)

print("Summary of DESeq2 results for Muscle Male_mel:")
cat("Muscle Male_mel - Number of significant DEGs (padj <= ", Wald_FDR, "): ", n_significant_Muscle_Male_mel, "\n")
summary(res_Muscle_Male_mel)

# Define a list of existing result objects and their corresponding names
result_list <- list(
  list(result = res_Brain_Male_parae, name = "Brain_Male_parae", n_significant = n_significant_Brain_Male_parae),
  list(result = res_Brain_Male_immac, name = "Brain_Male_immac", n_significant = n_significant_Brain_Male_immac),
  list(result = res_Brain_Male_mel, name = "Brain_Male_mel", n_significant = n_significant_Brain_Male_mel),
  list(result = res_Muscle_Male_parae, name = "Muscle_Male_parae", n_significant = n_significant_Muscle_Male_parae),
  list(result = res_Muscle_Male_immac, name = "Muscle_Male_immac", n_significant = n_significant_Muscle_Male_immac),
  list(result = res_Muscle_Male_mel, name = "Muscle_Male_mel", n_significant = n_significant_Muscle_Male_mel)
)

# Loop over the results and print summaries to files
for (res_obj in result_list) {
  result <- res_obj$result
  name <- res_obj$name
  n_significant <- res_obj$n_significant
  
  # Capture the output of the summary and number of significant DEGs
  summary_text <- capture.output({
    print(paste("Summary of DESeq2 results for", name, ":"))
    cat(name, "- Number of significant DEGs (padj <= Wald_FDR):", n_significant, "\n")
    summary(result)
  })
  
  # Save summary to the Summaries folder
  summary_file <- paste0("Summaries/Treatment/Wald/", name, "_DESeq2_Treatment_summary.txt")
  writeLines(summary_text, summary_file)
}

###################################
#   Print DESeq reuslts to files  #
###################################

# Perform apeglm log2fc shrinkage
{
  # Perform lfcShrink to get shrunken log2 fold changes
  print("Available coefficients for lfcShrink in Brain Male_parae:")
  print(resultsNames(dds_Brain_Male_parae))
  res_Brain_Male_parae_shrink <- lfcShrink(dds_Brain_Male_parae, coef = "Treatment_flx_vs_veh", type = "apeglm")
  res_Brain_Male_parae_shrink <- res_Brain_Male_parae_shrink[order(res_Brain_Male_parae_shrink$padj), ]
  
  # Perform lfcShrink to get shrunken log2 fold changes
  print("Available coefficients for lfcShrink in Brain Male_immac:")
  print(resultsNames(dds_Brain_Male_immac))
  res_Brain_Male_immac_shrink <- lfcShrink(dds_Brain_Male_immac, coef = "Treatment_flx_vs_veh", type = "apeglm")
  res_Brain_Male_immac_shrink <- res_Brain_Male_immac_shrink[order(res_Brain_Male_immac_shrink$padj), ]
  
  # Perform lfcShrink to get shrunken log2 fold changes
  print("Available coefficients for lfcShrink in Brain Male_mel:")
  print(resultsNames(dds_Brain_Male_mel))
  res_Brain_Male_mel_shrink <- lfcShrink(dds_Brain_Male_mel, coef = "Treatment_flx_vs_veh", type = "apeglm")
  res_Brain_Male_mel_shrink <- res_Brain_Male_mel_shrink[order(res_Brain_Male_mel_shrink$padj), ]
  
  # Perform lfcShrink to get shrunken log2 fold changes
  print("Available coefficients for lfcShrink in Muscle Male_parae:")
  print(resultsNames(dds_Muscle_Male_parae))
  res_Muscle_Male_parae_shrink <- lfcShrink(dds_Muscle_Male_parae, coef = "Treatment_flx_vs_veh", type = "apeglm")
  res_Muscle_Male_parae_shrink <- res_Muscle_Male_parae_shrink[order(res_Muscle_Male_parae_shrink$padj), ]
  
  # Perform lfcShrink to get shrunken log2 fold changes
  print("Available coefficients for lfcShrink in Muscle Male_immac:")
  print(resultsNames(dds_Muscle_Male_immac))
  res_Muscle_Male_immac_shrink <- lfcShrink(dds_Muscle_Male_immac, coef = "Treatment_flx_vs_veh", type = "apeglm")
  res_Muscle_Male_immac_shrink <- res_Muscle_Male_immac_shrink[order(res_Muscle_Male_immac_shrink$padj), ]
  
  # Perform lfcShrink to get shrunken log2 fold changes
  print("Available coefficients for lfcShrink in Muscle Male_mel:")
  print(resultsNames(dds_Muscle_Male_mel))
  res_Muscle_Male_mel_shrink <- lfcShrink(dds_Muscle_Male_mel, coef = "Treatment_flx_vs_veh", type = "apeglm")
  res_Muscle_Male_mel_shrink <- res_Muscle_Male_mel_shrink[order(res_Muscle_Male_mel_shrink$padj), ]
}

# Write results to files
{
# Save unshrunken results
write.table(as.data.frame(res_Brain_Male_parae),
            file = "Results_tables/Treatment/Brain/No_shrink/Brain_Male_parae_DE_results_unshrunken.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = NA)

# Save shrunken results
write.table(as.data.frame(res_Brain_Male_parae_shrink),
            file = "Results_tables/Treatment/Brain/l2fc_shrink/Brain_Male_parae_DE_results_shrunken.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = NA)

# Save unshrunken results
write.table(as.data.frame(res_Brain_Male_immac),
            file = "Results_tables/Treatment/Brain/No_shrink/Brain_Male_immac_DE_results_unshrunken.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = NA)

# Save shrunken results
write.table(as.data.frame(res_Brain_Male_immac_shrink),
            file = "Results_tables/Treatment/Brain/l2fc_shrink/Brain_Male_immac_DE_results_shrunken.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = NA)

# Save unshrunken results
write.table(as.data.frame(res_Brain_Male_mel),
            file = "Results_tables/Treatment/Brain/No_shrink/Brain_Male_mel_DE_results_unshrunken.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = NA)

# Save shrunken results
write.table(as.data.frame(res_Brain_Male_mel_shrink),
            file = "Results_tables/Treatment/Brain/l2fc_shrink/Brain_Male_mel_DE_results_shrunken.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = NA)

# Save unshrunken results
write.table(as.data.frame(res_Muscle_Male_parae),
            file = "Results_tables/Treatment/Muscle/No_shrink/Muscle_Male_parae_DE_results_unshrunken.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = NA)

# Save shrunken results
write.table(as.data.frame(res_Muscle_Male_parae_shrink),
            file = "Results_tables/Treatment/Muscle/l2fc_shrink/Muscle_Male_parae_DE_results_shrunken.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = NA)

# Save unshrunken results
write.table(as.data.frame(res_Muscle_Male_immac),
            file = "Results_tables/Treatment/Muscle/No_shrink/Muscle_Male_immac_DE_results_unshrunken.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = NA)

# Save shrunken results
write.table(as.data.frame(res_Muscle_Male_immac_shrink),
            file = "Results_tables/Treatment/Muscle/l2fc_shrink/Muscle_Male_immac_DE_results_shrunken.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = NA)

# Save unshrunken results
write.table(as.data.frame(res_Muscle_Male_mel),
            file = "Results_tables/Treatment/Muscle/No_shrink/Muscle_Male_mel_DE_results_unshrunken.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = NA)

# Save shrunken results
write.table(as.data.frame(res_Muscle_Male_mel_shrink),
            file = "Results_tables/Treatment/Muscle/l2fc_shrink/Muscle_Male_mel_DE_results_shrunken.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = NA)
}

#########################
###     PCA Plots     ###
#########################
# Function to create PCA plot and save as image
create_pca_plot <- function(dds, dkey_subset, filename, title) {
  # Get normalized counts
  counts_matrix <- counts(dds, normalized = TRUE)
  
  # Filter out genes with zero variance that cause PCA errors
  gene_vars <- apply(counts_matrix, 1, var)
  counts_matrix <- counts_matrix[gene_vars > 0, ]
  
  # Transpose for PCA (samples as rows)
  t_counts_matrix <- t(counts_matrix)
  
  # Run PCA
  prcomp_result <- prcomp(t_counts_matrix, scale. = TRUE)
  
  # Prepare PCA data frame
  pca_scores <- prcomp_result$x[, 1:2]
  pca_scores_df <- as.data.frame(pca_scores)
  pca_scores_df$Sample <- rownames(pca_scores_df)
  pca_df <- merge(pca_scores_df, dkey_subset, by = "Sample")
  pca_df$Treatment_label <- factor(pca_df$Treatment, levels = c("flx", "veh"),
                                   labels = c("Flutamide", "Vehicle"))
  
  # Variance explained for PC1 and PC2
  variance_explained <- prcomp_result$sdev[1:2]^2 / sum(prcomp_result$sdev^2) * 100
  
  # Plot and save to file
  png(filename = filename, units = "in", width = 10, height = 8, pointsize = 24, res = 600)
  print(
    ggplot(data = pca_df, aes(x = PC1, y = PC2, color = Treatment_label)) +
      geom_point(size = 3) +
      labs(
        title = title,
        x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
        y = paste0("PC2 (", round(variance_explained[2], 2), "%)"),
        color = "Treatment"
      ) +
      scale_color_manual(values = c("orange3", "blue")) +
      theme_classic() +
      guides(color = guide_legend(title = NULL)) +
      geom_encircle(
        aes(x = PC1, y = PC2),
        data = subset(pca_df, Treatment_label == "Flutamide"),
        color = "orange3", fill = "orange3", alpha = 0.2,
        s_shape = 1, expand = 0.01, show.legend = FALSE
      ) +
      geom_encircle(
        aes(x = PC1, y = PC2),
        data = subset(pca_df, Treatment_label == "Vehicle"),
        color = "blue", fill = "blue", alpha = 0.1,
        s_shape = 1, expand = 0.0, show.legend = FALSE
      ) +
      theme(
        text = element_text(size = 16, face = "bold"),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.box.background = element_rect(color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()
      )
  )
  dev.off()
  
  # Print information about filtering
  cat("PCA for", title, "- Filtered", sum(gene_vars == 0), "genes with zero variance\n")
}

# Create PCA_Plots directory
if (!dir.exists("PCA_Plots")) {
  dir.create("PCA_Plots")
}

# Generate PCA plots for all comparisons
create_pca_plot(dds_Brain_Male_parae, dkey_Brain_Male_parae, "PCA_Plots/Treatment/Brain_Male_parae_PCA.png", "PCA: Male Brain in parae morph")
create_pca_plot(dds_Brain_Male_immac, dkey_Brain_Male_immac, "PCA_Plots/Treatment/Brain_Male_immac_PCA.png", "PCA: Male Brain in immac morph")
create_pca_plot(dds_Brain_Male_mel, dkey_Brain_Male_mel, "PCA_Plots/Treatment/Brain_Male_mel_PCA.png", "PCA: Male Brain in mel morph")
create_pca_plot(dds_Muscle_Male_parae, dkey_Muscle_Male_parae, "PCA_Plots/Treatment/Muscle_Male_parae_PCA.png", "PCA: Male Muscle in parae morph")
create_pca_plot(dds_Muscle_Male_immac, dkey_Muscle_Male_immac, "PCA_Plots/Treatment/Muscle_Male_immac_PCA.png", "PCA: Male Muscle in immac morph")
create_pca_plot(dds_Muscle_Male_mel, dkey_Muscle_Male_mel, "PCA_Plots/Treatment/Muscle_Male_mel_PCA.png", "PCA: Male Muscle in mel morph")
create_pca_plot(dds_Brain_Female_parae, dkey_Brain_Female_parae, "PCA_Plots/Treatment/Brain_Female_parae_PCA.png", "PCA: Female Brain")
create_pca_plot(dds_Muscle_Female_parae, dkey_Muscle_Female_parae, "PCA_Plots/Treatment/Muscle_Female_parae_PCA.png", "PCA: Female Muscle")

#########################
###     Heatmaps      ###
#########################

# Function to create heatmap and save as image
create_heatmap <- function(dds, filename, title, pval_threshold = 0.05, colors = hmcol) {
  # Get DESeq2 results and filter significant genes by padj
  res <- results(dds, alpha = Wald_FDR)
  res_nona <- na.omit(res)
  DE_genes <- rownames(res_nona[res_nona$padj < pval_threshold, ])
  
  # Skip if there are no significant genes
  if (length(DE_genes) == 0) {
    message(paste("No significant genes for", title, "- skipping heatmap."))
    return(NULL)
  }
  
  # Get count matrix for significant genes
  matrix <- counts(dds, normalized = FALSE)[DE_genes, ]
  vst_matrix <- varianceStabilizingTransformation(dds[DE_genes, ])
  
  # Save heatmap as PNG
  png(filename = filename, units = "in", width = 8, height = 8, pointsize = 12, res = 300)
  heatmap3(assay(vst_matrix),     # Use the assay function to extract the matrix from VST object
           method = "complete",   # Hierarchical clustering method
           Rowv = TRUE,           # Enable row clustering
           Colv = NA,             # No column clustering
           col = colors,          # Color palette for heatmap
           scale = "row",         # Scale by row (Z-scores)
           labRow = NA,           # Remove row labels
           showRowDendro = TRUE)  # Show dendrogram for rows
  title(main = title)
  dev.off()
}

# Define heatmap color palette
hmcol <- colorRampPalette(c("blue3", "black", "orange3"))(15)

# Create Heatmaps directory
if (!dir.exists("Heatmaps/Treatment/Wald")) {
  dir.create("Heatmaps/Treatment/Wald")
}

# Generate heatmaps for all comparisons
create_heatmap(dds_Brain_Male_parae, "Heatmaps/Treatment/Wald/Brain_Male_parae_heatmap.png", "Heatmap: Male Brain in parae morph")
create_heatmap(dds_Brain_Male_immac, "Heatmaps/Treatment/Wald/Brain_Male_immac_heatmap.png", "Heatmap: Male Brain in immac morph")
create_heatmap(dds_Brain_Male_mel, "Heatmaps/Treatment/Wald/Brain_Male_mel_heatmap.png", "Heatmap: Male Brain in mel morph")
create_heatmap(dds_Muscle_Male_parae, "Heatmaps/Treatment/Wald/Muscle_Male_parae_heatmap.png", "Heatmap: Male Muscle in parae morph")
create_heatmap(dds_Muscle_Male_immac, "Heatmaps/Treatment/Wald/Muscle_Male_immac_heatmap.png", "Heatmap: Male Muscle in immac morph")
create_heatmap(dds_Muscle_Male_mel, "Heatmaps/Treatment/Wald/Muscle_Male_mel_heatmap.png", "Heatmap: Male Muscle in mel morph")
create_heatmap(dds_Brain_Female_parae, "Heatmaps/Treatment/Wald/Brain_Female_parae_heatmap.png", "Heatmap: Female Brain")
create_heatmap(dds_Muscle_Female_parae, "Heatmaps/Treatment/Wald/Muscle_Female_parae_heatmap.png", "Heatmap: Female Muscle")

###   Venn Diagrams   ###
#########################

############### Male Brain #####################
# Create lists of DE genes in each group
{
DE_Brain_Male_parae <- rownames(res_Brain_Male_parae)[!is.na(res_Brain_Male_parae$padj) & res_Brain_Male_parae$padj < 0.05]
DE_Brain_Male_mel <- rownames(res_Brain_Male_mel)[!is.na(res_Brain_Male_mel$padj) & res_Brain_Male_mel$padj < 0.05]
DE_Brain_Male_immac <- rownames(res_Brain_Male_immac)[!is.na(res_Brain_Male_immac$padj) & res_Brain_Male_immac$padj < 0.05]

length(DE_Brain_Male_immac)
length(DE_Brain_Male_mel)
length(DE_Brain_Male_parae)
}

# Create the sets for the Euler diagram
{
venn_sets_Brain_Male <- list(
  Parae_Brain = DE_Brain_Male_parae,
  Immaculata_Brain = DE_Brain_Male_immac,
  Melanzona_Brain = DE_Brain_Male_mel
)

# Define colors
venn_colors_Brain <- c("blue", "green4", "yellow")

# Generate an euler object
fit_Brain_Male <- euler(venn_sets_Brain_Male)
fit_grob_Brain_Male <- plot(
  fit_Brain_Male,
  fills = venn_colors_Brain,
  edges = "black",
  labels = FALSE,
  quantities = list(cex = 2),
  alpha = 0.4,
)
}

# Male Brain figure
{
png("EulerDiagrams/Male_Brain_DE_genes_euler.png", width = 9, height = 8.25, units = "in", res = 100)
grid.newpage()
pushViewport(viewport(width = 0.7, height = 0.7))
grid.draw(fit_grob_Brain_Male)
popViewport()
grid.text("DE genes in Male Brain",
          x = unit(0.5, "npc"), y = unit(0.96, "npc"),
          gp = gpar(cex = 2)
)
grid.text("Melanzona",
          x = unit(0.75, "npc"), y = unit(0.65, "npc"),
          gp = gpar(cex = 2)
)
grid.text("Parae",
          x = unit(0.12, "npc"), y = unit(0.65, "npc"),
          gp = gpar(cex = 2)
)
dev.off()
}
############### Male Muscle #####################
# Create lists of DE genes in each group
{
DE_Muscle_Male_parae <- rownames(res_Muscle_Male_parae)[!is.na(res_Muscle_Male_parae$padj) & res_Muscle_Male_parae$padj < 0.05]
DE_Muscle_Male_mel <- rownames(res_Muscle_Male_mel)[!is.na(res_Muscle_Male_mel$padj) & res_Muscle_Male_mel$padj < 0.05]
DE_Muscle_Male_immac <- rownames(res_Muscle_Male_immac)[!is.na(res_Muscle_Male_immac$padj) & res_Muscle_Male_immac$padj < 0.05]

length(DE_Muscle_Male_immac)
length(DE_Muscle_Male_mel)
length(DE_Muscle_Male_parae)
}

# Create the sets for the Euler diagram
{
venn_sets_Muscle_Male <- list(
  Parae_Muscle = DE_Muscle_Male_parae,
  Immaculata_Muscle = DE_Muscle_Male_immac,
  Melanzona_Muscle = DE_Muscle_Male_mel
)

# Define colors
venn_colors_Muscle <- c("blue", "green4", "yellow")

# Generate an euler object
fit_Muscle_Male <- euler(venn_sets_Muscle_Male)
fit_grob_Muscle_Male <- plot(
  fit_Muscle_Male,
  fills = venn_colors_Muscle,
  edges = "black",
  labels = FALSE,
  quantities = list(cex = 2),
  alpha = 0.4,
)
}

# Male Muscle figure
{
png("EulerDiagrams/Male_Muscle_DE_genes_euler.png", width = 9, height = 8.25, units = "in", res = 100)
grid.newpage()
pushViewport(viewport(width = 0.7, height = 0.7))
grid.draw(fit_grob_Muscle_Male)
popViewport()
grid.text("DE genes in Male Muscle",
          x = unit(0.5, "npc"), y = unit(0.96, "npc"),
          gp = gpar(cex = 2)
)
grid.text("Parae",
          x = unit(0.14, "npc"), y = unit(0.67, "npc"),
          gp = gpar(cex = 2)
)
grid.text("Melanzona",
          x = unit(0.86, "npc"), y = unit(0.79, "npc"),
          gp = gpar(cex = 2)
)
grid.text("Immaculata",
          x = unit(0.87, "npc"), y = unit(0.415, "npc"),
          gp = gpar(cex = 2)
)
dev.off()
}

#########################
###       UpSet       ###
#########################
# Create lists of DE genes in Brain=
{
DE_Brain_Male_parae <- rownames(res_Brain_Male_parae)[!is.na(res_Brain_Male_parae$padj) & res_Brain_Male_parae$padj < 0.05]
DE_Brain_Male_mel <- rownames(res_Brain_Male_mel)[!is.na(res_Brain_Male_mel$padj) & res_Brain_Male_mel$padj < 0.05]
DE_Brain_Male_immac <- rownames(res_Brain_Male_immac)[!is.na(res_Brain_Male_immac$padj) & res_Brain_Male_immac$padj < 0.05]

# Create lists of DE genes in Muscle
DE_Muscle_Male_parae <- rownames(res_Muscle_Male_parae)[!is.na(res_Muscle_Male_parae$padj) & res_Muscle_Male_parae$padj < 0.05]
DE_Muscle_Male_mel <- rownames(res_Muscle_Male_mel)[!is.na(res_Muscle_Male_mel$padj) & res_Muscle_Male_mel$padj < 0.05]
DE_Muscle_Male_immac <- rownames(res_Muscle_Male_immac)[!is.na(res_Muscle_Male_immac$padj) & res_Muscle_Male_immac$padj < 0.05]

# Brain
length(DE_Brain_Male_parae)
length(DE_Brain_Male_mel)
length(DE_Brain_Male_immac)

# Muscle
length(DE_Muscle_Male_immac)
length(DE_Muscle_Male_mel)
length(DE_Muscle_Male_parae)

# Create a binary matrix for Male Brain data
Male_Brain_genes <- unique(c(DE_Brain_Male_parae, DE_Brain_Male_immac, DE_Brain_Male_mel))
Male_Brain_data <- data.frame(
  Parae = ifelse(Male_Brain_genes %in% DE_Brain_Male_parae, 1, 0),
  Immaculata = ifelse(Male_Brain_genes %in% DE_Brain_Male_immac, 1, 0),
  Melanzona = ifelse(Male_Brain_genes %in% DE_Brain_Male_mel, 1, 0)
)
rownames(Male_Brain_data) <- Male_Brain_genes
# Create a binary matrix for Male Muscle data
Male_Muscle_genes <- unique(c(DE_Muscle_Male_parae, DE_Muscle_Male_immac, DE_Muscle_Male_mel))
Male_Muscle_data <- data.frame(
  Parae = ifelse(Male_Muscle_genes %in% DE_Muscle_Male_parae, 1, 0),
  Immaculata = ifelse(Male_Muscle_genes %in% DE_Muscle_Male_immac, 1, 0),
  Melanzona = ifelse(Male_Muscle_genes %in% DE_Muscle_Male_mel, 1, 0)
)
rownames(Male_Muscle_data) <- Male_Muscle_genes
}

# Generate UpSet plot for Male Brain
png(filename = "UpSet/Treatment/Wald/Brain_Male_UpSet.png", units = "in", width = 7 , height = 5, res = 600)
upset(
  Male_Brain_data,
  c("Parae",  "Immaculata", "Melanzona"),
  name = "DE Gene Overlap in Male Brain",
  width_ratio = 0.3,
  height_ratio = 0.3,
  stripes = 'grey95',
  set_sizes = (
    upset_set_size(
      geom = geom_bar(fill = "steelblue", stat = "count", width = 0.8),
      position = "right"
    ) +
      geom_text(
        aes(label = after_stat(count), y = after_stat(count)), # Use `after_stat(count)` for labels
        stat = "count",
        hjust = -0.2, # Adjust position of the text
        size = 3 # Adjust text size
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) + # Add extra space on the right
      theme(
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 9), # Adjust axis text size
        axis.title = element_text(size = 9) # Adjust axis title size
      )
  ),
  matrix = (
    intersection_matrix(
      geom = geom_point(
        size = 2.5, 
        color = "transparent"
      ),
      segment = geom_segment(
        linetype = 'solid',
        linewidth = 1,
        color = "steelblue" # Change intersection lines to orange3
      ),
      outline_color = list(
        active = "steelblue", # Set active intersection outline color to orange3
        inactive = "transparent" # Make the inactive outline transparent
      )
    )
  ),
  base_annotations = list(
    'Intersection size' =
      intersection_size(counts = FALSE, mapping = aes(fill = 'bars_color')) +
      scale_fill_manual(values = c('bars_color' = 'steelblue'), guide = 'none') +
      geom_text(
        aes(label = after_stat(count), y = after_stat(count)), # Use `after_stat(count)` for consistency
        stat = "count",
        vjust = -0.3, # Adjust this value to place the text above the bar
        size = 3.5 # Adjust text size if needed
      ) +
      theme(
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 9), # Intersection size axis text
        axis.title = element_text(size = 10) # Intersection size axis title
      ) +
      ylab('Number of DE genes in intersection')
  )
)
dev.off()
  
# Generate UpSet plot for Male Muscle
png(filename = "UpSet/Treatment/Wald/Muscle_Male_UpSet.png", units = "in", width = 7 , height = 5, res = 600)
upset(
  Male_Muscle_data,
  c("Parae", "Immaculata", "Melanzona"),
  name = "DE Gene Overlap in Male Muscle",
  width_ratio = 0.3,
  height_ratio = 0.3,
  stripes = 'grey95',
  set_sizes = (
    upset_set_size(
      geom = geom_bar(fill = "steelblue", stat = "count", width = 0.8),
      position = "right"
    ) +
      geom_text(
        aes(label = after_stat(count), y = after_stat(count)), # Use `after_stat(count)` for labels
        stat = "count",
        hjust = -0.2, # Adjust position of the text
        size = 3 # Adjust text size
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) + # Add extra space on the right
      theme(
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 9), # Adjust axis text size
        axis.title = element_text(size = 9) # Adjust axis title size
      )
  ),
  matrix = (
    intersection_matrix(
      geom = geom_point(
        size = 2.5, 
        color = "transparent"
      ),
      segment = geom_segment(
        linetype = 'solid',
        linewidth = 1,
        color = "steelblue" # Change intersection lines to blue
      ),
      outline_color = list(
        active = "steelblue", # Set active intersection outline color to blue
        inactive = "transparent" # Make the inactive outline transparent
      )
    )
  ),
  base_annotations = list(
    'Intersection size' =
      intersection_size(counts = FALSE, mapping = aes(fill = 'bars_color')) +
      scale_fill_manual(values = c('bars_color' = 'steelblue'), guide = 'none') +
      geom_text(
        aes(label = after_stat(count), y = after_stat(count)), # Use `after_stat(count)` for consistency
        stat = "count",
        vjust = -0.3, # Adjust this value to place the text above the bar
        size = 3.5 # Adjust text size if needed
      ) +
      theme(
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 9), # Intersection size axis text
        axis.title = element_text(size = 10) # Intersection size axis title
      ) +
      ylab('Number of DE genes in intersection')
  )
)
dev.off()


################################################################################
################################################################################

create_volcano_plot <- function(res_object, filename, title, pval_cutoff = 0.05, lfc_cutoff = 0, xlim_range = c(-15, 15), ylim_range = c(0, 15)) {
  # Convert results object to data frame and remove NA values
  res_table <- as.data.frame(na.omit(res_object))
  
  # Calculate -log10(padj)
  res_table$neglogpadj <- -log10(res_table$padj)
  
  # Assign significance based on padj and log2FoldChange thresholds
  res_table$Significance <- "Not Significant"
  res_table$Significance[res_table$log2FoldChange > lfc_cutoff & res_table$padj < pval_cutoff] <- "Flutamide"
  res_table$Significance[res_table$log2FoldChange < -lfc_cutoff & res_table$padj < pval_cutoff] <- "Vehicle"
  
  # Custom colors for different significance levels
  color_values <- c("Flutamide" = "red2", "Vehicle" = "blue3", "Not Significant" = "grey50")
  
  # Define y-axis limit
  ylim_max <- 6
  
  # Identify significant points that exceed y-axis limit
  out_of_range <- res_table[res_table$neglogpadj > (ylim_max + 0.25) & res_table$Significance != "Not Significant", ]
  
  # Save volcano plot as EPS using cairo_ps for transparency support
  cairo_ps(filename = filename, width = 4, height = 4)
  
  plot <- ggplot(res_table, aes(x = log2FoldChange, y = neglogpadj, color = Significance)) +
    geom_point(alpha = 0.8, size = 3) +
    scale_color_manual(values = color_values) +
    labs(
      x = expression(""),
      y = expression("")
    ) +
    theme_classic() +
    theme(
      text = element_text(size = 12),
      legend.title = element_blank(),
      legend.position = "",
      axis.line = element_line(linewidth = 0.8)
    ) +
    geom_segment(x = -8.15, xend = 8.15, y = -log10(pval_cutoff), yend = -log10(pval_cutoff), linetype = "dashed", color = "black", linewidth = 1) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "black", linewidth = 0.8) +
    coord_cartesian(xlim = c(-7.5, 7.5), ylim = c(0, 6))
  
  # Add arrows for out-of-range points
  if (nrow(out_of_range) > 0) {
    # Format p-values in scientific notation
    out_of_range$pval_label <- sapply(out_of_range$padj, function(p) {
      exponent <- floor(log10(p))
      mantissa <- round(p / 10^exponent, 1)
      paste0("(p=", mantissa, "e", exponent, ")")
    })
    
    plot <- plot +
      geom_segment(data = out_of_range, aes(x = log2FoldChange, xend = log2FoldChange, y = ylim_max - 0.15, yend = ylim_max + 0.2),
                   arrow = arrow(length = unit(0.15, "cm"), type = "closed"), linewidth = 2, linejoin = "mitre")
#      geom_text(data = out_of_range, aes(x = log2FoldChange, y = ylim_max - 0.4, label = pval_label), size = 3.5, hjust = 0.25)
  }
  
  print(plot)
  dev.off()
}
############################
# Generate volcano plots
{
  create_volcano_plot(res_Brain_Male_parae, "Volcano_Plots/Treatment/Wald/Brain_Male_parae_volcano.eps", "Volcano Plot: Male Brain in parae morph")
  create_volcano_plot(res_Brain_Male_immac, "Volcano_Plots/Treatment/Wald/Brain_Male_immac_volcano.eps", "Volcano Plot: Male Brain in immac morph")
  create_volcano_plot(res_Brain_Male_mel, "Volcano_Plots/Treatment/Wald/Brain_Male_mel_volcano.eps", "Volcano Plot: Male Brain in mel morph")
  create_volcano_plot(res_Muscle_Male_parae, "Volcano_Plots/Treatment/Wald/Muscle_Male_parae_volcano.eps", "Volcano Plot: Male Muscle in parae morph")
  create_volcano_plot(res_Muscle_Male_immac, "Volcano_Plots/Treatment/Wald/Muscle_Male_immac_volcano.eps", "Volcano Plot: Male Muscle in immac morph")
  create_volcano_plot(res_Muscle_Male_mel, "Volcano_Plots/Treatment/Wald/Muscle_Male_mel_volcano.eps", "Volcano Plot: Male Muscle in mel morph")
}
