setwd("")

library(WGCNA)
library(DESeq2)
library(RColorBrewer)

# Clear environment
rm(list = ls())
graphics.off()

# Enable WGCNA multi-threading
enableWGCNAThreads()

# Load previous session:
load("Brain_WGCNA_complete_analysis.RData")

###############################################################################
# Load Sample Information
###############################################################################
dkey_Brain <- read.table("Sample_Info.tsv", header=TRUE, sep="\t")
dkey_Brain <- dkey_Brain[dkey_Brain$Tissue == "Brain", ]
dkey_Brain <- dkey_Brain[dkey_Brain$Sex == "Male", ]
print(dkey_Brain)

###############################################################################
# Load Expression Data
###############################################################################
options(stringsAsFactors = FALSE)
df_Brain <- read.table("Combined_Counts.tsv", header=TRUE, row.names=1, sep="\t")

# Keep only columns matching the Brain samples in dkey_Brain
df_Brain <- df_Brain[, colnames(df_Brain) %in% dkey_Brain$Sample]
nrow(df_Brain)
length(df_Brain)

# Transpose: rows=Samples, columns=Genes
df_Brain <- t(df_Brain)

###############################################################################
# Load Traits Data
###############################################################################
datTraits_Brain <- read.table("WGCNA_Binary_Matrix_Brain.tsv", header = TRUE, 
                               sep = "\t", stringsAsFactors = FALSE)
rownames(datTraits_Brain) <- datTraits_Brain$Sample
datTraits_Brain$Sample <- NULL
datTraits_Brain$Treatment <- NULL
datTraits_Brain
###############################################################################
# Filter Low-Expression Genes
###############################################################################
coveragethrd_Brain <- 0.5  # min average reads per sample
cat("Number of genes before filtering:", ncol(df_Brain), "\n")

df_Brain <- df_Brain[, colSums(df_Brain) >= (nrow(df_Brain) * coveragethrd_Brain)]
cat("Number of genes after filtering:", ncol(df_Brain), "\n")

###############################################################################
# DESeq2 Normalization (VST)
###############################################################################
countData_Brain <- t(df_Brain)  # now rows=genes, columns=samples
colData_Brain   <- dkey_Brain
rownames(colData_Brain) <- colData_Brain$Sample

# Create a DESeqDataSet with a simple ~1 design
dds_Brain <- DESeqDataSetFromMatrix(countData = countData_Brain,
                                     colData   = colData_Brain,
                                     design    = ~ 1)

# Variance stabilizing transformation
vsd_Brain <- vst(dds_Brain, blind=TRUE)
dm_Brain  <- assay(vsd_Brain)

# Re-transpose so that rows=Samples, columns=Genes
df_Brain <- t(dm_Brain)

###############################################################################
# Check for Missing or Zero-Variance Genes
###############################################################################
gsg_Brain <- goodSamplesGenes(df_Brain, verbose=3)
if(!gsg_Brain$allOK){
  if(sum(!gsg_Brain$goodGenes) > 0) {
    cat("Removing genes:\n", paste(colnames(df_Brain)[!gsg_Brain$goodGenes], collapse=", "), "\n")
  }
  if(sum(!gsg_Brain$goodSamples) > 0){
    cat("Removing samples:\n", paste(rownames(df_Brain)[!gsg_Brain$goodSamples], collapse=", "), "\n")
  }
  df_Brain <- df_Brain[gsg_Brain$goodSamples, gsg_Brain$goodGenes]
}

###############################################################################
# Sample Clustering to Detect Outliers
###############################################################################
sampleTree_Brain <- hclust(dist(df_Brain), method="average")
plot(sampleTree_Brain, 
     main="Brain Sample Clustering", 
     xlab="", sub="", 
     cex.main=1.5)
abline(h=60, col="red")  
# If you want to remove outliers:
# clust_Brain <- cutreeStatic(sampleTree_Brain, cutHeight=40, minSize=10)
# keepSamples_Brain <- (clust_Brain==1)
# df_Brain <- df_Brain[keepSamples_Brain, ]
# cat("Samples kept:", sum(keepSamples_Brain), "\n")

nGenes_Brain   <- ncol(df_Brain)
nSamples_Brain <- nrow(df_Brain)
cat("Number of Genes:", nGenes_Brain, "\n")
cat("Number of Samples:", nSamples_Brain, "\n")

###############################################################################
# Pick Soft Threshold
###############################################################################
powers_Brain <- 1:20
sft_Brain <- pickSoftThreshold(df_Brain, powerVector=powers_Brain, verbose=5, networkType="signed")

par(mfrow=c(1,2))
# Scale-free topology fit
plot(sft_Brain$fitIndices[,1], 
     -sign(sft_Brain$fitIndices[,3])*sft_Brain$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, signed R^2",
     main="Scale independence", type="n")
text(sft_Brain$fitIndices[,1], 
     -sign(sft_Brain$fitIndices[,3])*sft_Brain$fitIndices[,2],
     labels=powers_Brain, col="red")
abline(h=0.90, col="red")

# Mean connectivity
plot(sft_Brain$fitIndices[,1], sft_Brain$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity",
     main="Mean connectivity", type="n")
text(sft_Brain$fitIndices[,1], sft_Brain$fitIndices[,5],
     labels=powers_Brain, col="red")

softPower_Brain <- 10  # Adjust based on the above plots
cat("Using softPower =", softPower_Brain, "\n")

###############################################################################
# Construct Network (Adjacency/TOM) and Gene Clustering
###############################################################################
adjacency_Brain <- adjacency(df_Brain, power=softPower_Brain, type="signed")
TOM_Brain       <- TOMsimilarity(adjacency_Brain)
dissTOM_Brain   <- 1 - TOM_Brain

geneTree_Brain <- hclust(as.dist(dissTOM_Brain), method="average")
plot(geneTree_Brain, main="Gene clustering on TOM-based dissimilarity in Brain tissue",
     labels=FALSE, hang=0.04)





###############################################################################
# Dynamic Tree Cutting & Module Merging
###############################################################################
minModuleSize_Brain <- 100

# Explore different deepSplit parameters
mColorh_Brain <- NULL
par(mfrow=c(1,1))
for (ds_Brain in 0:3) {
  tree_Brain <- cutreeHybrid(dendro=geneTree_Brain, distM=dissTOM_Brain,
                              deepSplit=ds_Brain, pamRespectsDendro=FALSE,
                              minClusterSize=minModuleSize_Brain, cutHeight=0.99)
  mColorh_Brain <- cbind(mColorh_Brain, labels2colors(tree_Brain$labels))
}
plotDendroAndColors(geneTree_Brain, mColorh_Brain, paste("deepSplit=", 0:3), dendroLabels=FALSE)

# Choose a deepSplit parameter
deepSplitParam_Brain <- 1
dynamicMods_Brain <- cutreeDynamic(dendro=geneTree_Brain,
                                    distM=dissTOM_Brain,
                                    deepSplit=deepSplitParam_Brain,
                                    pamRespectsDendro=FALSE,
                                    minClusterSize=minModuleSize_Brain)

dynamicColors_Brain <- labels2colors(dynamicMods_Brain)
table(dynamicColors_Brain)

plotDendroAndColors(geneTree_Brain, dynamicColors_Brain, "Dynamic Tree Cut",
                    dendroLabels=FALSE, hang=0.03,
                    addGuide=TRUE, guideHang=0.05)

# Calculate eigengenes, then merge modules if they are too similar
MEList_Brain <- moduleEigengenes(df_Brain, colors=dynamicColors_Brain)
MEs_Brain <- MEList_Brain$eigengenes
MEDiss_Brain <- 1 - cor(MEs_Brain)
METree_Brain <- hclust(as.dist(MEDiss_Brain), method="average")

plot(METree_Brain, main="Clustering of module eigengenes")
MEDissThres_Brain <- 0.2
abline(h=MEDissThres_Brain, col="red")

merge_Brain <- mergeCloseModules(df_Brain, dynamicColors_Brain, cutHeight=MEDissThres_Brain, verbose=3)
mergedColors_Brain <- merge_Brain$colors
mergedMEs_Brain <- merge_Brain$newMEs

plotDendroAndColors(geneTree_Brain, 
                    cbind(dynamicColors_Brain, mergedColors_Brain),
                    c("Dynamic Tree Cut","Merged"),
                    dendroLabels=FALSE,
                    hang=0.03)

png("Brain_geneTree_dendrograms.png", width=800, height=300)
plotDendroAndColors(geneTree_Brain, 
                    cbind(mergedColors_Brain),
                    c(""),
                    dendroLabels=FALSE,
                    hang=0.03)
dev.off()

moduleColors_Brain <- mergedColors_Brain
# Create numeric module labels
colorOrder_Brain <- c("grey", standardColors(50))
moduleLabels_Brain <- match(moduleColors_Brain, colorOrder_Brain) - 1
MEs_Brain <- mergedMEs_Brain

# Save final objects
save(MEs_Brain, moduleLabels_Brain, moduleColors_Brain, geneTree_Brain,
     file="Brain_WGCNA_network_construction.RData")


###############################################################################
# Create Full GeneTree dendrogram figure

colourMap <- c(
  "-1"   = "blue",
  "-0.5" = "skyblue",
  "0"    = "white",
  "0.5"  = "salmon",
  "1"    = "red"
)

# 1) Flutamide
geneTraitCor_Flx <- cor(df_Brain, datTraits_Brain$Flutamide, use="pairwise.complete.obs")
geneTraitDiscrete_Flx <- cut(geneTraitCor_Flx,
                             breaks = c(-Inf, -0.25, -0.10, 0.10, 0.25, Inf),
                             labels = c(-1, -0.5, 0, 0.5, 1))
traitColors_Flx <- colourMap[as.character(geneTraitDiscrete_Flx)]

# 2) Parae
geneTraitCor_Parae <- cor(df_Brain, datTraits_Brain$Parae, use="pairwise.complete.obs")
geneTraitDiscrete_Parae <- cut(geneTraitCor_Parae,
                               breaks = c(-Inf, -0.25, -0.10, 0.10, 0.25, Inf),
                               labels = c(-1, -0.5, 0, 0.5, 1))
traitColors_Parae <- colourMap[as.character(geneTraitDiscrete_Parae)]

# 3) Melanzona
geneTraitCor_Mel <- cor(df_Brain, datTraits_Brain$Melanzona, use="pairwise.complete.obs")
geneTraitDiscrete_Mel <- cut(geneTraitCor_Mel,
                             breaks = c(-Inf, -0.25, -0.10, 0.10, 0.25, Inf),
                             labels = c(-1, -0.5, 0, 0.5, 1))
traitColors_Mel <- colourMap[as.character(geneTraitDiscrete_Mel)]

# 4) Immaculata
geneTraitCor_Imm <- cor(df_Brain, datTraits_Brain$Immaculata, use="pairwise.complete.obs")
geneTraitDiscrete_Imm <- cut(geneTraitCor_Imm,
                             breaks = c(-Inf, -0.25, -0.10, 0.10, 0.25, Inf),
                             labels = c(-1, -0.5, 0, 0.5, 1))
traitColors_Imm <- colourMap[as.character(geneTraitDiscrete_Imm)]

# Combine module colours plus each trait
allColors <- cbind(
  Module     = mergedColors_Brain,
  Flutamide  = traitColors_Flx,
  Parae      = traitColors_Parae,
  Melanzona  = traitColors_Mel,
  Immaculata = traitColors_Imm
)

#png("Brain_Flutamide_geneTree_Cluster_dendrogram.png", width=3000, height=2400, res=300)
plotDendroAndColors(
  geneTree_Brain,
  allColors,
  groupLabels = c("Modules","All Morphs","Parae","Melanzona","Immaculata"),
  main = "Gene Tree Cluster Dendrogram in Brain Tissue",
  dendroLabels = FALSE,
  hang = 0.03
)
#dev.off()

###############################################################################
# Relate Modules to Traits
###############################################################################
nGenes_Brain <- ncol(df_Brain)
nSamples_Brain <- nrow(df_Brain)
# Recalculate MEs
MEs_Brain <- moduleEigengenes(df_Brain, colors=moduleColors_Brain)$eigengenes
MEs_Brain <- orderMEs(MEs_Brain)
moduleTraitCor_Brain <- cor(MEs_Brain, datTraits_Brain, use="p")
moduleTraitPvalue_Brain <- corPvalueStudent(moduleTraitCor_Brain, nSamples_Brain)

# Create labeled text matrix
textMatrix_Brain <- paste0(
  "Cor = ", signif(moduleTraitCor_Brain, 2),
  "\np = ~", signif(moduleTraitPvalue_Brain, 1)
)
dim(textMatrix_Brain) <- dim(moduleTraitCor_Brain)

png("Brain_Flutamide_by_Morph_Trait_relationships.png", width=2000, height=2000, res = 200)
{
par(mar=c(6, 8, 3, 3))
moduleNames  <- names(MEs_Brain)
moduleColors <- substring(moduleNames, 3)
moduleSizes  <- table(moduleColors_Brain)
moduleLabels <- sapply(moduleColors, function(col) {
  capName <- paste0(toupper(substr(col, 1, 1)), substr(col, 2, nchar(col)))
  paste0(capName, "\n (", moduleSizes[col], ")")
})
labeledHeatmap(
  Matrix      = moduleTraitCor_Brain,
  xLabels     = names(datTraits_Brain),
  yLabels     = moduleLabels,
  ySymbols    = moduleLabels,
  colorLabels = FALSE,
  colors      = blueWhiteRed(50),
  textMatrix  = textMatrix_Brain,
  setStdMargins = FALSE,
  cex.text    = 1.0,
  zlim        = c(-1,1),
  main        = "Module‑Trait Relationships:\nFlutamide Treatment by Morph in Brain tissue"
)
}
dev.off()

###############################################################################
# Signed KME and Gene-Module Membership
###############################################################################
signedKMEs_Brain <- signedKME(df_Brain, MEs_Brain)
write.table(signedKMEs_Brain, "Brain_Modules/Brain_signedKMEs.tsv", sep="\t", quote=FALSE)

geneModuleMembership_Brain <- as.data.frame(signedKMEs_Brain)
names(geneModuleMembership_Brain) <- paste("MM", names(MEs_Brain), sep="")
geneInfo_Brain <- data.frame(GeneID = colnames(df_Brain),
                              Module = moduleColors_Brain)
write.table(geneInfo_Brain, "Brain_Modules/Brain_gene_module_assignment.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)

# Save per-module gene lists
for (mod_Brain in unique(moduleColors_Brain)) {
  if (mod_Brain == "grey") next
  modGenes_Brain <- colnames(df_Brain)[moduleColors_Brain == mod_Brain]
  fname_Brain <- paste0("Brain_Modules/Brain_module_genes_", mod_Brain, ".txt")
  write.table(modGenes_Brain, fname_Brain, row.names=FALSE, col.names=FALSE, quote=FALSE)
}

###############################################################################
# Optional: TOM Heatmap of Subset of Genes
###############################################################################
nSelect_Brain <- 1000
if (nGenes_Brain < nSelect_Brain) nSelect_Brain <- nGenes_Brain
set.seed(10)
select_Brain <- sample(nGenes_Brain, nSelect_Brain)
selectTOM_Brain <- dissTOM_Brain[select_Brain, select_Brain]
selectTree_Brain <- hclust(as.dist(selectTOM_Brain), method="average")
selectColors_Brain <- moduleColors_Brain[select_Brain]

# Raise TOM to a power (e.g. 7) to highlight moderate connections
plotDiss_Brain <- selectTOM_Brain^7
diag(plotDiss_Brain) <- NA

TOMplot(plotDiss_Brain, selectTree_Brain, selectColors_Brain,
        main="Brain Network heatmap (subset)")

###############################################################################
# Eigengene Network Visualization
###############################################################################
par(cex=0.9)
plotEigengeneNetworks(MEs_Brain, "Brain Eigengene Network",
                      marDendro=c(0,4,2,0),
                      marHeatmap=c(3,4,2,2),
                      cex.lab=0.8,
                      xLabelsAngle=90)

###############################################################################
# Summary
###############################################################################
cat("\n===== Brain Network Analysis Summary =====\n")
cat("Number of genes: ", nGenes_Brain, "\n")
cat("Number of samples: ", nSamples_Brain, "\n")
moduleSize_Brain <- table(moduleColors_Brain)
cat("Number of modules (excluding grey): ", length(moduleSize_Brain) - 1, "\n")
print(moduleSize_Brain)

#save.image(file = "Brain_WGCNA_complete_analysis.RData")

