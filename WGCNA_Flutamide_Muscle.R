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
load("Muscle_WGCNA_complete_analysis.RData")

###############################################################################
# Load Sample Information
###############################################################################
dkey_Muscle <- read.table("Sample_Info.tsv", header=TRUE, sep="\t")
dkey_Muscle <- dkey_Muscle[dkey_Muscle$Tissue == "Muscle", ]
dkey_Muscle <- dkey_Muscle[dkey_Muscle$Sex == "Male", ]
print(dkey_Muscle)

###############################################################################
# Load Expression Data
###############################################################################
options(stringsAsFactors = FALSE)
df_Muscle <- read.table("Combined_Counts.tsv", header=TRUE, row.names=1, sep="\t")

# Keep only columns matching the Muscle samples in dkey_Muscle
df_Muscle <- df_Muscle[, colnames(df_Muscle) %in% dkey_Muscle$Sample]
nrow(df_Muscle)
length(df_Muscle)

# Transpose: rows=Samples, columns=Genes
df_Muscle <- t(df_Muscle)

###############################################################################
# Load Traits Data
###############################################################################
datTraits_Muscle <- read.table("WGCNA_Binary_Matrix_Muscle.tsv", header = TRUE, 
                               sep = "\t", stringsAsFactors = FALSE)
rownames(datTraits_Muscle) <- datTraits_Muscle$Sample
datTraits_Muscle$Sample <- NULL
datTraits_Muscle$Treatment <- NULL
datTraits_Muscle
###############################################################################
# Filter Low-Expression Genes
###############################################################################
coveragethrd_Muscle <- 0.5  # min average reads per sample
cat("Number of genes before filtering:", ncol(df_Muscle), "\n")

df_Muscle <- df_Muscle[, colSums(df_Muscle) >= (nrow(df_Muscle) * coveragethrd_Muscle)]
cat("Number of genes after filtering:", ncol(df_Muscle), "\n")

###############################################################################
# DESeq2 Normalization (VST)
###############################################################################
countData_Muscle <- t(df_Muscle)  # now rows=genes, columns=samples
colData_Muscle   <- dkey_Muscle
rownames(colData_Muscle) <- colData_Muscle$Sample

# Create a DESeqDataSet with a simple ~1 design
dds_Muscle <- DESeqDataSetFromMatrix(countData = countData_Muscle,
                                     colData   = colData_Muscle,
                                     design    = ~ 1)

# Variance stabilizing transformation
vsd_Muscle <- vst(dds_Muscle, blind=TRUE)
dm_Muscle  <- assay(vsd_Muscle)

# Re-transpose so that rows=Samples, columns=Genes
df_Muscle <- t(dm_Muscle)

###############################################################################
# Check for Missing or Zero-Variance Genes
###############################################################################
gsg_Muscle <- goodSamplesGenes(df_Muscle, verbose=3)
if(!gsg_Muscle$allOK){
  if(sum(!gsg_Muscle$goodGenes) > 0) {
    cat("Removing genes:\n", paste(colnames(df_Muscle)[!gsg_Muscle$goodGenes], collapse=", "), "\n")
  }
  if(sum(!gsg_Muscle$goodSamples) > 0){
    cat("Removing samples:\n", paste(rownames(df_Muscle)[!gsg_Muscle$goodSamples], collapse=", "), "\n")
  }
  df_Muscle <- df_Muscle[gsg_Muscle$goodSamples, gsg_Muscle$goodGenes]
}

###############################################################################
# Sample Clustering to Detect Outliers
###############################################################################
sampleTree_Muscle <- hclust(dist(df_Muscle), method="average")
plot(sampleTree_Muscle, 
     main="Muscle Sample Clustering", 
     xlab="", sub="", 
     cex.main=1.5)
abline(h=60, col="red")  
# If you want to remove outliers:
# clust_Muscle <- cutreeStatic(sampleTree_Muscle, cutHeight=40, minSize=10)
# keepSamples_Muscle <- (clust_Muscle==1)
# df_Muscle <- df_Muscle[keepSamples_Muscle, ]
# cat("Samples kept:", sum(keepSamples_Muscle), "\n")

nGenes_Muscle   <- ncol(df_Muscle)
nSamples_Muscle <- nrow(df_Muscle)
cat("Number of Genes:", nGenes_Muscle, "\n")
cat("Number of Samples:", nSamples_Muscle, "\n")

###############################################################################
# Pick Soft Threshold
###############################################################################
powers_Muscle <- 1:20
sft_Muscle <- pickSoftThreshold(df_Muscle, powerVector=powers_Muscle, verbose=5, networkType="signed")

par(mfrow=c(1,2))
# Scale-free topology fit
plot(sft_Muscle$fitIndices[,1], 
     -sign(sft_Muscle$fitIndices[,3])*sft_Muscle$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, signed R^2",
     main="Scale independence", type="n")
text(sft_Muscle$fitIndices[,1], 
     -sign(sft_Muscle$fitIndices[,3])*sft_Muscle$fitIndices[,2],
     labels=powers_Muscle, col="red")
abline(h=0.90, col="red")

# Mean connectivity
plot(sft_Muscle$fitIndices[,1], sft_Muscle$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity",
     main="Mean connectivity", type="n")
text(sft_Muscle$fitIndices[,1], sft_Muscle$fitIndices[,5],
     labels=powers_Muscle, col="red")

softPower_Muscle <- 8  # Adjust based on the above plots
cat("Using softPower =", softPower_Muscle, "\n")

###############################################################################
# Construct Network (Adjacency/TOM) and Gene Clustering
###############################################################################
adjacency_Muscle <- adjacency(df_Muscle, power=softPower_Muscle, type="signed")
TOM_Muscle       <- TOMsimilarity(adjacency_Muscle)
dissTOM_Muscle   <- 1 - TOM_Muscle

geneTree_Muscle <- hclust(as.dist(dissTOM_Muscle), method="average")
plot(geneTree_Muscle, main="Gene clustering on TOM-based dissimilarity in Muscle tissue",
     labels=FALSE, hang=0.04)

###############################################################################
# Dynamic Tree Cutting & Module Merging
###############################################################################
minModuleSize_Muscle <- 100

# Explore different deepSplit parameters
mColorh_Muscle <- NULL
par(mfrow=c(1,1))
for (ds_Muscle in 0:3) {
  tree_Muscle <- cutreeHybrid(dendro=geneTree_Muscle, distM=dissTOM_Muscle,
                              deepSplit=ds_Muscle, pamRespectsDendro=FALSE,
                              minClusterSize=minModuleSize_Muscle, cutHeight=0.99)
  mColorh_Muscle <- cbind(mColorh_Muscle, labels2colors(tree_Muscle$labels))
}
plotDendroAndColors(geneTree_Muscle, mColorh_Muscle, paste("deepSplit=", 0:3), dendroLabels=FALSE)

# Choose a deepSplit parameter
deepSplitParam_Muscle <- 1
dynamicMods_Muscle <- cutreeDynamic(dendro=geneTree_Muscle,
                                    distM=dissTOM_Muscle,
                                    deepSplit=deepSplitParam_Muscle,
                                    pamRespectsDendro=FALSE,
                                    minClusterSize=minModuleSize_Muscle)

dynamicColors_Muscle <- labels2colors(dynamicMods_Muscle)
table(dynamicColors_Muscle)

plotDendroAndColors(geneTree_Muscle, dynamicColors_Muscle, "Dynamic Tree Cut",
                    dendroLabels=FALSE, hang=0.03,
                    addGuide=TRUE, guideHang=0.05)

# Calculate eigengenes, then merge modules if they are too similar
MEList_Muscle <- moduleEigengenes(df_Muscle, colors=dynamicColors_Muscle)
MEs_Muscle <- MEList_Muscle$eigengenes
MEDiss_Muscle <- 1 - cor(MEs_Muscle)
METree_Muscle <- hclust(as.dist(MEDiss_Muscle), method="average")

plot(METree_Muscle, main="Clustering of module eigengenes")
MEDissThres_Muscle <- 0.2
abline(h=MEDissThres_Muscle, col="red")

merge_Muscle <- mergeCloseModules(df_Muscle, dynamicColors_Muscle, cutHeight=MEDissThres_Muscle, verbose=3)
mergedColors_Muscle <- merge_Muscle$colors
mergedMEs_Muscle <- merge_Muscle$newMEs

plotDendroAndColors(geneTree_Muscle, 
                    cbind(dynamicColors_Muscle, mergedColors_Muscle),
                    c("Dynamic Tree Cut","Merged"),
                    dendroLabels=FALSE,
                    hang=0.03)

png("Muscle_geneTree_dendrograms.png", width=800, height=300)
plotDendroAndColors(geneTree_Muscle, 
                    cbind(mergedColors_Muscle),
                    c(""),
                    dendroLabels=FALSE,
                    hang=0.03)
dev.off()

moduleColors_Muscle <- mergedColors_Muscle
# Create numeric module labels
colorOrder_Muscle <- c("grey", standardColors(50))
moduleLabels_Muscle <- match(moduleColors_Muscle, colorOrder_Muscle) - 1
MEs_Muscle <- mergedMEs_Muscle

# Save final objects
save(MEs_Muscle, moduleLabels_Muscle, moduleColors_Muscle, geneTree_Muscle,
     file="Muscle_WGCNA_network_construction.RData")


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
geneTraitCor_Flx <- cor(df_Muscle, datTraits_Muscle$Flutamide, use="pairwise.complete.obs")
geneTraitDiscrete_Flx <- cut(geneTraitCor_Flx,
                             breaks = c(-Inf, -0.25, -0.10, 0.10, 0.25, Inf),
                             labels = c(-1, -0.5, 0, 0.5, 1))
traitColors_Flx <- colourMap[as.character(geneTraitDiscrete_Flx)]

# 2) Parae
geneTraitCor_Parae <- cor(df_Muscle, datTraits_Muscle$Parae, use="pairwise.complete.obs")
geneTraitDiscrete_Parae <- cut(geneTraitCor_Parae,
                               breaks = c(-Inf, -0.25, -0.10, 0.10, 0.25, Inf),
                               labels = c(-1, -0.5, 0, 0.5, 1))
traitColors_Parae <- colourMap[as.character(geneTraitDiscrete_Parae)]

# 3) Melanzona
geneTraitCor_Mel <- cor(df_Muscle, datTraits_Muscle$Melanzona, use="pairwise.complete.obs")
geneTraitDiscrete_Mel <- cut(geneTraitCor_Mel,
                             breaks = c(-Inf, -0.25, -0.10, 0.10, 0.25, Inf),
                             labels = c(-1, -0.5, 0, 0.5, 1))
traitColors_Mel <- colourMap[as.character(geneTraitDiscrete_Mel)]

# 4) Immaculata
geneTraitCor_Imm <- cor(df_Muscle, datTraits_Muscle$Immaculata, use="pairwise.complete.obs")
geneTraitDiscrete_Imm <- cut(geneTraitCor_Imm,
                             breaks = c(-Inf, -0.25, -0.10, 0.10, 0.25, Inf),
                             labels = c(-1, -0.5, 0, 0.5, 1))
traitColors_Imm <- colourMap[as.character(geneTraitDiscrete_Imm)]

# Combine module colours plus each trait
allColors <- cbind(
  Module     = mergedColors_Muscle,
  Flutamide  = traitColors_Flx,
  Parae      = traitColors_Parae,
  Melanzona  = traitColors_Mel,
  Immaculata = traitColors_Imm
)

png("Muscle_Flutamide_geneTree_Cluster_dendrogram.png", width=3000, height=2400, res=300)
plotDendroAndColors(
  geneTree_Muscle,
  allColors,
  groupLabels = c("Modules","All Morphs","Parae","Melanzona","Immaculata"),
  main = "Gene Tree Cluster Dendrogram in Muscle Tissue",
  dendroLabels = FALSE,
  hang = 0.03
)
dev.off()

###############################################################################
# Relate Modules to Traits
###############################################################################
nGenes_Muscle <- ncol(df_Muscle)
nSamples_Muscle <- nrow(df_Muscle)
# Recalculate MEs
MEs_Muscle <- moduleEigengenes(df_Muscle, colors=moduleColors_Muscle)$eigengenes
MEs_Muscle <- orderMEs(MEs_Muscle)
moduleTraitCor_Muscle <- cor(MEs_Muscle, datTraits_Muscle, use="p")
moduleTraitPvalue_Muscle <- corPvalueStudent(moduleTraitCor_Muscle, nSamples_Muscle)

# Create labeled text matrix
textMatrix_Muscle <- paste0(
  "Cor = ", signif(moduleTraitCor_Muscle, 2),
  "\np = ~", signif(moduleTraitPvalue_Muscle, 1)
)
dim(textMatrix_Muscle) <- dim(moduleTraitCor_Muscle)

png("Muscle_Flutamide_by_Morph_Trait_relationships.png", width=2000, height=2000, res = 200)
{
  par(mar=c(6, 8, 3, 3))
  moduleNames  <- names(MEs_Muscle)
  moduleColors <- substring(moduleNames, 3)
  moduleSizes  <- table(moduleColors_Muscle)
  moduleLabels <- sapply(moduleColors, function(col) {
    capName <- paste0(toupper(substr(col, 1, 1)), substr(col, 2, nchar(col)))
    paste0(capName, "\n(", moduleSizes[col], ")")
  })
  labeledHeatmap(
    Matrix      = moduleTraitCor_Muscle,
    xLabels     = names(datTraits_Muscle),
    yLabels     = moduleLabels,
    ySymbols    = moduleLabels,
    colorLabels = FALSE,
    colors      = blueWhiteRed(50),
    textMatrix  = textMatrix_Muscle,
    setStdMargins = FALSE,
    cex.text    = 1.0,
    zlim        = c(-1,1),
    main        = "Module‑Trait Relationships:\nFlutamide Treatment by Morph in Muscle tissue"
  )
}
dev.off()

###############################################################################
# Signed KME and Gene-Module Membership
###############################################################################
signedKMEs_Muscle <- signedKME(df_Muscle, MEs_Muscle)
write.table(signedKMEs_Muscle, "Muscle_Modules/Muscle_signedKMEs.tsv", sep="\t", quote=FALSE)

geneModuleMembership_Muscle <- as.data.frame(signedKMEs_Muscle)
names(geneModuleMembership_Muscle) <- paste("MM", names(MEs_Muscle), sep="")
geneInfo_Muscle <- data.frame(GeneID = colnames(df_Muscle),
                              Module = moduleColors_Muscle)
write.table(geneInfo_Muscle, "Muscle_Modules/Muscle_gene_module_assignment.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)

# Save per-module gene lists
for (mod_Muscle in unique(moduleColors_Muscle)) {
  if (mod_Muscle == "grey") next
  modGenes_Muscle <- colnames(df_Muscle)[moduleColors_Muscle == mod_Muscle]
  fname_Muscle <- paste0("Muscle_Modules/Muscle_module_genes_", mod_Muscle, ".txt")
  write.table(modGenes_Muscle, fname_Muscle, row.names=FALSE, col.names=FALSE, quote=FALSE)
}

###############################################################################
# Optional: TOM Heatmap of Subset of Genes
###############################################################################
nSelect_Muscle <- 1000
if (nGenes_Muscle < nSelect_Muscle) nSelect_Muscle <- nGenes_Muscle
set.seed(10)
select_Muscle <- sample(nGenes_Muscle, nSelect_Muscle)
selectTOM_Muscle <- dissTOM_Muscle[select_Muscle, select_Muscle]
selectTree_Muscle <- hclust(as.dist(selectTOM_Muscle), method="average")
selectColors_Muscle <- moduleColors_Muscle[select_Muscle]

# Raise TOM to a power (e.g. 7) to highlight moderate connections
plotDiss_Muscle <- selectTOM_Muscle^7
diag(plotDiss_Muscle) <- NA

TOMplot(plotDiss_Muscle, selectTree_Muscle, selectColors_Muscle,
        main="Muscle Network heatmap (subset)")

###############################################################################
# Eigengene Network Visualization
###############################################################################
par(cex=0.9)
plotEigengeneNetworks(MEs_Muscle, "Muscle Eigengene Network",
                      marDendro=c(0,4,2,0),
                      marHeatmap=c(3,4,2,2),
                      cex.lab=0.8,
                      xLabelsAngle=90)

###############################################################################
# Summary
###############################################################################
cat("\n===== Muscle Network Analysis Summary =====\n")
cat("Number of genes: ", nGenes_Muscle, "\n")
cat("Number of samples: ", nSamples_Muscle, "\n")
moduleSize_Muscle <- table(moduleColors_Muscle)
cat("Number of modules (excluding grey): ", length(moduleSize_Muscle) - 1, "\n")
print(moduleSize_Muscle)

#save.image(file = "Muscle_WGCNA_complete_analysis.RData")
