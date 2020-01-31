## 3-way comparison of Mosquito samples ##

library(DESeq2)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(rebus)
library(readr)
library(tximport)
library(pheatmap)
library(tximportData)
library(RColorBrewer)
library(splitstackshape)
library(EnhancedVolcano)
library(PCAtools)
library(AnnotationHub)


# point to directory with count files
dir <- "~/OneDrive/Xie_lab/Mosquito_collab/Salmon_quants"
#list.files(dir)

# point to table with sample info, read Lean v control
samples <- as.data.frame(readxl::read_xlsx(file.path(dir, "Sample_Definitions.xlsx")))
files <- file.path(dir, samples$Sample_Name, "quant.sf")
names(files) <- paste0("sample", 1:110)

# be sure all files are present and/or labeled correctly
all(file.exists(files))

# create trancript to gene conversion table
Aedes.aegypti_GTF <- read.delim("~/OneDrive/Xie_lab/Mosquito_collab/Aedes-aegypti-LVP_AGWG_BASEFEATURES_AaegL5.2.gtf", header=FALSE, comment.char="#")
tx2gene <- strsplit(as.character(Aedes.aegypti_GTF$V9), split = ";")
tx2gene <- cSplit(data.frame(Aedes.aegypti_GTF), "V9", sep=";", drop=TRUE)[,9:10]
tx2gene$V9_1 <- gsub("gene_id ", "", tx2gene$V9_1)
tx2gene$V9_2 <- gsub("transcript_id ", "", tx2gene$V9_2)
names(tx2gene) <- c("gene_id", "transcript_id")
tx2gene <- tx2gene[complete.cases(tx2gene), ]
tx2gene <- tx2gene[,c(2,1)]
rm(Aedes.aegypti_GTF)
## Import transcript-level count data
# Can use tx counts using txOut = TRUE
txi <- tximport(files, type="salmon", tx2gene = tx2gene)
colnames(txi$counts) <- make.unique(samples$`Diet Group`)
# raw counts boxplot

ggplot(melt(as.data.frame(txi$counts), 
            variable.name = "Sample", value.name = "counts"), aes(Sample, counts)) +
  geom_boxplot() +
  xlab("sample") +
  ylab("counts") +
  ggtitle("Raw counts") +
  theme(axis.text.x=element_blank())

# use this to get back to genes
#txi.gene <- summarizeToGene(txi.tx, tx2gene)

# into DEseq
dds <- DESeqDataSetFromTximport(txi, samples, ~Condition)

# keep only non-zero genes/transcripts

# increase all counts by 1 later
keep <- rowSums(counts(ds)) > 1
dds <- ds[ keep, ]
# specify normal as control group
dds$Condition <- relevel(dds$Condition, ref = "C")

# do DE analysis
dds <- DESeq(dds, fitType = 'local')
dds.counts <- counts(dds)
plotDispEsts(dds)

ggplot(melt(as.data.frame(counts(dds, normalized = TRUE)), 
            variable.name = "Sample", value.name = "counts"), aes(Sample, counts)) +
  geom_boxplot() +
  xlab("sample") +
  ylab("counts") +
  ggtitle("Size normalized counts") +
  theme(axis.text.x=element_blank())

# vst normalization of counts (rld log-fold used for lower n, vst for high (>10) n)
#rld <- rlog(dds, blind = FALSE) # takes too long
vst <- varianceStabilizingTransformation(dds, fitType = 'local') #fitType=local for <8 counts, dispersion trend not linear
#head(assay(vst),5)
colnames(vst) <- colnames(dds)
rownames(vst) <- rownames(dds)
vst_df <- melt(as.data.frame(assay(vst)), variable.name = "Sample", value.name = "counts")

#visualize normalized counts
ggplot(vst_df, aes(Sample, counts)) +
  geom_boxplot() +
  xlab("Sample") +
  ylab("Counts (log2)") +
  ggtitle("Variance Stabilized Counts") +
  theme(legend.position = "none", 
    axis.text.x.bottom = element_text(angle = 315, hjust = 0, size = 6))

p <- pca(vst, intgroup="Diet Group")
plotPCA(vst, intgroup="Diet Group", color)

# Obese v Lean
plotPCA(vst[,-24:-54], intgroup='Diet Group')
# PEM v Lean
plotPCA(vst[,24:110], intgroup='Diet Group')


#res.OvC <- results(dds, contrast = c("Condition", "A", "B"))
res.PEMvCont <- results(dds, contrast = c("Condition", "B", "C"), alpha = 0.05, lfcThreshold = 0.5)
res.ObvCont <- results(dds, contrast = c("Condition", "A", "C"), alpha = 0.05, lfcThreshold = 0.5)

PEMvCont.shrink <- lfcShrink(dds, contrast = c("Condition", "B", "C"), type = 'ashr', alpha = 0.05)
ObvCont.shrink <- lfcShrink(dds, contrast = c("Condition", "A", "C"), type = 'ashr', alpha = 0.05)
ObvPEM.shrink <- lfcShrink(dds, contrast = c("Condition", "A", "B"), type = 'ashr', alpha = 0.05)

summary(ObvCont.shrink, alpha=0.05)
summary(PEMvCont.shrink, alpha=0.05)
summary(ObvPEM.shrink, alpha=0.05)


#Volcano plots of results
EnhancedVolcano(ObvCont.shrink, lab = rownames(ObvCont.shrink), x = 'log2FoldChange', y = 'padj',
                title='Obese v Lean Fed Volcano Plot',pointSize = 2.0,labSize = 4.0,
                xlim = c(-2, 4),pCutoff = 0.05, FCcutoff = 0.5, ylim = c(0,20),
                legend=c('NS','Log fold-change','p-value adj','p-value adj & Log fold-change'),
                legendPosition = 'right')

EnhancedVolcano(PEMvCont.shrink, lab = rownames(PEMvCont.shrink), x = 'log2FoldChange', y = 'padj',
                title='PEM v Lean Fed Volcano Plot',pointSize = 2.0,labSize = 4.0,
                xlim = c(-2, 4),pCutoff = 0.05, FCcutoff = 0.5, ylim = c(0,25),
                legend=c('NS','Log fold-change','p-value adj','p-value adj & Log fold-change'),
                legendPosition = 'right')

# Euclidean dist matrix
sampleDistMatrix <- as.matrix(dist(t(assay(vst))))
#make unique col names from sample list
colnames(sampleDistMatrix) <- samples$`Diet Group`
colnames(sampleDistMatrix) <- make.names(colnames(sampleDistMatrix), unique = TRUE)
rownames(sampleDistMatrix) <- colnames(sampleDistMatrix)
pheatmap(sampleDistMatrix,
         clustering_distance_rows="correlation",
         clustering_distance_cols="correlation",
         #annotation_col = mat_col,
         #annotation_colors = mat_colors,
         show_colnames = FALSE,
         show_rownames = FALSE)


# use alternative hyposthesis to find genes above or below difference thresholds
drawLines <- function() abline(h=c(-0.5,0.5),col="dodgerblue",lwd=2)
plotMA(ObvCont.shrink, ylim=c(-4,4), main="Obese v Control log-fold change distribution"); drawLines()
plotMA(PEMvCont.shrink, ylim=c(-4,4), main ="PEM v Control log-fold change distribution"); drawLines()


# Heatmap sorted by variance in genes, 
ObvCont_Sig <- subset(ObvCont.shrink, ObvCont.shrink$padj < 0.05)
PEMvCont_Sig <- subset(PEMvCont.shrink, PEMvCont.shrink$padj < 0.05)

# make lists of genes for later GO analysis
ObvCont_list <- cbind("p-adj" = ObvCont_Sig$padj, "LFC" = ObvCont_Sig$log2FoldChange)
rownames(ObvCont_list) <- rownames(ObvCont_Sig)

PEMvCont_list <- cbind("p-adj" = PEMvCont_Sig$padj, "LFC" = PEMvCont_Sig$log2FoldChange)
rownames(PEMvCont_list) <- rownames(PEMvCont_Sig)

rld_df.OvC <- as.matrix(assay(vst))[rownames(ObvCont_Sig),]
rld_df.PvC <- as.matrix(assay(vst))[rownames(PEMvCont_Sig),]
#highlight any genes that are in both significant gene datasets
dupgenes <- rld_df.OvC[which(rownames(rld_df.OvC) %in% rownames(rld_df.PvC)),]

## Heatmaps
# Obese v Control
topVarGenesOvC <- (order(-rowVars(rld_df.OvC)))
mat <- rld_df.OvC[topVarGenesOvC, ]
mat <- mat - rowMeans(mat)
colnames(mat) <- samples$`Diet Group`
#colnames(mat) <- make.names(colnames(mat), unique = TRUE)

mat <- mat[,-24:-54]
mat_cluster_cols <- hclust(dist(t(mat)))

mat_colors <- list(group = c("#377EB8","#E6AB02"))
col_groups <- c("Obese", "Lean")
names(mat_colors$group) <- unique(col_groups)
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- c("Obese", "Lean")

pheatmap(mat, show_colnames = FALSE,
         show_rownames     = FALSE,
         annotation_col    = mat_col,
         annotation_colors = mat_colors,
         treeheight_row = 0,
         #cluster_cols      = mat_cluster_cols,
         #clustering_callback = callback,
         clustering_method = "complete",
         fontsize_row = 8, fontsize_col = 6,
         main = "Obese v Lean Differentially expressed genes")

# PEM v Control
topVarGenesPvC <- (order(-rowVars(rld_df.PvC)))
mat <- rld_df.PvC[ topVarGenesPvC, ]
mat <- mat - rowMeans(mat)
colnames(mat) <- samples$`Diet Group`
#colnames(mat) <- make.names(colnames(mat), unique = TRUE)

mat <- mat[,24:110]
mat_cluster_cols <- hclust(dist(t(mat)))

mat_colors <- list(group = c("#A6761D","#E6AB02"))
col_groups <- c("PEM", "Lean")

names(mat_colors$group) <- unique(col_groups)
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- c("PEM", "Lean")
pheatmap(mat, show_colnames = F,
         show_rownames     = FALSE,
         annotation_col    = mat_col,
         annotation_colors = mat_colors,
         treeheight_row = 0,
         #cluster_cols      = mat_cluster_cols,
         #clustering_callback = callback,
         clustering_method = "complete",
         fontsize_row = 8, fontsize_col = 6,
         main = "PEM v Control Differentially expressed genes")


# DE genes common to all groups
topVarGenes <- (order(-rowVars(dupgenes)))
mat <- dupgenes[ topVarGenes, ]
mat <- mat - rowMeans(mat)
colnames(mat) <- samples$`Diet Group`
colnames(mat) <- make.names(colnames(mat), unique = TRUE)

mat_cluster_cols <- hclust(dist(t(mat)))

mat_colors <- list(group = c("#377EB8","#A6761D","#E6AB02"))
col_groups <- c("Obese", "PEM", "Lean")
names(mat_colors$group) <- unique(col_groups)
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- c("Obese", "PEM", "Lean")

pheatmap(mat, show_colnames = TRUE,
         show_rownames     = FALSE,
         annotation_col    = mat_col,
         annotation_colors = mat_colors,
         treeheight_row = 0,
         #cluster_cols    = mat_cluster_cols,
         #clustering_callback = callback,
         clustering_method = "complete",
         fontsize_row = 8, fontsize_col = 6,
         main = "Obese, PEM, and Control Differentially expressed genes")

## Two-column gene lists for up/down regulated genes compared to Control group (LEAN)
# contains Gene ID and LFC
OvC_upregulated <- subset(ObvCont_Sig, ObvCont_Sig$log2FoldChange > 0)
OvC_uplist <- do.call(rbind, Map(data.frame, A=rownames(OvC_upregulated), B=OvC_upregulated$log2FoldChange))
OvC_downregulated <- subset(ObvCont_Sig, ObvCont_Sig$log2FoldChange < 0)
OvC_downlist <- do.call(rbind, Map(data.frame, A=rownames(OvC_downregulated), B=OvC_downregulated$log2FoldChange))


PvC_upregulated <- subset(PEMvCont_Sig, PEMvCont_Sig$log2FoldChange > 0)
PvC_uplist <- do.call(rbind, Map(data.frame, A=rownames(PvC_upregulated), B=PvC_upregulated$log2FoldChange))
PvC_downregulated <- subset(PEMvCont_Sig, PEMvCont_Sig$log2FoldChange < 0)
PvC_downlist <- do.call(rbind, Map(data.frame, A=rownames(PvC_downregulated), B=PvC_downregulated$log2FoldChange))

# Plot counts of interesting genes:
plotCounts(dds, "AAEL021465", "Condition")
topVarGenesNames <- rownames(dupgenes[(order(-rowVars(dupgenes))),])
topVarGenesNamesOvC <- rownames(rld_df.OvC[(order(-rowVars(rld_df.OvC))),])
topVarGenesNamesPvC <- rownames(rld_df.PvC[(order(-rowVars(rld_df.PvC))),])


# create subset of interesting counts only with relevant info (gene name, norm counts etc)
tcounts <- t(assay(vst)[topVarGenesNames,]) %>%
  merge(colData(dds), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(topVarGenes)+1):ncol(.))

# plot
ggplot(tcounts, aes(Diet.Group, expression, fill=Condition))+
  geom_boxplot()+
  facet_wrap(~gene, scales="free_y")+
  theme(legend.position = "none", axis.text = element_text(size=12))+
  scale_fill_manual(values=c("#377EB8","#A6761D","#E6AB02"))+
  labs(x = "Condition", y="Expression (normalized log)", fill="Group",title="Top Results")

## GO analysis
library(clusterProfiler)
## ANNOTATIONHUB AND BIOMART DO NOT WORK FOR AEDES ##

## Did not use yet ##
write.table(topVarGenesNames, paste(dir,"/topVarGenes.txt",sep=""),sep="\t", row.names = F, col.names =F, quote = F)
write.table(ObvCont_list, paste(dir,"/OvL_topDEgenes.tsv",sep=""),sep="\t", row.names = T, col.names =NA, quote = F)
write.table(PEMvCont_list, paste(dir,"/PvL_topDEgenes.tsv",sep=""),sep="\t", row.names = T, col.names =NA, quote = F)
# Sample HL7J2AFXY_G3_3_E7 (B.15) (control) is minor outlier in PCA and HC plots

## Plot counts of GTPase genes
# Used online software DAVID for first list, need to filter to GO_ID and p-value for REVIGO
DAVID_out <- read.delim("~/OneDrive/Xie_lab/Mosquito_collab/DAVID_GOoutput_highestStringency.txt", header=T, comment.char="#", skip=1)
GTPgenes <- strsplit(as.character(DAVID_out$Genes[1]), split = ", ")[[1]]
GTPframe <- (data.frame(gene = GTPgenes))
rownames(GTPframe) <- strsplit(as.character(DAVID_out$Genes[1]), split = ", ")[[1]]
GTPgenes <-rownames(GTPframe)

GTPcounts <- t(assay(vst)[GTPgenes,]) %>%
  merge(colData(dds), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(GTPgenes)+1):ncol(.))

ggplot(GTPcounts, aes(Diet.Group, expression, fill=Condition))+
  geom_boxplot()+
  facet_wrap(~gene, scales="free_y")+
  theme(legend.position = "none", axis.text = element_text(size=12))+
  scale_fill_manual(values=c("#377EB8","#A6761D","#E6AB02"))+
  labs(x = "Condition", y="Expression (log2)", fill="Group",title="Small GTPase mediated signal transduction DE genes")
