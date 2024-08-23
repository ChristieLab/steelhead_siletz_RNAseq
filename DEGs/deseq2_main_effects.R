#---------------- Data Description -------------------------------------------------------------------------------------#
# script written by mark christie on 11/25/2019; contact at markchristie1500@gmail.com
# Rversion:
# DESeq2 version:
# script modified from Avril Harder, available at: https://github.com/ChristieLab/Salmo_salar_RNAseq
# covariates include sex, cross date
# relevant sample data can be found in "sample_info_all.csv"

#set working directory, load libraries, list files
setwd("~/local/go/siletz/analyses/differential_expression/main_effects")
library('DESeq2')
library("rtracklayer")
library("ggplot2")
library("gplots")
library("beyonce")
list.files()

#---------------- Data loading and prep --------------------------------------------------------------------------------#
# read in count data generated with featureCounts
gene.counts <- read.table("featureCount_counts.txt", sep="\t", header=T, row.names=1, check.names=F)
drop.vars   <- names(gene.counts) %in% c("Chr","Start", "End", "Strand", "Length")
gene.counts <- gene.counts[!drop.vars]
head(gene.counts)

samples <- read.csv("sample_info_all_including_reciprocal.csv")
samples
#---------------- Sample processing ------------------------------------------------------------------------------------#
# remove HxW and WxH offspring to test for main effects (HxH vs WxW)
m1 <- which(samples[, 2] == "HxW")
m2 <- which(samples[, 2] == "WxH")
crosses <- c(m1, m2)
samples[crosses, ] # should all be reciprocal

# with sample 22 excluded (outlier) and all HxW and WxH offspring excluded
gene.counts <- gene.counts[, -crosses]


# read in CSV with sample information
col.data <- read.csv("sample_info_no_reciprocal.csv", row.names = 1)

# check to make sure that sample names are in the same order in the gene count table and the sample info table; must be TRUE or do not proceed
all(rownames(col.data) %in% colnames(gene.counts))
all(rownames(col.data) == colnames(gene.counts))


#----------------------- Construction of DESeqDataSet object and DE analysis ---------------------------------------------------
# test for effects of treatment while controlling for the effect of family; dds = DESeqDataSet
# variable order matters! Covariates go first, variable of interest goes last in design 


# testing for effect of HxH vs WxW after controlling for cross date , sex of fry, and sequencing run (seq_group)
# sex has little effect
#dds <- DESeqDataSetFromMatrix(countData = gene.counts,
#                              colData = col.data,
#                              design = ~ as.factor(cross_date.2) + as.factor(sex) + as.factor(seq_group) + cross_type)


dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = col.data,
                              design = ~ as.factor(cross_date.2) + as.factor(sex)  + cross_type)

dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = col.data,
                              design = ~ as.factor(cross_date) + as.factor(sex)  + cross_type)




dds$treatment <- relevel(dds$cross_type, "WxW") ## set the untreated group as the control
dds

# renames genes in a DESeqDataSet to match gene names in the reference annotation (i.e., replaces MSTRG wherever gene name is known)
assembly <- readGFF("stringtie_all_merged.gtf")             ## read in merged GTF
gene_idx <- match(dds@rowRanges@partitioning@NAMES, assembly$gene_id)

# create gene_names, which is a table linking gene name, transcript ID, and xloc information for later distinguishing different isoforms of the same gene (same gene start and end) and transcripts generated from different paralogs of the same gene (by xloc)
gene.name      <- assembly$gene_name[gene_idx]
transcript.id  <- assembly$transcript_id[gene_idx]
xloc           <- assembly$exon_number[gene_idx]
gene_names     <- cbind(gene.name, transcript.id, xloc)


dds@rowRanges@partitioning@NAMES <- paste(gene_names[,1],gene_names[,2],gene_names[,3], sep="|") # adds unique gene names to dds
which(duplicated(dds@rowRanges@partitioning@NAMES)) # check to make sure that no gene names in dds are duplicates

# DESeq() = differential expression analysis based on the negative binomial distribution
dds <- DESeq(dds)

# examine distribution of dispersion values
plotDispEsts(dds, ylim=c(1e-6, 1e1))
plotDispEsts(dds)

save.image("steelhead_all.RData")
#load("steelhead_all.RData")

#-----------------------------------------------------------------------------------------------------------------------#

#----------------- Extract and examine results from DESeq analysis -----------------------------------------------------#
res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
length(which(res[, 6] < 0.05))
summary(res)

## examine MA plot (normalized counts vs. log2fold changes)
plotMA(res, ylim=c(-3,3), xlim=c(1e-2,1e6))
plotMA(res, ylim = c(-5, 5))

## examine histograms of p-values and adjusted p-values
hist(res$pvalue, breaks=20, col="grey")
hist(res$padj, breaks=20, col="grey")
#-----------------------------------------------------------------------------------------------------------------------#

################### DATA EXPLORATION ############################
#--<>--<>--<>--<>-- Log transformation and distance calculation --<>--<>--<>--<>--
## rlog() = transforms count data to the log2 scale in a way that minimizes differences between
## samples for genes with small counts, and which normalizes with respect to library size
#rld <- rlog(dds)
# vst is faster alternative
rld <- vst(dds)

head(assay(rld))

## calculate Euclidean distances between all samples to examine overall similarity
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)

## construct heatmaps to examine sample patterns with different clustering methods
rownames(sampleDistMatrix) <- paste(rld$family, rld$treatment, sep="-")
colnames(sampleDistMatrix) <- paste(rld$family, rld$treatment, sep="-")
dists.ward <- hclust(sampleDists, method="ward.D2")
heatmap.2(sampleDistMatrix, trace="none", Rowv=as.dendrogram(dists.ward),
          margins=c(6,7))
dists.avg <- hclust(sampleDists, method="average")
heatmap.2(sampleDistMatrix, trace="none", Rowv=as.dendrogram(dists.avg),
          margins=c(6,7))
heatmap.2(sampleDistMatrix, trace="none", dendrogram="column",
          scale="row", margins=c(6,7))

#-----------------------------------------------------------------------------------------------------------------------#
## for plotPCA(), >getMethod("plotPCA","DESeqTransform") to see source code

# --<>--<>--<>--<>-- PCA to examine effect of treatment --<>--<>--<>--<>--<>--<>--
plot.treat.data <- plotPCA(rld, intgroup = c("treatment"), returnData=TRUE, n=500) ## uses n most variable genes -- not necessarily sig DEGs!
percentVar <- round(100*attr(plot.treat.data, "percentVar"))
my.colors <- beyonce_palette(94,2,type=c("discrete"))
my.colors <- c("dodgerblue", "orange")
ggplot(plot.treat.data, aes(PC1, PC2, color=treatment)) +
  scale_color_manual(values=c(my.colors)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()
#-----------------------------------------------------------------------------------------------------------------------#


################### ID OF DEGs ############################################
#--<>--<>--<>--<>-- How many and which genes are sig DE? --<>--<>--<>--<>--
sum(res$padj < 0.05, na.rm=T)            
sum(res$padj < 0.05 & abs(res$log2FoldChange) >= 0, na.rm=T)
resSig       <- res[which(res$padj < 0.05),]   ## put all genes w/ padj < 0.05 in resSig
plotMA(resSig, ylim=c(-3,3), xlim=c(1e-2,1e6)) ## examine MA plot 
head(resSig[order(resSig$log2FoldChange),])    ## sort by log2foldchange and show most down-regulated
tail(resSig[order(resSig$log2FoldChange),])    ## sort by log2foldchange and show most up-regulated
write.csv(as.data.frame(resSig), file="p_05_degs.csv")     ## write CSV with all DEGs, padj < 0.05
#-----------------------------------------------------------------------------------------------------------------------#

################### VISUALIZING DEGs #######################################
#--<>--<>--<>--<>-- Heatmap w/ top significant DEGs --<>--<>--<>--<>--

topSigGenes <- head(resSig[order(resSig$padj),],n=10)      ## get 10 genes with lowest padj values
colours <- beyonce_palette(64,300,type=c("continuous"))

#!# genes on y, samples on x
heatmap.2(assay(rld)[rownames(topSigGenes),], 
          scale="row", ## plot heatmap using these genes
          trace="none", dendrogram="column",margins=c(4,10), key.title=NA,
          cexRow=1, cexCol=1,
          col = colours)
#-----------------------------------------------------------------------------------------------------------------------#

#--<>--<>--<>--<>-- PCA plots - DEGs by cut-offs --<>--<>--<>--<>--
p.cutoff <- 0.05
fc.cutoff <- 0

topSigGenes <- res[which(res$padj < p.cutoff & abs(res$log2FoldChange) >= fc.cutoff),]
num.genes <- topSigGenes@nrows
rownames(rld)
rownames(topSigGenes)
rld.Sig <- rld[rownames(rld) %in% rownames(topSigGenes)]

plot.all.data <- plotPCA(rld.Sig, intgroup = c("family", "treatment"), returnData=TRUE, ntop=num.genes)
percentVar <- round(100*attr(plot.all.data, "percentVar"))
my.colors <- beyonce_palette(66,9,type=c("continuous"))



## color PCA by treatment
ggplot(plot.all.data, aes(PC1, PC2, color=treatment, shape=treatment)) +
  ggtitle(paste("Top",num.genes,"DEGs, padj < ",p.cutoff,"and log2FC >",fc.cutoff)) +       
  scale_color_manual(values=c("firebrick4","darkgoldenrod3")) + 
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()
#-----------------------------------------------------------------------------------------------------------------------#

dat  <- plot.all.data
#cols <-colorRampPalette(c("blue", "orange"))(length(unique(dat[, 4])) )
cols <-c("blue", "orange", "purple", "red", "black", "white", "gray", "brown", "yellow", "green", "violet")
cols <- grDevices::hcl.colors(18, "BrBG")
cols <- sample(colorspace::sequential_hcl(18, "viridis"), 18, replace = FALSE)

plot(dat[, 1], dat[, 2], cex = 0.001)
for(n in 1:length(unique(dat[, 4]))){
  f <- unique(dat[, 4])[n] 
  points(dat[which(dat[, 4] == f), 1], dat[which(dat[, 4] == f), 2], pch = 21, bg = cols[n], cex = 2)
  #text(dat[which(dat[, 4] == f), 1], dat[which(dat[, 4] == f), 2], f, cex = 1.2)
}

