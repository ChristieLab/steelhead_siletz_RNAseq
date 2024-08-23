#---------------- Data Description -------------------------------------------------------------------------------------#
# script written by mark christie on 11/25/2019; contact at markchristie1500@gmail.com

#set working directory, load libraries, list files
#setwd("C:/Users/fishf/Dropbox/manuscripts/steelhead_rna-seq/analyses/siletz/differential_expression/limma")
library("DESeq2");
library("rtracklayer");
#library("ggplot2");
#library("gplots");
#library("beyonce");
library(limma); 
library(edgeR); 
#library(Glimma)
list.files()

#---------------- Data loading and prep --------------------------------------------------------------------------------#
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

gene.counts <- gene.counts[, -crosses]


# read in CSV with sample information
col.data <- read.csv("sample_info_no_reciprocal.csv", row.names = 1)

# check to make sure that sample names are in the same order in the gene count table and the sample info table; must be TRUE or do not proceed
all(rownames(col.data) %in% colnames(gene.counts))
all(rownames(col.data) == colnames(gene.counts))



#----------------------- Construction of DESeqDataSet object and DE analysis ---------------------------------------------------


dge <- DGEList(counts=gene.counts)

assembly <- readGFF("stringtie_all_merged.gtf")             ## read in merged GTF
gene_idx <- match(row.names(dge$counts), assembly$gene_id)
gene.name      <- assembly$gene_name[gene_idx]
transcript.id  <- assembly$transcript_id[gene_idx]
xloc           <- assembly$exon_number[gene_idx]
gene_names     <- cbind(gene.name, transcript.id, xloc)
paste(gene_names[,1],gene_names[,2],gene_names[,3], sep="|")
rownames(dge$counts) <- paste(gene_names[,1],gene_names[,2],gene_names[,3], sep="|") # adds unique gene names to dds


#L = mean and M = median cpm values
L <- mean(dge$samples$lib.size) * 1e-6
M <- median(dge$samples$lib.size) * 1e-6
c(L, M)

# how many genes have a count of 0 across all 121 samples
table(rowSums(dge$counts==0)==121)

x <- dge
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)


#filter out lowly expressed genes
#filtering is uneccsearily stringent (see online posts)
x <- dge
keep.exprs <- filterByExpr(x, group=NULL)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

dim(x)#normalization
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

# create design matrix
treatment = col.data[, 1]
machine =  col.data[, 5]
cross.date = col.data[, 3]
sex = col.data[, 4]
family = col.data[, 2]

#design <- model.matrix(~treatment+machine+cross.date+sex)
design <- model.matrix(~treatment+cross.date+sex)
#design <- model.matrix(~treatment+sex)
#design <- model.matrix(~treatment)

colnames(design) <- gsub("treatment", "", colnames(design))
design


# estimate correlation among siblings
v <- voom(x, design)
cor <- duplicateCorrelation(v, design, block = family)
cor$consensus

#final model  - note that random effect (family) and seq machine are confounded due to experimental design

y      <- voom(x,design)
corfit <- duplicateCorrelation(y,design,block=family)
y      <- voom(x,design,block=family,correlation=corfit$consensus)
fit    <- lmFit(y,design,block=family,correlation=corfit$consensus)
fit    <- eBayes(fit)
summary(decideTests(fit))
topTable(fit)
topTable(fit,coef="WxW",number=400,sort.by="p")
degenes <- topTable(fit,coef="WxW",number=400,sort.by="p")
write.csv(as.data.frame(degenes), file="p_05_degs_limma.csv")     ## write CSV with all DEGs, padj < 0.05

glMDPlot(fit, status=dt, counts=vm, groups=groups, side.main="Symbols")


y      <- voom(x,design)
fit    <- lmFit(y,)
fit    <- eBayes(fit)
summary(decideTests(fit))
topTable(fit,coef="WxW",number=1500,sort.by="p")
degenes <- topTable(fit,coef="WxW",number=1500,sort.by="p")
write.csv(as.data.frame(degenes), file="p_05_degs_no_random_effect.csv")     ## write CSV with all DEGs, padj < 0.05


##alternative models
#fit <- lmFit(v, design, block = family, correlation = cor$consensus)
#fit <- lmFit(v, design, correlation = cor$consensus)
#fit <- lmFit(v, design)
#fit <- eBayes(fit)
#summary(decideTests(fit))
#topTable(fit)
#topTable(fit,coef="WxW",number=170,sort.by="p")

