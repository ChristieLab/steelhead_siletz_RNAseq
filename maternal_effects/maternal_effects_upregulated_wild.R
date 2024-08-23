#---------------- Data Description -------------------------------------------------------------------------------------#
# script written by Mark Christie on June 5, 2024; contact at markchristie1500@gmail.com
# Rversion: 4.4.0
# all plots exported as pdf at 5.98 x 5.33 dimensions


#set working directory, load libraries, list files
#setwd("C:/Users/fishf/Dropbox/manuscripts/steelhead_gene_expression/differential_expression/hood_2010_all")
library("DESeq2");library("rtracklayer");library("ggplot2");library("gplots");library("beyonce")
list.files()

#---------------- Data loading and prep --------------------------------------------------------------------------------#
# read in count data generated with featureCounts
gene.counts <- read.table("featureCount_counts.txt", sep="\t", header=T, row.names=1, check.names=F)
drop.vars   <- names(gene.counts) %in% c("Chr","Start", "End", "Strand", "Length")
gene.counts <- gene.counts[!drop.vars]
head(gene.counts)


#---------------- Sample processing ------------------------------------------------------------------------------------#
# remove HxW and WxH offspring to test for main effects (HxH vs WxW)


# read in CSV with sample information
col.data <- read.csv("sample_info_all_including_reciprocal.csv", row.names = 1)
samples  <- read.csv("sample_info_all_including_reciprocal.csv")


m1 <- which(samples[, 2] == "HxW")
m2 <- which(samples[, 2] == "WxH")
crosses <- c(m1, m2)
samples[crosses, ] # should all be reciprocal

# with sample 22 excluded (outlier) and all HxW and WxH offspring excluded
gene.counts <- gene.counts[, -crosses]
col.data    <- col.data[-crosses, ]
samples     <- samples[-crosses, ]


# check to make sure that sample names are in the same order in the gene count table and the sample info table; must be TRUE or do not proceed
all(samples[, 1] %in% colnames(gene.counts))
all(samples[, 1] == colnames(gene.counts))


#----------------------- Construction of DESeqDataSet object and DE analysis ---------------------------------------------------
# test for effects of treatment while controlling for the effect of family; dds = DESeqDataSet
# variable order matters! Covariates go first, variable of interest goes last in design 


# testing for effect of HxH vs WxW after controlling for cross date , sex of fry, and sequencing run (seq_group)
# sex has little effect


dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = col.data,
                              design = ~ as.factor(cross_date) + as.factor(sex)  + cross_type)


dds$treatment <- relevel(dds$cross_type, "WxW") ## set the untreated group as the control
dds
dds <- DESeq(dds)

#----------------- Extract and examine results from DESeq analysis -----------------------------------------------------#
# isolate DEGs based on counts

res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
length(which(res[, 6] < 0.05))
summary(res)

genes <- res[which(res[, 6] < 0.05), ] # genes and their names and p-vals
cts   <- as.data.frame(counts( dds, normalized=TRUE )) # normalized counts

gnames <- rownames(genes)
m1     <- match(rownames(genes), rownames(cts)) 

decounts <- cts[m1, ]  # counts of DEGs!

# checks
rownames(decounts) == rownames(genes)

# split counts into HH and WW

# only from same matrices (else comment out) ===================================================#
# need them in exact order, hence lack of better indexing
s1  <- which(col.data$family =="1A")
s2  <- which(col.data$family =="1C")

s3  <- which(col.data$family =="2A")
s4  <- which(col.data$family =="2C")

s5  <- which(col.data$family =="3A")
s6  <- which(col.data$family =="3C")

s7  <- which(col.data$family =="4A")
s8  <- which(col.data$family =="4C")

s9  <- which(col.data$family =="5A")
s10  <- which(col.data$family =="5C")

s11 <- which(col.data$family =="6A")
s12  <- which(col.data$family =="6C")

s13  <- which(col.data$family =="7A")
s14  <- which(col.data$family =="7C")

s15  <- which(col.data$family =="8A")
s16  <- which(col.data$family =="8C")

s17  <- which(col.data$family =="9A")
s18  <- which(col.data$family =="9C")

s19  <- which(col.data$family =="10A")
s20  <- which(col.data$family =="10C")

s21  <- which(col.data$family =="11A")
s22  <- which(col.data$family =="11C")

s23  <- which(col.data$family =="12A")
s24  <- which(col.data$family =="12C")

s25  <- which(col.data$family =="13A")
s26  <- which(col.data$family =="13C")

s27  <- which(col.data$family =="14A")
s28  <- which(col.data$family =="14C")

s29  <- which(col.data$family =="15A")
s30  <- which(col.data$family =="15C")


wilds <- rownames(col.data[c(s1, s3, s5, s7, s9, s11, s13, s15, s17, s19, s21, s23, s25, s27, s29), ])
hatch <- rownames(col.data[c(s2, s4, s6, s8, s10, s12, s14, s16, s18, s20, s22, s24, s26, s28, s30), ])
#===========================================================================#

#UNCOMMENT IF WANT GENES FROM ALL HH and WWW============================#
#samps  <- c(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10)

#wilds <- which(col.data$condition == "WxW")
#wilds <- rownames(col.data[wilds, ])

#hatch <- which(col.data$condition == "HxH")
#hatch <- rownames(col.data[hatch, ])
#=======================================================================#

m1    <- match(wilds, colnames(decounts))
wilds <- decounts[, m1] 

m1     <- match(hatch, colnames(decounts))
hatch  <- decounts[, m1] 

# begin counts per gene
nrow(hatch)
nrow(wilds)

hmeans <- apply(hatch, 1, mean)
wmeans <- apply(wilds, 1, mean)

# upregulated in wild =======================================================
#up   <-  which(genes$log2FoldChange > 0)
up   <-  which(genes$log2FoldChange > 0)
hdwn <-  hmeans[up]
wup  <-  wmeans[up]

MEANUPWILD <- mean(wup)  # used for reciprocal crosses
# figure
diffs    <- cbind(wup, hdwn)
#meandiff <- mean(wup - hdwn) # mean difference between wild and hatchery genes for those that are upregulated in wild fish
#MEANDIFF <- mean(wup - hdwn)
#mdiff    <- (wup-hdwn)/(meandiff) # standardized mean difference

#y1       <- meandiff + mdiff
#y2       <- mdiff

#Y1 <- y1


#x1 <- rep(1, length(y1))
#x2 <- rep(4, length(y2))
#library('scales')
#plot(-100, -100, xlim = c(0, 5), ylim = c(0, 700))

#upl  <- length(y1) - (length(y1)*0.975)
#pts <- y1
#down <-  sort(y1)[upl]  
#up   <-  sort(y1, decreasing = TRUE)[upl]  
#pts <- y1[(y1 > down) & (y1 < up)]
#points(jitter(rep(1, length(pts))), pts, pch=19, cex=0.8, col=alpha("deepskyblue", 0.4))
#points(jitter(rep(1, length(pts)), factor = 18), pts, pch=19, cex=1.8, col=alpha("deepskyblue", 0.4))


#upl  <- length(y2) - (length(y2)*0.975)
#down <-  sort(y2)[upl]  
#up   <-  sort(y2, decreasing = TRUE)[upl] 
#pts <- y2[(y2 > down) & (y2 < up)]
##pts <- y2
#points(jitter(rep(4, length(pts)), factor = 4), pts, pch=19, cex=1.8, col=alpha("deepskyblue", 0.4))

#points(1, mean(y1), cex=2, pch=19, col="dodgerblue")
#points(4, mean(y2), cex=2, pch=19, col="dodgerblue")
#meany1 <- mean(y1)


#==load in reciprocal DEGS=====================================================================#
# read in count data generated with featureCounts
gene.counts <- read.table("featureCount_counts.txt", sep="\t", header=T, row.names=1, check.names=F)
drop.vars   <- names(gene.counts) %in% c("Chr","Start", "End", "Strand", "Length")
gene.counts <- gene.counts[!drop.vars]
head(gene.counts)

#---------------- Sample processing ------------------------------------------------------------------------------------#

# read in CSV with sample information
col.data <- read.csv("sample_info_all_including_reciprocal.csv", row.names = 1)
# check to make sure that sample names are in the same order in the gene count table and the sample info table; must be TRUE or do not proceed
all(rownames(col.data) %in% colnames(gene.counts))
all(rownames(col.data) == colnames(gene.counts))

#==================================================================#
samps1  <- which(col.data$cross_type =="WxH") #WFxHM
samps2  <- which(col.data$cross_type =="HxW")
samps  <- c(samps1, samps2)


col.data    <- col.data[samps, ]
table(col.data$cross_type)

gene.counts <- gene.counts[, samps]

all(rownames(col.data) %in% colnames(gene.counts))
all(rownames(col.data) == colnames(gene.counts))

#----------------------- Construction of DESeqDataSet object and DE analysis ---------------------------------------------------
# test for effects of treatment while controlling for the effect of family; dds = DESeqDataSet
# variable order matters! Covariates go first, variable of interest goes last in design 


dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = col.data,
                              design = ~ as.factor(cross_date) + as.factor(sex)  + cross_type)

dds <- DESeq(dds)
#----------------- Extract and examine results from DESeq analysis -----------------------------------------------------#
res2    <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
ndegs2  <- length(which(res2[, 6] < 0.05))

#genes <- res[which(res[, 6] < 0.05), ] # genes and their names and p-vals
cts      <- as.data.frame(counts( dds, normalized=TRUE )) # normalized counts
gnames   <- rownames(genes)
m1       <- match(rownames(genes), rownames(cts)) 
decounts <- cts[m1, ]  # counts of DEGs!

# checks
rownames(decounts) == rownames(genes)

# split counts into HH and WW
wh <- which(col.data$cross_type == "WxH")
wh <- rownames(col.data[wh, ])

hw <- which(col.data$cross_type == "HxW")
hw <- rownames(col.data[hw, ])

m1    <- match(wh, colnames(decounts))
wh    <- decounts[, m1] 

m1     <- match(hw, colnames(decounts))
hw     <- decounts[, m1] 

# begin counts per gene
nrow(wh)
nrow(hw)

whmeans <- apply(wh, 1, mean)
hwmeans <- apply(hw, 1, mean)

# upregulated =======================================================
up    <-  which(genes$log2FoldChange > 0)
hwdwn <-  hwmeans[up]
whup  <-  whmeans[up]
mean(hwdwn)

diffs <- cbind(diffs, whup, hwdwn)

# plot each gene separately; useful don't delete code! commented out till needed.
#for(n in 1:nrow(diffs)){
#  dat <- diffs[n, ]
#  plot(c(1,4,3,2), dat, main = n)
#  
#}


# standardize to log scale
diffs2 <- log(diffs)
HH <- mean(diffs2[, 1])
HW <- mean(diffs2[, 3])
WH <- mean(diffs2[, 4])
WW <- mean(diffs2[, 2])

#HH <- median(diffs2[, 1])
#HW <- median(diffs2[, 3])
#WH <- median(diffs2[, 4])
#WW <- median(diffs2[, 2])

#####==========Hypo prediciton plot; not real data ========================================#

plot(c(1,2,3,4), c(HH, HH, WW, WW), bty = "l", xaxt='n', ylab = "Standardized log gene counts", xlab="", xlim = c(0.75, 4.25), ylim = c(5.25, 5.73))
perfect <- (HH+WW)/2
lines(c(1,2,3,4), c(HH, perfect, perfect, WW), lwd=2, col="grey", lty=2)
points(c(1,2,3,4),  c(HH, perfect, perfect, WW), cex=2.5, pch=21, bg=c("grey70"))
axis(side=1, at=c(1,2,3,4), labels = c("HH", "HW", "WH", "WW"))
lines(c(1,2,3,4), c(HH, HH, WW, WW), lwd=2, col="grey40", lty=2)
points(c(1,2,3,4),  c(HH, HH, WW, WW), cex=2.5, pch=21, bg=c("grey40"))

#####=======================================================================================#


d1 <- mean(diffs2[, 1] - diffs2[, 3])
plot(c(1,2,3,4), c(HH, HW, WH, WW), bty = "l", xaxt='n', ylab = "Standardized log gene counts", xlab="", xlim = c(0.75, 4.25), ylim = c(5.25, 5.73))
lines(c(1,2,3,4), c(HH, HW, WH, WW), lwd=2, col="grey", lty=2)
points(c(1,2,3,4), c(HH, HW, WH, WW), cex=2.5, pch=21, bg=c("orange", "forestgreen", "forestgreen", "dodgerblue"))
axis(side=1, at=c(1,2,3,4), labels = c("HH", "HW", "WH", "WW"))

# standardize to mean of HH gene count
xbar   <- mean(diffs2[, 1]) # mean of first column
diffs3 <- NULL

x      <- xbar - diffs2[, 1]
diffs3 <- cbind(diffs3, diffs2[, 1] + x)
diffs3 <- cbind(diffs3, diffs2[, 2] + x)
diffs3 <- cbind(diffs3, diffs2[, 3] + x)
diffs3 <- cbind(diffs3, diffs2[, 4] + x)

#reorder by HH, HW, WH and WW and add col names
# original 1=HH,2=WW,  3=WH 4=HW
diffs4 <- cbind(diffs3[, 1], diffs3[, 4], diffs3[, 3], diffs3[, 2])
colnames(diffs4) <- c("HH", "HW", "WH", "WW")
head(diffs4)

plot(rep(1, nrow(diffs4)), diffs4[, 1], bty="l", ylab="Standardized log gene counts", pch=19, col = alpha("blue", 0.1), cex=1.2,  xaxt='n', xlab="", xlim = c(0.75, 4.25), ylim = c(5.25, 5.73))
#plot(rep(1, nrow(diffs4)), diffs4[, 1], bty="l", ylab="Standardized log gene counts", pch=19, col = alpha("blue", 0.1), cex=1.2,  xaxt='n', xlab="", xlim = c(0.75, 4.25), ylim = c(5.25, 6.73))
for(n in 1:nrow(diffs4)){
  dat <- diffs4[n, ]
  lines(c(1,2), c(dat[1], dat[2]), col=alpha("blue", 0.2))
  lines(c(2,3), c(dat[2], dat[3]), col = alpha("gray30", 0.2))
  #lines(c(2,3), c(dat[2], dat[3]), col = alpha("blue", 0.2))
  lines(c(3,4), c(dat[3], dat[4]), col=alpha("blue", 0.2))
  
}

points(rep(2, nrow(diffs4)), diffs4[, 2], pch=19, col = alpha("blue", 0.1), cex=1.2)
points(rep(3, nrow(diffs4)), diffs4[, 3], pch=19, col = alpha("blue", 0.1), cex=1.2)
points(rep(4, nrow(diffs4)), diffs4[, 4], pch=19, col = alpha("blue", 0.1), cex=1.2)
#abline(h=median(diffs4[, 2]), col = "red", lwd=2)
axis(side=1, at=c(1,2,3,4), labels = c("HH", "HW", "WH", "WW"))


#boxplot2(diffs4, ylim = c(5.8, 6.5))

# stastical test 
d1 <- diffs4[, 1]-diffs4[, 2]
d2 <- diffs4[, 2]-diffs4[, 3]
d3 <- diffs4[, 3]-diffs4[, 4]

d1 <- diffs4[, 2]-diffs4[, 1]
d2 <- diffs4[, 3]-diffs4[, 2]
d3 <- diffs4[, 4]-diffs4[, 3]
dall <- cbind(d1, d2, d3)

## to conservative, what if increases between HW and WH - not accounted for yet!
OUT <- NULL
for(n in 1:nrow(dall)){
  rown <- dall[n, ]
  out  <- order(rown)  
  OUT <- rbind(OUT, out)  
}

table(OUT[, 1]) # what is lowest value; for 2 = smallest change between HW and WH
table(OUT[, 2]) 
table(OUT[, 3]) # what was the highest change

stat <- 69/95

# start permutation

OUT <- NULL
for(n in 1:10000){ # number of replicate data sets
    
  OUT1 <- NULL
  for(n in 1:95){ # number of genes
    rown <- sample(dall, 3, replace = FALSE)
    out  <- order(rown)  
    OUT1 <- rbind(OUT1, out)  
  }
  
  out2 <- length(which(OUT1[, 1] == 2))
  OUT  <- rbind(OUT, out2)
}

range(OUT)
hist(OUT, breaks = 30, xlim = c(0, 80), col="plum2", xlab="Smallest difference in gene counts HW vs WH",border="plum3", main="")
abline(v=69, col="blue", lty=2, lwd=2)
