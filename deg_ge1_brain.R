# set working directory
setwd("E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-brain/")

# install requisite packages and add to library
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager", repos="https://repo.miserver.it.umich.edu/cran/")
BiocManager::install(c('edgeR','limma','DESeq2','ctc','Biobase','Glimma'), ask=FALSE)

install.packages('gplots', repos="https://repo.miserver.it.umich.edu/cran/");
install.packages('ape', repos="https://repo.miserver.it.umich.edu/cran/"); 
install.packages('heatmap3', repos="https://repo.miserver.it.umich.edu/cran/")

library(edgeR); 
library(limma); 
library(DESeq2); 
library(ctc); 
library(Biobase); 
library(gplots); 
library(ape); 
library(Glimma); 
library(heatmap3); 
library(colorspace)

# read in data matrix
data = read.table("Z:/edgeR-1/gcmatrix/meta1_rsem_gene.counts.matrix", header=T, row.names=1, com='')

data1 = data[,c("lmc_brn_1","lmc_brn_2","lmc_brn_3","lmc_brn_4","lmt_brn_1","lmt_brn_2","lmt_brn_3","lmt_brn_4")]

# round counts; remove lowly expressed transcripts, and convert to CPM
# remove zeros after decimal point
rnaseqMatrix = round(data1)

# removal of lowly expressed transcripts AFTER first converting to CPM to account for differences in library size
# removes any genes with fewer than 5 cpm (6 to 7 regular counts) OR those expressed in 3 or fewer libraries

rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 3) >= 4,] 

# reduce matricx to columns associated with samples from brain_ge1
names(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[,1:8]

# plot counts per sample
colnames(rnaseqMatrix) <- gsub("_brn","", colnames(rnaseqMatrix))
barplot(colSums(rnaseqMatrix), main = "Counts of RNAseq fragments aligned to Trinity transcripts per sample", col = "blue")

# establish Conditions (treatments + populations)
sample_id <- colnames(rnaseqMatrix)

# create vector for treatment type T=treatment; C=control
treatment = factor(c(rep("Con", 4), rep("Trt", 4)))
group = factor(c(rep("lm_c", 4), rep("lm_t", 4)))
targets <- data.frame(sample_id, treatment, group)

# design
design <- model.matrix(~0+group)

# rename column names to make typing out contrasts less tedious
colnames(design) <- gsub("group", "",colnames(design))

# create DGEList object
# turn counts matrix into an object of class DGElist; samples grouped based on condition vector established earlier
y = DGEList(counts=rnaseqMatrix, group=group)

# b/c highly expressed genes can consume substantial proportion of lib size for sample other genes can falsely appear to be down regulated. This fctn normalizes samples based on RNA comp. by finding set of scaling factors for the library sizes that minimizes log-FC btwn samples for most genes (see edgeR manual)
y = calcNormFactors(y)
rownames(design) <- colnames(y)
# rownames(design) <- gsub("_mus_ge2","", rownames(design))

# plotMDS samples to view how these cluster in multidimensional space
logcpm <- cpm(y, prior.count=2, log=TRUE)
sym = as.numeric(treatment) + 14
colz = rep("blue",8)
pdf("E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-brain/mds_sample_clust-brain.pdf", width=10, height=10, useDingbats=FALSE, onefile=TRUE)
# using the top 500 differentially expressed genes by default; change 'top' param to change this
plotMDS(logcpm, pch= sym, col=colz, main = "GE1 Brain MDS Plot using top 500 DEGs", cex = 1)
legend("bottomright", legend = c("lm-c","lm-t"), col = colz, pch = unique(sym), ncol =1, cex = 1)
dev.off()

# Glimma interactive MDS plot to see how samples cluster
glMDSPlot(y, groups=group)

# generate pseudo-counts to speed computation of likelihood; estimate dispersion BTWN samples due to biological variation and potential technical variation
y <- estimateGLMRobustDisp(y, design)
# y <- estimateDisp(y, design, robust = T)
# these three lines generate three different types of dispersion generated in previous line in separate steps
# y <- estimateGLMCommonDisp(y,design)
# y <- estimateGLMTrendedDisp(y,design)
# y <- estimateGLMTagwiseDisp(y,design)

# view Dispersion estimates
# biological coefficient of biological variation; invalid is estimateGLMRobustDisp used instead of estimateDisp
bcv <- sqrt(y$common.dispersion)
# gives a plot showing bcv (red line) and tagwise dispersion estimates
plotBCV(y, main = "Dispersion Estimates Across All Populations")

# fit negative binomial generalized log-linear model (GLM) for each read count (i.e. each tag) /gene
fit <- glmFit(y, design)
# use quasi-likelihood F-test instead of GLM; more conservative and better control of type I error
# fit <- glmQLFit(y,design, robust = TRUE)

# make contrasts
# contrasts for comparing c vs t within each population
my.contrasts <- makeContrasts(lm_cvt=lm_t-lm_c, levels=design)

### Comparison: Compare DEGs between treated and untreated lamprey in Lake Mich Pop.
# p-values determined by combination of log FC (how different was expression btwn treatment groups) and logCPM (a proxy of how highly expressed gene was)
lrt_lm <- glmLRT(fit, contrast = my.contrasts[,"lm_cvt"])
# adjust p-val is FDR adjusted p-value using Benjamini-Hochberg
tTags_lm <- topTags(lrt_lm, n=NULL, p.value = 0.05, adjust.method = "BH")
result_table = tTags_lm$table
result_table = data.frame(sampleA="control", sampleB="treatment", result_table)
# multiply by negative 1 to ensure that change 
# result_table_lm$logFC = -1 * result_table_lm$logFC
top <- rownames(topTags(lrt_lm))
cpm(y)[top,]
summary(decideTests(lrt_lm))

pdf("E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-brain/MDplot_lm.pdf", width=2.8, height=2.8, useDingbats=FALSE, onefile=TRUE)
# jpeg("E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-brain/MDplot_lm-brain.jpeg")
# plot genes significantly DE at FDR = 0.05 level
par(mgp=c(0.1,0.05,0),mar=c(0.9,0.9,0.9,0.9))
plotMD(lrt_lm, main=NA, xlab=NA, ylab=NA, bg.cex=0.1,legend=FALSE,cex=0.5,cex.axis=0.5,tck=-0.01)
#axis(2, mgp=c(0,0.2,0),cex.axis=0.5,tck=-0.01)
#axis(1, mgp=c(0,0,0),cex.axis=0.5,tck=-0.01)
abline(h=c(-1, 1), col="blue")
dev.off()

write.csv(result_table, file='E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-brain/edgeR_CvsT_brain_LM.csv', quote=F, row.names=T)
