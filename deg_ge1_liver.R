# set working directory
setwd("E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-liver/")

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

data1 = data[,c("ctc_liv_1","ctc_liv_2","ctc_liv_3","ctt_liv_1","ctt_liv_2","ctt_liv_3","lcc_liv_1","lcc_liv_2","lcc_liv_3","lct_liv_1","lct_liv_2","lct_liv_3","lmc_liv_1","lmc_liv_2","lmc_liv_3","lmt_liv_1","lmt_liv_2","lmt_liv_3","lmt_liv_4")]

# round counts; remove lowly expressed transcripts, and convert to CPM
# remove zeros after decimal point
rnaseqMatrix = round(data1)

# removal of lowly expressed transcripts AFTER first converting to CPM to account for differences in library size
# removes any genes with fewer than 5 cpm (6 to 7 regular counts) OR those expressed in 3 or fewer libraries
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 3) >= 3,]

# reduce matrix to columns associated with samples from liver_ge1
names(rnaseqMatrix)

# establish Conditions (treatments + populations)
sample_id <- colnames(rnaseqMatrix)

# create vector for treatment type T=treatment; C=control
treatment = factor(c(rep("Con", 3), rep("Trt", 3), rep("Con", 3), rep("Trt", 3), rep("Con", 3), rep("Trt", 4)))

# create vector for population
population = factor(c(rep("ct",6), rep("lc",6), rep("lm",7)))

group = factor(c(rep("ct_c", 3), rep("ct_t", 3), rep("lc_c", 3), rep("lc_t", 3), rep("lm_c", 3), rep("lm_t", 4)))

targets <- data.frame(sample_id, treatment, population, group)

# alter command if you want to exclude samples (rows) from comparison
rows.keep <- c(1:19)
comp <- targets[rows.keep,]

comp.group <- comp$group
comp.group = droplevels(comp.group)

barplot(colSums(rnaseqMatrix), main = "Counts of RNAseq fragments aligned to Trinity transcripts per sample", col = c(rep("darkgoldenrod1",6), rep("chocolate1", 6), rep("firebrick3",7)))

# readjust RNA matrix based on samples kept, but they will be columns in rna matrix
rnaseqMatrix <- rnaseqMatrix[,rows.keep]

# design
design <- model.matrix(~0+comp.group)

# rename column names to make typing out contrasts less tedious
colnames(design) <- gsub("comp.group", "",colnames(design))

# rename column names to make typing out contrasts less tedious
# colnames(design) <- gsub("comp.", "",colnames(design))

# create DGEList object
# turn counts matrix into an object of class DGElist; samples grouped based on condition vector established earlier
y = DGEList(counts=rnaseqMatrix, group=comp$group)

# b/c highly expressed genes can consume substantial proportion of lib size for sample other genes can falsely appear to be down regulated. This fctn normalizes samples based on RNA comp. by finding set of scaling factors for the library sizes that minimizes log-FC btwn samples for most genes (see edgeR manual pg. 12)
y = calcNormFactors(y)
rownames(design) <- colnames(y)
# rownames(design) <- gsub("_mus_ge2","", rownames(design))

# plotMDS samples to view how these cluster in multidimensional space
logcpm <- cpm(y, prior.count=2, log=TRUE)
sym = as.numeric(comp$treatment) + 14
colz = c(rep("green", 6), rep("orange",6), rep("blue",7))

pdf("E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-liver/mds_sample_clust-liver.pdf", width=10, height=10, useDingbats=FALSE, onefile=TRUE)
m <- plotMDS(logcpm, pch= sym, col=colz, main = "GE1 Liver MDS Plot using top 500 DEGs")
legend("bottomleft", legend = c("ct-c","ct-t","lc-c","lc-t","lm-c","lm-t"), col = c(rep("green",2), rep("orange", 2), rep("blue",2)), pch = rep(unique(sym), 3), ncol = 3, cex = 0.5)
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
# biological coefficient of biological variation; coefficient of variation with which the (unknown) true abundance of gene varies btwn replicate RNA samples due to biologival causes
bcv <- sqrt(y$common.dispersion)
# gives a plot showing bcv (red line) and tagwise dispersion estimates
plotBCV(y, main = "Dispersion Estimates Across All Populations")

# fit negative binomial generalized log-linear model (GLM) for each read count (i.e. each tag) /gene
fit <- glmFit(y, design)

# use quasi-likelihood F-test instead of GLM; more conservative and better control of type I error
# fit <- glmQLFit(y,design, robust = TRUE)

# make contrasts
# contrasts for comparing c vs t within each population
my.contrasts <- makeContrasts(lm_cvt=(lm_t)-lm_c, lc_cvt=(lc_t)-lc_c, ct_cvt=(ct_t)-ct_c, levels=design)

# contrast for all pops comparing c vs t; divide by three to get average treatment log-fold change across 3 populations
allpops_contrast <- makeContrasts((lm_t + lc_t + ct_t)/3 - (lm_c + lc_c + ct_c)/3,levels=design)

# Comparison1: Compare DEGs between treated and untreated lamprey across ALL populations treatment
# p-values determined by combination of log FC (how different was expression btwn treatment groups) and logCPM (a proxy of how highly expressed gene was)
# lrt <- glmLRT(fit, contrast = allpops_contrast)
lrt <- glmLRT(fit, contrast = allpops_contrast)
# adjust p-val is FDR adjusted p-value using Benjamini-Hochberg
tTags <- topTags(lrt, n=NULL, p.value = 0.05, adjust.method = "BH")
result_table = tTags$table
result_table = data.frame(sampleA="control", sampleB="treatment", result_table)
# multiply by negative 1 to ensure that change 
# result_table$logFC = -1 * result_table$logFC
top <- rownames(topTags(lrt))
cpm(y)[top,]
summary(decideTests(lrt))

pdf("E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-muscle-vf/MDplot_allpops-liver.pdf", width=2.8, height=2.8, useDingbats=FALSE, onefile=TRUE)
# plot genes significantly DE at FDR = 0.05 level
par(mgp=c(0.1,0.05,0),mar=c(0.9,0.9,0.9,0.9))
plotMD(lrt, main=NA, xlab=NA, ylab=NA, bg.cex=0.1,legend=FALSE,cex=0.5,cex.axis=0.5,tck=-0.01)
#axis(2, mgp=c(0,0.2,0),cex.axis=0.5,tck=-0.01)
#axis(1, mgp=c(0,0,0),cex.axis=0.5,tck=-0.01)
abline(h=c(-1, 1), col="blue")
dev.off()

write.csv(result_table, file='E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-muscle-vf/edgeR_CvT_liv_allpops.csv', quote=F, row.names=T)

# perform DGE and give edgeR a specific log Fold change below which we ignore DEGs -- ONLY applicable if GLM test used
# tr <- glmTreat(fit, contrast = allpops_contrast, lfc=1)
# lfc_tTags <- topTags(tr, n=NULL,p.value = 0.05)
# lfc_table <- lfc_tTags$table

# Comparison2: Compare DEGs between treated and untreated lamprey in Lake Mich Pop.
lrt_lm <- glmLRT(fit, contrast = my.contrasts[,"lm_cvt"])
# adjust p-val is FDR adjusted p-value using Benjamini-Hochberg
tTags_lm <- topTags(lrt_lm, n=NULL, p.value = 0.05, adjust.method = "BH")
result_table_lm = tTags_lm$table
result_table_lm = data.frame(sampleA="control", sampleB="treatment", result_table_lm)
top_lm <- rownames(topTags(lrt_lm))
cpm(y)[top_lm,]
summary(decideTests(lrt_lm))

pdf("E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-liver/MDplot_lm-liver.pdf", width=2.8, height=2.8, useDingbats=FALSE, onefile=TRUE)
par(mgp=c(0.1,0.05,0),mar=c(0.9,0.9,0.9,0.9))
# plot genes significantly DE at FDR = 0.05 level
plotMD(lrt_lm, main=NA, xlab=NA, ylab=NA, bg.cex=0.1,legend=FALSE,cex=0.5,cex.axis=0.5,tck=-0.01)
abline(h=c(-1, 1), col="blue")
dev.off()

write.csv(result_table_lm, file='E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-liver/edgeR_CvsT_liv_LM.csv', quote=F, row.names=T)

# Comparison3: Compare DEGs between treated and untreated lamprey in Lake Champlain Pop.
lrt_lc <- glmLRT(fit, contrast = my.contrasts[,"lc_cvt"])
# adjust p-val is FDR adjusted p-value using Benjamini-Hochberg
tTags_lc <- topTags(lrt_lc, n=NULL, p.value = 0.05, adjust.method = "BH")
result_table_lc = tTags_lc$table
result_table_lc = data.frame(sampleA="control", sampleB="treatment", result_table_lc)
top_lc <- rownames(topTags(lrt_lc))
cpm(y)[top_lc,]
summary(decideTests(lrt_lc))

pdf("E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-liver/MDplot_lc-liver.pdf", width=2.8, height=2.8, useDingbats=FALSE, onefile=TRUE)
par(mgp=c(0.1,0.05,0),mar=c(0.9,0.9,0.9,0.9))
# plot genes significantly DE at FDR = 0.05 level
plotMD(lrt_lc, main=NA, xlab=NA, ylab=NA, bg.cex=0.1,legend=FALSE,cex=0.5,cex.axis=0.5,tck=-0.01)
abline(h=c(-1, 1), col="blue")
dev.off()

write.csv(result_table_lc, file='E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-liver/edgeR_CvsT_liv_LC.csv', quote=F, row.names=T)

# Comparison4: Compare DEGs between treated and untreated lamprey in Connecticut Pop.
lrt_ct <- glmLRT(fit, contrast = my.contrasts[,"ct_cvt"])
# adjust p-val is FDR adjusted p-value using Benjamini-Hochberg
tTags_ct <- topTags(lrt_ct, n=NULL, p.value = 0.05, adjust.method = "BH")
result_table_ct = tTags_ct$table
result_table_ct = data.frame(sampleA="control", sampleB="treatment", result_table_ct)
top_ct <- rownames(topTags(lrt_ct))
cpm(y)[top_ct,]
summary(decideTests(lrt_ct))

pdf("E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-liver/MDplot_ct-liver.pdf", width=2.8, height=2.8, useDingbats=FALSE, onefile=TRUE)
par(mgp=c(0.1,0.05,0),mar=c(0.9,0.9,0.9,0.9))
# plot genes significantly DE at FDR = 0.05 level 
plotMD(lrt_ct, main=NA, xlab=NA, ylab=NA, bg.cex=0.1,legend=FALSE,cex=0.5,cex.axis=0.5,tck=-0.01)
abline(h=c(-1, 1), col="blue")
dev.off()

write.csv(result_table_ct, file='E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-liver/edgeR_CvsT_liv_CT.csv', quote=F, row.names=T)
