# set working directory
setwd("E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-vf/")

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

# read in data matrix (matrix created from mapped reads using RSEM for met_ge1 transcriptome, assembled from brain, liver, and muscle tissue)
data = read.table("Z:/edgeR-1/gcmatrix/meta1_rsem_gene.counts.matrix", header=T, row.names = 1, com = '')

rnaseqMatrix = data[,c("ct_c_1","ct_c_2","ct_c_3","ct_th_1","ct_th_2","ct_th_3","lc_c_1","lc_c_2","lc_th_1","lc_th_2","lc_th_3","lm_c_1","lm_c_2","lm_c_4","lm_c_5","lm_th_1","lm_th_2")]

# round counts; remove lowly expressed transcripts, and convert to CPM
# remove zeros after decimal point
rnaseqMatrix = round(rnaseqMatrix)
# removal of lowly expressed transcripts AFTER first converting to CPM to account for differences in library size
# removes any genes with fewer than 3 cpm (6 to 7 regular counts) OR those expressed in 3 or fewer libraries
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 3) >= 2,]

# reduce matrix to columns associated with samples from muscle_ge1
names(rnaseqMatrix)


# establish conditions (treatments + populations)
sample_id <- colnames(rnaseqMatrix)

# create vector for treatment type T=treatment; C=control
treatment = factor(c(rep("Con", 3), rep("Trt-H", 3), rep("Con", 2), rep("Trt-H", 3), rep("Con", 4), rep("Trt-H", 2)))

group = factor(c(rep("ct_c", 3), rep("ct_th", 3), rep("lc_c", 2) , rep("lc_th", 3) , rep("lm_c", 4), rep("lm_th", 2)))

seq<-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1)

targets <- data.frame(sample_id, treatment, group, seq)
comp<-targets

comp.group <- comp$group
comp.seq <- comp$seq

# design
design <- model.matrix(~0+comp.seq+comp.group)

# rename column names to make typing out contrasts less tedious
colnames(design) <- gsub("comp.group", "",colnames(design))
# rename column names to make typing out contrasts less tedious
colnames(design) <- gsub("comp.", "",colnames(design))

# create DGEList object
# turn counts matrix into an object of class DGElist; samples grouped based on condition vector established earlier
y = DGEList(counts=rnaseqMatrix, group=comp$group)

# b/c highly expressed genes can consume substantial proportion of lib size for sample
# other genes can falsely appear to be down regulated. This fctn normalizes samples based on RNA comp. by finding set of scaling factors for the library sizes that minimizes log-FC btwn samples for most genes (see edgeR manual pg. 12)
y = calcNormFactors(y)

rownames(design) <- colnames(y)

# plotMDS samples to view how these cluster in multidimensional space
logcpm <- cpm(y, prior.count=2, log=TRUE)
sym = as.numeric(comp$treatment) + 14
colz = c(rep("green", 6), rep("orange", 5), rep("blue", 6))

pdf("E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-vf/plotmds-ge1-hvc-w415-hiseq.pdf", width=10, height=10, useDingbats=FALSE, onefile=TRUE)
plotMDS(logcpm, pch= sym, col=colz)
legend("bottomright", legend = c("ct-c","ct-h","lc-c","lc-h","lm-c","lm-h"),col = c(rep("green",2), rep("orange",2), rep("blue",2)), pch = rep(unique(sym), 3), ncol = 3, cex = 0.5)
dev.off()

# generate pseudo-counts to speed computation of likelihood; estimate dispersion BTWN samples due to biological variation and potential technical variation
y <- estimateGLMRobustDisp(y, design)

# view dispersion estimates
# biological coefficient of biological variation; coefficient of variation with which the (unknown) true abundance of gene varies btwn replicate RNA samples due to biologival causes
bcv <- sqrt(y$common.dispersion)
# gives a plot showing bcv (red line) and tagwise dispersion estimates
plotBCV(y, main = "Dispersion Estimates Across All Populations")

# fit negative binomial generalized log-linear model (GLM) for each read count (i.e. each tag)/gene
fit <- glmFit(y, design)

# use quasi-likelihood F-test instead of GLM; more conservative and better control of type I error
# fit <- glmQLFit(y,design, robust = TRUE)

# make contrasts
# contrasts for comparing c vs t within each population
my.contrasts <- makeContrasts(lm_cvt=lm_th-lm_c, lc_cvt=lc_th-lc_c, ct_cvt=ct_th-ct_c, levels=design)

# contrast for all pops comparing c vs t; divide by three to get average treatment log-fold change across 3 populations
allpops_contrast <- makeContrasts((lm_th + lc_th + ct_th)/3 - (lm_c + lc_c + ct_c)/3, levels=design)

pdf("E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-vf/MDplot-ge1-hvc-w415-hiseq.pdf")

par(mfrow=c(2,2))

# Comparison1: Compare DEGs between treated and untreated lamprey across ALL populations treatment
# p-values determined by combination of log FC (how different was expression btwn treatment groups) and logCPM (a proxy of how highly expressed gene was)
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

# pdf("E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-vf/MDplot_allpops-ge1-hvc-w415-hiseq.pdf", width=10, height=10, useDingbats=FALSE, onefile=TRUE)
# plot genes significantly DE at FDR = 0.05 level
plotMD(lrt, main = "DEGs: All Pops")
abline(h=c(-1, 1), col="blue")
# dev.off()

write.csv(result_table, file='E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-vf/edgeR_CvsT_muscle_allpops-ge1-hvc-w415-hiseq.csv', quote=F, row.names=T)

# perform DGE and give edgeR a specific log Fold change below which we ignore DEGs -- ONLY applicable if GLM test used
# tr <- glmTreat(fit, contrast = allpops_contrast, lfc=1)
# lfc_tTags <- topTags(tr, n=NULL,p.value = 0.05)
# lfc_table <- lfc_tTags$table

# Comparison2: Compare DEGs between treated and untreated lamprey in Lake Mich Pop.
# p-values determined by combination of log FC (how different was expression btwn treatment groups) and logCPM (a proxy of how highly expressed gene was)
lrt_lm <- glmLRT(fit, contrast = my.contrasts[,"lm_cvt"])
# adjust p-val is FDR adjusted p-value using Benjamini-Hochberg
tTags_lm <- topTags(lrt_lm, n=NULL, p.value = 0.05, adjust.method = "BH")
result_table_lm = tTags_lm$table
result_table_lm = data.frame(sampleA="control", sampleB="treatment", result_table_lm)
# multiply by negative 1 to ensure that change 
# result_table_lm$logFC = -1 * result_table_lm$logFC
top_lm <- rownames(topTags(lrt_lm))
cpm(y)[top_lm,]
summary(decideTests(lrt_lm))

# pdf("E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-vf/MDplot_lm-ge1-hvc-w415-hiseq.pdf", width=10, height=10, useDingbats=FALSE, onefile=TRUE)
plotMD(lrt_lm, main = "DEGs: LM") # plot genes significantly DE at FDR = 0.05 level 
abline(h=c(-1, 1), col="blue")
# dev.off()

write.csv(result_table_lm, file='E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-vf/edgeR_CvsT_muscle_LM-ge1-hvc-w415-hiseq.csv', quote=F, row.names=T)

# perform DGE and give edgeR a specific log Fold change below which we ignore DEGs
# tr_lm <- glmTreat(fit, contrast = my.contrasts[,"lm_cvt"], lfc=0.5)
# lfc_tTags_lm <- topTags(tr_lm, n=NULL,p.value = 0.05)
# lfc_table_lm <- lfc_tTags_lm$table

# Comparison3: Compare DEGs between treated and untreated lamprey in Lake Champlain Pop.
# p-values determined by combination of log FC (how different was expression btwn treatment groups) and logCPM (a proxy of how highly expressed gene was)
lrt_lc <- glmLRT(fit, contrast = my.contrasts[,"lc_cvt"])
# adjust p-val is FDR adjusted p-value using Benjamini-Hochberg
tTags_lc <- topTags(lrt_lc, n=NULL, p.value = 0.05, adjust.method = "BH")
result_table_lc = tTags_lc$table
result_table_lc = data.frame(sampleA="control", sampleB="treatment", result_table_lc)
# multiply by negative 1 to ensure that change
# result_table_lc$logFC = -1 * result_table_lc$logFC
top_lc <- rownames(topTags(lrt_lc))
cpm(y)[top_lc,]
summary(decideTests(lrt_lc))

# pdf("E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-vf/MDplot_lc-ge1-hvc-w415-hiseq.pdf", width=10, height=10, useDingbats=FALSE, onefile=TRUE)
plotMD(lrt_lc, main = "DEGs: LC") # plot genes significantly DE at FDR = 0.05 level 
abline(h=c(-1, 1), col="blue")
# dev.off()

write.csv(result_table_lc, file='E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-vf/edgeR_CvsT_muscle_LC-ge1-hvc-w415-hiseq.csv', quote=F, row.names=T)

# perform DGE and give edgeR a specific log Fold change below which we ignore DEGs
# tr_lc <- glmTreat(fit, contrast = my.contrasts[,"lc_cvt"], lfc=0.5)
# lfc_tTags_lc <- topTags(tr_lc, n=NULL,p.value = 0.05)
# lfc_table_lc <- lfc_tTags_lc$table

# Comparison4: Compare DEGs between treated and untreated lamprey in Connecticut Pop.
# p-values determined by combination of log FC (how different was expression btwn treatment groups) and logCPM (a proxy of how highly expressed gene was)
lrt_ct <- glmLRT(fit, contrast = my.contrasts[,"ct_cvt"])
# adjust p-val is FDR adjusted p-value using Benjamini-Hochberg
tTags_ct <- topTags(lrt_ct, n=NULL, p.value = 0.05, adjust.method = "BH")
result_table_ct = tTags_ct$table
result_table_ct = data.frame(sampleA="control", sampleB="treatment", result_table_ct)
# multiply by negative 1 to ensure that change
# result_table_ct$logFC = -1 * result_table_ct$logFC
top_ct <- rownames(topTags(lrt_ct))
cpm(y)[top_ct,]
summary(decideTests(lrt_ct))

# pdf("E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-vf/MDplot_ct-ge1-hvc-w415-hiseq.pdf", width=10, height=10, useDingbats=FALSE, onefile=TRUE)
plotMD(lrt_ct, main = "DEGs: CT") # plot genes significantly DE at FDR = 0.05 level 
abline(h=c(-1, 1), col="blue")
# dev.off()

write.csv(result_table_ct, file='E:/Purdue/pmarinus/TFM/DEGs/todolist-v1/edgeR-vf/edgeR_CvsT_muscle_CT-ge1-hvc-w415-hiseq.csv', quote=F, row.names=T)

# perform DGE and give edgeR a specific log Fold change below which we ignore DEGs
# tr_ct <- glmTreat(fit, contrast = my.contrasts[,"ct_cvt"], lfc=0.5)
# lfc_tTags_ct <- topTags(tr_ct, n=NULL,p.value = 0.05)
# lfc_table_ct <- lfc_tTags_ct$table

dev.off()
