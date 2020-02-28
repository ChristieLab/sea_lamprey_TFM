###################################################
## check prop of loci genotyped for each individual
###################################################

data0<-read.table(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/GG-muscle-drcom-mindp5-site80indv80.GT.FORMAT", head=TRUE)
data0$locus<-paste(data0$CHROM,"-",data0$POS)

nrow<-dim(data0)[1]
ncol<-dim(data0)[2]
nsample<-42

data_sum<-matrix(,nrow=4,ncol=ncol)
colnames(data_sum)<-colnames(data0)
rownames(data_sum)<-c("missing","others","total","geno_prop")

for (i in (ncol-nsample):(ncol-1)){
data_sum[1,i]<-sum(data0[1:nrow,i]=="./.")
data_sum[2,i]<-sum(data0[1:nrow,i]=="0/0" | data0[1:nrow,i]=="0/1" | data0[1:nrow,i]=="0/2" | data0[1:nrow,i]=="0/3" | data0[1:nrow,i]=="1/0" | data0[1:nrow,i]=="1/1" | data0[1:nrow,i]=="1/2" | data0[1:nrow,i]=="1/3" | data0[1:nrow,i]=="2/0" | data0[1:nrow,i]=="2/1" | data0[1:nrow,i]=="2/2" | data0[1:nrow,i]=="2/3" | data0[1:nrow,i]=="3/0" | data0[1:nrow,i]=="3/1" | data0[1:nrow,i]=="3/2" | data0[1:nrow,i]=="3/3")
data_sum[3,i]<-data_sum[1,i]+data_sum[2,i]
data_sum[4,i]<-(1-data_sum[1,i]/data_sum[3,i])
}

min(data_sum[4,3:(ncol-1)])
max(data_sum[4,3:(ncol-1)])






#################################
## separate samples by population
#################################

data_all<-read.table(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/GG-muscle-drcom-mindp5-site80indv80.GT.FORMAT",head=TRUE)
data_all$locus<-paste(data_all$CHROM,"-",data_all$POS)

data_lm<-data.frame(data_all$CHROM,data_all$POS,data_all$locus,data_all$X403M,data_all$X407M,data_all$X409,data_all$X415M,data_all$X418,data_all$X419M,data_all$X422M,data_all$X423,data_all$X427M,data_all$X567M,data_all$X573M,data_all$X576M,data_all$X577M,data_all$X581M,data_all$X591M)

colnames(data_lm)<-c("CHROM","POS","locus","X403M","X407M","X409","X415M","X418","X419M","X422M","X423","X427M","X567M","X573M","X576M","X577M","X581M","X591M")

data_lc<-data.frame(data_all$CHROM,data_all$POS,data_all$locus,data_all$X400,data_all$X404M,data_all$X405,data_all$X408,data_all$X410M,data_all$X417M,data_all$X420,data_all$X566M,data_all$X571M,data_all$X572M,data_all$X580M,data_all$X587M,data_all$X588M)

colnames(data_lc)<-c("CHROM","POS","locus","X400","X404M","X405","X408","X410M","X417M","X420","X566M","X571M","X572M","X580M","X587M","X588M")

data_ct<-data.frame(data_all$CHROM,data_all$POS,data_all$locus,data_all$X406,data_all$X411M,data_all$X412,data_all$X413,data_all$X414M,data_all$X416M,data_all$X421,data_all$X424,data_all$X429,data_all$X565M,data_all$X569M,data_all$X579M,data_all$X582M,data_all$X583M)

colnames(data_ct)<-c("CHROM","POS","locus","X406","X411M","X412","X413","X414M","X416M","X421","X424","X429","X565M","X569M","X579M","X582M","X583M")

write.csv(data_lm,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data_lm.csv")
write.csv(data_lc,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data_lc.csv")
write.csv(data_ct,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data_ct.csv")
