############################################################
## test Hardy-Windberg equilibrium for the Connecticut River
############################################################

install.packages("HWxtest",repos="https://repo.miserver.it.umich.edu/cran/",INSTALL_opts=c('--no-lock'))
library(HWxtest)

data1<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data_lc.csv")

ngeno<-dim(data1)[2]

data1$missing<-rowSums(data1[,5:ngeno]=="./.")
data1$g00<-rowSums(data1[,5:ngeno]=="0/0")
data1$g01<-rowSums(data1[,5:ngeno]=="0/1")
data1$g02<-rowSums(data1[,5:ngeno]=="0/2")
data1$g03<-rowSums(data1[,5:ngeno]=="0/3")
data1$g10<-rowSums(data1[,5:ngeno]=="1/0")
data1$g11<-rowSums(data1[,5:ngeno]=="1/1")
data1$g12<-rowSums(data1[,5:ngeno]=="1/2")
data1$g13<-rowSums(data1[,5:ngeno]=="1/3")
data1$g20<-rowSums(data1[,5:ngeno]=="2/0")
data1$g21<-rowSums(data1[,5:ngeno]=="2/1")
data1$g22<-rowSums(data1[,5:ngeno]=="2/2")
data1$g23<-rowSums(data1[,5:ngeno]=="2/3")
data1$g30<-rowSums(data1[,5:ngeno]=="3/0")
data1$g31<-rowSums(data1[,5:ngeno]=="3/1")
data1$g32<-rowSums(data1[,5:ngeno]=="3/2")
data1$g33<-rowSums(data1[,5:ngeno]=="3/3")
data1$total<-(data1$missing+data1$g00+data1$g01+data1$g02+data1$g03+data1$g10+data1$g11+data1$g12+data1$g13+data1$g20+data1$g21+data1$g22+data1$g23+data1$g30+data1$g31+data1$g32+data1$g33)
data1$G11<-data1$g00
data1$G21<-data1$g10+data1$g01
data1$G22<-data1$g11
data1$G31<-data1$g20+data1$g02
data1$G32<-data1$g21+data1$g12
data1$G33<-data1$g22
data1$G41<-data1$g30+data1$g03
data1$G42<-data1$g31+data1$g13
data1$G43<-data1$g32+data1$g23
data1$G44<-data1$g33
data1$TOTAL<-data1$missing+data1$G11+data1$G21+data1$G22+data1$G31+data1$G32+data1$G33+data1$G41+data1$G42+data1$G43+data1$G44

write.csv(data1,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data1_lc.csv")

####################################################################

data1<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data1_lc.csv")

data1[,c("X","X.1")]<- list(NULL)

data2<-data.frame(data1$locus,data1$G11,data1$G21,data1$G22,data1$G31,data1$G32,data1$G33,data1$G41,data1$G42,data1$G43,data1$G44)
colnames(data2)<-c("locus","G11","G21","G22","G31","G32","G33","G41","G42","G43","G44")

data3<-data.frame(data2,data1$missing)
colnames(data3)<-c("locus","G11","G21","G22","G31","G32","G33","G41","G42","G43","G44","missing")
data3$total<-(data3$G11+data3$G21+data3$G22+data3$G31+data3$G32+data3$G33+data3$G41+data3$G42+data3$G43+data3$G44+data3$missing)

unique(data3$total)

data2$total<-(data2$G11+data2$G21+data2$G22+data2$G31+data2$G32+data2$G33+data2$G41+data2$G42+data2$G43+data2$G44)

data2$check<-rowSums(data2[,c("G11","G21","G22","G31","G32","G33","G41","G42","G43","G44")]==data2$total)

write.csv(data2,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data2-lc.csv")

##############################################################

data2<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data2-lc.csv")

data2[,c("X")]<- list(NULL)

data4<-data2[data2$check==0,]

for (i in 1:dim(data4)[1]){
gc<-c(data4$G11[i],data4$G21[i],data4$G22[i],data4$G31[i],data4$G32[i],data4$G33[i],data4$G41[i],data4$G42[i],data4$G43[i],data4$G44[i])
hwe<-hwx.test(gc)
data4$LLR[i]<-hwe$Pvalues[1]
data4$Prob[i]<-hwe$Pvalues[2]
data4$U[i]<-hwe$Pvalues[3]
data4$Chisq[i]<-hwe$Pvalues[4]
}

write.csv(data4, file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data4-lc.csv")

#####################################################################

data4<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data4-lc.csv")

data4[,c("X")]<- list(NULL)

data2<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data2-lc.csv")
data2[,c("X")]<- list(NULL)
data51<-data2[data2$check==TRUE & (data2$G11!=0 | data2$G22!=0 | data2$G33!=0 | data2$G44!=0),]
data51$LLR<-"na"
data51$Prob<-"na"
data51$U<-"na"
data51$Chisq<-"na"
data52<-data4[!(data4$LLR<0.05),]
data5<-rbind(data51,data52)

## data5<-data4[!(data4$LLR<0.05),]

## write.csv(data51, file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data51-lc-homofix.csv")

write.csv(data5, file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data4-lc-hwe.csv")

data61<-data2[((data2$check==TRUE)&(data2$G21==data2$total | data2$G31==data2$total | data2$G32==data2$total | data2$G41==data2$total | data2$G42==data2$total | data2$G43==data2$total)),]
data61[,c("X")]<-list(NULL)

data62<-data4[(data4$LLR<0.05),]
data62[,c("LLR","Prob","U","Chisq")]<-list(NULL)

data6<-rbind(data61,data62)

write.csv(data6, file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data4-lc-nonhwe.csv")






#################################
## cactulate Hind for ct with hwe
#################################

ct<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data_lc.csv")

cthwe<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data4-lc-hwe.csv")

hweloci_lc<-merge(ct,cthwe,by.x="locus",by.y="locus")

hweloci_lc[,c("X.x","X.y","G11","G21","G22","G31","G32","G33","G41","G42","G43","G44","total","check","LLR","Prob","U","Chisq")] <- list(NULL)

nrow<-dim(hweloci_lc)[1]
ncol<-dim(hweloci_lc)[2]

hwesum_lc<-matrix(,nrow=19,ncol=ncol)

colnames(hwesum_lc)<-colnames(hweloci_lc)
rownames(hwesum_lc)<-c("missing","g00","g01","g02","g03","g10","g11","g12","g13","g20","g21","g22","g23","g30","g31","g32","g33","total","Hind")

for (i in 4:ncol){
hwesum_lc[1,i]<-sum(hweloci_lc[,i]=="./.")
hwesum_lc[2,i]<-sum(hweloci_lc[,i]=="0/0")
hwesum_lc[3,i]<-sum(hweloci_lc[,i]=="0/1")
hwesum_lc[4,i]<-sum(hweloci_lc[,i]=="0/2")
hwesum_lc[5,i]<-sum(hweloci_lc[,i]=="0/3")
hwesum_lc[6,i]<-sum(hweloci_lc[,i]=="1/0")
hwesum_lc[7,i]<-sum(hweloci_lc[,i]=="1/1")
hwesum_lc[8,i]<-sum(hweloci_lc[,i]=="1/2")
hwesum_lc[9,i]<-sum(hweloci_lc[,i]=="1/3")
hwesum_lc[10,i]<-sum(hweloci_lc[,i]=="2/0")
hwesum_lc[11,i]<-sum(hweloci_lc[,i]=="2/1")
hwesum_lc[12,i]<-sum(hweloci_lc[,i]=="2/2")
hwesum_lc[13,i]<-sum(hweloci_lc[,i]=="2/3")
hwesum_lc[14,i]<-sum(hweloci_lc[,i]=="3/0")
hwesum_lc[15,i]<-sum(hweloci_lc[,i]=="3/1")
hwesum_lc[16,i]<-sum(hweloci_lc[,i]=="3/2")
hwesum_lc[17,i]<-sum(hweloci_lc[,i]=="3/3")
hwesum_lc[18,i]<-(hwesum_lc[1,i]+hwesum_lc[2,i]+hwesum_lc[3,i]+hwesum_lc[4,i]+hwesum_lc[5,i]+hwesum_lc[6,i]+hwesum_lc[7,i]+hwesum_lc[8,i]+hwesum_lc[9,i]+hwesum_lc[10,i]+hwesum_lc[11,i]+hwesum_lc[12,i]+hwesum_lc[13,i]+hwesum_lc[14,i]+hwesum_lc[15,i]+hwesum_lc[16,i]+hwesum_lc[17,i])
hwesum_lc[19,i]<-(hwesum_lc[3,i]+hwesum_lc[4,i]+hwesum_lc[5,i]+hwesum_lc[6,i]+hwesum_lc[8,i]+hwesum_lc[9,i]+hwesum_lc[10,i]+hwesum_lc[11,i]+hwesum_lc[13,i]+hwesum_lc[14,i]+hwesum_lc[15,i]+hwesum_lc[16,i])/(hwesum_lc[2,i]+hwesum_lc[3,i]+hwesum_lc[4,i]+hwesum_lc[5,i]+hwesum_lc[6,i]+hwesum_lc[7,i]+hwesum_lc[8,i]+hwesum_lc[9,i]+hwesum_lc[10,i]+hwesum_lc[11,i]+hwesum_lc[12,i]+hwesum_lc[13,i]+hwesum_lc[14,i]+hwesum_lc[15,i]+hwesum_lc[16,i]+hwesum_lc[17,i])
}

write.csv(hwesum_lc,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/sum_hwe_lc.csv")






####################################
## cactulate Hind for ct with nonhwe
####################################

ct<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data_lc.csv")

ctnonhwe<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data4-lc-nonhwe.csv")

nonhweloci_lc<-merge(ct,ctnonhwe,by.x="locus",by.y="locus")

nonhweloci_lc[,c("X.x","X.y","G11","G21","G22","G31","G32","G33","G41","G42","G43","G44","total","check","LLR","Prob","U","Chisq")] <- list(NULL)

nrow<-dim(nonhweloci_lc)[1]
ncol<-dim(nonhweloci_lc)[2]

nonhwesum_lc<-matrix(,nrow=19,ncol=ncol)

colnames(nonhwesum_lc)<-colnames(nonhweloci_lc)
rownames(nonhwesum_lc)<-c("missing","g00","g01","g02","g03","g10","g11","g12","g13","g20","g21","g22","g23","g30","g31","g32","g33","total","Hind")

for (i in 4:ncol){
nonhwesum_lc[1,i]<-sum(nonhweloci_lc[,i]=="./.")
nonhwesum_lc[2,i]<-sum(nonhweloci_lc[,i]=="0/0")
nonhwesum_lc[3,i]<-sum(nonhweloci_lc[,i]=="0/1")
nonhwesum_lc[4,i]<-sum(nonhweloci_lc[,i]=="0/2")
nonhwesum_lc[5,i]<-sum(nonhweloci_lc[,i]=="0/3")
nonhwesum_lc[6,i]<-sum(nonhweloci_lc[,i]=="1/0")
nonhwesum_lc[7,i]<-sum(nonhweloci_lc[,i]=="1/1")
nonhwesum_lc[8,i]<-sum(nonhweloci_lc[,i]=="1/2")
nonhwesum_lc[9,i]<-sum(nonhweloci_lc[,i]=="1/3")
nonhwesum_lc[10,i]<-sum(nonhweloci_lc[,i]=="2/0")
nonhwesum_lc[11,i]<-sum(nonhweloci_lc[,i]=="2/1")
nonhwesum_lc[12,i]<-sum(nonhweloci_lc[,i]=="2/2")
nonhwesum_lc[13,i]<-sum(nonhweloci_lc[,i]=="2/3")
nonhwesum_lc[14,i]<-sum(nonhweloci_lc[,i]=="3/0")
nonhwesum_lc[15,i]<-sum(nonhweloci_lc[,i]=="3/1")
nonhwesum_lc[16,i]<-sum(nonhweloci_lc[,i]=="3/2")
nonhwesum_lc[17,i]<-sum(nonhweloci_lc[,i]=="3/3")
nonhwesum_lc[18,i]<-(nonhwesum_lc[1,i]+nonhwesum_lc[2,i]+nonhwesum_lc[3,i]+nonhwesum_lc[4,i]+nonhwesum_lc[5,i]+nonhwesum_lc[6,i]+nonhwesum_lc[7,i]+nonhwesum_lc[8,i]+nonhwesum_lc[9,i]+nonhwesum_lc[10,i]+nonhwesum_lc[11,i]+nonhwesum_lc[12,i]+nonhwesum_lc[13,i]+nonhwesum_lc[14,i]+nonhwesum_lc[15,i]+nonhwesum_lc[16,i]+nonhwesum_lc[17,i])
nonhwesum_lc[19,i]<-(nonhwesum_lc[3,i]+nonhwesum_lc[4,i]+nonhwesum_lc[5,i]+nonhwesum_lc[6,i]+nonhwesum_lc[8,i]+nonhwesum_lc[9,i]+nonhwesum_lc[10,i]+nonhwesum_lc[11,i]+nonhwesum_lc[13,i]+nonhwesum_lc[14,i]+nonhwesum_lc[15,i]+nonhwesum_lc[16,i])/(nonhwesum_lc[2,i]+nonhwesum_lc[3,i]+nonhwesum_lc[4,i]+nonhwesum_lc[5,i]+nonhwesum_lc[6,i]+nonhwesum_lc[7,i]+nonhwesum_lc[8,i]+nonhwesum_lc[9,i]+nonhwesum_lc[10,i]+nonhwesum_lc[11,i]+nonhwesum_lc[12,i]+nonhwesum_lc[13,i]+nonhwesum_lc[14,i]+nonhwesum_lc[15,i]+nonhwesum_lc[16,i]+nonhwesum_lc[17,i])
}

write.csv(nonhwesum_lc,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/sum_nonhwe_lc.csv")






######################################################################
## calculate allele frequency, Hexp and Hobs for the Connecticut River
######################################################################

data1<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data1_lc.csv")

data1[,c("X","X.1")]<- list(NULL)

data7<-data.frame(data1$CHROM,data1$POS,data1$locus,data1$missing,data1$G11,data1$G21,data1$G22,data1$G31,data1$G32,data1$G33,data1$G41,data1$G42,data1$G43,data1$G44)
colnames(data7)<-c("CHROM","POS","locus","missing","G11","G21","G22","G31","G32","G33","G41","G42","G43","G44")

data7$allele1<-data7$G11+data7$G11+data7$G21+data7$G31+data7$G41
data7$allele2<-data7$G21+data7$G22+data7$G22+data7$G32+data7$G42
data7$allele3<-data7$G31+data7$G32+data7$G33+data7$G33+data7$G43
data7$allele4<-data7$G41+data7$G42+data7$G43+data7$G44+data7$G44

data7$atotal<-data7$allele1+data7$allele2+data7$allele3+data7$allele4+data7$missing+data7$missing

data7$freq1<-(data7$G11+data7$G11+data7$G21+data7$G31+data7$G41)/(2*(data7$G11+data7$G21+data7$G22+data7$G31+data7$G32+data7$G33+data7$G41+data7$G42+data7$G43+data7$G44))
data7$freq2<-(data7$G21+data7$G22+data7$G22+data7$G32+data7$G42)/(2*(data7$G11+data7$G21+data7$G22+data7$G31+data7$G32+data7$G33+data7$G41+data7$G42+data7$G43+data7$G44))
data7$freq3<-(data7$G31+data7$G32+data7$G33+data7$G33+data7$G43)/(2*(data7$G11+data7$G21+data7$G22+data7$G31+data7$G32+data7$G33+data7$G41+data7$G42+data7$G43+data7$G44))
data7$freq4<-(data7$G41+data7$G42+data7$G43+data7$G44+data7$G44)/(2*(data7$G11+data7$G21+data7$G22+data7$G31+data7$G32+data7$G33+data7$G41+data7$G42+data7$G43+data7$G44))

data7$ftotal<-data7$freq1+data7$freq2+data7$freq3+data7$freq4

data7$Hexp<-(1-data7$freq1^2-data7$freq2^2-data7$freq3^2-data7$freq4^2)
data7$Hobs<-(data7$G21+data7$G31+data7$G32+data7$G41+data7$G42+data7$G43)/(data7$G11+data7$G21+data7$G22+data7$G31+data7$G32+data7$G33+data7$G41+data7$G42+data7$G43+data7$G44)

unique(data7$atotal)
unique(data7$ftotal)

write.csv(data7,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fh-lc.csv")
