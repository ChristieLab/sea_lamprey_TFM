################################################################
## calculate Fst between Lake Michigan and the Connecticut River
################################################################

install.packages("splitstackshape",repos="https://repo.miserver.it.umich.edu/cran/",INSTALL_opts=c('--no-lock'))
library(splitstackshape)

flm<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fh-lm.csv")
fct<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fh-ct.csv")

flmtemp<-data.frame(flm$locus,flm$G11,flm$G21,flm$G22,flm$G31,flm$G32,flm$G33,flm$G41,flm$G42,flm$G43,flm$G44,flm$freq1,flm$freq2,flm$freq3,flm$freq4)
fcttemp<-data.frame(fct$locus,fct$G11,fct$G21,fct$G22,fct$G31,fct$G32,fct$G33,fct$G41,fct$G42,fct$G43,fct$G44,fct$freq1,fct$freq2,fct$freq3,fct$freq4)

flmct<-merge(flmtemp,fcttemp,by.x="flm.locus",by.y="fct.locus")



flmct$nlm<-(flmct$flm.G11+flmct$flm.G21+flmct$flm.G22+flmct$flm.G31+flmct$flm.G32+flmct$flm.G33+flmct$flm.G41+flmct$flm.G42+flmct$flm.G43+flmct$flm.G44)
flmct$nct<-(flmct$fct.G11+flmct$fct.G21+flmct$fct.G22+flmct$fct.G31+flmct$fct.G32+flmct$fct.G33+flmct$fct.G41+flmct$fct.G42+flmct$fct.G43+flmct$fct.G44)

r<-2

flmct$nbar<-(flmct$nlm+flmct$nct)/r
flmct$nc<-(r*flmct$nbar-((flmct$nlm^2+flmct$nct^2)/(r*flmct$nbar)))/(r-1)

flmct$pbar1<-(flmct$nlm*flmct$flm.freq1+flmct$nct*flmct$fct.freq1)/(r*flmct$nbar)
flmct$pbar2<-(flmct$nlm*flmct$flm.freq2+flmct$nct*flmct$fct.freq2)/(r*flmct$nbar)
flmct$pbar3<-(flmct$nlm*flmct$flm.freq3+flmct$nct*flmct$fct.freq3)/(r*flmct$nbar)
flmct$pbar4<-(flmct$nlm*flmct$flm.freq4+flmct$nct*flmct$fct.freq4)/(r*flmct$nbar)

flmct$ssq1<-(flmct$nlm*(flmct$flm.freq1-flmct$pbar1)^2+flmct$nct*(flmct$fct.freq1-flmct$pbar1)^2)/((r-1)*flmct$nbar)
flmct$ssq2<-(flmct$nlm*(flmct$flm.freq2-flmct$pbar2)^2+flmct$nct*(flmct$fct.freq2-flmct$pbar2)^2)/((r-1)*flmct$nbar)
flmct$ssq3<-(flmct$nlm*(flmct$flm.freq3-flmct$pbar3)^2+flmct$nct*(flmct$fct.freq3-flmct$pbar3)^2)/((r-1)*flmct$nbar)
flmct$ssq4<-(flmct$nlm*(flmct$flm.freq4-flmct$pbar4)^2+flmct$nct*(flmct$fct.freq4-flmct$pbar4)^2)/((r-1)*flmct$nbar)

flmct$flm.h1<-(flmct$flm.G21+flmct$flm.G31+flmct$flm.G41)/(flmct$flm.G11+flmct$flm.G21+flmct$flm.G22+flmct$flm.G31+flmct$flm.G32+flmct$flm.G33+flmct$flm.G41+flmct$flm.G42+flmct$flm.G43+flmct$flm.G44)
flmct$fct.h1<-(flmct$fct.G21+flmct$fct.G31+flmct$fct.G41)/(flmct$fct.G11+flmct$fct.G21+flmct$fct.G22+flmct$fct.G31+flmct$fct.G32+flmct$fct.G33+flmct$fct.G41+flmct$fct.G42+flmct$fct.G43+flmct$fct.G44)

flmct$flm.h2<-(flmct$flm.G21+flmct$flm.G32+flmct$flm.G42)/(flmct$flm.G11+flmct$flm.G21+flmct$flm.G22+flmct$flm.G31+flmct$flm.G32+flmct$flm.G33+flmct$flm.G41+flmct$flm.G42+flmct$flm.G43+flmct$flm.G44)
flmct$fct.h2<-(flmct$fct.G21+flmct$fct.G32+flmct$fct.G42)/(flmct$fct.G11+flmct$fct.G21+flmct$fct.G22+flmct$fct.G31+flmct$fct.G32+flmct$fct.G33+flmct$fct.G41+flmct$fct.G42+flmct$fct.G43+flmct$fct.G44)

flmct$flm.h3<-(flmct$flm.G31+flmct$flm.G32+flmct$flm.G43)/(flmct$flm.G11+flmct$flm.G21+flmct$flm.G22+flmct$flm.G31+flmct$flm.G32+flmct$flm.G33+flmct$flm.G41+flmct$flm.G42+flmct$flm.G43+flmct$flm.G44)
flmct$fct.h3<-(flmct$fct.G31+flmct$fct.G32+flmct$fct.G43)/(flmct$fct.G11+flmct$fct.G21+flmct$fct.G22+flmct$fct.G31+flmct$fct.G32+flmct$fct.G33+flmct$fct.G41+flmct$fct.G42+flmct$fct.G43+flmct$fct.G44)

flmct$flm.h4<-(flmct$flm.G41+flmct$flm.G42+flmct$flm.G43)/(flmct$flm.G11+flmct$flm.G21+flmct$flm.G22+flmct$flm.G31+flmct$flm.G32+flmct$flm.G33+flmct$flm.G41+flmct$flm.G42+flmct$flm.G43+flmct$flm.G44)
flmct$fct.h4<-(flmct$fct.G41+flmct$fct.G42+flmct$fct.G43)/(flmct$fct.G11+flmct$fct.G21+flmct$fct.G22+flmct$fct.G31+flmct$fct.G32+flmct$fct.G33+flmct$fct.G41+flmct$fct.G42+flmct$fct.G43+flmct$fct.G44)

flmct$hbar1<-(flmct$nlm*flmct$flm.h1+flmct$nct*flmct$fct.h1)/(r*flmct$nbar)
flmct$hbar2<-(flmct$nlm*flmct$flm.h2+flmct$nct*flmct$fct.h2)/(r*flmct$nbar)
flmct$hbar3<-(flmct$nlm*flmct$flm.h3+flmct$nct*flmct$fct.h3)/(r*flmct$nbar)
flmct$hbar4<-(flmct$nlm*flmct$flm.h4+flmct$nct*flmct$fct.h4)/(r*flmct$nbar)

flmct$a1<-(flmct$nbar/flmct$nc)*(flmct$ssq1-(1/(flmct$nbar-1))*(flmct$pbar1*(1-flmct$pbar1)-(r-1)*flmct$ssq1/r-1/4*flmct$hbar1))
flmct$a2<-(flmct$nbar/flmct$nc)*(flmct$ssq2-(1/(flmct$nbar-1))*(flmct$pbar2*(1-flmct$pbar2)-(r-1)*flmct$ssq2/r-1/4*flmct$hbar2))
flmct$a3<-(flmct$nbar/flmct$nc)*(flmct$ssq3-(1/(flmct$nbar-1))*(flmct$pbar3*(1-flmct$pbar3)-(r-1)*flmct$ssq3/r-1/4*flmct$hbar3))
flmct$a4<-(flmct$nbar/flmct$nc)*(flmct$ssq4-(1/(flmct$nbar-1))*(flmct$pbar4*(1-flmct$pbar4)-(r-1)*flmct$ssq4/r-1/4*flmct$hbar4))

flmct$b1<-(flmct$nbar/(flmct$nbar-1))*(flmct$pbar1*(1-flmct$pbar1)-(r-1)*flmct$ssq1/r-(2*flmct$nbar-1)*flmct$hbar1/(4*flmct$nbar))
flmct$b2<-(flmct$nbar/(flmct$nbar-1))*(flmct$pbar2*(1-flmct$pbar2)-(r-1)*flmct$ssq2/r-(2*flmct$nbar-1)*flmct$hbar2/(4*flmct$nbar))
flmct$b3<-(flmct$nbar/(flmct$nbar-1))*(flmct$pbar3*(1-flmct$pbar3)-(r-1)*flmct$ssq3/r-(2*flmct$nbar-1)*flmct$hbar3/(4*flmct$nbar))
flmct$b4<-(flmct$nbar/(flmct$nbar-1))*(flmct$pbar4*(1-flmct$pbar4)-(r-1)*flmct$ssq4/r-(2*flmct$nbar-1)*flmct$hbar4/(4*flmct$nbar))

flmct$c1<-1/2*flmct$hbar1
flmct$c2<-1/2*flmct$hbar2
flmct$c3<-1/2*flmct$hbar3
flmct$c4<-1/2*flmct$hbar4

flmct$FST<-((flmct$a1+flmct$a2+flmct$a3+flmct$a4)/((flmct$a1+flmct$b1+flmct$c1)+(flmct$a2+flmct$b2+flmct$c2)+(flmct$a3+flmct$b3+flmct$c3)+(flmct$a4+flmct$b4+flmct$c4)))

flmct$FIT<-(((flmct$a1+flmct$b1)+(flmct$a2+flmct$b2)+(flmct$a3+flmct$b3)+(flmct$a4+flmct$b4))/((flmct$a1+flmct$b1+flmct$c1)+(flmct$a2+flmct$b2+flmct$c2)+(flmct$a3+flmct$b3+flmct$c3)+(flmct$a4+flmct$b4+flmct$c4)))

flmct$FIS<-((flmct$b1+flmct$b2+flmct$b3+flmct$b4)/((flmct$b1+flmct$c1)+(flmct$b2+flmct$c2)+(flmct$b3+flmct$c3)+(flmct$b4+flmct$c4)))

write.csv(flmct,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/flmct.csv")

fst_lmct_vcf<-read.table(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fst_lmct.weir.fst",head=TRUE)
fst_lmct_vcf$locus<-paste(fst_lmct_vcf$CHROM,"-",fst_lmct_vcf$POS)

fst_lmct_vs<-merge(flmct,fst_lmct_vcf,by.x="flm.locus",by.y="locus")

pdf("/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fst_lmct_vs.pdf")
plot(fst_lmct_vs$FST,fst_lmct_vs$WEIR_AND_COCKERHAM_FST,xlab="fst_lmct_vcf",ylab="fst_lmct_calc")
dev.off()

#####################################################################################

data4_lm_hwe<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data4-lm-hwe.csv")
data4_ct_hwe<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data4-ct-hwe.csv")

hwe_lmct<-merge(data4_lm_hwe,data4_ct_hwe,by.x="locus",by.y="locus")

flmct_hwe<-merge(flmct,hwe_lmct,by.x="flm.locus",by.y="locus")
fst_lmct_hwe<-data.frame(flmct_hwe$flm.locus,flmct_hwe$FST)
colnames(fst_lmct_hwe)<-c("locus","FST")

#################################################################################################
# for (i in 1:dim(fst_lmct_hwe)[1]){
# fst_lmct_hwe$CHROM[i]<-strsplit(as.character(fst_lmct_hwe$locus[i]),"-",fixed=TRUE)[[1]][1]
# fst_lmct_hwe$SCAF[i]<-strsplit(as.character(fst_lmct_hwe$CHROM[i]),"_",fixed=TRUE)[[1]][2]
# fst_lmct_hwe$POS[i]<-strsplit(as.character(fst_lmct_hwe$locus[i]),"-",fixed=TRUE)[[1]][2]
# }
#################################################################################################

fst_lmct_hwe_1<-cSplit(fst_lmct_hwe,"locus",sep = "-")
fst_lmct_hwe_2<-cSplit(fst_lmct_hwe_1,"locus_1",sep = "_")
fst_lmct_hwe_2$CHROM<-sprintf("scaf_%05d",fst_lmct_hwe_2$locus_1_2)
fst_lmct_hwe_2$locus<-paste(fst_lmct_hwe_2$CHROM,"-",fst_lmct_hwe_2$locus_2)
fst_lmct_hwe_2$SCAF<-sprintf("%05d",fst_lmct_hwe_2$locus_1_2)
fst_lmct_hwe_3<-data.frame(fst_lmct_hwe_2$locus,fst_lmct_hwe_2$FST,fst_lmct_hwe_2$CHROM,fst_lmct_hwe_2$SCAF,fst_lmct_hwe_2$locus_2)
colnames(fst_lmct_hwe_3)<-c("locus","FST","CHROM","SCAF","POS")

fstlmcthwe<-fst_lmct_hwe_3[order(fst_lmct_hwe_3$SCAF,fst_lmct_hwe_3$POS),]

fstlmcthwe<-fstlmcthwe[fstlmcthwe$FST!="NaN",]

write.csv(fstlmcthwe,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fstlmcthwe.csv")






#######################
## plot Fst by scaffold
#######################

## ZFST

fstlmcthwe<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fstlmcthwe.csv")

avgFst_lmct_hwe<-colMeans(as.matrix(fstlmcthwe$FST))
sdFst_lmct_hwe<-sd(fstlmcthwe$FST)
fstlmcthwe$ZFST<-(fstlmcthwe$FST-avgFst_lmct_hwe)/sdFst_lmct_hwe

fstlmcthwe90<-fstlmcthwe[fstlmcthwe$SCAF<91,]

max.pos<-sapply(1:89,function(i){(max(fstlmcthwe90$POS[fstlmcthwe90$SCAF==i]))})
add1<-c(0,cumsum(max.pos))
add2<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90)

added.length<-(add1+add2)

fstlmcthwe90$gp<-fstlmcthwe90$POS+added.length[fstlmcthwe90$SCAF]

#########################################
# for (i in 1:dim(fstlmcthwe90)[1]){
# if((fstlmcthwe90$SCAF[i] %% 2) == 0) {
# fstlmcthwe90$color[i]<-"darkorchid3"
# } else {
# fstlmcthwe90$color[i]<-"black"
# }}
#########################################

fstlmcthwe90$color<-((fstlmcthwe90$SCAF%%2)==0)
fstlmcthwe90$color[fstlmcthwe90$color==TRUE]<-"darkorchid3"
fstlmcthwe90$color[fstlmcthwe90$color==FALSE]<-"black"

## plot(fstlmcthwe90$gp,fstlmcthwe90$ZFST,xlab="Genomic Position",ylab="Z(FST)",xlim=c(0,900000000),ylim=c(-2,10),col=fstlmcthwe90$color,pch=20)
## abline(h=5,col="red",lty=2)

write.csv(fstlmcthwe90,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fstlmcthwe90.csv")
write.csv(fstlmcthwe,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fstlmcthwe_zfst.csv")






##############################
## plot ZFST and ZFST-window
##############################

fstlmcthwe90<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fstlmcthwe90.csv")

pdf("/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fstlmcthwe90.pdf")
plot(fstlmcthwe90$gp,fstlmcthwe90$ZFST,xlab="Genomic Position",ylab="Z(FST)",xlim=c(0,900000000),ylim=c(-2,10),col=fstlmcthwe90$color,pch=20)
abline(h=5,col="red",lty=2)
dev.off()
