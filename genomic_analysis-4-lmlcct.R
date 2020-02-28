################################################################################
## calculate Fst between Lake Michigan, Lake Champlain and the Connecticut River
################################################################################

install.packages("splitstackshape",repos="https://repo.miserver.it.umich.edu/cran/",INSTALL_opts=c('--no-lock'))
library(splitstackshape)

flm<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fh-lm.csv")
flc<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fh-lc.csv")
fct<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fh-ct.csv")

flmtemp<-data.frame(flm$locus,flm$G11,flm$G21,flm$G22,flm$G31,flm$G32,flm$G33,flm$G41,flm$G42,flm$G43,flm$G44,flm$freq1,flm$freq2,flm$freq3,flm$freq4)
flctemp<-data.frame(flc$locus,flc$G11,flc$G21,flc$G22,flc$G31,flc$G32,flc$G33,flc$G41,flc$G42,flc$G43,flc$G44,flc$freq1,flc$freq2,flc$freq3,flc$freq4)
fcttemp<-data.frame(fct$locus,fct$G11,fct$G21,fct$G22,fct$G31,fct$G32,fct$G33,fct$G41,fct$G42,fct$G43,fct$G44,fct$freq1,fct$freq2,fct$freq3,fct$freq4)

flcct<-merge(flctemp,fcttemp,by.x="flc.locus",by.y="fct.locus")
flmlcct<-merge(flcct,flmtemp,by.x="flc.locus",by.y="flm.locus")

flmlcct$nlm<-(flmlcct$flm.G11+flmlcct$flm.G21+flmlcct$flm.G22+flmlcct$flm.G31+flmlcct$flm.G32+flmlcct$flm.G33+flmlcct$flm.G41+flmlcct$flm.G42+flmlcct$flm.G43+flmlcct$flm.G44)
flmlcct$nlc<-(flmlcct$flc.G11+flmlcct$flc.G21+flmlcct$flc.G22+flmlcct$flc.G31+flmlcct$flc.G32+flmlcct$flc.G33+flmlcct$flc.G41+flmlcct$flc.G42+flmlcct$flc.G43+flmlcct$flc.G44)
flmlcct$nct<-(flmlcct$fct.G11+flmlcct$fct.G21+flmlcct$fct.G22+flmlcct$fct.G31+flmlcct$fct.G32+flmlcct$fct.G33+flmlcct$fct.G41+flmlcct$fct.G42+flmlcct$fct.G43+flmlcct$fct.G44)

r<-3

flmlcct$nbar<-(flmlcct$nlm+flmlcct$nlc+flmlcct$nct)/r
flmlcct$nc<-(r*flmlcct$nbar-((flmlcct$nlm^2+flmlcct$nlc^2+flmlcct$nct^2)/(r*flmlcct$nbar)))/(r-1)

flmlcct$pbar1<-(flmlcct$nlm*flmlcct$flm.freq1+flmlcct$nlc*flmlcct$flc.freq1+flmlcct$nct*flmlcct$fct.freq1)/(r*flmlcct$nbar)
flmlcct$pbar2<-(flmlcct$nlm*flmlcct$flm.freq2+flmlcct$nlc*flmlcct$flc.freq2+flmlcct$nct*flmlcct$fct.freq2)/(r*flmlcct$nbar)
flmlcct$pbar3<-(flmlcct$nlm*flmlcct$flm.freq3+flmlcct$nlc*flmlcct$flc.freq3+flmlcct$nct*flmlcct$fct.freq3)/(r*flmlcct$nbar)
flmlcct$pbar4<-(flmlcct$nlm*flmlcct$flm.freq4+flmlcct$nlc*flmlcct$flc.freq4+flmlcct$nct*flmlcct$fct.freq4)/(r*flmlcct$nbar)

flmlcct$ssq1<-(flmlcct$nlm*(flmlcct$flm.freq1-flmlcct$pbar1)^2+flmlcct$nlc*(flmlcct$flc.freq1-flmlcct$pbar1)^2+flmlcct$nct*(flmlcct$fct.freq1-flmlcct$pbar1)^2)/((r-1)*flmlcct$nbar)
flmlcct$ssq2<-(flmlcct$nlm*(flmlcct$flm.freq2-flmlcct$pbar2)^2+flmlcct$nlc*(flmlcct$flc.freq2-flmlcct$pbar2)^2+flmlcct$nct*(flmlcct$fct.freq2-flmlcct$pbar2)^2)/((r-1)*flmlcct$nbar)
flmlcct$ssq3<-(flmlcct$nlm*(flmlcct$flm.freq3-flmlcct$pbar3)^2+flmlcct$nlc*(flmlcct$flc.freq3-flmlcct$pbar3)^2+flmlcct$nct*(flmlcct$fct.freq3-flmlcct$pbar3)^2)/((r-1)*flmlcct$nbar)
flmlcct$ssq4<-(flmlcct$nlm*(flmlcct$flm.freq4-flmlcct$pbar4)^2+flmlcct$nlc*(flmlcct$flc.freq4-flmlcct$pbar4)^2+flmlcct$nct*(flmlcct$fct.freq4-flmlcct$pbar4)^2)/((r-1)*flmlcct$nbar)

flmlcct$flm.h1<-(flmlcct$flm.G21+flmlcct$flm.G31+flmlcct$flm.G41)/(flmlcct$flm.G11+flmlcct$flm.G21+flmlcct$flm.G22+flmlcct$flm.G31+flmlcct$flm.G32+flmlcct$flm.G33+flmlcct$flm.G41+flmlcct$flm.G42+flmlcct$flm.G43+flmlcct$flm.G44)
flmlcct$flc.h1<-(flmlcct$flc.G21+flmlcct$flc.G31+flmlcct$flc.G41)/(flmlcct$flc.G11+flmlcct$flc.G21+flmlcct$flc.G22+flmlcct$flc.G31+flmlcct$flc.G32+flmlcct$flc.G33+flmlcct$flc.G41+flmlcct$flc.G42+flmlcct$flc.G43+flmlcct$flc.G44)
flmlcct$fct.h1<-(flmlcct$fct.G21+flmlcct$fct.G31+flmlcct$fct.G41)/(flmlcct$fct.G11+flmlcct$fct.G21+flmlcct$fct.G22+flmlcct$fct.G31+flmlcct$fct.G32+flmlcct$fct.G33+flmlcct$fct.G41+flmlcct$fct.G42+flmlcct$fct.G43+flmlcct$fct.G44)

flmlcct$flm.h2<-(flmlcct$flm.G21+flmlcct$flm.G32+flmlcct$flm.G42)/(flmlcct$flm.G11+flmlcct$flm.G21+flmlcct$flm.G22+flmlcct$flm.G31+flmlcct$flm.G32+flmlcct$flm.G33+flmlcct$flm.G41+flmlcct$flm.G42+flmlcct$flm.G43+flmlcct$flm.G44)
flmlcct$flc.h2<-(flmlcct$flc.G21+flmlcct$flc.G32+flmlcct$flc.G42)/(flmlcct$flc.G11+flmlcct$flc.G21+flmlcct$flc.G22+flmlcct$flc.G31+flmlcct$flc.G32+flmlcct$flc.G33+flmlcct$flc.G41+flmlcct$flc.G42+flmlcct$flc.G43+flmlcct$flc.G44)
flmlcct$fct.h2<-(flmlcct$fct.G21+flmlcct$fct.G32+flmlcct$fct.G42)/(flmlcct$fct.G11+flmlcct$fct.G21+flmlcct$fct.G22+flmlcct$fct.G31+flmlcct$fct.G32+flmlcct$fct.G33+flmlcct$fct.G41+flmlcct$fct.G42+flmlcct$fct.G43+flmlcct$fct.G44)

flmlcct$flm.h3<-(flmlcct$flm.G31+flmlcct$flm.G32+flmlcct$flm.G43)/(flmlcct$flm.G11+flmlcct$flm.G21+flmlcct$flm.G22+flmlcct$flm.G31+flmlcct$flm.G32+flmlcct$flm.G33+flmlcct$flm.G41+flmlcct$flm.G42+flmlcct$flm.G43+flmlcct$flm.G44)
flmlcct$flc.h3<-(flmlcct$flc.G31+flmlcct$flc.G32+flmlcct$flc.G43)/(flmlcct$flc.G11+flmlcct$flc.G21+flmlcct$flc.G22+flmlcct$flc.G31+flmlcct$flc.G32+flmlcct$flc.G33+flmlcct$flc.G41+flmlcct$flc.G42+flmlcct$flc.G43+flmlcct$flc.G44)
flmlcct$fct.h3<-(flmlcct$fct.G31+flmlcct$fct.G32+flmlcct$fct.G43)/(flmlcct$fct.G11+flmlcct$fct.G21+flmlcct$fct.G22+flmlcct$fct.G31+flmlcct$fct.G32+flmlcct$fct.G33+flmlcct$fct.G41+flmlcct$fct.G42+flmlcct$fct.G43+flmlcct$fct.G44)

flmlcct$flm.h4<-(flmlcct$flm.G41+flmlcct$flm.G42+flmlcct$flm.G43)/(flmlcct$flm.G11+flmlcct$flm.G21+flmlcct$flm.G22+flmlcct$flm.G31+flmlcct$flm.G32+flmlcct$flm.G33+flmlcct$flm.G41+flmlcct$flm.G42+flmlcct$flm.G43+flmlcct$flm.G44)
flmlcct$flc.h4<-(flmlcct$flc.G41+flmlcct$flc.G42+flmlcct$flc.G43)/(flmlcct$flc.G11+flmlcct$flc.G21+flmlcct$flc.G22+flmlcct$flc.G31+flmlcct$flc.G32+flmlcct$flc.G33+flmlcct$flc.G41+flmlcct$flc.G42+flmlcct$flc.G43+flmlcct$flc.G44)
flmlcct$fct.h4<-(flmlcct$fct.G41+flmlcct$fct.G42+flmlcct$fct.G43)/(flmlcct$fct.G11+flmlcct$fct.G21+flmlcct$fct.G22+flmlcct$fct.G31+flmlcct$fct.G32+flmlcct$fct.G33+flmlcct$fct.G41+flmlcct$fct.G42+flmlcct$fct.G43+flmlcct$fct.G44)

flmlcct$hbar1<-(flmlcct$nlm*flmlcct$flm.h1+flmlcct$nlc*flmlcct$flc.h1+flmlcct$nct*flmlcct$fct.h1)/(r*flmlcct$nbar)
flmlcct$hbar2<-(flmlcct$nlm*flmlcct$flm.h2+flmlcct$nlc*flmlcct$flc.h2+flmlcct$nct*flmlcct$fct.h2)/(r*flmlcct$nbar)
flmlcct$hbar3<-(flmlcct$nlm*flmlcct$flm.h3+flmlcct$nlc*flmlcct$flc.h3+flmlcct$nct*flmlcct$fct.h3)/(r*flmlcct$nbar)
flmlcct$hbar4<-(flmlcct$nlm*flmlcct$flm.h4+flmlcct$nlc*flmlcct$flc.h4+flmlcct$nct*flmlcct$fct.h4)/(r*flmlcct$nbar)

flmlcct$a1<-(flmlcct$nbar/flmlcct$nc)*(flmlcct$ssq1-(1/(flmlcct$nbar-1))*(flmlcct$pbar1*(1-flmlcct$pbar1)-(r-1)*flmlcct$ssq1/r-1/4*flmlcct$hbar1))
flmlcct$a2<-(flmlcct$nbar/flmlcct$nc)*(flmlcct$ssq2-(1/(flmlcct$nbar-1))*(flmlcct$pbar2*(1-flmlcct$pbar2)-(r-1)*flmlcct$ssq2/r-1/4*flmlcct$hbar2))
flmlcct$a3<-(flmlcct$nbar/flmlcct$nc)*(flmlcct$ssq3-(1/(flmlcct$nbar-1))*(flmlcct$pbar3*(1-flmlcct$pbar3)-(r-1)*flmlcct$ssq3/r-1/4*flmlcct$hbar3))
flmlcct$a4<-(flmlcct$nbar/flmlcct$nc)*(flmlcct$ssq4-(1/(flmlcct$nbar-1))*(flmlcct$pbar4*(1-flmlcct$pbar4)-(r-1)*flmlcct$ssq4/r-1/4*flmlcct$hbar4))

flmlcct$b1<-(flmlcct$nbar/(flmlcct$nbar-1))*(flmlcct$pbar1*(1-flmlcct$pbar1)-(r-1)*flmlcct$ssq1/r-(2*flmlcct$nbar-1)*flmlcct$hbar1/(4*flmlcct$nbar))
flmlcct$b2<-(flmlcct$nbar/(flmlcct$nbar-1))*(flmlcct$pbar2*(1-flmlcct$pbar2)-(r-1)*flmlcct$ssq2/r-(2*flmlcct$nbar-1)*flmlcct$hbar2/(4*flmlcct$nbar))
flmlcct$b3<-(flmlcct$nbar/(flmlcct$nbar-1))*(flmlcct$pbar3*(1-flmlcct$pbar3)-(r-1)*flmlcct$ssq3/r-(2*flmlcct$nbar-1)*flmlcct$hbar3/(4*flmlcct$nbar))
flmlcct$b4<-(flmlcct$nbar/(flmlcct$nbar-1))*(flmlcct$pbar4*(1-flmlcct$pbar4)-(r-1)*flmlcct$ssq4/r-(2*flmlcct$nbar-1)*flmlcct$hbar4/(4*flmlcct$nbar))

flmlcct$c1<-1/2*flmlcct$hbar1
flmlcct$c2<-1/2*flmlcct$hbar2
flmlcct$c3<-1/2*flmlcct$hbar3
flmlcct$c4<-1/2*flmlcct$hbar4

flmlcct$FST<-((flmlcct$a1+flmlcct$a2+flmlcct$a3+flmlcct$a4)/((flmlcct$a1+flmlcct$b1+flmlcct$c1)+(flmlcct$a2+flmlcct$b2+flmlcct$c2)+(flmlcct$a3+flmlcct$b3+flmlcct$c3)+(flmlcct$a4+flmlcct$b4+flmlcct$c4)))

flmlcct$FIT<-(((flmlcct$a1+flmlcct$b1)+(flmlcct$a2+flmlcct$b2)+(flmlcct$a3+flmlcct$b3)+(flmlcct$a4+flmlcct$b4))/((flmlcct$a1+flmlcct$b1+flmlcct$c1)+(flmlcct$a2+flmlcct$b2+flmlcct$c2)+(flmlcct$a3+flmlcct$b3+flmlcct$c3)+(flmlcct$a4+flmlcct$b4+flmlcct$c4)))

flmlcct$FIS<-((flmlcct$b1+flmlcct$b2+flmlcct$b3+flmlcct$b4)/((flmlcct$b1+flmlcct$c1)+(flmlcct$b2+flmlcct$c2)+(flmlcct$b3+flmlcct$c3)+(flmlcct$b4+flmlcct$c4)))

write.csv(flmlcct,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/flmlcct.csv")

fst_lmlcct_vcf<-read.table(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fst_lmlcct.weir.fst",head=TRUE)
fst_lmlcct_vcf$locus<-paste(fst_lmlcct_vcf$CHROM,"-",fst_lmlcct_vcf$POS)

fst_lmlcct_vs<-merge(flmlcct,fst_lmlcct_vcf,by.x="flc.locus",by.y="locus")

pdf("/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fst_lmlcct_vs.pdf")
plot(fst_lmlcct_vs$FST,fst_lmlcct_vs$WEIR_AND_COCKERHAM_FST,xlab="fst_lmlcct_vcf",ylab="fst_lmlcct_calc")
dev.off()

#####################################################################################

data4_lm_hwe<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data4-lm-hwe.csv")
data4_lc_hwe<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data4-lc-hwe.csv")
data4_ct_hwe<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data4-ct-hwe.csv")

hwe_lcct<-merge(data4_lc_hwe,data4_ct_hwe,by.x="locus",by.y="locus")
hwe_lmlcct<-merge(hwe_lcct,data4_lm_hwe,by.x="locus",by.y="locus")

flmlcct_hwe<-merge(flmlcct,hwe_lmlcct,by.x="flc.locus",by.y="locus")
fst_lmlcct_hwe<-data.frame(flmlcct_hwe$flc.locus,flmlcct_hwe$FST)
colnames(fst_lmlcct_hwe)<-c("locus","FST")

#################################################################################################
# for (i in 1:dim(fst_lmlcct_hwe)[1]){
# fst_lmlcct_hwe$CHROM[i]<-strsplit(as.character(fst_lmlcct_hwe$locus[i]),"-",fixed=TRUE)[[1]][1]
# fst_lmlcct_hwe$SCAF[i]<-strsplit(as.character(fst_lmlcct_hwe$CHROM[i]),"_",fixed=TRUE)[[1]][2]
# fst_lmlcct_hwe$POS[i]<-strsplit(as.character(fst_lmlcct_hwe$locus[i]),"-",fixed=TRUE)[[1]][2]
# }
#################################################################################################

fst_lmlcct_hwe_1<-cSplit(fst_lmlcct_hwe,"locus",sep = "-")
fst_lmlcct_hwe_2<-cSplit(fst_lmlcct_hwe_1,"locus_1",sep = "_")
fst_lmlcct_hwe_2$CHROM<-sprintf("scaf_%05d",fst_lmlcct_hwe_2$locus_1_2)
fst_lmlcct_hwe_2$locus<-paste(fst_lmlcct_hwe_2$CHROM,"-",fst_lmlcct_hwe_2$locus_2)
fst_lmlcct_hwe_2$SCAF<-sprintf("%05d",fst_lmlcct_hwe_2$locus_1_2)
fst_lmlcct_hwe_3<-data.frame(fst_lmlcct_hwe_2$locus,fst_lmlcct_hwe_2$FST,fst_lmlcct_hwe_2$CHROM,fst_lmlcct_hwe_2$SCAF,fst_lmlcct_hwe_2$locus_2)
colnames(fst_lmlcct_hwe_3)<-c("locus","FST","CHROM","SCAF","POS")

fstlmlccthwe<-fst_lmlcct_hwe_3[order(fst_lmlcct_hwe_3$SCAF,fst_lmlcct_hwe_3$POS),]

fstlmlccthwe<-fstlmlccthwe[fstlmlccthwe$FST!="NaN",]

write.csv(fstlmlccthwe,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fstlmlccthwe.csv")






#######################
## plot Fst by scaffold
#######################

## ZFST

fstlmlccthwe<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fstlmlccthwe.csv")

avgFst_lmlcct_hwe<-colMeans(as.matrix(fstlmlccthwe$FST))
sdFst_lmlcct_hwe<-sd(fstlmlccthwe$FST)
fstlmlccthwe$ZFST<-(fstlmlccthwe$FST-avgFst_lmlcct_hwe)/sdFst_lmlcct_hwe

fstlmlccthwe90<-fstlmlccthwe[fstlmlccthwe$SCAF<91,]

max.pos<-sapply(1:89,function(i){(max(fstlmlccthwe90$POS[fstlmlccthwe90$SCAF==i]))})
add1<-c(0,cumsum(max.pos))
add2<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90)

added.length<-(add1+add2)

fstlmlccthwe90$gp<-fstlmlccthwe90$POS+added.length[fstlmlccthwe90$SCAF]

#########################################
# for (i in 1:dim(fstlmlccthwe90)[1]){
# if((fstlmlccthwe90$SCAF[i] %% 2) == 0) {
# fstlmlccthwe90$color[i]<-"darkorchid3"
# } else {
# fstlmlccthwe90$color[i]<-"black"
# }}
#########################################

fstlmlccthwe90$color<-((fstlmlccthwe90$SCAF%%2)==0)
fstlmlccthwe90$color[fstlmlccthwe90$color==TRUE]<-"darkorchid3"
fstlmlccthwe90$color[fstlmlccthwe90$color==FALSE]<-"black"

## plot(fstlmlccthwe90$gp,fstlmlccthwe90$ZFST,xlab="Genomic Position",ylab="Z(FST)",xlim=c(0,900000000),ylim=c(-2,10),col=fstlmlccthwe90$color,pch=20)
## abline(h=5,col="red",lty=2)

write.csv(fstlmlccthwe90,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fstlmlccthwe90.csv")
write.csv(fstlmlccthwe,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fstlmlccthwe_zfst.csv")






##############################
## plot ZFST and ZFST-window
##############################

fstlmlccthwe90<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fstlmlccthwe90.csv")

pdf("/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fstlmlccthwe90.pdf")
plot(fstlmlccthwe90$gp,fstlmlccthwe90$ZFST,xlab="Genomic Position",ylab="Z(FST)",xlim=c(0,900000000),ylim=c(-2,10),col=fstlmlccthwe90$color,pch=20)
abline(h=5,col="red",lty=2)
dev.off()
