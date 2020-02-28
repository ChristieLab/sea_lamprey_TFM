#################################################################
## calculate Fst between Lake Champlain and the Connecticut River
#################################################################

install.packages("splitstackshape",repos="https://repo.miserver.it.umich.edu/cran/",INSTALL_opts=c('--no-lock'))
library(splitstackshape)

flc<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fh-lc.csv")
fct<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fh-ct.csv")

flctemp<-data.frame(flc$locus,flc$G11,flc$G21,flc$G22,flc$G31,flc$G32,flc$G33,flc$G41,flc$G42,flc$G43,flc$G44,flc$freq1,flc$freq2,flc$freq3,flc$freq4)
fcttemp<-data.frame(fct$locus,fct$G11,fct$G21,fct$G22,fct$G31,fct$G32,fct$G33,fct$G41,fct$G42,fct$G43,fct$G44,fct$freq1,fct$freq2,fct$freq3,fct$freq4)

flcct<-merge(flctemp,fcttemp,by.x="flc.locus",by.y="fct.locus")



flcct$nlc<-(flcct$flc.G11+flcct$flc.G21+flcct$flc.G22+flcct$flc.G31+flcct$flc.G32+flcct$flc.G33+flcct$flc.G41+flcct$flc.G42+flcct$flc.G43+flcct$flc.G44)
flcct$nct<-(flcct$fct.G11+flcct$fct.G21+flcct$fct.G22+flcct$fct.G31+flcct$fct.G32+flcct$fct.G33+flcct$fct.G41+flcct$fct.G42+flcct$fct.G43+flcct$fct.G44)

r<-2

flcct$nbar<-(flcct$nlc+flcct$nct)/r
flcct$nc<-(r*flcct$nbar-((flcct$nlc^2+flcct$nct^2)/(r*flcct$nbar)))/(r-1)

flcct$pbar1<-(flcct$nlc*flcct$flc.freq1+flcct$nct*flcct$fct.freq1)/(r*flcct$nbar)
flcct$pbar2<-(flcct$nlc*flcct$flc.freq2+flcct$nct*flcct$fct.freq2)/(r*flcct$nbar)
flcct$pbar3<-(flcct$nlc*flcct$flc.freq3+flcct$nct*flcct$fct.freq3)/(r*flcct$nbar)
flcct$pbar4<-(flcct$nlc*flcct$flc.freq4+flcct$nct*flcct$fct.freq4)/(r*flcct$nbar)

flcct$ssq1<-(flcct$nlc*(flcct$flc.freq1-flcct$pbar1)^2+flcct$nct*(flcct$fct.freq1-flcct$pbar1)^2)/((r-1)*flcct$nbar)
flcct$ssq2<-(flcct$nlc*(flcct$flc.freq2-flcct$pbar2)^2+flcct$nct*(flcct$fct.freq2-flcct$pbar2)^2)/((r-1)*flcct$nbar)
flcct$ssq3<-(flcct$nlc*(flcct$flc.freq3-flcct$pbar3)^2+flcct$nct*(flcct$fct.freq3-flcct$pbar3)^2)/((r-1)*flcct$nbar)
flcct$ssq4<-(flcct$nlc*(flcct$flc.freq4-flcct$pbar4)^2+flcct$nct*(flcct$fct.freq4-flcct$pbar4)^2)/((r-1)*flcct$nbar)

flcct$flc.h1<-(flcct$flc.G21+flcct$flc.G31+flcct$flc.G41)/(flcct$flc.G11+flcct$flc.G21+flcct$flc.G22+flcct$flc.G31+flcct$flc.G32+flcct$flc.G33+flcct$flc.G41+flcct$flc.G42+flcct$flc.G43+flcct$flc.G44)
flcct$fct.h1<-(flcct$fct.G21+flcct$fct.G31+flcct$fct.G41)/(flcct$fct.G11+flcct$fct.G21+flcct$fct.G22+flcct$fct.G31+flcct$fct.G32+flcct$fct.G33+flcct$fct.G41+flcct$fct.G42+flcct$fct.G43+flcct$fct.G44)

flcct$flc.h2<-(flcct$flc.G21+flcct$flc.G32+flcct$flc.G42)/(flcct$flc.G11+flcct$flc.G21+flcct$flc.G22+flcct$flc.G31+flcct$flc.G32+flcct$flc.G33+flcct$flc.G41+flcct$flc.G42+flcct$flc.G43+flcct$flc.G44)
flcct$fct.h2<-(flcct$fct.G21+flcct$fct.G32+flcct$fct.G42)/(flcct$fct.G11+flcct$fct.G21+flcct$fct.G22+flcct$fct.G31+flcct$fct.G32+flcct$fct.G33+flcct$fct.G41+flcct$fct.G42+flcct$fct.G43+flcct$fct.G44)

flcct$flc.h3<-(flcct$flc.G31+flcct$flc.G32+flcct$flc.G43)/(flcct$flc.G11+flcct$flc.G21+flcct$flc.G22+flcct$flc.G31+flcct$flc.G32+flcct$flc.G33+flcct$flc.G41+flcct$flc.G42+flcct$flc.G43+flcct$flc.G44)
flcct$fct.h3<-(flcct$fct.G31+flcct$fct.G32+flcct$fct.G43)/(flcct$fct.G11+flcct$fct.G21+flcct$fct.G22+flcct$fct.G31+flcct$fct.G32+flcct$fct.G33+flcct$fct.G41+flcct$fct.G42+flcct$fct.G43+flcct$fct.G44)

flcct$flc.h4<-(flcct$flc.G41+flcct$flc.G42+flcct$flc.G43)/(flcct$flc.G11+flcct$flc.G21+flcct$flc.G22+flcct$flc.G31+flcct$flc.G32+flcct$flc.G33+flcct$flc.G41+flcct$flc.G42+flcct$flc.G43+flcct$flc.G44)
flcct$fct.h4<-(flcct$fct.G41+flcct$fct.G42+flcct$fct.G43)/(flcct$fct.G11+flcct$fct.G21+flcct$fct.G22+flcct$fct.G31+flcct$fct.G32+flcct$fct.G33+flcct$fct.G41+flcct$fct.G42+flcct$fct.G43+flcct$fct.G44)

flcct$hbar1<-(flcct$nlc*flcct$flc.h1+flcct$nct*flcct$fct.h1)/(r*flcct$nbar)
flcct$hbar2<-(flcct$nlc*flcct$flc.h2+flcct$nct*flcct$fct.h2)/(r*flcct$nbar)
flcct$hbar3<-(flcct$nlc*flcct$flc.h3+flcct$nct*flcct$fct.h3)/(r*flcct$nbar)
flcct$hbar4<-(flcct$nlc*flcct$flc.h4+flcct$nct*flcct$fct.h4)/(r*flcct$nbar)

flcct$a1<-(flcct$nbar/flcct$nc)*(flcct$ssq1-(1/(flcct$nbar-1))*(flcct$pbar1*(1-flcct$pbar1)-(r-1)*flcct$ssq1/r-1/4*flcct$hbar1))
flcct$a2<-(flcct$nbar/flcct$nc)*(flcct$ssq2-(1/(flcct$nbar-1))*(flcct$pbar2*(1-flcct$pbar2)-(r-1)*flcct$ssq2/r-1/4*flcct$hbar2))
flcct$a3<-(flcct$nbar/flcct$nc)*(flcct$ssq3-(1/(flcct$nbar-1))*(flcct$pbar3*(1-flcct$pbar3)-(r-1)*flcct$ssq3/r-1/4*flcct$hbar3))
flcct$a4<-(flcct$nbar/flcct$nc)*(flcct$ssq4-(1/(flcct$nbar-1))*(flcct$pbar4*(1-flcct$pbar4)-(r-1)*flcct$ssq4/r-1/4*flcct$hbar4))

flcct$b1<-(flcct$nbar/(flcct$nbar-1))*(flcct$pbar1*(1-flcct$pbar1)-(r-1)*flcct$ssq1/r-(2*flcct$nbar-1)*flcct$hbar1/(4*flcct$nbar))
flcct$b2<-(flcct$nbar/(flcct$nbar-1))*(flcct$pbar2*(1-flcct$pbar2)-(r-1)*flcct$ssq2/r-(2*flcct$nbar-1)*flcct$hbar2/(4*flcct$nbar))
flcct$b3<-(flcct$nbar/(flcct$nbar-1))*(flcct$pbar3*(1-flcct$pbar3)-(r-1)*flcct$ssq3/r-(2*flcct$nbar-1)*flcct$hbar3/(4*flcct$nbar))
flcct$b4<-(flcct$nbar/(flcct$nbar-1))*(flcct$pbar4*(1-flcct$pbar4)-(r-1)*flcct$ssq4/r-(2*flcct$nbar-1)*flcct$hbar4/(4*flcct$nbar))

flcct$c1<-1/2*flcct$hbar1
flcct$c2<-1/2*flcct$hbar2
flcct$c3<-1/2*flcct$hbar3
flcct$c4<-1/2*flcct$hbar4

flcct$FST<-((flcct$a1+flcct$a2+flcct$a3+flcct$a4)/((flcct$a1+flcct$b1+flcct$c1)+(flcct$a2+flcct$b2+flcct$c2)+(flcct$a3+flcct$b3+flcct$c3)+(flcct$a4+flcct$b4+flcct$c4)))

flcct$FIT<-(((flcct$a1+flcct$b1)+(flcct$a2+flcct$b2)+(flcct$a3+flcct$b3)+(flcct$a4+flcct$b4))/((flcct$a1+flcct$b1+flcct$c1)+(flcct$a2+flcct$b2+flcct$c2)+(flcct$a3+flcct$b3+flcct$c3)+(flcct$a4+flcct$b4+flcct$c4)))

flcct$FIS<-((flcct$b1+flcct$b2+flcct$b3+flcct$b4)/((flcct$b1+flcct$c1)+(flcct$b2+flcct$c2)+(flcct$b3+flcct$c3)+(flcct$b4+flcct$c4)))

write.csv(flcct,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/flcct.csv")

fst_lcct_vcf<-read.table(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fst_lcct.weir.fst",head=TRUE)
fst_lcct_vcf$locus<-paste(fst_lcct_vcf$CHROM,"-",fst_lcct_vcf$POS)

fst_lcct_vs<-merge(flcct,fst_lcct_vcf,by.x="flc.locus",by.y="locus")

pdf("/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fst_lcct_vs.pdf")
plot(fst_lcct_vs$FST,fst_lcct_vs$WEIR_AND_COCKERHAM_FST,xlab="fst_lcct_vcf",ylab="fst_lcct_calc")
dev.off()

#####################################################################################

data4_lc_hwe<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data4-lc-hwe.csv")
data4_ct_hwe<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data4-ct-hwe.csv")

hwe_lcct<-merge(data4_lc_hwe,data4_ct_hwe,by.x="locus",by.y="locus")

flcct_hwe<-merge(flcct,hwe_lcct,by.x="flc.locus",by.y="locus")
fst_lcct_hwe<-data.frame(flcct_hwe$flc.locus,flcct_hwe$FST)
colnames(fst_lcct_hwe)<-c("locus","FST")

#################################################################################################
# for (i in 1:dim(fst_lcct_hwe)[1]){
# fst_lcct_hwe$CHROM[i]<-strsplit(as.character(fst_lcct_hwe$locus[i]),"-",fixed=TRUE)[[1]][1]
# fst_lcct_hwe$SCAF[i]<-strsplit(as.character(fst_lcct_hwe$CHROM[i]),"_",fixed=TRUE)[[1]][2]
# fst_lcct_hwe$POS[i]<-strsplit(as.character(fst_lcct_hwe$locus[i]),"-",fixed=TRUE)[[1]][2]
# }
#################################################################################################

fst_lcct_hwe_1<-cSplit(fst_lcct_hwe,"locus",sep = "-")
fst_lcct_hwe_2<-cSplit(fst_lcct_hwe_1,"locus_1",sep = "_")
fst_lcct_hwe_2$CHROM<-sprintf("scaf_%05d",fst_lcct_hwe_2$locus_1_2)
fst_lcct_hwe_2$locus<-paste(fst_lcct_hwe_2$CHROM,"-",fst_lcct_hwe_2$locus_2)
fst_lcct_hwe_2$SCAF<-sprintf("%05d",fst_lcct_hwe_2$locus_1_2)
fst_lcct_hwe_3<-data.frame(fst_lcct_hwe_2$locus,fst_lcct_hwe_2$FST,fst_lcct_hwe_2$CHROM,fst_lcct_hwe_2$SCAF,fst_lcct_hwe_2$locus_2)
colnames(fst_lcct_hwe_3)<-c("locus","FST","CHROM","SCAF","POS")

fstlccthwe<-fst_lcct_hwe_3[order(fst_lcct_hwe_3$SCAF,fst_lcct_hwe_3$POS),]

fstlccthwe<-fstlccthwe[fstlccthwe$FST!="NaN",]

write.csv(fstlccthwe,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fstlccthwe.csv")






#######################
## plot Fst by scaffold
#######################

## ZFST

fstlccthwe<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fstlccthwe.csv")

avgFst_lcct_hwe<-colMeans(as.matrix(fstlccthwe$FST))
sdFst_lcct_hwe<-sd(fstlccthwe$FST)
fstlccthwe$ZFST<-(fstlccthwe$FST-avgFst_lcct_hwe)/sdFst_lcct_hwe

fstlccthwe90<-fstlccthwe[fstlccthwe$SCAF<91,]

max.pos<-sapply(1:89,function(i){(max(fstlccthwe90$POS[fstlccthwe90$SCAF==i]))})
add1<-c(0,cumsum(max.pos))
add2<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90)

added.length<-(add1+add2)

fstlccthwe90$gp<-fstlccthwe90$POS+added.length[fstlccthwe90$SCAF]

#########################################
# for (i in 1:dim(fstlccthwe90)[1]){
# if((fstlccthwe90$SCAF[i] %% 2) == 0) {
# fstlccthwe90$color[i]<-"darkorchid3"
# } else {
# fstlccthwe90$color[i]<-"black"
# }}
#########################################

fstlccthwe90$color<-((fstlccthwe90$SCAF%%2)==0)
fstlccthwe90$color[fstlccthwe90$color==TRUE]<-"darkorchid3"
fstlccthwe90$color[fstlccthwe90$color==FALSE]<-"black"

## plot(fstlccthwe90$gp,fstlccthwe90$ZFST,xlab="Genomic Position",ylab="Z(FST)",xlim=c(0,900000000),ylim=c(-2,10),col=fstlccthwe90$color,pch=20)
## abline(h=5,col="red",lty=2)

write.csv(fstlccthwe90,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fstlccthwe90.csv")
write.csv(fstlccthwe,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fstlccthwe_zfst.csv")






##############################
## plot ZFST and ZFST-window
##############################

fstlccthwe90<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fstlccthwe90.csv")

pdf("/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fstlccthwe90.pdf")
plot(fstlccthwe90$gp,fstlccthwe90$ZFST,xlab="Genomic Position",ylab="Z(FST)",xlim=c(0,900000000),ylim=c(-2,10),col=fstlccthwe90$color,pch=20)
abline(h=5,col="red",lty=2)
dev.off()
