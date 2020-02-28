#########################################################
## calculate Fst between Lake Michigan and Lake Champlain
#########################################################

install.packages("splitstackshape",repos="https://repo.miserver.it.umich.edu/cran/",INSTALL_opts=c('--no-lock'))
library(splitstackshape)

flm<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fh-lm.csv")
flc<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fh-lc.csv")

flmtemp<-data.frame(flm$locus,flm$G11,flm$G21,flm$G22,flm$G31,flm$G32,flm$G33,flm$G41,flm$G42,flm$G43,flm$G44,flm$freq1,flm$freq2,flm$freq3,flm$freq4)
flctemp<-data.frame(flc$locus,flc$G11,flc$G21,flc$G22,flc$G31,flc$G32,flc$G33,flc$G41,flc$G42,flc$G43,flc$G44,flc$freq1,flc$freq2,flc$freq3,flc$freq4)

flmlc<-merge(flmtemp,flctemp,by.x="flm.locus",by.y="flc.locus")



flmlc$nlm<-(flmlc$flm.G11+flmlc$flm.G21+flmlc$flm.G22+flmlc$flm.G31+flmlc$flm.G32+flmlc$flm.G33+flmlc$flm.G41+flmlc$flm.G42+flmlc$flm.G43+flmlc$flm.G44)
flmlc$nlc<-(flmlc$flc.G11+flmlc$flc.G21+flmlc$flc.G22+flmlc$flc.G31+flmlc$flc.G32+flmlc$flc.G33+flmlc$flc.G41+flmlc$flc.G42+flmlc$flc.G43+flmlc$flc.G44)

r<-2

flmlc$nbar<-(flmlc$nlm+flmlc$nlc)/r
flmlc$nc<-(r*flmlc$nbar-((flmlc$nlm^2+flmlc$nlc^2)/(r*flmlc$nbar)))/(r-1)

flmlc$pbar1<-(flmlc$nlm*flmlc$flm.freq1+flmlc$nlc*flmlc$flc.freq1)/(r*flmlc$nbar)
flmlc$pbar2<-(flmlc$nlm*flmlc$flm.freq2+flmlc$nlc*flmlc$flc.freq2)/(r*flmlc$nbar)
flmlc$pbar3<-(flmlc$nlm*flmlc$flm.freq3+flmlc$nlc*flmlc$flc.freq3)/(r*flmlc$nbar)
flmlc$pbar4<-(flmlc$nlm*flmlc$flm.freq4+flmlc$nlc*flmlc$flc.freq4)/(r*flmlc$nbar)

flmlc$ssq1<-(flmlc$nlm*(flmlc$flm.freq1-flmlc$pbar1)^2+flmlc$nlc*(flmlc$flc.freq1-flmlc$pbar1)^2)/((r-1)*flmlc$nbar)
flmlc$ssq2<-(flmlc$nlm*(flmlc$flm.freq2-flmlc$pbar2)^2+flmlc$nlc*(flmlc$flc.freq2-flmlc$pbar2)^2)/((r-1)*flmlc$nbar)
flmlc$ssq3<-(flmlc$nlm*(flmlc$flm.freq3-flmlc$pbar3)^2+flmlc$nlc*(flmlc$flc.freq3-flmlc$pbar3)^2)/((r-1)*flmlc$nbar)
flmlc$ssq4<-(flmlc$nlm*(flmlc$flm.freq4-flmlc$pbar4)^2+flmlc$nlc*(flmlc$flc.freq4-flmlc$pbar4)^2)/((r-1)*flmlc$nbar)

flmlc$flm.h1<-(flmlc$flm.G21+flmlc$flm.G31+flmlc$flm.G41)/(flmlc$flm.G11+flmlc$flm.G21+flmlc$flm.G22+flmlc$flm.G31+flmlc$flm.G32+flmlc$flm.G33+flmlc$flm.G41+flmlc$flm.G42+flmlc$flm.G43+flmlc$flm.G44)
flmlc$flc.h1<-(flmlc$flc.G21+flmlc$flc.G31+flmlc$flc.G41)/(flmlc$flc.G11+flmlc$flc.G21+flmlc$flc.G22+flmlc$flc.G31+flmlc$flc.G32+flmlc$flc.G33+flmlc$flc.G41+flmlc$flc.G42+flmlc$flc.G43+flmlc$flc.G44)

flmlc$flm.h2<-(flmlc$flm.G21+flmlc$flm.G32+flmlc$flm.G42)/(flmlc$flm.G11+flmlc$flm.G21+flmlc$flm.G22+flmlc$flm.G31+flmlc$flm.G32+flmlc$flm.G33+flmlc$flm.G41+flmlc$flm.G42+flmlc$flm.G43+flmlc$flm.G44)
flmlc$flc.h2<-(flmlc$flc.G21+flmlc$flc.G32+flmlc$flc.G42)/(flmlc$flc.G11+flmlc$flc.G21+flmlc$flc.G22+flmlc$flc.G31+flmlc$flc.G32+flmlc$flc.G33+flmlc$flc.G41+flmlc$flc.G42+flmlc$flc.G43+flmlc$flc.G44)

flmlc$flm.h3<-(flmlc$flm.G31+flmlc$flm.G32+flmlc$flm.G43)/(flmlc$flm.G11+flmlc$flm.G21+flmlc$flm.G22+flmlc$flm.G31+flmlc$flm.G32+flmlc$flm.G33+flmlc$flm.G41+flmlc$flm.G42+flmlc$flm.G43+flmlc$flm.G44)
flmlc$flc.h3<-(flmlc$flc.G31+flmlc$flc.G32+flmlc$flc.G43)/(flmlc$flc.G11+flmlc$flc.G21+flmlc$flc.G22+flmlc$flc.G31+flmlc$flc.G32+flmlc$flc.G33+flmlc$flc.G41+flmlc$flc.G42+flmlc$flc.G43+flmlc$flc.G44)

flmlc$flm.h4<-(flmlc$flm.G41+flmlc$flm.G42+flmlc$flm.G43)/(flmlc$flm.G11+flmlc$flm.G21+flmlc$flm.G22+flmlc$flm.G31+flmlc$flm.G32+flmlc$flm.G33+flmlc$flm.G41+flmlc$flm.G42+flmlc$flm.G43+flmlc$flm.G44)
flmlc$flc.h4<-(flmlc$flc.G41+flmlc$flc.G42+flmlc$flc.G43)/(flmlc$flc.G11+flmlc$flc.G21+flmlc$flc.G22+flmlc$flc.G31+flmlc$flc.G32+flmlc$flc.G33+flmlc$flc.G41+flmlc$flc.G42+flmlc$flc.G43+flmlc$flc.G44)

flmlc$hbar1<-(flmlc$nlm*flmlc$flm.h1+flmlc$nlc*flmlc$flc.h1)/(r*flmlc$nbar)
flmlc$hbar2<-(flmlc$nlm*flmlc$flm.h2+flmlc$nlc*flmlc$flc.h2)/(r*flmlc$nbar)
flmlc$hbar3<-(flmlc$nlm*flmlc$flm.h3+flmlc$nlc*flmlc$flc.h3)/(r*flmlc$nbar)
flmlc$hbar4<-(flmlc$nlm*flmlc$flm.h4+flmlc$nlc*flmlc$flc.h4)/(r*flmlc$nbar)

flmlc$a1<-(flmlc$nbar/flmlc$nc)*(flmlc$ssq1-(1/(flmlc$nbar-1))*(flmlc$pbar1*(1-flmlc$pbar1)-(r-1)*flmlc$ssq1/r-1/4*flmlc$hbar1))
flmlc$a2<-(flmlc$nbar/flmlc$nc)*(flmlc$ssq2-(1/(flmlc$nbar-1))*(flmlc$pbar2*(1-flmlc$pbar2)-(r-1)*flmlc$ssq2/r-1/4*flmlc$hbar2))
flmlc$a3<-(flmlc$nbar/flmlc$nc)*(flmlc$ssq3-(1/(flmlc$nbar-1))*(flmlc$pbar3*(1-flmlc$pbar3)-(r-1)*flmlc$ssq3/r-1/4*flmlc$hbar3))
flmlc$a4<-(flmlc$nbar/flmlc$nc)*(flmlc$ssq4-(1/(flmlc$nbar-1))*(flmlc$pbar4*(1-flmlc$pbar4)-(r-1)*flmlc$ssq4/r-1/4*flmlc$hbar4))

flmlc$b1<-(flmlc$nbar/(flmlc$nbar-1))*(flmlc$pbar1*(1-flmlc$pbar1)-(r-1)*flmlc$ssq1/r-(2*flmlc$nbar-1)*flmlc$hbar1/(4*flmlc$nbar))
flmlc$b2<-(flmlc$nbar/(flmlc$nbar-1))*(flmlc$pbar2*(1-flmlc$pbar2)-(r-1)*flmlc$ssq2/r-(2*flmlc$nbar-1)*flmlc$hbar2/(4*flmlc$nbar))
flmlc$b3<-(flmlc$nbar/(flmlc$nbar-1))*(flmlc$pbar3*(1-flmlc$pbar3)-(r-1)*flmlc$ssq3/r-(2*flmlc$nbar-1)*flmlc$hbar3/(4*flmlc$nbar))
flmlc$b4<-(flmlc$nbar/(flmlc$nbar-1))*(flmlc$pbar4*(1-flmlc$pbar4)-(r-1)*flmlc$ssq4/r-(2*flmlc$nbar-1)*flmlc$hbar4/(4*flmlc$nbar))

flmlc$c1<-1/2*flmlc$hbar1
flmlc$c2<-1/2*flmlc$hbar2
flmlc$c3<-1/2*flmlc$hbar3
flmlc$c4<-1/2*flmlc$hbar4

flmlc$FST<-((flmlc$a1+flmlc$a2+flmlc$a3+flmlc$a4)/((flmlc$a1+flmlc$b1+flmlc$c1)+(flmlc$a2+flmlc$b2+flmlc$c2)+(flmlc$a3+flmlc$b3+flmlc$c3)+(flmlc$a4+flmlc$b4+flmlc$c4)))

flmlc$FIT<-(((flmlc$a1+flmlc$b1)+(flmlc$a2+flmlc$b2)+(flmlc$a3+flmlc$b3)+(flmlc$a4+flmlc$b4))/((flmlc$a1+flmlc$b1+flmlc$c1)+(flmlc$a2+flmlc$b2+flmlc$c2)+(flmlc$a3+flmlc$b3+flmlc$c3)+(flmlc$a4+flmlc$b4+flmlc$c4)))

flmlc$FIS<-((flmlc$b1+flmlc$b2+flmlc$b3+flmlc$b4)/((flmlc$b1+flmlc$c1)+(flmlc$b2+flmlc$c2)+(flmlc$b3+flmlc$c3)+(flmlc$b4+flmlc$c4)))

write.csv(flmlc,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/flmlc.csv")

fst_lmlc_vcf<-read.table(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fst_lmlc.weir.fst",head=TRUE)
fst_lmlc_vcf$locus<-paste(fst_lmlc_vcf$CHROM,"-",fst_lmlc_vcf$POS)

fst_lmlc_vs<-merge(flmlc,fst_lmlc_vcf,by.x="flm.locus",by.y="locus")

pdf("/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fst_lmlc_vs.pdf")
plot(fst_lmlc_vs$FST,fst_lmlc_vs$WEIR_AND_COCKERHAM_FST,xlab="fst_lmlc_vcf",ylab="fst_lmlc_calc")
dev.off()

#####################################################################################

data4_lm_hwe<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data4-lm-hwe.csv")
data4_lc_hwe<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data4-lc-hwe.csv")

hwe_lmlc<-merge(data4_lm_hwe,data4_lc_hwe,by.x="locus",by.y="locus")

flmlc_hwe<-merge(flmlc,hwe_lmlc,by.x="flm.locus",by.y="locus")
fst_lmlc_hwe<-data.frame(flmlc_hwe$flm.locus,flmlc_hwe$FST)
colnames(fst_lmlc_hwe)<-c("locus","FST")

#################################################################################################
# for (i in 1:dim(fst_lmlc_hwe)[1]){
# fst_lmlc_hwe$CHROM[i]<-strsplit(as.character(fst_lmlc_hwe$locus[i]),"-",fixed=TRUE)[[1]][1]
# fst_lmlc_hwe$SCAF[i]<-strsplit(as.character(fst_lmlc_hwe$CHROM[i]),"_",fixed=TRUE)[[1]][2]
# fst_lmlc_hwe$POS[i]<-strsplit(as.character(fst_lmlc_hwe$locus[i]),"-",fixed=TRUE)[[1]][2]
# }
#################################################################################################

fst_lmlc_hwe_1<-cSplit(fst_lmlc_hwe,"locus",sep = "-")
fst_lmlc_hwe_2<-cSplit(fst_lmlc_hwe_1,"locus_1",sep = "_")
fst_lmlc_hwe_2$CHROM<-sprintf("scaf_%05d",fst_lmlc_hwe_2$locus_1_2)
fst_lmlc_hwe_2$locus<-paste(fst_lmlc_hwe_2$CHROM,"-",fst_lmlc_hwe_2$locus_2)
fst_lmlc_hwe_2$SCAF<-sprintf("%05d",fst_lmlc_hwe_2$locus_1_2)
fst_lmlc_hwe_3<-data.frame(fst_lmlc_hwe_2$locus,fst_lmlc_hwe_2$FST,fst_lmlc_hwe_2$CHROM,fst_lmlc_hwe_2$SCAF,fst_lmlc_hwe_2$locus_2)
colnames(fst_lmlc_hwe_3)<-c("locus","FST","CHROM","SCAF","POS")

fstlmlchwe<-fst_lmlc_hwe_3[order(fst_lmlc_hwe_3$SCAF,fst_lmlc_hwe_3$POS),]

fstlmlchwe<-fstlmlchwe[fstlmlchwe$FST!="NaN",]

write.csv(fstlmlchwe,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fstlmlchwe.csv")






#######################
## plot Fst by scaffold
#######################

## ZFST

fstlmlchwe<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fstlmlchwe.csv")

avgFst_lmlc_hwe<-colMeans(as.matrix(fstlmlchwe$FST))
sdFst_lmlc_hwe<-sd(fstlmlchwe$FST)
fstlmlchwe$ZFST<-(fstlmlchwe$FST-avgFst_lmlc_hwe)/sdFst_lmlc_hwe

fstlmlchwe90<-fstlmlchwe[fstlmlchwe$SCAF<91,]

max.pos<-sapply(1:89,function(i){(max(fstlmlchwe90$POS[fstlmlchwe90$SCAF==i]))})
add1<-c(0,cumsum(max.pos))
add2<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90)

added.length<-(add1+add2)

fstlmlchwe90$gp<-fstlmlchwe90$POS+added.length[fstlmlchwe90$SCAF]

#########################################
# for (i in 1:dim(fstlmlchwe90)[1]){
# if((fstlmlchwe90$SCAF[i] %% 2) == 0) {
# fstlmlchwe90$color[i]<-"darkorchid3"
# } else {
# fstlmlchwe90$color[i]<-"black"
# }}
#########################################

fstlmlchwe90$color<-((fstlmlchwe90$SCAF%%2)==0)
fstlmlchwe90$color[fstlmlchwe90$color==TRUE]<-"darkorchid3"
fstlmlchwe90$color[fstlmlchwe90$color==FALSE]<-"black"

## plot(fstlmlchwe90$gp,fstlmlchwe90$ZFST,xlab="Genomic Position",ylab="Z(FST)",xlim=c(0,900000000),ylim=c(-2,10),col=fstlmlchwe90$color,pch=20)
## abline(h=5,col="red",lty=2)

write.csv(fstlmlchwe90,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fstlmlchwe90.csv")
write.csv(fstlmlchwe,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fstlmlchwe_zfst.csv")






##############################
## plot ZFST and ZFST-window
##############################

fstlmlchwe90<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fstlmlchwe90.csv")

pdf("/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fstlmlchwe90.pdf")
plot(fstlmlchwe90$gp,fstlmlchwe90$ZFST,xlab="Genomic Position",ylab="Z(FST)",xlim=c(0,900000000),ylim=c(-2,10),col=fstlmlchwe90$color,pch=20)
abline(h=5,col="red",lty=2)
dev.off()
