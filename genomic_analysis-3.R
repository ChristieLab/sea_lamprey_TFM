hwe_lm<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data4-lm-hwe.csv")
hwe_lc<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data4-lc-hwe.csv")
hwe_ct<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data4-ct-hwe.csv")

nonhwe_lm<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data4-lm-nonhwe.csv")
nonhwe_lc<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data4-lc-nonhwe.csv")
nonhwe_ct<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/data4-ct-nonhwe.csv")

fh_lm<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fh-lm.csv")
fh_lc<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fh-lc.csv")
fh_ct<-read.csv(file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/fh-ct.csv")

fh_lm_hwe<-merge(fh_lm,hwe_lm,by.x="locus",by.y="locus")
fh_lc_hwe<-merge(fh_lc,hwe_lc,by.x="locus",by.y="locus")
fh_ct_hwe<-merge(fh_ct,hwe_ct,by.x="locus",by.y="locus")

fh_lm_nonhwe<-merge(fh_lm,nonhwe_lm,by.x="locus",by.y="locus")
fh_lc_nonhwe<-merge(fh_lc,nonhwe_lc,by.x="locus",by.y="locus")
fh_ct_nonhwe<-merge(fh_ct,nonhwe_ct,by.x="locus",by.y="locus")

het_lm_hwe<-data.frame(fh_lm_hwe$locus,fh_lm_hwe$Hexp,fh_lm_hwe$Hobs)
het_lc_hwe<-data.frame(fh_lc_hwe$locus,fh_lc_hwe$Hexp,fh_lc_hwe$Hobs)
het_ct_hwe<-data.frame(fh_ct_hwe$locus,fh_ct_hwe$Hexp,fh_ct_hwe$Hobs)

het_lm_nonhwe<-data.frame(fh_lm_nonhwe$locus,fh_lm_nonhwe$Hexp,fh_lm_nonhwe$Hobs)
het_lc_nonhwe<-data.frame(fh_lc_nonhwe$locus,fh_lc_nonhwe$Hexp,fh_lc_nonhwe$Hobs)
het_ct_nonhwe<-data.frame(fh_ct_nonhwe$locus,fh_ct_nonhwe$Hexp,fh_ct_nonhwe$Hobs)

write.csv(het_lm_hwe,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/het_lm_hwe.csv")
write.csv(het_lc_hwe,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/het_lc_hwe.csv")
write.csv(het_ct_hwe,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/het_ct_hwe.csv")

write.csv(het_lm_nonhwe,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/het_lm_nonhwe.csv")
write.csv(het_lc_nonhwe,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/het_lc_nonhwe.csv")
write.csv(het_ct_nonhwe,file="/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/het_ct_nonhwe.csv")

pdf("/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/het_lm_hwe.pdf")
plot(het_lm_hwe$fh_lm_hwe.Hexp,het_lm_hwe$fh_lm_hwe.Hobs,xlab="Hexp",ylab="Hobs",xlim=c(0,1),ylim=c(0,1),main="LM HWE")
dev.off()

pdf("/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/het_lc_hwe.pdf")
plot(het_lc_hwe$fh_lc_hwe.Hexp,het_lc_hwe$fh_lc_hwe.Hobs,xlab="Hexp",ylab="Hobs",xlim=c(0,1),ylim=c(0,1),main="LC HWE")
dev.off()

pdf("/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/het_ct_hwe.pdf")
plot(het_ct_hwe$fh_ct_hwe.Hexp,het_ct_hwe$fh_ct_hwe.Hobs,xlab="Hexp",ylab="Hobs",xlim=c(0,1),ylim=c(0,1),main="CT HWE")
dev.off()

pdf("/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/het_lm_nonhwe.pdf")
plot(het_lm_nonhwe$fh_lm_nonhwe.Hexp,het_lm_nonhwe$fh_lm_nonhwe.Hobs,xlab="Hexp",ylab="Hobs",xlim=c(0,1),ylim=c(0,1),main="LM nonHWE")
dev.off()

pdf("/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/het_lc_nonhwe.pdf")
plot(het_lc_nonhwe$fh_lc_nonhwe.Hexp,het_lc_nonhwe$fh_lc_nonhwe.Hobs,xlab="Hexp",ylab="Hobs",xlim=c(0,1),ylim=c(0,1),main="LC nonHWE")
dev.off()

pdf("/scratch/snyder/y/yin168/gatksep/rnavdna/pmarinus1/test/dp5/popmuscle-2/het_ct_nonhwe.pdf")
plot(het_ct_nonhwe$fh_ct_nonhwe.Hexp,het_ct_nonhwe$fh_ct_nonhwe.Hobs,xlab="Hexp",ylab="Hobs",xlim=c(0,1),ylim=c(0,1),main="CT nonHWE")
dev.off()
