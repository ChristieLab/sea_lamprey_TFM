## install packages needed in survival analysis

install.packages("survival")
install.packages("KMsurv")
install.packages("survsim")
install.packages("broom")
install.packages("ggfortify")


library(survival)
library(KMsurv)
library(survsim)
library(broom)
library(plyr)
library(dplyr)
library(ggfortify)
library(ggplot2)
library(survminer)

## working dir
setwd("E:/Purdue/pmarinus/TFM/Fig-v4/Fig2/")

## load and visualize data
data = read.table("tox3_surv_analysis.csv", header=TRUE, as.is=TRUE, sep=",")
attach(data)
names(data)
data <- transform(data,pop=revalue(pop, c("lm" = "Michigan", "lc" = "Champlain", "ct" = "Connecticut"))) #rename locations for figure purposes
head(data)
dim(data)
data$pop = as.factor(data$pop)
data$pop <- relevel(data$pop, "Michigan")

## surv. analysis for all three populations
sur.obj2 <- Surv(time = time , event = delta) #Create survival object, event = status (died =1  or censored =0) 

srvfit2 = survfit(sur.obj2 ~ pop, data = data, conf.type = "plain") #Compare survival distributions of populations
str(srvfit2)
summary(srvfit2)

pdf("E:/Purdue/pmarinus/TFM/Fig-v4/Fig2/survival_plot.pdf", width=3.6, height=4, useDingbats=FALSE, onefile=FALSE)
ggsurvplot(srvfit2,data=data, conf.int=T, censor=T, surv.median.line="hv", ggtheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour = "black")), palette=c("navy","dodgerblue2","cadetblue2"), xlab="", ylab="", main="", font.tickslab=c(12,"plain","black"), legend=c(0.26,0.17), legend.labs=c("Lake Michigan","Lake Champlain","Connecticut River"), legend.title="")
dev.off()
