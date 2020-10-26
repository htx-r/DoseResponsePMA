# Run sperate PMA dose-response for each drug (on the origional dose-level)

library(rms) # for rcs()
library(MASS) # for truehist()
library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponsePMA",force=TRUE)
library(DoseResponseNMA)
library(meta)
library(dplyr)
source('Functions needed for dosres MA antidep.R')
source('functions for antidep plots.R')

########################################
#     load data an prepare

# load and exclude single arm studies
mydata <-  read.csv('DOSEmainanalysis.csv')
antidep=mydata[mydata$exc==F,]

# compute sample per drug
drug_name <- unique(antidep$Drug)[-1]

sample_drug <- cbind( c(sample_per_drug(data_per_drug(antidep,drug_name[1]),drug_name[1]),
                   sample_per_drug(data_per_drug(antidep,drug_name[2]),drug_name[2]),
                   sample_per_drug(data_per_drug(antidep,drug_name[3]),drug_name[3]),
                   sample_per_drug(data_per_drug(antidep,drug_name[4]),drug_name[4]),
                   sample_per_drug(data_per_drug(antidep,drug_name[5]),drug_name[5])
                   ))
rownames(sample_drug) <-drug_name

# Generate separate dataset for each drug

# additional plots:  run the whole functions for antidep plots.R to get two plots
res <- sapply(unique(antidep$Drug)[-1], run3modelsDRmeta,simplify = FALSE)
names(res) <- drug_name

length(res[['fluoxetine']]$doseresORsplineBin$BUGSoutput$mean$p.drug)
data_per_drug(antidep,'sertraline')$Dose_delivered_mean

data1 <- plotdata.fun(dataDrug = data_per_drug(antidep,drug_name[1]),drug_name=drug_name[1],res[['fluoxetine']]$doseresORsplineBin)
data2 <- plotdata.fun(dataDrug = data_per_drug(antidep,drug_name[2]),drug_name=drug_name[2],res[['paroxetine']]$doseresORsplineBin)
data3 <- plotdata.fun(dataDrug = data_per_drug(antidep,drug_name[3]),drug_name=drug_name[3],res[['citalopram']]$doseresORsplineBin)
data4 <- plotdata.fun(dataDrug = data_per_drug(antidep,drug_name[4]),drug_name=drug_name[4],res[['escitalopram']]$doseresORsplineBin)
data5 <- plotdata.fun(dataDrug = data_per_drug(antidep,drug_name[5]),drug_name=drug_name[5],res[['sertraline']]$doseresORsplineBin)

plotdata <- rbind.data.frame(data1$plotdata,data2$plotdata,data3$plotdata,data4$plotdata,data5$plotdata)
df <- rbind.data.frame(data1$df,data2$df,data3$df,data4$df,data5$df)

# Absolute effect plot
myggplot(plotdata=list(plotdata=plotdata,df=df))#+annotation_custom(ggplotGrob(k2), xmin=120, xmax=180, ymin=-4,ymax=0)

#
data1 <- plotdata.fun2(dataDrug = data_per_drug(antidep,drug_name[1]),drug_name=drug_name[1],res[['fluoxetine']])
data2 <- plotdata.fun2(dataDrug = data_per_drug(antidep,drug_name[2]),drug_name=drug_name[2],res[['paroxetine']])
data3 <- plotdata.fun2(dataDrug = data_per_drug(antidep,drug_name[3]),drug_name=drug_name[3],res[['citalopram']])
data4 <- plotdata.fun2(dataDrug = data_per_drug(antidep,drug_name[4]),drug_name=drug_name[4],res[['escitalopram']])
data5 <- plotdata.fun2(dataDrug = data_per_drug(antidep,drug_name[5]),drug_name=drug_name[5],res[['sertraline']])

plotdata <- rbind.data.frame(data1$plotdata,data2$plotdata,data3$plotdata,data4$plotdata,data5$plotdata)
df <- rbind.data.frame(data1$df,data2$df,data3$df,data4$df,data5$df)

myggplot2(plotdata=list(plotdata=plotdata,df=df))
#




# sapply(unique(antidep$Drug)[-1], function(i) plotdata(dataDrug = data_per_drug(antidep,i),drug_name=i,res[[i]]$doseresORsplineBin))
# dataplot <- list()
#
# for (i in 1:length(unique(antidep$Drug)[-1])) {
#   drug_name <- unique(antidep$Drug)[-1][i]
#   dataplot[[i]] <- plotdata(dataDrug = data_per_drug(antidep,drug_name),drug_name=drug_name,res[[drug_name]]$doseresORsplineBin)
# drug_name
#   }

#** Run the analysis for antidepressant with knot on 10%, 50%, 90%

# add OR

antidep$studyid <- as.numeric(as.factor(antidep$Study_No))
antidep$nonResponders <- antidep$No_randomised- antidep$Responders
logORmat <- sapply(unique(antidep$studyid),function(i) createORreference.fun(antidep$Responders[antidep$studyid==i],antidep$No_randomised[antidep$studyid==i]),simplify = FALSE)
logORmat <- do.call(rbind,logORmat)
antidep$logOR <- c(logORmat[,1])
antidep$selogOR <- c(logORmat[,2])

# transform doses, rcs
knots= quantile(antidep$hayasaka_ddd[antidep$hayasaka_ddd!=0],c(0.10,0.50,0.90))
antidep$dose1 <- as.matrix(rcs(antidep$hayasaka_ddd,knots))[,1]
antidep$dose2 <- as.matrix(rcs(antidep$hayasaka_ddd,knots))[,2]

# jags data
jagsdataORspline<- makejagsDRmeta(studyid=studyid,logOR,dose1=dose1,dose2=dose2,cases=Responders,noncases=nonResponders,se=selogOR,type=type,data=antidep,splines=T)

# additional arguments into jagsdata to compute the absolute response for the placebo and drug arms
jagsdataORspline$np <- 58 # sum(jagsdataORspline$X1[,1]==0)
jagsdataORspline$nn <- jagsdataORspline$n[-56,]
jagsdataORspline$rr <- jagsdataORspline$r[-56,]
jagsdataORspline$new.dose <-  seq(1,80,1)
jagsdataORspline$f.new.dose <- rcspline.eval(jagsdataORspline$new.dose,knots,inclx = T)[,2]
jagsdataORspline$nd.new <- length(jagsdataORspline$new.dose)


########## ##### ##### #########################
# 1. univariate normal priors for beta1 and beta2

## Frequentist: one-stage model using dosresmeta
doseresORsplineFreq <- dosresmeta(formula=logOR~dose1+dose2, proc="1stage",id=Study_No, type=type,cases=Responders,n=No_randomised,se=selogOR,data=antidep,method = 'reml')
summary(doseresORsplineFreq)

# Bayes with normal likelihood
doseresORsplineNor <- jags.parallel(data = jagsdataORspline,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelNorSplineDRmeta,
                                    n.chains=3,n.iter = 10000,n.burnin = 4000,DIC=F,n.thin = 1)
doseresORsplineNor$BUGSoutput$summary
# Bayes with binomial likelihood
doseresORsplineBin <- jags.parallel(data = jagsdataORspline,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau','Z','p.drug','p.drug3020','p.drug4030'),model.file = modelBinSplineDRmetaOR,
                                    n.chains=3,n.iter = 10000,n.burnin = 4000,DIC=F,n.thin = 1)
doseresORsplineBin$BUGSoutput$summary[c('beta1.pooled','beta2.pooled','tau'),]
# compute the probabilities of





