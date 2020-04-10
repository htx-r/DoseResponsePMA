##################################################################################
#         Master analysis for dose-response meta-analysis of Antidepressant
##################################################################################

# load libraries
library(rms) # for rcs()
library(MASS) # for truehist()
library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponsePMA",force=TRUE)
library(DoseResponseNMA)
library(meta)
source('Functions needed for dosres MA antidep.R')

########################################
#     load data an prepare

# load and exclude single arm studies
mydata <-  read.csv('DOSEmainanalysis.csv')
antidep=mydata[mydata$exc==F,]


(as.matrix(table(antidep$Study_No,antidep$Drug)))
#
antidep$studyid <- as.numeric(as.factor(antidep$Study_No))
antidep$nonResponders <- antidep$No_randomised- antidep$Responders


# Response:  odds ratio
logORmat <- sapply(unique(antidep$studyid),function(i) createORreference.fun(antidep$Responders[antidep$studyid==i],antidep$No_randomised[antidep$studyid==i]),simplify = FALSE)
logORmat <- do.call(rbind,logORmat)
antidep$logOR <- c(logORmat[,1])
antidep$selogOR <- c(logORmat[,2])

# Dose: cubic spline transformation
# antidep$hayasaka_ddd <- antidep$hayasaka_ddd/100
knots = c(10,20,50)
antidep$dose1 <- as.matrix(rcs(antidep$hayasaka_ddd,knots))[,1]
antidep$dose2 <- as.matrix(rcs(antidep$hayasaka_ddd,knots))[,2]

# transform data to jags format
jagsdataORspline<- makejagsDRmeta(studyid=studyid,logOR,dose1=dose1,dose2=dose2,cases=Responders,noncases=nonResponders,se=selogOR,type=type,data=antidep,splines=T)


# jagsdataORspline$nstudy <- 60-sapply(1:6, function(i) sum(as.matrix(table(antidep$Study_No,antidep$Drug))[,i]==0))
# names(nstudy) <- levels(antidep$Drug)

jagsdataORspline$drug <- as.numeric(unlist(sapply(1:60, function(i) unique(antidep$Drug[antidep$studyid==i][antidep$Drug[antidep$studyid==i]!='placebo']))))



# Bayes with binomial likelihood
doseresORsplineBindrugcluster <- jags.parallel(data = jagsdataORspline,inits=NULL,parameters.to.save = c('b1','b2','tau','sigma'),model.file = modelBinSplineDRmetaORdrugcluster,
                                    n.chains=3,n.iter = 100000,n.burnin = 20000,DIC=F,n.thin = 1)
doseresORsplineBindrugcluster$BUGSoutput$summary[c('b1','b2','tau','sigma'),]

doseresORsplineBindrugclusterBiv <- jags.parallel(data = jagsdataORspline,inits=NULL,parameters.to.save = c('b1','b2','tau','sigma','rho','sigma','rhoD'),model.file = modelBinSplineDRmetaORdrugclusterBiv,
                                               n.chains=3,n.iter = 100000,n.burnin = 20000,DIC=F,n.thin = 1)
doseresORsplineBindrugclusterBiv$BUGSoutput$summary


