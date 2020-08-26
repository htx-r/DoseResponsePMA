# Run sperate PMA dose-response for each drug (on the origional dose-level)

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

# Generate separate dataset for each drug
data_per_drug<-function(data, name_of_drug){
  data2<-antidep
  data2$count<-1
  data_drug<-data2%>%filter(data2$Drug==name_of_drug)%>%select(Study_No)
  data2 <-data2 %>% filter(Study_No%in%data_drug$Study_No)%>%filter(Drug=='placebo'||Drug==name_of_drug)
  # data_count<-data2 %>% group_by(Study_No, Study_No) %>% summarise(arms=sum(count)) %>% filter(arms>2)
  # data2<-data2 %>% filter(Study_No%in%data_count$Study_No)
data2
}

# Dose: cubic spline transformation
run3modelsDRmeta <- function(drug_name){
  drug_name <- 'escitalopram'
dataDrug <- data_per_drug(antidep,drug_name)
knots = quantile(dataDrug$Dose_delivered_mean,probs=c(0.10,0.50,0.90))
dataDrug$dose1 <- as.matrix(rcs(dataDrug$Dose_delivered_mean,knots))[,1]
dataDrug$dose2 <- as.matrix(rcs(dataDrug$Dose_delivered_mean,knots))[,2]

# data_per_drug(antidep,'citalopram')
# data_per_drug(antidep,'escitalopram')
# data_per_drug(antidep,'sertraline')

#

#
dataDrug$studyid <- as.numeric(as.factor(dataDrug$Study_No))
dataDrug$nonResponders <- dataDrug$No_randomised- dataDrug$Responders

# Response:  odds ratio
logORmat <- sapply(unique(dataDrug$studyid),function(i) createORreference.fun(dataDrug$Responders[dataDrug$studyid==i],dataDrug$No_randomised[dataDrug$studyid==i]),simplify = FALSE)
logORmat <- do.call(rbind,logORmat)
dataDrug$logOR <- c(logORmat[,1])
dataDrug$selogOR <- c(logORmat[,2])




# transform data to jags format
jagsdataORspline<- makejagsDRmeta(studyid=studyid,logOR,dose1=dose1,dose2=dose2,cases=Responders,noncases=nonResponders,se=selogOR,type=type,data=dataDrug,splines=T)

# additional arguments into jagsdata to compute the absolute response for the placebo and drug arms
jagsdataORspline$np <- sum(jagsdataORspline$X1[,1]==0)#length(unique(dataDrug$Study_No))#58 # sum(jagsdataORspline$X1[,1]==0)
jagsdataORspline$nn <- jagsdataORspline$n
jagsdataORspline$rr <- jagsdataORspline$r
jagsdataORspline$new.dose <-1:max(dataDrug$Dose_delivered_mean)
jagsdataORspline$f.new.dose <- rcspline.eval(jagsdataORspline$new.dose,knots,inclx = T)[,2]
jagsdataORspline$nd.new <- length(jagsdataORspline$new.dose)


## Frequentist: one-stage model using dosresmeta
doseresORsplineFreq <- dosresmeta(formula=logOR~dose1+dose2, proc="1stage",id=Study_No, type=type,cases=Responders,n=No_randomised,se=selogOR,data=dataDrug,method = 'reml')

# Bayes with normal likelihood
doseresORsplineNor <- jags.parallel(data = jagsdataORspline,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelNorSplineDRmeta,
                                    n.chains=3,n.iter = 1000,n.burnin = 400,DIC=F,n.thin = 1)
# Bayes with binomial likelihood
doseresORsplineBin <- jags.parallel(data = jagsdataORspline,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau','Z','p.drug','p.drug3020','p.drug4030'),model.file = modelBinSplineDRmetaOR,
                                    n.chains=3,n.iter = 1000,n.burnin = 400,DIC=F,n.thin = 1)

return(list(doseresORsplineBin=doseresORsplineBin,doseresORsplineNor=doseresORsplineNor, doseresORsplineFreq=doseresORsplineFreq))
}

# additional plots:  run the whole functions for antidep plots.R to get two plots
res_paro <- run3modelsDRmeta('sertraline')
res <- sapply(unique(antidep$Drug)[-1], run3modelsDRmeta,simplify = FALSE)
names(res) <- unique(antidep$Drug)[-1]

