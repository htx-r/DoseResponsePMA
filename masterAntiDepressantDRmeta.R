# Libraries
library(rms)
library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force=TRUE)
library(DoseResponseNMA)

##################################################################
#       Data
##################################################################


antidep <-  read.csv('~/Desktop/TasnimPhD/DoseResponseNMA/DoseResponseNMA/DOSEmainanalysis.csv')
mymoredata=antidep[antidep$exc==F,]

max.nd <- max(as.numeric(table(mymoredata$Study_No)))
ns <- length(unique(mymoredata$Study_No))
nd <- as.numeric(table(mymoredata$Study_No))-1

mymoredata$studyid <- as.numeric(as.factor(mymoredata$Study_No))
mymoredata$nonResponders <- mymoredata$No_randomised- mymoredata$Responders

# Spline transformation
knots = c(10,20,50)
mymoredata$dose1 <- as.matrix(rcs(mymoredata$hayasaka_ddd,knots))[,1]
mymoredata$dose2 <- as.matrix(rcs(mymoredata$hayasaka_ddd,knots))[,2]

# transform data to make jags data
jagsdata<- makejagsDRmeta(studyid=studyid,logRR,dose1=dose1,dose2=dose2,cases=Responders,noncases=nonResponders,se=selogRR,type=type,data=mymoredata,Splines=T,new.dose.range = c(5,10))

##################################################################
#     ANALYSIS
#################################################################

## 1.Frequentist
doseresRRsplineFreq <- dosresmeta(formula=logRR~rcs(hayasaka_ddd,knots), proc="1stage",id=Study_No, type='ci',cases=Responders,n=No_randomised,se=selogRR,data=mymoredata)


# 2. Bayes with normal likelihood
doseresRRsplineNor <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelNorSplineDRmeta,
                              n.chains=2,n.iter = 100000,n.burnin = 20000,DIC=F,n.thin = 10)


# 3. Bayes with binomial likelihood
doseresRRsplineBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','beta3.pooled'),model.file = modelBinSplineDRmetaRR,
                                          n.chains=2,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 1)

#%% combine the three results
cbind(bayesBin=doseresRRsplineBin$BUGSoutput$mean[1:2],bayesNor=doseresRRsplineNor$BUGSoutput$mean[1:2],Freq=coef(doseresRRsplineFreq))

###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Linear: 1. odds ratio (OR) 2. risk ratio (RR)
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
newdata=data.frame(hayasaka_ddd=seq(0,80,1))
xref=min(mymoredata$hayasaka_ddd)
with(predict(doseresRR, newdata,xref, exp = TRUE), {
  plot(get("rcs(hayasaka_ddd, knots)hayasaka_ddd"),pred, log = "y", type = "l",
       xlim = c(0, 80), ylim = c(.5, 5),xlab="Dose",ylab="RR",main=c("Response"))
  matlines(get("rcs(hayasaka_ddd, knots)hayasaka_ddd"),cbind(ci.ub,ci.lb),col=1,lty="dashed")})
with(mymoredata,rug(hayasaka_ddd, quiet = TRUE))
