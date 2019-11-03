# Libraries
library(rms) # for rcs()
library(MASS) # for truehist()
library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force=TRUE)
library(DoseResponseNMA)
library(meta)
devAskNewPage(ask=F)

###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Splines: risk ratio (RR) + odds ratio (OR)
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##################################################################
#       Data
##################################################################

# read data
mydata <-  read.csv('~/Google Drive/DoseResponseNMA/DoseResponseNMA/DOSEmainanalysis.csv')
antidep=mydata[mydata$exc==F,]

#
antidep$studyid <- as.numeric(as.factor(antidep$Study_No))
antidep$nonResponders <- antidep$No_randomised- antidep$Responders

# compute relative odds ratio for all studies
createORreference.fun=function(r,n)
{

  logOR=c(0)
  selogOR=c(NA)

  for(i in 2:c(length(n)))
  {
    calculate=metabin(r[i],n[i],r[1],n[1],sm="OR")
    logOR=c(logOR,calculate$TE)
    selogOR=c(selogOR,calculate$seTE)

  }
  return(cbind(logOR=logOR,selogOR=selogOR))
}
logORmat <- sapply(unique(antidep$studyid),function(i) createORreference.fun(antidep$Responders[antidep$studyid==i],antidep$No_randomised[antidep$studyid==i]),simplify = FALSE)
logORmat <- do.call(rbind,logORmat)
antidep$logOR <- c(logORmat[,1])
antidep$selogOR <- c(logORmat[,2])

# Spline transformation
knots = c(10,20,50)
antidep$dose1 <- as.matrix(rcs(antidep$hayasaka_ddd,knots))[,1]
antidep$dose2 <- as.matrix(rcs(antidep$hayasaka_ddd,knots))[,2]

# transform data into jags format
jagsdataRRspline<- makejagsDRmeta(studyid=studyid,logRR,dose1=dose1,dose2=dose2,cases=Responders,noncases=nonResponders,se=selogRR,type=type,data=antidep,Splines=T,new.dose.range = c(5,10))
jagsdataORspline<- makejagsDRmeta(studyid=studyid,logOR,dose1=dose1,dose2=dose2,cases=Responders,noncases=nonResponders,se=selogOR,type=type,data=antidep,Splines=T,new.dose.range = c(5,10))
jagsdataRRlinear<- makejagsDRmeta(studyid=studyid,logRR,dose1=hayasaka_ddd,dose2=NULL,cases=Responders,noncases=nonResponders,se=selogRR,type=type,data=antidep,Splines=F,new.dose.range = c(5,10))
jagsdataORlinear<- makejagsDRmeta(studyid=studyid,logOR,dose1=hayasaka_ddd,dose2=NULL,cases=Responders,noncases=nonResponders,se=selogOR,type=type,data=antidep,Splines=F,new.dose.range = c(5,10))

# new dose for predictions and plot
new.dose <- seq(0,80,1)
new.dose1 <- c(rcs(new.dose,knots)[,1])
new.dose2 <- c(rcs(new.dose,knots)[,2])

##################################################################
#     ANALYSIS: spline RR
#################################################################

## 1.Frequentist
doseresRRsplineFreq <- dosresmeta(formula=logRR~rcs(hayasaka_ddd,knots), proc="1stage",id=Study_No, type='ci',cases=Responders,n=No_randomised,se=selogRR,data=antidep,method = 'reml')


# 2. Bayes with normal likelihood
doseresRRsplineNor <- jags.parallel(data = jagsdataRRspline,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelNorSplineDRmeta,
                              n.chains=3,n.iter = 1000000,n.burnin = 20000,DIC=F,n.thin = 10)


# 3. Bayes with binomial likelihood
doseresRRsplineBin <- jags.parallel(data = jagsdataRRspline,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaRR,
                                    n.chains=3,n.iter = 1000000,n.burnin = 20000,DIC=F,n.thin = 10)

#%% combine the three results
beta1fRR <- coef(doseresRRsplineFreq)[1]
beta2fRR <- coef(doseresRRsplineFreq)[2]

beta1nRR <- doseresRRsplineNor$BUGSoutput$mean$beta1.pooled
beta2nRR <- doseresRRsplineNor$BUGSoutput$mean$beta2.pooled

beta1bRR <- doseresRRsplineBin$BUGSoutput$mean$beta1.pooled
beta2bRR <- doseresRRsplineBin$BUGSoutput$mean$beta2.pooled

taubRR <- doseresRRsplineBin$BUGSoutput$mean$tau
taunRR <- doseresRRsplineNor$BUGSoutput$mean$tau
taufRR <- NA

cbind(bayesBin=c(beta1bRR,beta2bRR,taubRR),bayesNor=c(beta1nRR,beta2nRR,taunRR),Freq=c(beta1fRR,beta2fRR,taufRR))

 ##  check autocorraltion
 # normal
acf(doseresRRsplineNor$BUGSoutput$sims.array[,1,'beta1.pooled'],main='Normal: beta1.pooled')
acf(doseresRRsplineNor$BUGSoutput$sims.array[,1,'beta2.pooled'],main='Normal: beta2.pooled')

 # binomial
acf(doseresRRsplineBin$BUGSoutput$sims.array[,1,'beta1.pooled'],main='Binomial: beta1.pooled')
acf(doseresRRsplineBin$BUGSoutput$sims.array[,1,'beta2.pooled'],main='Binomial: beta2.pooled')

###%%%% it is highly autocorrelated

## check convergance


# normal: it is not converged!!!!!!!!!!!!!!!!!!!!!!
traceplot(doseresRRsplineNor$BUGSoutput,varname='beta1.pooled')
traceplot(doseresRRsplineNor$BUGSoutput,varname='beta2.pooled')

# binomial: it is converged
traceplot(doseresRRsplineBin$BUGSoutput,varname='beta1.pooled')
traceplot(doseresRRsplineBin$BUGSoutput,varname='beta2.pooled')

## beta's distributions

# normal
beta1.pooled.sim.norRR <- c(doseresRRsplineNor$BUGSoutput$sims.array[,,'beta1.pooled']) ## chain 1+2+3 for beta1.pooled
beta2.pooled.sim.norRR <- c(doseresRRsplineNor$BUGSoutput$sims.array[,,'beta2.pooled']) ## chain 1+2+3 for beta2.pooled

truehist(beta1.pooled.sim.norRR)
truehist(beta2.pooled.sim.norRR)

# binomial
beta1.pooled.sim.binRR <- c(doseresRRsplineBin$BUGSoutput$sims.array[,,'beta1.pooled']) ## chain 1+2+3 for beta1.pooled
beta2.pooled.sim.binRR <- c(doseresRRsplineBin$BUGSoutput$sims.array[,,'beta2.pooled']) ## chain 1+2+3 for beta2.pooled

truehist(beta1.pooled.sim.binRR)
truehist(beta2.pooled.sim.binRR)

# plot the model based on the three apporaches: freq, bayes normal and bayes binomial

plot(new.dose1,exp(beta1fRR*new.dose1+beta2fRR*new.dose2),col=1,type='l',ylim = c(0.5,2)
     ,las=1,ylab='RR',xlab='dose',lwd=3,cex.axis=1.4,cex.lab=1.4) #  freq
lines(exp(beta1nRR*new.dose1+beta2nRR*new.dose2),col=2,lwd=3) # bayes normal
lines(exp(beta1bRR*new.dose1+beta2bRR*new.dose2),col=3,lwd=3) # bayes binomial
legend(10,2,legend=c('Freq', 'normal Bayes', 'binomial Bayes'),col=1:3,horiz = T,lty=1,
       bty='n',xjust = 0.1,cex = 1,lwd=3)

##################################################################
#     ANALYSIS: spline OR
#################################################################


## 1.Frequentist
doseresORsplineFreq <- dosresmeta(formula=logOR~rcs(hayasaka_ddd,knots), proc="1stage",id=Study_No, type=type,cases=Responders,n=No_randomised,se=selogOR,data=antidep,method = 'reml')


# 2. Bayes with normal likelihood
doseresORsplineNor <- jags.parallel(data = jagsdataORspline,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelNorSplineDRmeta,
                                    n.chains=3,n.iter = 1000000,n.burnin = 20000,DIC=F,n.thin = 10)


# 3. Bayes with binomial likelihood
doseresORsplineBin <- jags.parallel(data = jagsdataORspline,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinSplineDRmetaOR,
                                    n.chains=3,n.iter = 1000000,n.burnin = 20000,DIC=F,n.thin = 10)

#%% combine the three results
beta1fOR <- coef(doseresORsplineFreq)[1]
beta2fOR <- coef(doseresORsplineFreq)[2]

beta1nOR <- doseresORsplineNor$BUGSoutput$mean$beta1.pooled
beta2nOR <- doseresORsplineNor$BUGSoutput$mean$beta2.pooled

beta1bOR <- doseresORsplineBin$BUGSoutput$mean$beta1.pooled
beta2bOR <- doseresORsplineBin$BUGSoutput$mean$beta2.pooled

taubOR <- doseresORsplineBin$BUGSoutput$mean$tau
taunOR <- doseresORsplineNor$BUGSoutput$mean$tau
taufOR <- NA


cbind(bayesBin=c(beta1bOR,beta2bOR,taubOR),bayesNor=c(beta1nOR,beta2nOR,taunOR),Freq=c(beta1fOR,beta2fOR,taufOR))

## check convergance
# binomial: converged
traceplot(doseresORsplineBin$BUGSoutput,varname='beta1.pooled')

traceplot(doseresORsplineBin$BUGSoutput,varname='beta2.pooled')
traceplot(doseresORsplineBin$BUGSoutput,varname='tau')

# normal: not converged!!!!!!!!!!!!!!!!!!!!!!
traceplot(doseresORsplineNor$BUGSoutput,varname='beta1.pooled')
traceplot(doseresORsplineNor$BUGSoutput,varname='beta2.pooled')
traceplot(doseresORsplineNor$BUGSoutput,varname='tau')

## beta's distributions
beta1.pooled.sim.norOR <- c(doseresORsplineNor$BUGSoutput$sims.array[,,'beta1.pooled']) ## chain 1+2+3 for beta1.pooled
beta2.pooled.sim.norOR <- c(doseresORsplineNor$BUGSoutput$sims.array[,,'beta2.pooled']) ## chain 1+2+3 for beta2.pooled

truehist(beta1.pooled.sim.norOR)
truehist(beta2.pooled.sim.norOR)

beta1.pooled.sim.binOR <- c(doseresORsplineBin$BUGSoutput$sims.array[,,'beta1.pooled']) ## chain 1+2+3 for beta1.pooled
beta2.pooled.sim.binOR <- c(doseresORsplineBin$BUGSoutput$sims.array[,,'beta2.pooled']) ## chain 1+2+3 for beta2.pooled

truehist(beta1.pooled.sim.binOR)
truehist(beta2.pooled.sim.binOR)

# plot the model based on the three apporaches: freq, bayes normal and bayes binomial
plot(new.dose1,exp(beta1fOR*new.dose1+beta2fOR*new.dose2),col=1,type='l',ylim = c(0.5,3),
     las=1,ylab='OR',xlab='dose',lwd=2) #  freq
lines(exp(beta1nOR*new.dose1+beta2nOR*new.dose2),col=2,lwd=2) # bayes normal
lines(exp(beta1bOR*new.dose1+beta2bOR*new.dose2),col=3,lwd=2) # bayes binomial
legend('topleft',legend=c('Freq', 'normalBayes', 'binomialBayes'),col=1:3,horiz = T,lty=1,
       bty='n',xjust = 0,cex = 0.8,lwd=2)


#






###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Linear: risk ratio (RR) + odds ratio (OR)
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


##################################################################
#     ANALYSIS: RR linear
#################################################################

## 1.Frequentist
doseresRRlinearFreq <- dosresmeta(formula=logRR~hayasaka_ddd, proc="1stage",id=Study_No, type='ci',cases=Responders,n=No_randomised,se=selogRR,data=antidep,method = 'reml')


# 2. Bayes with normal likelihood
doseresRRlinearNor <- jags.parallel(data = jagsdataRRlinear,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelNorLinearDRmeta,
                                    n.chains=3,n.iter = 1000000,n.burnin = 20000,DIC=F,n.thin = 10)


# 3. Bayes with binomial likelihood
doseresRRlinearBin <- jags.parallel(data = jagsdataRRlinear,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelBinLinearDRmetaRR,
                                    n.chains=3,n.iter = 1000000,n.burnin = 20000,DIC=F,n.thin = 10)

#%% combine the three results
betafRR <- coef(doseresRRlinearFreq)[1]

betanRR <- doseresRRlinearNor$BUGSoutput$mean$beta.pooled

betabRR <- doseresRRlinearBin$BUGSoutput$mean$beta.pooled

taunRRL <- doseresRRlinearNor$BUGSoutput$mean$tau

taubRRL <- doseresRRlinearBin$BUGSoutput$mean$tau

taufRRL <- NA
cbind(bayesBin=c(betabRR,taubRRL),bayesNor=c(betanRR,taunRRL),Freq=c(betafRR,taufRRL))

## check convergance

# normal:
traceplot(doseresRRlinearNor$BUGSoutput,varname='beta.pooled')

# binomial: it is converged
traceplot(doseresRRlinearBin$BUGSoutput,varname='beta.pooled')

## beta's distributions
beta.pooled.sim.norRR <- c(doseresRRlinearNor$BUGSoutput$sims.array[,,'beta.pooled']) ## chain 1+2+3 for beta1.pooled

truehist(beta.pooled.sim.norRR)

beta.pooled.sim.binRR <- c(doseresRRlinearBin$BUGSoutput$sims.array[,,'beta.pooled']) ## chain 1+2+3 for beta1.pooled

truehist(beta.pooled.sim.binRR)

# plot the model based on the three apporaches: freq, bayes normal and bayes binomial

plot(new.dose1,exp(betafRR*new.dose1),col=1,type='l') #  freq
lines(exp(betanRR*new.dose1),col=2) # bayes normal
lines(exp(betabRR*new.dose1),col=3) # bayes binomial

#


##################################################################
#     ANALYSIS: OR linear
#################################################################

## 1.Frequentist
doseresORlinearFreq <- dosresmeta(formula=logOR~hayasaka_ddd, proc="1stage",id=Study_No, type='cc',cases=Responders,n=No_randomised,se=selogOR,data=antidep,method = 'reml')


# 2. Bayes with normal likelihood
doseresORlinearNor <- jags.parallel(data = jagsdataORlinear,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelNorLinearDRmeta,
                                    n.chains=3,n.iter = 1000000,n.burnin = 20000,DIC=F,n.thin = 10)


# 3. Bayes with binomial likelihood
doseresORlinearBin <- jags.parallel(data = jagsdataORlinear,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelBinLinearDRmetaOR,
                                    n.chains=3,n.iter = 1000000,n.burnin = 20000,DIC=F,n.thin = 10)

#%% combine the three results
betafOR <- coef(doseresORlinearFreq)[1]

betanOR <- doseresORlinearNor$BUGSoutput$mean$beta.pooled

betabOR <- doseresORlinearBin$BUGSoutput$mean$beta.pooled

taunORL <- doseresORlinearNor$BUGSoutput$mean$tau

taubORL <- doseresORlinearBin$BUGSoutput$mean$tau

taufORL <- NA

cbind(bayesBin=c(betabOR,taubORL),bayesNor=c(betanOR,taunORL),Freq=c(betafOR,taufORL))

## check convergance
# normal: converge
traceplot(doseresORlinearNor$BUGSoutput,varname='beta.pooled')

# binomial: converge
traceplot(doseresORlinearBin$BUGSoutput,varname='beta.pooled')

## beta's distributions
beta.pooled.sim.norOR <- c(doseresORlinearNor$BUGSoutput$sims.array[,,'beta.pooled']) ## chain 1+2+3 for beta1.pooled

truehist(beta.pooled.sim.norOR)

beta.pooled.sim.binOR <- c(doseresORlinearBin$BUGSoutput$sims.array[,,'beta.pooled']) ## chain 1+2+3 for beta1.pooled

truehist(beta.pooled.sim.binOR)

# plot the model based on the three apporaches: freq, bayes normal and bayes binomial

plot(new.dose1,exp(betafOR*new.dose1),col=1,type='l',ylim = c(1,3)) #  freq
lines(new.dose1,exp(betanOR*new.dose1),col=2) # bayes normal
lines(new.dose1,exp(betabOR*new.dose1),col=3) # bayes binomial
#
#
#
#
#
#
# newdata=data.frame(hayasaka_ddd=seq(0,80,1))
# xref=min(antidep$hayasaka_ddd)
# with(predict(doseresRRsplineFreq, newdata,xref, exp = TRUE), {
#   plot(get("rcs(hayasaka_ddd, knots)hayasaka_ddd"),pred, log = "y", type = "l",
#        xlim = c(0, 80), ylim = c(.5, 5),xlab="Dose",ylab="RR",main=c("Response"))
# matlines(get("rcs(hayasaka_ddd, knots)hayasaka_ddd"),cbind(ci.ub,ci.lb),col=1,lty="dashed")
# })
# with(antidep,rug(hayasaka_ddd, quiet = TRUE))
# dose1 <- c(rcs(newdata$hayasaka_ddd,knots)[,1])
# dose2 <- c(rcs(newdata$hayasaka_ddd,knots)[,2])
#
# lines(exp(beta1fRR*dose1+beta2fRR*dose2),col=2)
