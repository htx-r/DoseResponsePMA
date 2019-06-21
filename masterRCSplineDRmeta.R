# Libraries
library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force = T)

library(DoseResponseNMA)

# Data: exclude the 9 studies with NA log RR.
antiDep <-  read.csv('~/Desktop/TasnimPhD/DoseResponseNMA1/DOSEmainanalysis.csv')
#antiDep <-  read.csv("~/Google Drive/_mydrive/HTx/HTx-R/DoseResponseNMA/DOSEmainanalysis.csv")
NAstudyid <- antiDep$Study_No[is.na(antiDep$logRR)]
antiDep1 <-antiDep[!antiDep$Study_No %in% NAstudyid,]

# JAGS data
antiDep1$nonResponders <- antiDep1$No_randomised- antiDep1$Responders
jagsdataRCS<- makeJAGSDRmeta(Study_No,logRR,hayasaka_ddd,Responders,nonResponders,data=antiDep1,Splines=T,knots=c(10,20,50))

# Anaylsis:
# a. Bayes: JAGS model
rcsplineDRmetaJAGSmodel <- jags(data = jagsdataRCS,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau1','tau2'),model.file = modelRCSplineDRmeta,
                                 n.chains=2,n.iter = 1000000,n.burnin = 20000,DIC=F,n.thin = 10)
traceplot(rcsplineDRmetaJAGSmodel,varname='beta1.pooled') ## looks good YOU NEEDED MORE ITERATIONS
traceplot(rcsplineDRmetaJAGSmodel,varname='beta2.pooled') ## looks good
traceplot(rcsplineDRmetaJAGSmodel,varname='tau1') ## looks good
traceplot(rcsplineDRmetaJAGSmodel,varname='tau2') ## looks good

# b. Frequentist
knots=c(10,20,50)
rcsplineDRmetaFreq=dosresmeta(formula = logRR~rcs(hayasaka_ddd,knots), type = type, id = Study_No,
                                            se = selogRR, cases = Responders, n = No_randomised, data = antiDep1, proc='1stage')#!!!!!!!!!!!!!!

# c. Compare Bayesian and Frequentist
Bayescoeff<-unlist(rcsplineDRmetaJAGSmodel$BUGSoutput$mean[1:2])

FreqBayesCoef <- cbind(coef(rcsplineDRmetaFreq)[1:2],rcsplineDRmetaJAGSmodel$BUGSoutput$mean[1:2])
colnames(FreqBayesCoef) <- c('Freq', 'Bayes')
FreqBayesCoef


##I GOT PRETTY DIFFERENT ESTIMATES, FOR THE FIRST 2 COEFFICIENTS
#
#
#                                               Freq        Bayes
#rcs(hayasaka_ddd, c(10, 20, 50))hayasaka_ddd  0.01119403  -5.680354e-05
#rcs(hayasaka_ddd, c(10, 20, 50))hayasaka_ddd' -0.02019236 0.005018968

### d. Plot the Splines model from JAGs
nd=(rcs(seq(0,80,1),c(10,20,50)))
new.data=cbind.data.frame(d=as.vector(nd[,1]),trans.d=as.vector(nd[,2]))
new.data$RR<-exp(apply(t(t(new.data)*Bayescoeff),1,sum))
#BAYESIAN LINE IN BLUE
with(new.data,
  plot(d,RR, type = "l",xlim = c(0, 80), c(.5, 1.5),xlab="Dose",ylab="RR",main=c("Response"),col="deepskyblue", lwd=2) )
xref <- 0
#ADD frequentist line in red
with(predict(rcsplineDRmetaFreq,data.frame(hayasaka_ddd=seq(0,80,1)),xref, exp = TRUE),
 lines(get("rcs(hayasaka_ddd, knots)hayasaka_ddd"),pred,
       xlim = c(0, 80), ylim = c(.5, 1.5),xlab="Dose",ylab="RR",main=c("Response"),col="red", lwd=2)
)




