# Libraries
library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA")
library(DoseResponseNMA)

# Data: exclude the 9 studies with NA log RR.
antiDep <-  read.csv('~/Desktop/TasnimPhD/DoseResponseNMA1/DOSEmainanalysis.csv')
NAstudyid <- antiDep$Study_No[is.na(antiDep$logRR)]
antiDep1 <-antiDep[!antiDep$Study_No %in% NAstudyid,]

# JAGS data
antiDep1$nonResponders <- antiDep1$No_randomised- antiDep1$Responders
jagsdataRCS <- makeJAGSDRmeta(Study_No,logRR,hayasaka_ddd,Responders,nonResponders,data=antiDep1,LQ=F)

# Anaylsis:
# a. Bayes: JAGS model
rcsplineDRmetaJAGSmodel <- jags(data = jagsdataRCS,inits=NULL,parameters.to.save = c('beta3.pooled','beta2.pooled','beta1.pooled','tau'),model.file = modelRCSplineDRmeta,
                                 n.chains=2,n.iter = 10000,n.burnin = 200,DIC=F,n.thin = 10)
traceplot(rcsplineDRmetaJAGSmodel,varname='beta1.pooled') ## looks good

# b. Frequntist
rcsplineDRmetaFreq <- dosresmeta::dosresmeta(formula = logRR ~ rcs(hayasaka_ddd,knots=c(10,20,30)), type = type, id = Study_No,
                                            se = selogRR, cases = Responders, n = No_randomised  , data = antiDep1,covariance = 'gl',proc = '1stage')#!!!!!!!!!!!!!!
# c. Compare Bayesian and Frequentist
FreqBayesCoef <- cbind(coef(rcsplineDRmetaFreq),rcsplineDRmetaJAGSmodel$BUGSoutput$mean[1:3])
colnames(FreqBayesCoef) <- c('Freq', 'Bayes')
FreqBayesCoef

