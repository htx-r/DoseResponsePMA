# Libraries
library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA")
library(DoseResponseNMA)

# Data
antiDep <-  read.csv('~/Desktop/TasnimPhD/DoseResponseNMA1/DOSEmainanalysis.csv')
#antiDep <-  read.csv("~/Google Drive/_mydrive/HTx/HTx-R/DoseResponseNMA/DOSEmainanalysis.csv")
NAstudyid <- antiDep$Study_No[is.na(antiDep$logRR)]
antiDep1 <-antiDep[!antiDep$Study_No %in% NAstudyid,]

# JAGS data
antiDep1$nonResponders <- antiDep1$No_randomised- antiDep1$Responders
jagsdataLinear <- makeJAGSDRmeta(Study_No,logRR,hayasaka_ddd,Responders,nonResponders,data=antiDep1,Splines=F)


# Anaylsis:
# a. Bayes: JAGS model
linearDRmetaJAGSmodel <- jags(data = jagsdataLinear,inits=NULL,parameters.to.save = c('beta','beta.pooled','tau'),model.file = modelLinearDRmeta,
                              n.chains=2,n.iter = 100000,n.burnin = 20000,DIC=F,n.thin = 10)
traceplot(linearDRmetaJAGSmodel,varname='beta.pooled') ## looks good
traceplot(linearDRmetaJAGSmodel,varname='tau') ## ALWAYS check the tau for convergence
# b.Frequentist: dosresmeta

linearDRmetaFreq <- dosresmeta::dosresmeta(formula = logRR ~ hayasaka_ddd, type = type, id = Study_No,
                                               se = selogRR, cases = Responders, n = No_randomised  , data = antiDep1,covariance = 'gl',proc = '1stage')#!!!!!!!!!!!!!!

# c.Compare Bayesian and Frequentist
linearDRmetaJAGSmodel$BUGSoutput$mean$beta.pooled - coef(linearDRmetaFreq)

