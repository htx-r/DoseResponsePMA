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
jagsdataQuadratic <- makeJAGSDRmeta(Study_No,logRR,hayasaka_ddd,Responders,nonResponders,data=antiDep1,LQ=T)

# Anaylsis:
# a. Bayes: JAGS model
quadraticDRmetaJAGSmodel <- jags(data = jagsdataQuadratic,inits=NULL,parameters.to.save = c('beta2.pooled','beta1.pooled','tau'),model.file = modelQuadraticDRmeta,
                              n.chains=2,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 10)
traceplot(quadraticDRmetaJAGSmodel,varname='beta1.pooled') ## looks good

# b. Frequntist
quadraticDRmetaFreq <- dosresmeta::dosresmeta(formula = logRR ~ hayasaka_ddd+ I(hayasaka_ddd^2), type = type, id = Study_No,
                                                  se = selogRR, cases = Responders, n = No_randomised  , data = antiDep1,covariance = 'gl',proc = '1stage')#!!!!!!!!!!!!!!
# c. Compare Bayes vs Freq.

quadraticDRmetaJAGSmodel$BUGSoutput$mean$beta1.pooled - coef(quadraticDRmetaFreq)[1]
quadraticDRmetaJAGSmodel$BUGSoutput$mean$beta2.pooled - coef(quadraticDRmetaFreq)[2]
