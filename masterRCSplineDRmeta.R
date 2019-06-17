# Libraries
library(R2jags)
library(devtools)
install_github("htx-r/DoseResponseNMA")
library(DoseResponseNMA)

# Data: exclude the 9 studies with NA log RR.
antiDep <-  read.csv('~/Desktop/TasnimPhD/DoseResponseNMA1/DOSEmainanalysis.csv')
NAstudyid <- antiDep$Study_No[is.na(antiDep$logRR)]
antiDep1 <-antiDep[!antiDep$Study_No %in% NAstudyid,]

# Anaylsis: JAGS model
rcsplineDRmetaJAGSmodel <- jags(data = makeJAGSRCSplineDRmeta(antiDep1),inits=NULL,parameters.to.save = c('beta2.pooled','beta1.pooled','tau'),model.file = modelRCSplineDRmeta,
                                 n.chains=2,n.iter = 10000,n.burnin = 200,DIC=F,n.thin = 10)
traceplot(rcsplineDRmetaJAGSmodel,varname='beta1.pooled') ## looks good


