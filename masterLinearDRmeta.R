# Libraries
library(R2jags)
library(devtools)
install_github("htx-r/DoseResponseNMA")
library(DoseResponseNMA)

# Data
antiDep <-  read.csv('~/Desktop/TasnimPhD/DoseResponseNMA/DOSEmainanalysis.csv')
NAstudyid <- antiDep$Study_No[is.na(antiDep$logRR)]
antiDep1 <-antiDep[!antiDep$Study_No %in% NAstudyid,]

# Anaylsis: JAGS model
linearDRmetaJAGSmodel <- jags(data = makeJAGSLinearQuadraticDRmeta(antiDep1),inits=NULL,parameters.to.save = c('beta','beta.pooled','tau'),model.file = modelLinearDRmeta,
                              n.chains=2,n.iter = 10000,n.burnin = 200,DIC=F,n.thin = 10)
traceplot(linearDRmetaJAGSmodel,varname='beta.pooled') ## looks good


