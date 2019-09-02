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
antiDep1$studyid <- as.numeric(as.factor(antiDep1$Study_No))
dose  <- rcs(antiDep1$hayasaka_ddd)

antiDep1$dose1 <- dose[1,]
antiDep1$dose2 <- dose[2,]
jagsdataLinear <- makeJAGSDRmeta(studyid=studyid,logrr =  logRR,dose1 =  dose1,dose2=dose[2,],cases=Responders,noncases=nonResponders,data=antiDep1,Splines=T)
jagsdataLinear$new.dose <- c(5,10,15)
jagsdataLinear$new.n <- length(jagsdataLinear$new.dose)

# Anaylsis:
# a.Frequentist: dosresmeta

linearDRmetaFreq <-dosresmeta(formula = logRR ~ rcs(hayasaka_ddd,knots=c(10,20,30)), type = type, id = Study_No,
                                           se = selogRR, cases = Responders, n = No_randomised  , data = antiDep1,covariance = 'gl',proc = '1stage')#!!!!!!!!!!!!!!

# b. Bayes: JAGS model

prec <-  sapply(linearDRmetaFreq$Slist,solve,simplify = F)

tncomp <- sum(as.numeric(table(antiDep1$Study_No))-1) ## total number of non-zero comparisons
max.nd <- max(as.numeric(table(antiDep1$Study_No)))

precmat <- matrix(NA,tncomp,max.nd-1)
b <- vector()

for (i in 1:ns) {
      b[1] <- 0
      nd[i] <- as.numeric(table(antiDep1$Study_No)[i])-1
      precmat[(b[i]+1):(b[i]+nd[i]),1:(nd[i])] <- prec[[i]]
      b[i+1] <- b[i]+ nd[i]
      precmat
    }
jagsdataLinear$prec <- precmat
jagsdataLinear$nd <- jagsdataLinear$nd-1
linearDRmetaJAGSmodel <- jags(data = jagsdataLinear,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelLinearDRmeta,
                              n.chains=2,n.iter = 100000,n.burnin = 20000,DIC=F,n.thin = 10)
traceplot(linearDRmetaJAGSmodel,varname='beta.pooled') ## looks good
traceplot(linearDRmetaJAGSmodel,varname='tau') ## ALWAYS check the tau for convergence

# c.Compare Bayesian and Frequentist
linearDRmetaJAGSmodel$BUGSoutput$mean$beta.pooled - coef(linearDRmetaFreq)
jagsdataLinearBin <- makeBinomialJAGSDRmeta(studyid=Study_No,dose1 = hayasaka_ddd,dose2=NULL,cases=Responders,noncases=nonResponders,data=antiDep1,Splines=F)

linearDRmetaJAGSmodelBin <- jags.parallel(data = jagsdataLinearBin,inits=NULL,parameters.to.save = c('beta.pooled','beta','tau'),model.file = modelBinomialLinearDRmeta,
                                          n.chains=2,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 1)
linearDRmetaJAGSmodelBin$BUGSoutput$mean$beta.pooled



