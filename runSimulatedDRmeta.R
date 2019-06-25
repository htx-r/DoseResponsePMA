### In this file I run jags and dosresmeta for the simulated data assuming: cubic spline, linear and quadratic.

library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force = T)

library(DoseResponseNMA)

## Cubic Splines
# compare dosresmeta vs Bayes

res1 <- res2 <- list()
b1 <- b2<-f1<-f2 <- vector()
n.sim.data <- 5
ns <- 20
start <- Sys.time()

for (j in 1:n.sim.data) {

sim <- simulateDRsplinedata.fun(beta1.pooled=0.01,beta2.pooled=0.01,tau=0.001,doserange = c(1,10))
sim.data <- sim$simulatedDRdata
knots <- sim$knots
# Bayes
jagsdataRCS<- makeJAGSDRmeta(Study_No,logRR,dose,cases,noncases,data=sim.data,Splines=T,knots=knots)
## Precision
precmat <- matrix(NA,2*ns,2)
b <- nd <- vector()

for (i in 1:ns) {
  b[1] <- 0
  nd[i] <- as.numeric(table(sim.data$Study_No)[i])-1
  precmat[(b[i]+1):(b[i]+nd[i]),1:(nd[i])] <- diag(sim.data[sim.data$Study_No==i,]$selogRR[2:3])
  b[i+1] <- b[i]+ nd[i]
  precmat
}
jagsdataRCS$prec <- precmat
# jagsdataRCS$X <- jagsdataRCS$X1
# jagsdataRCS$Xref <- jagsdataRCS$X1ref


rcsplineDRmetaJAGSmodel <- jags.parallel(data = jagsdataRCS,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau1','tau2'),model.file = modelRCSplineDRmeta,
                                n.chains=2,n.iter = 100000,n.burnin = 2000,DIC=F,n.thin = 10)
# Freq
rcsplineDRmetaFreq <- dosresmeta(formula = logRR~rcs(sim.data$dose,knots), id = Study_No,type=type,
                              se = selogRR, cases = cases, n = cases+noncases, data = sim.data, proc='1stage')#!!!!!!!!!!!!!!

f1[j] <-coef(rcsplineDRmetaFreq)[1]
b1[j] <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta1.pooled

f2[j] <-coef(rcsplineDRmetaFreq)[2]
b2[j] <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta2.pooled

res1[[j]] <- c(b1[j], f1[j])
res2[[j]] <- c(b2[j], f2[j])

}
end <- Sys.time()
(end-start)
resultDRsplinemeta1 <- apply(matrix(unlist(res1),n.sim.data,2,byrow = T),2,mean)
resultDRsplinemeta2 <- apply(matrix(unlist(res2),n.sim.data,2,byrow = T),2,mean)
bias1<- 0.01-resultDRsplinemeta1
bias2<- 0.01-resultDRsplinemeta2

#  0.004700924 0.006625585: 1000 simulations, 0.01
#

## Linear

# compare dosresmeta vs Bayes

res1 <- list()
b1 <-f1 <- vector()
n.sim.data <- 2
start <- Sys.time()
for (j in 1:n.sim.data) {


  sim.data <- simulateDRlineardata.fun()

  # Bayes
  jagsdatalinear<- makeJAGSDRmeta(Study_No,logRR,dose,cases,noncases,data=sim.data,Splines=F,knots=knots)

  ## Precision
  precmat <- matrix(NA,2*20,2)
  b <- nd <- vector()

  for (i in 1:ns) {
    b[1] <- 0
    nd[i] <- as.numeric(table(sim.data$Study_No)[i])-1
    precmat[(b[i]+1):(b[i]+nd[i]),1:(nd[i])] <- diag(sim.data[sim.data$Study_No==i,]$selogRR[2:3])
    b[i+1] <- b[i]+ nd[i]
    precmat
  }
  jagsdatalinear$prec <- precmat


  linearDRmetaJAGSmodel <- jags.parallel(data = jagsdatalinear,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelLinearDRmeta,
                                           n.chains=2,n.iter = 100000,n.burnin = 20000,DIC=F,n.thin = 10)
  # Freq
  linearDRmetaFreq=dosresmeta(formula = logRR~dose, id = Study_No,type=type,
                                se = selogRR, cases = cases, n = cases+noncases, data = sim.data, proc='1stage')#!!!!!!!!!!!!!!

  f1[j] <-coef(linearDRmetaFreq)[1]
  b1[j] <- linearDRmetaJAGSmodel$BUGSoutput$mean$beta.pooled

  res1[[j]] <- c(b1[j], f1[j])

}
end <- Sys.time()
(end-start)  # it
resultDRlinearmeta <- apply(matrix(unlist(res1),n.sim.data,2,byrow = T),2,mean) ##  -0.4206411  0.9744341
bias<- 0.02-resultDRlinearmeta
#  0.01023117 0.01435658: 50 simulations, 0.01 (true)

## Quadratic
res1 <- res2 <- list()
b1 <- b2<-f1<-f2 <- vector()
n.sim.data <- 10
start <- Sys.time()
for (i in 1:n.sim.data) {

  sim.data <- simulateDRquadraticdata.fun(beta1.pooled=0.01,beta2.pooled=0.02,tau=0.001,doserange = c(1,10))

  # Bayes
  jagsdataquadratic<- makeJAGSDRmeta(Study_No,logRR,dose,cases,noncases,data=sim.data,Splines=F,knots=knots)

  ## Precision
  precmat <- matrix(NA,2*20,2)
  b <- nd <- vector()

  for (i in 1:ns) {
    b[1] <- 0
    nd[i] <- as.numeric(table(sim.data$Study_No)[i])-1
    precmat[(b[i]+1):(b[i]+nd[i]),1:(nd[i])] <- diag(sim.data[sim.data$Study_No==i,]$selogRR[2:3])
    b[i+1] <- b[i]+ nd[i]
    precmat
  }
  jagsdataquadratic$prec <- precmat


  quadraticDRmetaJAGSmodel <- jags.parallel(data = jagsdataquadratic,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelQuadraticDRmeta,
                                           n.chains=2,n.iter = 100000,n.burnin = 20000,DIC=F,n.thin = 10)
  # Freq
quadraticDRmetaFreq=dosresmeta(formula = logRR~dose+I(dose^2), id = Study_No,type=type,
                                se = selogRR, cases = cases, n = cases+noncases, data = sim.data, proc='1stage',method = 'fixed')#!!!!!!!!!!!!!!
lm(logRR~dose+I(dose^2)-1,data=sim.data)

f1[j] <-coef(quadraticDRmetaFreq)[1]
  b1[j] <- quadraticDRmetaJAGSmodel$BUGSoutput$mean$beta1.pooled

  f2[j] <-coef(quadraticDRmetaFreq)[2]
  b2[j] <- quadraticDRmetaJAGSmodel$BUGSoutput$mean$beta2.pooled

  res1[[j]] <- c(b1[j], f1[j])
  res2[[j]] <- c(b2[j], f2[j])

}
end <- Sys.time()
(end-start)*10/60  # it
resultDRquadraticmeta1 <- apply(matrix(unlist(res1),n.sim.data,2,byrow = T),2,mean) ##  -0.4206411  0.9744341
resultDRquadraticmeta2 <- apply(matrix(unlist(res2),n.sim.data,2,byrow = T),2,mean) ##  -0.4206411  0.9744341

#



























# # Data: exclude the 9 studies with NA log RR.
# #antiDep <-  read.csv('~/Desktop/TasnimPhD/DoseResponseNMA1/DOSEmainanalysis.csv')
# antiDep <-  read.csv('C:/Users/th19g983/Desktop/DoseResponseNMA/DOSEmainanalysis.csv')
# #antiDep <-  read.csv("~/Google Drive/_mydrive/HTx/HTx-R/DoseResponseNMA/DOSEmainanalysis.csv")
# NAstudyid <- antiDep$Study_No[is.na(antiDep$logRR)]
# antiDep1 <-antiDep[!antiDep$Study_No %in% NAstudyid,]
