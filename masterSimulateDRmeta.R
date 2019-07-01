### In this file I run jags and dosresmeta for the simulated data assuming: cubic spline, linear and quadratic.

library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force = T)

library(DoseResponseNMA)

## Cubic Splines
# compare dosresmeta vs Bayes

b1 <- b2<-f1<-f2 <- vector()
n.sim.data <- 100
ns <- 20
start <- Sys.time()
beta1.pooled=0.03
beta2.pooled=0.05
for (j in 1:n.sim.data) {

sim.data <- simulateDRsplinedata.fun(beta1.pooled,beta2.pooled,tau=0.001,doserange = c(1,10))

# Freq
rcsplineDRmetaFreq <- dosresmeta(formula = logRR~dose1+dose2, id = Study_No,type=type,
                                 se = selogRR, cases = cases, n = cases+noncases, data = sim.data, proc='1stage',covariance = 'gl')#!!!!!!!!!!!!!!

# Bayes
jagsdataRCS<- makeJAGSDRmeta(Study_No,logRR,dose1,dose2,cases,noncases,data=sim.data,Splines=T,new.dose.range = c(5,10))
#jagsdataRCS$prec.beta <- 1/(0.001)^2 I tried a fixed precision but it gives similar result as if I do not
jagsdataRCS$prec <-  matrix(unlist(sapply(rcsplineDRmetaFreq$Slist,solve,simplify = F)),40,2,byrow = T)
rcsplineDRmetaJAGSmodel <- jags.parallel(data = jagsdataRCS,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau1','tau2','newRR'),model.file = modelRCSplineDRmeta,
                                n.chains=2,n.iter = 10000,n.burnin = 200,DIC=F,n.thin = 1)

f1[j] <-coef(rcsplineDRmetaFreq)[1]
b1[j] <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta1.pooled

f2[j] <-coef(rcsplineDRmetaFreq)[2]
b2[j] <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta2.pooled

}
end <- Sys.time()
(end-start)
biasF1<- mean(beta1.pooled-f1)
biasF2<-mean( beta2.pooled-f2)
biasB1<- mean(beta1.pooled-b1)
biasB2<- mean(beta2.pooled-b2)
RCSres.dist.tau.100 <- cbind(freq=c(biasF1,biasF2),bayes=c(biasB1,biasB2))


# Linear
# compare dosresmeta vs. Bayes

b1 <-f1 <- vector()
n.sim.data <- 100
beta.pooled = 0.03
start <- Sys.time()
for (j in 1:n.sim.data) {
sim.data <- simulateDRlineardata.fun(beta.pooled)



  linearDRmetaFreq=dosresmeta(formula = logRR~dose, id = Study_No,type=type,
                              se = selogRR, cases = cases, n = cases+noncases, data = sim.data, proc='1stage')
  #jagsdatalinear$prec <-  matrix(unlist(sapply(linearDRmetaFreq$Slist,solve,simplify = F)),40,2,byrow = T)

  # Bayes
  jagsdatalinear<- makeJAGSDRmeta(Study_No,logRR,dose,dose2=NULL,cases,noncases,data=sim.data,Splines=F,new.dose.range = c(5,10))

  linearDRmetaJAGSmodel <- jags.parallel(data = jagsdatalinear,inits=NULL,parameters.to.save = c('beta.pooled','tau','newRR'),model.file = modelLinearDRmeta,
                                           n.chains=2,n.iter = 100000,n.burnin = 20000,DIC=F,n.thin = 10)
  # Freq

  f1[j] <-coef(linearDRmetaFreq)[1]
  b1[j] <- linearDRmetaJAGSmodel$BUGSoutput$mean$beta.pooled


}
end <- Sys.time()
(end-start)
biasb1<- quantile(beta.pooled-b1) #the Bayesian is biased- wrong?
biasf1<- quantile(beta.pooled-f1) #the frequentist is unbiased
Linearres.100 <- cbind(freq=biasf1,bayes=biasb1)


## Quadratic
b1 <- b2<-f1<-f2 <- vector()
n.sim.data <- 100
beta1.pooled=0.01
beta2.pooled=0.02
start <- Sys.time()
for (i in 1:n.sim.data) {

  sim.data <- simulateDRquadraticdata.fun(beta1.pooled,beta2.pooled,tau=0.001,doserange = c(1,10))

  # Bayes
  jagsdataquadratic<- makeJAGSDRmeta(Study_No,logRR,dose,dose2=NULL,cases,noncases,data=sim.data,Splines=F,new.dose.range = c(5,10))



  # Freq
  quadraticDRmetaFreq=dosresmeta(formula = logRR~dose+I(dose^2), id = Study_No,type=type,
                                 se = selogRR, cases = cases, n = cases+noncases, data = sim.data, proc='1stage',method = 'fixed')#!!!!!!!!!!!!!!


  quadraticDRmetaJAGSmodel <- jags.parallel(data = jagsdataquadratic,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau','newRR'),model.file = modelQuadraticDRmeta,
                                           n.chains=2,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 10)


  f1[j] <-coef(quadraticDRmetaFreq)[1]
  b1[j] <- quadraticDRmetaJAGSmodel$BUGSoutput$mean$beta1.pooled

  f2[j] <-coef(quadraticDRmetaFreq)[2]
  b2[j] <- quadraticDRmetaJAGSmodel$BUGSoutput$mean$beta2.pooled


}
end <- Sys.time()

biasF1<- mean(beta1.pooled-f1)
biasF2<-mean( beta2.pooled-f2)
biasB1<- mean(beta1.pooled-b1)
biasB2<- mean(beta2.pooled-b2)
Quadres.dist.tau.100 <- cbind(freq=c(biasF1,biasF2),bayes=c(biasB1,biasB2))
#


























#res.fix.tau.200 <-cbind(freq=c(biasF1,biasF2),bayes=c(biasB1,biasB2))
## I assumed a distribution for tau, run it for 200 simulations and I modified the precision matrix
# TRUE : beta1 = 0.03, beta2 = 0.05
#          distribution for tau
# biasF1  -2.823532e-05           -2.237929e-05
# biasF2  6.61811e-06             -1.206093e-05
# biasB1  0.003293799              0.001281889
# biasB2  0.01202379               0.013902880

## Try another time with fixed tau, it is quicker
## Linear




## Precision
# precmat <- matrix(NA,2*20,2)
# b <- nd <- vector()
#
# for (i in 1:ns) {
#   b[1] <- 0
#   nd[i] <- as.numeric(table(sim.data$Study_No)[i])-1
#   precmat[(b[i]+1):(b[i]+nd[i]),1:(nd[i])] <- diag(sim.data[sim.data$Study_No==i,]$selogRR[2:3])
#   b[i+1] <- b[i]+ nd[i]
#   precmat
# }
# jagsdataquadratic$prec <- precmat

# m <- matrix(unlist(rcsplineDRmetaFreq$Slist),40,2,byrow = T)
# m.inv <- solve(m[1:2,1:2])
#
# jagsdataRCS$prec



# # Data: exclude the 9 studies with NA log RR.
# #antiDep <-  read.csv('~/Desktop/TasnimPhD/DoseResponseNMA1/DOSEmainanalysis.csv')
# antiDep <-  read.csv('C:/Users/th19g983/Desktop/DoseResponseNMA/DOSEmainanalysis.csv')
# #antiDep <-  read.csv("~/Google Drive/_mydrive/HTx/HTx-R/DoseResponseNMA/DOSEmainanalysis.csv")
# NAstudyid <- antiDep$Study_No[is.na(antiDep$logRR)]
# antiDep1 <-antiDep[!antiDep$Study_No %in% NAstudyid,]
