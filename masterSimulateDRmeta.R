### In this file I run jags and dosresmeta for the simulated data assuming: cubic spline, linear and quadratic.

library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force = T)

library(DoseResponseNMA)

## Cubic Splines
# compare dosresmeta vs Bayes


n.sim.data <- 100
ns <- 20
start <- Sys.time()
beta1.pooled=0.03
beta2.pooled=0.05
bayesCoef1RR<-bayesCoef2RR<-freqCoef1RR<-freqCoef2RR<-c()

for (j in 1:n.sim.data) {
  sim.data <- simulateDRsplinedata.fun(beta1.pooled,beta2.pooled,tau=0.001,doserange = c(1,10))

  # Freq
  rcsplineDRmetaFreq <- dosresmeta(formula = logRR~dose1+dose2, id = Study_No,type=type,
                                   se = selogRR, cases = cases, n = cases+noncases, data = sim.data$simulatedDRdata, proc='1stage',covariance = 'gl')#!!!!!!!!!!!!!!
  #rcsplineDRmetaFreq2 <- dosresmeta(formula = logRR~rcs(dose1,sim.data$knots), id = Study_No,type=type,
   #                                 se = selogRR, cases = cases, n = cases+noncases, data = sim.data$simulatedDRdata, proc='1stage',covariance = 'gl')#!!!!!!!!!!!!!!
  #coef2<-cbind(coef2,c(coef(rcsplineDRmetaFreq2)[1],coef(rcsplineDRmetaFreq2)[2]))

  # Bayes normal
  jagsdataRCS<- makeJAGSDRmeta(Study_No,logRR,dose1,dose2,cases,noncases,data=sim.data$simulatedDRdata,Splines=T,new.dose.range = c(1,10))
  ## ----!!!!!---- PLEASE PUT THE COMMAND OF THE PRECISION IN THE makeJAGSDRmeta FUNCTION. And as we agreed you need to find out why you calculated it wrong-------------
  jagsdataRCS$prec <-  matrix(unlist(sapply(rcsplineDRmetaFreq$Slist,solve,simplify = F)),40,2,byrow = T)
  rcsplineDRmetaJAGSmodel <- jags.parallel(data = jagsdataRCS,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau1','tau2','newRR'),model.file = modelRCSplineDRmeta,
                                           n.chains=2,n.iter = 10000,n.burnin = 200,DIC=F,n.thin = 1)


  #Bayes binomial
  #jagsdataSplineBin <- makeBinomialJAGSDRmeta(studyid=Study_No,dose1 = dose1,dose2=dose2,cases=cases,noncases=noncases,data=sim.data$simulatedDRdata,Splines=T)

  #splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdataSplineBin,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelBinomialRCSsplineDRmeta,
   #                                         n.chains=1,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 5)

  bayesCoef1RR <- c(bayesCoef1RR,rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta1.pooled)
  bayesCoef2RR <- c(bayesCoef2RR,rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta2.pooled)

  freqCoef1RR<-c(freqCoef1RR,coef(rcsplineDRmetaFreq)[1])
  freqCoef2RR<-c(freqCoef2RR,coef(rcsplineDRmetaFreq)[2])
}

(mean(bayesCoef1RR)-beta1.pooled) ## 20/07 bias = -0.0002586919
(mean(bayesCoef2RR)-beta2.pooled) ## 20/07 bias = 0.0004295099


(mean(freqCoef1RR)-beta1.pooled) # 20/07 bias = -7.182011e-06
(mean(freqCoef2RR)-beta2.pooled) # 20/07 bias = 2.164394e-05


cbind(bayes1=quantile(bayesCoef1RR), freq1=quantile(freqCoef1RR))

# bayes1      freq1
# 0%   0.02606195 0.02912352
# 25%  0.02890433 0.02978318
# 50%  0.02976993 0.02999954
# 75%  0.03061950 0.03019181
# 100% 0.03360143 0.03073798

cbind(bayes2=quantile(bayesCoef2RR), freq2=quantile(freqCoef2RR))

# bayes2      freq2
# 0%   0.04200586 0.04818960
# 25%  0.04876260 0.04947566
# 50%  0.05053595 0.05009610
# 75%  0.05221470 0.05057717
# 100% 0.05873851 0.05233345






# Linear
# compare dosresmeta vs. Bayes

freqCoefRR <-bayesCoefRR <- vector()
n.sim.data <- 100
beta.pooled = 0.02
tau <- 0.01
start <- Sys.time()
for (j in 1:n.sim.data) {
  sim.data <- simulateDRlineardata.fun(beta.pooled,tau)

  # Freq
  linearDRmetaFreq=dosresmeta(formula = logRR~dose, id = Study_No,type=type,
                              se = selogRR, cases = cases, n = cases+noncases, data = sim.data, proc='1stage')

  # Bayes
  jagsdatalinear<- makeJAGSDRmeta(Study_No,logRR,dose,dose2=NULL,cases,noncases,data=sim.data,Splines=F,new.dose.range = c(5,10))

  #####!!!!!!!! I THOUGHT YOU SAID YOU HAVE TO USE THE CORRECT MATRIX FROM THE FREQ MODEL!!!################

  jagsdatalinear$prec <-  matrix(unlist(sapply(linearDRmetaFreq$Slist,solve,simplify = F)),40,2,byrow = T)
  linearDRmetaJAGSmodel <- jags.parallel(data = jagsdatalinear,inits=NULL,parameters.to.save = c('beta.pooled','tau','newRR'),model.file = modelLinearDRmeta,
                                         n.chains=2,n.iter = 100000,n.burnin = 20000,DIC=F,n.thin = 10)

  freqCoefRR<-c(freqCoefRR,coef(linearDRmetaFreq))
  bayesCoefRR <- c(bayesCoefRR,linearDRmetaJAGSmodel$BUGSoutput$mean$beta.pooled)


}

mean(bayesCoefRR)-beta.pooled ## 20/07: bias=0.000409433, true=0.02
mean(freqCoefRR)-beta.pooled # 20/07: bias = 0.0005221882, true =0.02
cbind(bayes=quantile(bayesCoefRR), freq=quantile(freqCoefRR))

# bayes       freq
# 0%   0.01331337 0.01354443
# 25%  0.01801244 0.01808566
# 50%  0.02033683 0.02046003
# 75%  0.02243847 0.02251887
# 100% 0.03094213 0.03130403





## Quadratic
b1 <- b2<-f1<-f2 <- vector()
n.sim.data <- 100
beta1.pooled=0.01
beta2.pooled=0.01
start <- Sys.time()
for (i in 1:n.sim.data) {

  sim.data <- simulateDRquadraticdata.fun(beta1.pooled,beta2.pooled,tau=0.001,doserange = c(1,10))

  # Freq
  quadraticDRmetaFreq=dosresmeta(formula = logRR~dose+I(dose^2), id = Study_No,type=type,
                                 se = selogRR, cases = cases, n = cases+noncases, data = sim.data, proc='1stage',method = 'fixed')#!!!!!!!!!!!!!!

  # Bayes normal
  jagsdataquadratic<- makeJAGSDRmeta(Study_No,logRR,dose,dose2=NULL,cases,noncases,data=sim.data,Splines=F,new.dose.range = c(5,10))

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































# biasF1<- (c(beta1.pooled, beta2.pooled)-coef1)
# biasF2<- (c(beta1.pooled, beta2.pooled)-coef2)
# biasBNorm<- (c(beta1.pooled, beta2.pooled)-coef3)
# biasBBin<- (c(beta1.pooled, beta2.pooled)-coef4)
# apply(coef1-coef2,1,mean)#model F2 with direct splines within dosres gives *exactly* the same results as the F1 that takes the two transformations
# ## Tasnim: Yes thats true! But what I meant that rcs withing makejags function does not give the same result as F.
# BIAS=rbind(Freq=apply(biasF1,1,mean),
#            BNorm=apply(biasBNorm,1,mean),
#            BBinom=apply(biasBBin,1,mean))
# relativeBIAS<-round(BIAS/c(beta1.pooled, beta2.pooled)*100,2)
# # the bayesian binomial model has the largest bias, which is weird... we need to investigate this better


# the results from 500 simulations are
#> relativeBIAS
#         dose1 dose2
#Freq   -0.06 -0.07
#BNorm  -0.01 -0.08
#BBinom -1.89  3.99
#> BIAS
#         dose1         dose2
#Freq   -1.684855e-05 -3.526103e-05
#BNorm  -5.643893e-06 -2.285033e-05
#BBinom -5.673424e-04  1.996080e-03









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
