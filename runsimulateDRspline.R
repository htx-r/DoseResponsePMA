###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Restricted Cubic Spline
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# There are two functions:
# OneRunSimulateDRspline: run the simulation for one set of arguments and for one simulation
# MultiRunSimulateDRspline: run the simulation for one set of arguments and for 'nrep' simulations

# And to run the simulation for different arguments, I apply the MultiRunSimulateDRspline() for the different scenarios
     # using mapply (we can make it as a function if Georgia would like)
library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force = T)
library(DoseResponseNMA)

OneRunSimulateDRspline <- function(beta1.pooled=0.03,beta2.pooled=0.05,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){
  sim.data <- simulateDRsplinedata.fun(beta1.pooled=beta1.pooled,beta2.pooled=beta2.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize)

  # Freq
  rcsplineDRmetaFreq <- dosresmeta(formula = logRR~dose1+dose2, id = Study_No,type=type,
                                   se = selogRR, cases = cases, n = cases+noncases, data = sim.data$simulatedDRdata, proc='1stage',covariance = 'gl')#!!!!!!!!!!!!!!

  # Bayes
  jagsdataRCS<- makeJAGSDRmeta(Study_No,logRR,dose1,dose2,cases,noncases,data=sim.data$simulatedDRdata,Splines=T,new.dose.range = c(5,10))
  #jagsdataRCS$prec.beta <- 1/(0.001)^2 I tried a fixed precision but it gives similar result
  jagsdataRCS$prec <-  matrix(unlist(sapply(rcsplineDRmetaFreq$Slist,solve,simplify = F)),40,2,byrow = T)
  rcsplineDRmetaJAGSmodel <- jags.parallel(data = jagsdataRCS,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau1','tau2','newRR'),model.file = modelRCSplineDRmeta,
                                           n.chains=2,n.iter = 10000,n.burnin = 200,DIC=F,n.thin = 1)

  f1 <-coef(rcsplineDRmetaFreq)[1]
  b1 <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta1.pooled

  f2 <-coef(rcsplineDRmetaFreq)[2]
  b2 <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta2.pooled

  t1 <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$tau1
  t2 <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$tau2

  return(cbind(Bayes=c(b1,b2),Freq=c(f1,f2),tau.hat=c(t1,t2)))

}



MultiRunSimulateDRspline <- function(nrep=3,beta1.pooled=0.02,beta2.pooled=0.05,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){

  # Repeat the simulation nrep times
  res <- replicate(nrep,OneRunSimulateDRspline(beta1.pooled=beta1.pooled,beta2.pooled = beta2.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize),simplify = F)
  res.mat1 <- t(sapply(1:nrep, function(i) res[[i]][1,]))
  res.mat2 <- t(sapply(1:nrep, function(i) res[[i]][2,]))


  # Biases
  biasB1 <- beta1.pooled-colMeans(res.mat1)[1]
  biasF1 <- beta1.pooled-colMeans(res.mat1)[2]

  biasB2 <- beta2.pooled-colMeans(res.mat2)[1]
  biasF2 <- beta2.pooled-colMeans(res.mat2)[2]


  # standard error
  tau1.hat <- colMeans(res.mat1)[3]
  se1 <- tau1.hat/nrep

  tau2.hat <- colMeans(res.mat2)[3]
  se2 <- tau2.hat/nrep
  # Mean square error
  mseB1 <- se1^2 + biasB1^2
  mseF1 <- se1^2 + biasF1^2

  mseB2 <- se2^2 + biasB2^2
  mseF2 <- se2^2 + biasF2^2

  ret.obj <- list(sum.coef1=cbind(beta1.pooled=beta1.pooled,beta2.pooled=beta2.pooled,tau=tau,beta1.pooled.hatB=colMeans(res.mat1)[1],beta1.pooled.hatF=colMeans(res.mat1)[2],tau1.hat=tau1.hat,biasB1=biasB1,biasF1=biasF1,mseB1=mseB1,mseF1=mseF1),
                  sum.coef2=cbind(beta1.pooled=beta1.pooled,beta2.pooled=beta2.pooled,tau=tau,beta2.pooled.hatB=colMeans(res.mat2)[1],beta2.pooled.hatF=colMeans(res.mat2)[2],tau2.hat=tau2.hat,biasB2=biasB2,biasF2=biasF2,mseB2=mseB2,mseF2=mseF2))
  row.names(ret.obj[[1]]) <-row.names(ret.obj[[2]]) <- NULL
  return(ret.obj)
}
beta1.pooled <- c(0,0.05,0.02,0.03,0.2,0.03)
beta2.pooled <- c(0,0,0.02,0.05,0.03,0.2)
tau <- c(0.001,0.05)

# Scenario 1
MultiRunSimulateDRspline(nrep=100,beta1.pooled = beta1.pooled[1],beta2.pooled = beta2.pooled[1],tau=tau[1])

# $sum.coef1
# beta1.pooled beta2.pooled   tau beta1.pooled.hatB beta1.pooled.hatF    tau1.hat       biasB1        biasF1
# [1,]            0            0 0.001     -0.0002493762      3.204149e-05 0.003038812 0.0002493762 -3.204149e-05
# mseB1        mseF1
# [1,] 6.311191e-08 1.950095e-09
#
# $sum.coef2
# beta1.pooled beta2.pooled   tau beta2.pooled.hatB beta2.pooled.hatF    tau2.hat        biasB2       biasF2
# [1,]            0            0 0.001      0.0001307393     -0.0001365649 0.007214545 -0.0001307393 0.0001365649
# mseB2        mseF2
# [1,] 2.229772e-08 2.385494e-08

# Scenario 2
MultiRunSimulateDRspline(nrep=100,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[2],tau=tau[1])

# $sum.coef1
# beta1.pooled beta2.pooled   tau beta1.pooled.hatB beta1.pooled.hatF    tau1.hat       biasB1       biasF1
# [1,]         0.05            0 0.001        0.04956291        0.04999326 0.003903633 0.0004370915 6.737608e-06
# mseB1       mseF1
# [1,] 1.925728e-07 1.56923e-09
#
# $sum.coef2
# beta1.pooled beta2.pooled   tau beta2.pooled.hatB beta2.pooled.hatF    tau2.hat        biasB2        biasF2
# [1,]         0.05            0 0.001      0.0006700803      5.909034e-05 0.008654569 -0.0006700803 -5.909034e-05
# mseB2        mseF2
# [1,] 4.564978e-07 1.098182e-08

# Scenario 3
MultiRunSimulateDRspline(nrep=100,beta1.pooled = beta1.pooled[3],beta2.pooled = beta2.pooled[3],tau=tau[1])

# $sum.coef1
# beta1.pooled beta2.pooled   tau beta1.pooled.hatB beta1.pooled.hatF    tau1.hat       biasB1       biasF1
# [1,]         0.02         0.02 0.001        0.01951247        0.01995735 0.003696527 0.0004875335 4.264764e-05
# mseB1        mseF1
# [1,] 2.390553e-07 3.185252e-09
#
# $sum.coef2
# beta1.pooled beta2.pooled   tau beta2.pooled.hatB beta2.pooled.hatF    tau2.hat        biasB2       biasF2
# [1,]         0.02         0.02 0.001        0.02052635        0.01988887 0.008198946 -0.0005263478 0.0001111275
# mseB2        mseF2
# [1,] 2.837643e-07 1.907159e-08
# Scenario 4
MultiRunSimulateDRspline(nrep=100,beta1.pooled = beta1.pooled[4],beta2.pooled = beta2.pooled[4],tau=tau[1])

# $sum.coef1
# beta1.pooled beta2.pooled   tau beta1.pooled.hatB beta1.pooled.hatF    tau1.hat       biasB1       biasF1
# [1,]         0.03         0.05 0.001        0.02991095        0.02999295 0.004383463 8.905421e-05 7.046931e-06
# mseB1        mseF1
# [1,] 9.852128e-09 1.971134e-09
#
# $sum.coef2
# beta1.pooled beta2.pooled   tau beta2.pooled.hatB beta2.pooled.hatF    tau2.hat        biasB2       biasF2
# [1,]         0.03         0.05 0.001        0.05000051        0.04995082 0.008958395 -5.133388e-07 4.917665e-05
# mseB2        mseF2
# [1,] 8.025548e-09 1.044363e-08

# Scenario 5
MultiRunSimulateDRspline(nrep=100,beta1.pooled = beta1.pooled[5],beta2.pooled = beta2.pooled[5],tau=tau[1])

# $sum.coef1
# beta1.pooled beta2.pooled   tau beta1.pooled.hatB beta1.pooled.hatF    tau1.hat       biasB1       biasF1
# [1,]          0.2         0.03 0.001         0.1995374         0.1999708 0.009117236 0.0004625822 2.919596e-05
# mseB1        mseF1
# [1,] 2.222947e-07 9.164803e-09
#
# $sum.coef2
# beta1.pooled beta2.pooled   tau beta2.pooled.hatB beta2.pooled.hatF   tau2.hat        biasB2       biasF2
# [1,]          0.2         0.03 0.001         0.0303164        0.02999838 0.01732281 -0.0003164031 1.617229e-06
# mseB2        mseF2
# [1,] 1.301189e-07 3.001057e-08

# Scenario 6
MultiRunSimulateDRspline(nrep=100,beta1.pooled = beta1.pooled[6],beta2.pooled = beta2.pooled[6],tau=tau[1])

# $sum.coef1
# beta1.pooled beta2.pooled   tau beta1.pooled.hatB beta1.pooled.hatF    tau1.hat       biasB1       biasF1
# [1,]         0.03          0.2 0.001        0.02955996        0.02997596 0.007166366 0.0004400372 2.404132e-05
# mseB1        mseF1
# [1,] 1.987684e-07 5.713666e-09
#
# $sum.coef2
# beta1.pooled beta2.pooled   tau beta2.pooled.hatB beta2.pooled.hatF  tau2.hat        biasB2        biasF2
# [1,]         0.03          0.2 0.001         0.2006544         0.2001472 0.0144096 -0.0006543923 -0.0001472431
# mseB2       mseF2
# [1,] 4.489929e-07 4.24442e-08

# Scenario 7:
MultiRunSimulateDRspline(nrep=100,beta1.pooled = beta1.pooled[1],beta2.pooled = beta2.pooled[1],tau=tau[2])

# $sum.coef1
# beta1.pooled beta2.pooled  tau beta1.pooled.hatB beta1.pooled.hatF   tau1.hat       biasB1       biasF1
# [1,]            0            0 0.05       0.005059108       0.005414216 0.02872869 -0.005059108 -0.005414216
# mseB1        mseF1
# [1,] 2.567711e-05 2.939626e-05
#
# $sum.coef2
# beta1.pooled beta2.pooled  tau beta2.pooled.hatB beta2.pooled.hatF   tau2.hat       biasB2       biasF2
# [1,]            0            0 0.05       0.004047606       0.004526225 0.05908485 -0.004047606 -0.004526225
# mseB2        mseF2
# [1,] 1.673221e-05 2.083581e-05

# Scenario 8
MultiRunSimulateDRspline(nrep=100,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[2],tau=tau[2])

# $sum.coef1
# beta1.pooled beta2.pooled  tau beta1.pooled.hatB beta1.pooled.hatF   tau1.hat       biasB1       biasF1
# [1,]         0.05            0 0.05        0.05434472        0.05564241 0.02771095 -0.004344716 -0.005642407
# mseB1        mseF1
# [1,] 1.895335e-05 3.191355e-05
#
# $sum.coef2
# beta1.pooled beta2.pooled  tau beta2.pooled.hatB beta2.pooled.hatF   tau2.hat       biasB2       biasF2
# [1,]         0.05            0 0.05       0.006353568       0.005674696 0.05805583 -0.006353568 -0.005674696
# mseB2        mseF2
# [1,] 4.070487e-05 3.253922e-05
# Scenario 9
MultiRunSimulateDRspline(nrep=100,beta1.pooled = beta1.pooled[3],beta2.pooled = beta2.pooled[3],tau=tau[2])

# $sum.coef1
# beta1.pooled beta2.pooled  tau beta1.pooled.hatB beta1.pooled.hatF  tau1.hat        biasB1        biasF1
# [1,]         0.02         0.02 0.05        0.02049488        0.02097067 0.0275847 -0.0004948843 -0.0009706729
# mseB1        mseF1
# [1,] 3.21002e-07 1.018297e-06
#
# $sum.coef2
# beta1.pooled beta2.pooled  tau beta2.pooled.hatB beta2.pooled.hatF   tau2.hat       biasB2       biasF2
# [1,]         0.02         0.02 0.05        0.02385763        0.02436841 0.05537632 -0.003857629 -0.004368412
# mseB2        mseF2
# [1,] 1.518795e-05 1.938968e-05

#Scenario 10
MultiRunSimulateDRspline(nrep=100,beta1.pooled = beta1.pooled[4],beta2.pooled = beta2.pooled[4],tau=tau[2])

# $sum.coef1
# beta1.pooled beta2.pooled  tau beta1.pooled.hatB beta1.pooled.hatF   tau1.hat        biasB1       biasF1
# [1,]         0.03         0.05 0.05        0.03063292          0.031895 0.02491122 -0.0006329198 -0.001895003
# mseB1        mseF1
# [1,] 4.626444e-07 3.653093e-06
#
# $sum.coef2
# beta1.pooled beta2.pooled  tau beta2.pooled.hatB beta2.pooled.hatF   tau2.hat       biasB2       biasF2
# [1,]         0.03         0.05 0.05        0.05482805        0.05407731 0.05121914 -0.004828049 -0.004077311
# mseB2       mseF2
# [1,] 2.35724e-05 1.68868e-05

# Scenario 11
MultiRunSimulateDRspline(nrep=100,beta1.pooled = beta1.pooled[5],beta2.pooled = beta2.pooled[5],tau=tau[2]-0.04)

# $sum.coef1
# beta1.pooled beta2.pooled  tau beta1.pooled.hatB beta1.pooled.hatF   tau1.hat        biasB1        biasF1
# [1,]          0.2         0.03 0.01          0.200102         0.2004506 0.01084626 -0.0001020014 -0.0004505848
# mseB1        mseF1
# [1,] 2.216842e-08 2.147908e-07
#
# $sum.coef2
# beta1.pooled beta2.pooled  tau beta2.pooled.hatB beta2.pooled.hatF   tau2.hat      biasB2       biasF2
# [1,]          0.2         0.03 0.01        0.02986585         0.0297256 0.01952418 0.000134147 0.0002744017
# mseB2        mseF2
# [1,] 5.61148e-08 1.134157e-07

# Scenario 12
MultiRunSimulateDRspline(nrep=100,beta1.pooled = beta1.pooled[6],beta2.pooled = beta2.pooled[6],tau=tau[2]-0.04)

# $sum.coef1
# beta1.pooled beta2.pooled  tau beta1.pooled.hatB beta1.pooled.hatF    tau1.hat      biasB1       biasF1
# [1,]         0.03          0.2 0.01        0.02893096        0.02982251 0.008884698 0.001069036 0.0001774922
# mseB1        mseF1
# [1,] 1.150731e-06 3.939728e-08
#
# $sum.coef2
# beta1.pooled beta2.pooled  tau beta2.pooled.hatB beta2.pooled.hatF   tau2.hat       biasB2        biasF2
# [1,]         0.03          0.2 0.01         0.2021122          0.200677 0.01849595 -0.002112178 -0.0006769527
# mseB2       mseF2
# [1,] 4.495506e-06 4.92475e-07


# Arguments
# beta1.pooled <- rep(c(0,0.02,0.03,0.05,0.03),2)
# beta2.pooled <- rep(c(0,0,0.02,0.03,0.05),2)
# tau <- rep(c(0.001,0.01),each=6)
# ns <- 20
# doserange <- c(0,10)
# samplesize<- 200
# nrep <- 100
#
# # Results
# #MultiArgSimulateDRlinear <- function(nrep=3,beta.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){
# res <- mapply(MultiRunSimulateDRspline, beta1.pooled=beta1.pooled,beta2.pooled=beta2.pooled,tau=tau,MoreArgs=list(ns=ns,doserange=doserange,samplesize=samplesize,nrep=nrep),SIMPLIFY =F)
# res.mat1 <- t(sapply(1:length(beta1.pooled), function(i) res[[i]]$sum.coef1[1,]))
# res.mat2 <- t(sapply(1:length(beta1.pooled), function(i) res[[i]]$sum.coef2[1,]))
# rownames(res.mat1) <- paste0('Scenario ',1:nrow(res.mat1))
# rownames(res.mat2) <- paste0('Scenario ',1:nrow(res.mat2))
# res.mat1
# res.mat2
#return(list(sum.coef1=res.mat1,sum.coef2=res.mat2))
#}
