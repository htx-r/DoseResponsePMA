library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force = T)
library(DoseResponseNMA)

OneRunSimulateDRlinear <- function(beta.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){
  sim.data <- simulateDRlineardata.fun(beta.pooled=beta.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize)

  # Freq
  linearDRmetaFreq<-dosresmeta(formula = logRR~dose, id = Study_No,type=type,
                              se = selogRR, cases = cases, n = cases+noncases, data = sim.data, proc='2stage',method = 'reml',covariance = 'gl')

  # Bayes
  jagsdatalinear<- makeJAGSDRmeta(Study_No,logRR,dose,dose2=NULL,cases,noncases,data=sim.data,Splines=F,new.dose.range = c(5,10))
  jagsdatalinear$prec <-  matrix(unlist(sapply(linearDRmetaFreq$Slist,solve,simplify = F)),40,2,byrow = T)

  linearDRmetaJAGSmodel <- jags.parallel(data = jagsdatalinear,inits=NULL,parameters.to.save = c('beta.pooled','tau','newRR'),model.file = modelLinearDRmeta,
                                         n.chains=2,n.iter = 100000,n.burnin = 20000,DIC=F,n.thin = 10)

  f <-coef(linearDRmetaFreq)[1]
  b <- linearDRmetaJAGSmodel$BUGSoutput$mean$beta.pooled
  t <- linearDRmetaJAGSmodel$BUGSoutput$mean$tau
return(cbind(Bayes=b,Freq=f,tau=t))

}


MultiRunSimulateDRlinear <- function(nrep=3,beta.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){

# Repeat the simulation nrep times
res <- replicate(nrep,OneRunSimulateDRlinear(beta.pooled=beta.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize),simplify = F)
res.mat <- do.call(rbind,res)

# Biases
biasB <- beta.pooled-colMeans(res.mat)[1]
biasF <- beta.pooled-colMeans(res.mat)[2]

# standard error
tau.hat <- colMeans(res.mat)[3]
se <- tau.hat/nrep

# Mean square error
mseB <- se^2 + biasB^2
mseF <- se^2 + biasF^2

ret.obj <- cbind(beta.pooled=beta.pooled,tau=tau,beta.pooled.hatB=colMeans(res.mat)[1],beta.pooled.hatF=colMeans(res.mat)[2],tau.hat=tau.hat,biasB=biasB,biasF=biasF,mseB=mseB,mseF=mseF)
row.names(ret.obj) <- NULL
return(ret.obj)
  }


#MultiArgSimulateDRlinear <- function(nrep=3,beta.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200){
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Linear
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Scenario 1
MultiRunSimulateDRlinear(100,beta.pooled=0,tau=0.001)
# beta.pooled   tau beta.pooled.hatB beta.pooled.hatF     tau.hat        biasB         biasF         mseB
# [1,]           0 0.001    -1.929476e-05     4.335371e-05 0.002752483 1.929476e-05 -4.335371e-05 1.129904e-09
# mseF
# [1,] 2.63716e-09

# Scenario 2
MultiRunSimulateDRlinear(100,beta.pooled=0.05,tau=0.001)
# beta.pooled   tau beta.pooled.hatB beta.pooled.hatF     tau.hat        biasB        biasF         mseB
# [1,]        0.05 0.001       0.04995209       0.05003014 0.003565083 4.790848e-05 -3.01432e-05 3.566204e-09
# mseF
# [1,] 2.179594e-09

# Scenario 3
MultiRunSimulateDRlinear(100,beta.pooled=0.02,tau=0.001)

# beta.pooled   tau beta.pooled.hatB beta.pooled.hatF     tau.hat        biasB         biasF         mseB
# [1,]        0.02 0.001       0.01995331        0.0200119 0.003119205 4.668803e-05 -1.189812e-05 3.152716e-09
# mseF
# [1,] 1.114509e-09

# Scenario 4
MultiRunSimulateDRlinear(100,beta.pooled=0.03,tau=0.001)

# beta.pooled   tau beta.pooled.hatB beta.pooled.hatF     tau.hat        biasB        biasF        mseB
# [1,]        0.03 0.001       0.02985735       0.02992357 0.003219418 0.0001426528 7.643393e-05 2.13863e-08
# mseF
# [1,] 6.878611e-09

# Scenario 5
MultiRunSimulateDRlinear(100,beta.pooled=0.2,tau=0.001)

# beta.pooled   tau beta.pooled.hatB beta.pooled.hatF     tau.hat        biasB         biasF         mseB
# [1,]         0.2 0.001        0.1999413        0.2000793 0.007511882 5.873492e-05 -7.934617e-05 9.092628e-09
# mseF
# [1,] 1.193865e-08

# Scenario 6: tau= 0.1 and 0.05 does not work
MultiRunSimulateDRlinear(100,beta.pooled=0,tau=0.01)

# beta.pooled  tau beta.pooled.hatB beta.pooled.hatF     tau.hat         biasB         biasF         mseB
# [1,]           0 0.01     0.0008028447     0.0009025927 0.004211825 -0.0008028447 -0.0009025927 6.463336e-07
# mseF
# [1,] 8.164476e-07

# Scenario 7: tau=0.1 does not work
MultiRunSimulateDRlinear(100,beta.pooled=0.05,tau=0.01)
# beta.pooled  tau beta.pooled.hatB beta.pooled.hatF     tau.hat         biasB         biasF         mseB
# [1,]        0.05 0.01        0.0506469       0.05077669 0.004732937 -0.0006468989 -0.0007766859 4.207182e-07
# mseF
# [1,] 6.054811e-07

# Scenario 8: tau=0.1 and 0.05 do not work
MultiRunSimulateDRlinear(100,beta.pooled=0.02,tau=0.01)

# beta.pooled  tau beta.pooled.hatB beta.pooled.hatF     tau.hat         biasB         biasF         mseB
# [1,]        0.02 0.01       0.02057344       0.02068407 0.004289026 -0.0005734367 -0.0006840749 3.306692e-07
# mseF
# [1,] 4.69798e-07

# Scenario 9:tau = 0.1 does not work
MultiRunSimulateDRlinear(100,beta.pooled=0.03,tau=0.05)

# beta.pooled  tau beta.pooled.hatB beta.pooled.hatF    tau.hat        biasB        biasF         mseB
# [1,]        0.03 0.05       0.03348436       0.03366172 0.03367269 -0.003484356 -0.003661719 1.225412e-05
# mseF
# [1,] 1.352157e-05

# Scenario 10: for tau=0.1 and 0.05 do not work
MultiRunSimulateDRlinear(100,beta.pooled=0.2,tau=0.01)

# beta.pooled  tau beta.pooled.hatB beta.pooled.hatF    tau.hat        biasB       biasF         mseB
# [1,]         0.2 0.01        0.1996736        0.1997921 0.00853372 0.0003263693 0.000207876 1.137994e-07
# mseF
# [1,] 5.049487e-08


# beta.pooled <- rep(c(0.02)) ## 0,0.05
# tau <- rep(c(0.001)) # ,0.01
# ns <- 20
# doserange <- c(0,10)
# samplesize<- 200
# nrep <- 100
# res <- mapply(MultiRunSimulateDRlinear, beta.pooled=beta.pooled,tau=tau,MoreArgs=list(ns=ns,doserange=doserange,samplesize=samplesize,nrep=nrep),SIMPLIFY =F)
# result <- do.call(rbind,res)
# rownames(result) <- paste0('Scenario ',1:nrow(result))
# result

#}
## Run all scenarios
