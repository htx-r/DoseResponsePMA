source('FunctionsForSimulations.R')
library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force = T)
library(DoseResponseNMA)
library('rsimsum')
library(tidyr)


# I did the analysis in these settings
nsim <- 100
beta.pooled <- 0.1
tau <- c(0.001, 0.03, 0.05)
alpha1 <- tau-tau/4
alpha2 <- tau+tau/4
alpha1
alpha2

## 1. Unif~(truevalue-truevalue/5, truevalue+truevalue/5)
modelBinLinearDRmetaOR <- function(){

  for (i in 1:ns) { ## for each study
    r[i,1] ~ dbinom(p[i,1],n[i,1])
    logit(p[i,1])<- u[i]
    for (j in 2:(nd[i])) { ## for each dose
      # Likelihood
      r[i,j] ~ dbinom(p[i,j],n[i,j])
      logit(p[i,j])<- u[i] + delta[i,j]
      delta[i,j] <-  beta[i]*(X[i,(j)]-X[i,1])

    }
  }
  # Random effect

  for(i in 1:ns) {
    beta[i] ~dnorm(beta.pooled,prec.tau)
    u[i]~dnorm(0,0.1)
  }

  # Priors
  prec.tau<-1/variance
  variance<-tau*tau
  tau~dunif(0.03750,0.06250)#(0.03750,0.06250)#(0.02250,0.03750)#(0.00075,0.00125)
  beta.pooled ~ dnorm(0,16)
}



# Scenario 1
S1ORlineartau <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[1],OR=TRUE,splines=FALSE)

# informative tau prior and informative beta prior
# true.beta  BayesBbias   BayesNbias    Freqbias BayesBcover BayesNcover Freqcover   BayesBmse
# 1       0.1 -0.06642342 0.0003934013 0.001262958        0.06           1      0.95 0.004603701
# BayesNmse      Freqmse BayesBpower BayesNpower Freqpower BayesBse2mean BayesNse2mean  Freqse2mean
# 1 0.0001454643 7.529782e-05        0.55           0         1   0.000306583    0.00550949 7.130285e-05
# mcse.biasB  mcse.biasN   mcse.biasF    RhatB    RhatN true.tau    tau.hatB   tau.hatN    tau.hatF
# 1 0.001391281 0.001211517 0.0008628281 2.132221 1.001985    0.001 0.001000626 0.06497772 0.009872852


# informative tau prior and Non-informative beta prior
# true.beta  BayesBbias  BayesNbias     Freqbias BayesBcover BayesNcover Freqcover   BayesBmse
# 1       0.1 -0.06526795 0.001652694 0.0006945844        0.04           1         1 0.004438145
# BayesNmse      Freqmse BayesBpower BayesNpower Freqpower BayesBse2mean BayesNse2mean  Freqse2mean
# 1 0.0001160161 4.754157e-05        0.64           0         1  0.0002898252   0.005515213 6.686179e-05
# mcse.biasB  mcse.biasN   mcse.biasF    RhatB    RhatN true.tau    tau.hatB   tau.hatN    tau.hatF
# 1 0.001341786 0.001069715 0.0006894524 2.043359 1.001677    0.001 0.001001033 0.06490197 0.007337904

# Scenario 2
S2ORlineartau <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[2],OR=TRUE,splines=FALSE)

# informative tau prior and Non-informative beta prior
# true.beta  BayesBbias   BayesNbias    Freqbias BayesBcover BayesNcover Freqcover    BayesBmse
# 1       0.1 -0.00271466 0.0002297896 0.001670818        0.96           1      0.92 0.0005458064
# BayesNmse    Freqmse BayesBpower BayesNpower Freqpower BayesBse2mean BayesNse2mean  Freqse2mean
# 1 0.0001730421 0.00011884        0.99           0         1  0.0005503226   0.005474095 0.0001039017
# mcse.biasB mcse.biasN  mcse.biasF   RhatB    RhatN true.tau   tau.hatB   tau.hatN   tau.hatF
# 1 0.002332114 0.00132188 0.001082685 1.00676 1.001737     0.03 0.03000682 0.06506521 0.02563138

# informative tau prior and informative beta prior
# true.beta   BayesBbias  BayesNbias    Freqbias BayesBcover BayesNcover Freqcover    BayesBmse
# 1       0.1 0.0003824406 0.001861967 0.003280239        0.96           1      0.95 0.0004807285
# BayesNmse      Freqmse BayesBpower BayesNpower Freqpower BayesBse2mean BayesNse2mean  Freqse2mean
# 1 0.0001634588 0.0001031294        0.98           0         1  0.0005548605   0.005505381 0.0001015051
# mcse.biasB  mcse.biasN   mcse.biasF    RhatB    RhatN true.tau   tau.hatB   tau.hatN   tau.hatF
# 1 0.002203263 0.001271251 0.0009659319 1.005156 1.001892     0.03 0.03001913 0.06516375 0.02441259

# Scenario 3

S3ORlineartau <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[3],OR=TRUE,splines=FALSE)

# informative tau prior and Non-informative beta prior
# true.beta  BayesBbias  BayesNbias    Freqbias BayesBcover BayesNcover Freqcover    BayesBmse
# 1       0.1 0.006169978 0.001214205 0.003498081        0.96           1      0.98 0.0006738484
# BayesNmse      Freqmse BayesBpower BayesNpower Freqpower BayesBse2mean BayesNse2mean  Freqse2mean
# 1 0.0002897547 0.0001843148        0.97           0         1  0.0006815896   0.005442607 0.0001974957
# mcse.biasB  mcse.biasN  mcse.biasF    RhatB    RhatN true.tau   tau.hatB   tau.hatN   tau.hatF
# 1 0.00253417 0.001706436 0.001318394 1.003183 1.001628     0.05 0.05020189 0.06557535 0.04938595

# informative tau prior and informative beta prior
# true.beta    BayesBbias   BayesNbias    Freqbias BayesBcover BayesNcover Freqcover   BayesBmse
# 1       0.1 -0.0009567751 -0.001921319 0.001558428        0.97           1      0.92 0.000548357
# BayesNmse      Freqmse BayesBpower BayesNpower Freqpower BayesBse2mean BayesNse2mean  Freqse2mean
# 1 0.0002686271 0.0002252963        0.97           0         1  0.0006718406   0.005417716 0.0001946698
# mcse.biasB  mcse.biasN  mcse.biasF    RhatB    RhatN true.tau   tau.hatB  tau.hatN  tau.hatF
# 1 0.002351534 0.001635884 0.001500396 1.002956 1.001826     0.05 0.05026539 0.0653721 0.0486875


# 2. normal prior for all taus
modelBinLinearDRmetaOR <- function(){

  for (i in 1:ns) { ## for each study
    r[i,1] ~ dbinom(p[i,1],n[i,1])
    logit(p[i,1])<- u[i]
    for (j in 2:(nd[i])) { ## for each dose
      # Likelihood
      r[i,j] ~ dbinom(p[i,j],n[i,j])
      logit(p[i,j])<- u[i] + delta[i,j]
      delta[i,j] <-  beta[i]*(X[i,(j)]-X[i,1])

    }
  }
  # Random effect

  for(i in 1:ns) {
    beta[i] ~dnorm(beta.pooled,prec.tau)
    u[i]~dnorm(0,0.1)
  }

  # Priors


  prec.tau<-1/variance
  variance<-tau*tau
  tau~dnorm(0,400)%_%T(0,)
  beta.pooled ~ dnorm(0,10)
}



# Scenario 1
S1ORlineartauN <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[1],OR=TRUE,splines=FALSE)

# Scenario 2
S2ORlineartauN <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=0.01,OR=TRUE,splines=FALSE)


# Scenario 3
S3ORlineartauN <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[2],OR=TRUE,splines=FALSE)


# Scenario 4

S4ORlineartauN <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[3],OR=TRUE,splines=FALSE)






# 3. uniform prior for all taus
modelBinLinearDRmetaOR <- function(){

  for (i in 1:ns) { ## for each study
    r[i,1] ~ dbinom(p[i,1],n[i,1])
    logit(p[i,1])<- u[i]
    for (j in 2:(nd[i])) { ## for each dose
      # Likelihood
      r[i,j] ~ dbinom(p[i,j],n[i,j])
      logit(p[i,j])<- u[i] + delta[i,j]
      delta[i,j] <-  beta[i]*(X[i,(j)]-X[i,1])

    }
  }
  # Random effect

  for(i in 1:ns) {
    beta[i] ~dnorm(beta.pooled,prec.tau)
    u[i]~dnorm(0,0.1)
  }

  # Priors


  prec.tau<-1/variance
  variance<-tau*tau
  #tau~dunif(0,0.1)
  log.tau~dunif(-10,-3)
  tau <- exp(log.tau)
  beta.pooled ~ dnorm(0,20)
}



# Scenario 1
S1ORlineartauU <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[1],OR=TRUE,splines=FALSE)

# Scenario 2
S2ORlineartauU <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=0.01,OR=TRUE,splines=FALSE)


# Scenario 3
S3ORlineartauU <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[2],OR=TRUE,splines=FALSE)


# Scenario 4

S4ORlineartauU <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[3],OR=TRUE,splines=FALSE)

# bind all results and save them
rval <- rbind(unifInftau1 =S1ORlineartau,unifInftau2=S2ORlineartau,unifInftau3=S3ORlineartau,
              nortau1=S1ORlineartauN,nortau2=S2ORlineartauN,nortau3=S3ORlineartauN,nortau4=S4ORlineartauN,
              uniftau1=S1ORlineartauU,uniftau2=S2ORlineartauU,uniftau3=S3ORlineartauU,uniftau4=S4ORlineartauU)
write.csv(rval,file=paste0(Sys.Date(),'varytaupriorORlinear.csv'))

rvalBin <- rval[,c('true.beta','BayesBbias','true.tau',   'tau.hatB')]

write.csv(rvalBin,file=paste0(Sys.Date(),'varytaupriorORlinearBin.csv'))



## Normal DR linear model with strongly informative uniform prior

modelNorLinearDRmeta <- function(){
  b[1] <-0
  for (i in 1:ns) { ## for each study

    # Likelihood

    Y[i,1:(nd[i])]  ~ dmnorm(mean[i,1:(nd[i])], prec[(b[i]+1):(b[i]+nd[i]),1:(nd[i])])

    mean[i,1:(nd[i])] <-  beta[i]*(X[i, 2:(nd[i]+1)]-X[i, 1])

    b[i+1] <- b[i]+ nd[i]
  }

  # Random effect

  for(i in 1:ns) {

    beta[i] ~dnorm(beta.pooled,prec.tau)

  }

  # Priors
  prec.tau<-1/variance
  variance<-tau*tau
  tau~dunif(0.00075,0.00125)#(0.03750,0.06250)#(0.02250,0.03750)#(0.00075,0.00125)
  beta.pooled ~ dnorm(0,0.1)

}


# Scenario 1
S1ORlineartau <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[1],OR=TRUE,splines=FALSE)

# Scenario 2
S2ORlineartau <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[2],OR=TRUE,splines=FALSE)

# Scenario 3
S3ORlineartau <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[3],OR=TRUE,splines=FALSE)

# the estimates of tau are unbiased now


###############################################################
# Based on our 29/10 meeting:
# Binomial model vs Freq, only beta=0.5, 2 different vague priors and OR
#1.	Generate data where dose<-dose/5, use the two vague priors for tau (RE model)
#2.	Fixed effect the    beta[i] ~dnorm(beta.pooled,prec.tau) will be just beta[i]<-beta.pooled and no priors for tau only for beta.pooled
###############################################################

# I did the analysis in these settings
nsim <- 100
beta.pooled <- 0.5
tau <- c(0.001, 0.01,0.03, 0.05)



# 1. RE model normal prior for all taus
modelBinLinearDRmetaOR <- function(){

  for (i in 1:ns) { ## for each study
    r[i,1] ~ dbinom(p[i,1],n[i,1])
    logit(p[i,1])<- u[i]
    for (j in 2:(nd[i])) { ## for each dose
      # Likelihood
      r[i,j] ~ dbinom(p[i,j],n[i,j])
      logit(p[i,j])<- u[i] + delta[i,j]
      delta[i,j] <-  beta[i]*(X[i,(j)]-X[i,1])

    }
  }
  # Random effect

  for(i in 1:ns) {
    beta[i] ~ dnorm(beta.pooled,prec.tau)
    u[i]~dnorm(0,0.001)
  }

  # Priors


  prec.tau<-1/variance
  variance<-tau*tau
  tau~dnorm(0,400)%_%T(0,)
  beta.pooled ~ dnorm(0,0.001)
}



# Scenario 1
S1ORlineartauREN <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[1],dose=c(0,5),OR=TRUE,splines=FALSE)

# Scenario 2
S2ORlineartauREN <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[2],dose=c(0,5),OR=TRUE,splines=FALSE)


# Scenario 3
S3ORlineartauREN <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[3],dose=c(0,5),OR=TRUE,splines=FALSE)


# Scenario 4

S4ORlineartauREN <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[4],dose=c(0,5),OR=TRUE,splines=FALSE)






# 3. RE model with uniform prior for all taus
modelBinLinearDRmetaOR <- function(){

  for (i in 1:ns) { ## for each study
    r[i,1] ~ dbinom(p[i,1],n[i,1])
    logit(p[i,1])<- u[i]
    for (j in 2:(nd[i])) { ## for each dose
      # Likelihood
      r[i,j] ~ dbinom(p[i,j],n[i,j])
      logit(p[i,j])<- u[i] + delta[i,j]
      delta[i,j] <-  beta[i]*(X[i,(j)]-X[i,1])

    }
  }
  # Random effect

  for(i in 1:ns) {
    beta[i] ~ dnorm(beta.pooled,prec.tau)
    u[i]~dnorm(0,0.001)
  }

  # Priors


  prec.tau<-1/variance
  variance<-tau*tau
  tau~dunif(0,0.1)
  # log.tau~dunif(-10,-3)
  # tau <- exp(log.tau)
  beta.pooled ~ dnorm(0,0.001)
}



# Scenario 1
S1ORlineartauREU <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[1],dose=c(0,5),OR=TRUE,splines=FALSE)

# Scenario 2
S2ORlineartauREU <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[2],dose=c(0,5),OR=TRUE,splines=FALSE)


# Scenario 3
S3ORlineartauREU <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[3],dose=c(0,5),OR=TRUE,splines=FALSE)


# Scenario 4

S4ORlineartauREU <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[4],dose=c(0,5),OR=TRUE,splines=FALSE)

# bind all results and save them
rvalRE <- rbind(nortau1=S1ORlineartauREN,nortau2=S2ORlineartauREN,nortau3=S3ORlineartauREN,nortau4=S4ORlineartauREN,
              uniftau1=S1ORlineartauREU,uniftau2=S2ORlineartauREU,uniftau3=S3ORlineartauREU,uniftau4=S4ORlineartauREU)
write.csv(rvalRE,file=paste0(Sys.Date(),'varytauRE.csv'))




simulateDRmeta.fun=function(beta1.pooled=0.01,beta2.pooled=0.02,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,p0=0.1,OR=TRUE,splines=TRUE){ #

  # This function generate a dataset based on linear dose-response model.
  # Arguments:
  # beta.pooled: (vector) the underlying slope of the linear dose-response model and rcs trnasformation
  # tau: the heterogenity across studies.
  # ns: number of studies
  # doserange: the range of the doses in the generated dataset.
  # sample size: it is not the actual sample size but it is used to draw a sample size for each dose from a uinform distribution around this value(200), i.e. Unif(samplesize-20,samplesize+20)
  #  Within each study, each dose assumed to have the same drawn sample size
  library(Hmisc)
  # 1. Create the doses and its transformation
  d<-cbind(rep(0,ns),matrix(round(c(runif(ns,doserange[1],doserange[2]/2),runif(ns,doserange[2]/2+1,doserange[2])),2),nrow=ns))##
  dose <- c(t(d))

  #nr of observations: I assume each study has 3 levels of doses
  nobs<-ns*3


  # 2. Create the dose-specific logOR,
  if(splines==TRUE){

    knots<-unlist(round(quantile(d[,2:3],c(0.25,0.5,0.75))))
    trans.d<-rcspline.eval(c(t(d)),knots,inclx = T)
    dose1 <- c(trans.d[,1])
    dose2 <- c(trans.d[,2])

    maxlogRR<- (beta1.pooled)*max(trans.d[,1]) +(beta2.pooled)*max(trans.d[,2])#
    beta1<-beta1.pooled #random effects of the slopes
    beta2<-beta2.pooled #random effects of the slopes
    logrr<- beta1*trans.d[,1]+beta2*trans.d[,2]   #derive study-specific logOR using regression

  }else{
    maxlogRR<- (beta1.pooled)*max(dose) #the maximum possible value of logRR
    beta  <- beta1.pooled #random effects of the slopes
    logrr <- beta*dose  #derive study-specific logOR using regression
    dose1 <- dose2 <- dose
  }


  # 3.
  if(OR==TRUE){ ## odds ratio
    odds0 <- p0/(1-p0)
    odds1 <- exp(logrr)*odds0
    p1<- odds1/(1+odds1)#estimate p1 the probability of a case at any dose level

    uniquess<-round(runif(ns,samplesize-20,samplesize+20))#sample size of each study arm
    ss<-c(sapply(uniquess,rep,3))#sample size for all study arms
    cases<- matrix(rbinom(ns*3,ss,p1),nrow=3) #drow the cases for all levels
    noncases<-matrix(c(ss-cases),nrow=3)
    hatodds<-cases/noncases
    hatlogOR<-t(log(apply(hatodds,1,"/",hatodds[1,]))) #sample logOR
    SEhatlogOR<-sqrt(1/cases[1,]+1/cases[c(2,3),]+1/(noncases[1,]) + 1/(noncases[c(2,3),])) #SE of the sample logOR
    SEhatlogOR<-rbind(rep(NA,ns),SEhatlogOR)

    hatlogrr <- hatlogOR
    SEhatlogrr <- SEhatlogOR
    type=rep('cc',3*ns)
  }else{ ## Risk ratio
    maxRR<-exp(maxlogRR)
    p0<-ifelse(maxRR>10,0.05,0.5/maxRR)#set p0 to be half the maximum allowed, just to be on the safe side!
    p1 <-ifelse(exp(logrr)*p0>1,0.95,exp(logrr)*p0)

    uniquess<-round(runif(ns,samplesize-20,samplesize+20))#sample size in study arm of zero dose
    ss<-c(sapply(uniquess,rep,3)) #sample size per study arm
    cases<-matrix(rbinom(ns*3,ss,p1),nrow = 3)     #events per study at zero dose
    noncases<-matrix(c(ss-cases),nrow = 3)     #events per study at zero dose
    hatRR <- cases[2:3,]/cases[1,]
    hatlogRR <- log(rbind(rep(1,ns),hatRR))

    #a much easier way to calculate the SE
    SEhatlogRR<-sqrt(1/cases[2:3,]+1/cases[1,]-2/uniquess)
    selogRR<-c(rbind(NA,SEhatlogRR))

    hatlogrr <- hatlogRR
    SEhatlogrr <- selogRR
    type=rep('ci',3*ns)
  }


  Study_No<-rep(1:ns,each=3)


  simulatedDRdata<-cbind.data.frame(Study_No=Study_No,logrr=c(hatlogrr),dose1=dose1,dose2=dose2,cases=c(cases),noncases=c(noncases),
                                    selogrr =c(SEhatlogrr), type=type)

  # To avoid errors that occur due to having zero or n cases I replace them by 1 or n-1, respectively.
  simulatedDRdata$cases[simulatedDRdata$cases==0]  <-1
  simulatedDRdata$noncases[simulatedDRdata$noncases==0] <-1

  return(simulatedDRdata=simulatedDRdata)
}
#

# 1. FE model normal prior for all taus
modelBinLinearDRmetaOR <- function(){

  for (i in 1:ns) { ## for each study
    r[i,1] ~ dbinom(p[i,1],n[i,1])
    logit(p[i,1])<- u[i]
    for (j in 2:(nd[i])) { ## for each dose
      # Likelihood
      r[i,j] ~ dbinom(p[i,j],n[i,j])
      logit(p[i,j])<- u[i] + delta[i,j]
      delta[i,j] <-  beta[i]*(X[i,(j)]-X[i,1])

    }
  }
  # Random effect

  for(i in 1:ns) {
    beta[i] <- beta.pooled
    u[i]~dnorm(0,0.001)
  }

  # Priors


  prec.tau<-1/variance
  variance<-tau*tau
  tau~dnorm(0,400)%_%T(0,)
  beta.pooled ~ dnorm(0,0.001)
}



# Scenario 1
S1ORlineartauFEN <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[1],OR=TRUE,splines=FALSE)

# Scenario 2
S2ORlineartauFEN <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[2],OR=TRUE,splines=FALSE)


# Scenario 3
S3ORlineartauFEN <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[3],OR=TRUE,splines=FALSE)


# Scenario 4

S4ORlineartauFEN <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[4],OR=TRUE,splines=FALSE)






# 3. FE model with uniform prior for all taus
modelBinLinearDRmetaOR <- function(){

  for (i in 1:ns) { ## for each study
    r[i,1] ~ dbinom(p[i,1],n[i,1])
    logit(p[i,1])<- u[i]
    for (j in 2:(nd[i])) { ## for each dose
      # Likelihood
      r[i,j] ~ dbinom(p[i,j],n[i,j])
      logit(p[i,j])<- u[i] + delta[i,j]
      delta[i,j] <-  beta[i]*(X[i,(j)]-X[i,1])

    }
  }
  # Random effect

  for(i in 1:ns) {
    beta[i] <- beta.pooled
    u[i]~dnorm(0,0.001)
  }

  # Priors


  prec.tau<-1/variance
  variance<-tau*tau
  tau~dunif(0,0.1)
  # log.tau~dunif(-10,-3)
  # tau <- exp(log.tau)
  beta.pooled ~ dnorm(0,0.001)
}



# Scenario 1
S1ORlineartauFEU <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[1],OR=TRUE,splines=FALSE)

# Scenario 2
S2ORlineartauFEU <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[2],OR=TRUE,splines=FALSE)


# Scenario 3
S3ORlineartauFEU <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[3],OR=TRUE,splines=FALSE)


# Scenario 4

S4ORlineartauFEU <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[4],OR=TRUE,splines=FALSE)

# bind all results and save them
rvalFE <- rbind(nortau1=S1ORlineartauFEN,nortau2=S2ORlineartauFEN,nortau3=S3ORlineartauFEN,nortau4=S4ORlineartauFEN,
                uniftau1=S1ORlineartauFEU,uniftau2=S2ORlineartauFEU,uniftau3=S3ORlineartauFEU,uniftau4=S4ORlineartauFEU)
write.csv(rvalFE,file=paste0(Sys.Date(),'varytauFE.csv'))


### orsini's scenarios: ORsplines beta1=-0.32 and beta2 =0.02, tau small, same dose=c(0,10)
# priors for beta1.pooled, beta2.pooled and u ~ dnorm(0,0.001)
# prior for tau ~ dnorm(0,400)%_%T(0,)
nsim <- 500
beta1 <- -0.32
beta2 <- 0.02
tau <- 0.001
ns <- 20
OR=FALSE
#splineDRmetaFreq<-dosresmeta(formula = logrr~dose1+dose2, id = Study_No,type=type,
#                             se = selogrr, cases = cases, n = cases+noncases, data = sim.data, proc='1stage',method = 'reml',covariance = 'gl')

orsiniORspline <- simpower(nsim=nsim,beta1.pooled=beta1,beta2.pooled=beta2,tau=tau,ns=ns,OR=FALSE,splines=TRUE)
simulateDRmeta.fun(beta1.pooled=beta1,beta2.pooled=beta2,tau=tau,ns=20,doserange=c(1, 10),samplesize=200,OR=OR,splines = TRUE)
# run the simulation 10'000 times only for freq with 1stage or 2stage models
freqfun <- function(beta1.pooled=beta1.pooled,beta2.pooled=beta2.pooled,tau,OR=FALSE){

  sim.data <- simulateDRmeta.fun(beta1.pooled=beta1.pooled,beta2.pooled=beta2.pooled,tau=tau,ns=20,doserange=c(1, 10),samplesize=200,OR=OR,splines = TRUE)

  freq2<-dosresmeta(formula = logrr~dose1+dose2, id = Study_No,type=type,
                         se = selogrr, cases = cases, n = cases+noncases, data = sim.data, proc='2stage',method = 'reml',covariance = 'gl')
  freq1<-dosresmeta(formula = logrr~dose1+dose2, id = Study_No,type=type,
                         se = selogrr, cases = cases, n = cases+noncases, data = sim.data, proc='1stage',method = 'reml',covariance = 'gl')

    return(cbind(coef(freq1),coef(freq2)))
}
nrep <- 10000
tauhat1 <- replicate(nrep=3,freqfun(beta1.pooled=beta1,beta2.pooled=beta2,tau = tau))


## simulate the antidepressant scenario for OR: beta1=0.25

nsim <- 1000
beta1 <- 0.02
beta2 <- -0.03
tau <- 0.007
antidepSim20 <-simpower(nsim=nsim,beta1.pooled=beta1,beta2.pooled=beta2,tau=tau,OR=TRUE,splines=TRUE)
antidepSim40 <-simpower(nsim=nsim,beta1.pooled=beta1,beta2.pooled=beta2,tau=tau,OR=TRUE,ns=40,splines=TRUE)


