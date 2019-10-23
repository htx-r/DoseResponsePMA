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
  tau~dunif(0.03750,0.06250)#(0.02250,0.03750)#(0.00075,0.00125)
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
S2ORlineartauN <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[2],OR=TRUE,splines=FALSE)


# Scenario 3

S3ORlineartauN <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[3],OR=TRUE,splines=FALSE)







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
  tau~dunif(0,0.05)
  beta.pooled ~ dnorm(0,20)
}



# Scenario 1
S1ORlineartauU <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[1],OR=TRUE,splines=FALSE)

# Scenario 2
S2ORlineartauU <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=0.01,OR=TRUE,splines=FALSE)


# Scenario 3
S3ORlineartauU <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[2],OR=TRUE,splines=FALSE)


# Scenario 4

S3ORlineartauU <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[3],OR=TRUE,splines=FALSE)
mm <- rbind(r1 =S3ORlineartauU,r2=S2ORlineartauU)
mm[,c('true.beta','BayesBbias','BayesNbias','Freqbias','true.tau',   'tau.hatB',   'tau.hatN',   'tau.hatF')]





