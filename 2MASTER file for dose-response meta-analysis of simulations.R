#load libraries and functions needed
source('Functions needed for dosres MA simulations.R')
source('functions for additional simulations.R')
library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponsePMA",force = T)
library(DoseResponseNMA)
library('rsimsum')
library(tidyr)
#library(rsimsum)
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Restricted Cubic Splines (RCS): 1. odds ratio (OR) 2. risk ratio (RR)
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# spline simulation settings for OR and RR
nsim <- 1000
beta1.pooled <- c(0,0.04,0.1,0.2)
beta2.pooled <- c(0,0,0.03,-0.2 )
tau <- c(0.001,0.01)
ns <- c(4,8)


### 1. odds ratio (OR): spline
## %% small tau

# Scenario 1: ns=4
start <- Sys.time()
set.seed('122')
S1ORspline2 <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[3],tau=tau[1],OR=TRUE,ns=ns[1],splines = TRUE)
end <- Sys.time()
save(S1ORspline2,file='S1ORspline2')
load('S1ORspline2')

# Scenario 2: ns=8
start <- Sys.time()
set.seed('123')
S2ORspline2 <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[3],tau=tau[1],OR=TRUE,ns=ns[2],splines = TRUE)
end <- Sys.time()
save(S2ORspline2,file='S2ORspline2')


# Scenario 3: d~ unif(1,7) U unif(2,10)
ns <- 20
start <- Sys.time()
set.seed('123')
S3ORspline2 <- simpowerDose(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[3],tau=tau[1],OR=TRUE,ns=ns,splines = TRUE)
end <- Sys.time()
save(S3ORspline2,file='S3ORspline2')

# Scenario 4: samplesize~ unif(20,100) ''instead of (180,220)''
ns <- 20
start <- Sys.time()
set.seed('124')
S4ORspline2 <- simpowerSample(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[3],tau=tau[1],OR=TRUE,ns=ns,splines = TRUE)
end <- Sys.time()
save(S4ORspline2,file='S4ORspline2')


# Scenario 5: RCS ---> quadratic
start <- Sys.time()
set.seed('125')
S5ORspline2 <- simpowerQuad(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[3],tau=tau[1],OR=TRUE,ns=ns,splines = TRUE)
end <- Sys.time()
save(S5ORspline2,file='S5ORspline2')


# Scenario 6: discrete doses
# see the doses distribution in antidep dataset
#
mydata <-  read.csv('DOSEmainanalysis.csv')
antidep=mydata[mydata$exc==F,]
library(MASS)
truehist(antidep$hayasaka_ddd,h=1,prob = F,las=1,xlab='hayasaka dose')

# change the dose generation from dose~U(0, 10) to dose~smaple(1,10)

start <- Sys.time()
set.seed('126')
S6ORspline2 <- simpowerDiscDose(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[3],tau=tau[1],OR=TRUE,ns=20,splines = TRUE)
end <- Sys.time()
save(S6ORspline2,file='S6ORspline2')


# Scenario 7: ns=15
start <- Sys.time()
set.seed('1234')
S7ORspline2 <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[3],tau=tau[1],OR=TRUE,ns=16,splines = TRUE)
end <- Sys.time()
save(S7ORspline2,file='S7ORspline2')


# Scenario 8: samplesize~(20,100)
start <- Sys.time()
set.seed('1281')
S8ORspline2.1 <- simpowerSample(nsim=nsim,beta1.pooled = beta1.pooled[1],beta2.pooled = beta2.pooled[1],tau=tau[1],OR=TRUE,ns=20,splines = TRUE)
end <- Sys.time()
save(S8ORspline2.1,file='S8ORspline2.1')

set.seed('1282')
S8ORspline2.2 <- simpowerSample(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[2],tau=tau[1],OR=TRUE,ns=20,splines = TRUE)
save(S8ORspline2.2,file='S8ORspline2.2')

set.seed('1283')
S8ORspline2.3 <- simpowerSample(nsim=nsim,beta1.pooled = beta1.pooled[3],beta2.pooled = beta2.pooled[3],tau=tau[1],OR=TRUE,ns=20,splines = TRUE)
save(S8ORspline2.3,file='S8ORspline2.3')

set.seed('1284')
S8ORspline2.4 <- simpowerSample(nsim=nsim,beta1.pooled = beta1.pooled[4],beta2.pooled = beta2.pooled[4],tau=tau[1],OR=TRUE,ns=20,splines = TRUE)
save(S8ORspline2.4,file='S8ORspline2.4')

set.seed('1285')
start <- Sys.time()
S8ORspline2.5 <- simpowerSample(nsim=nsim,beta1.pooled = beta1.pooled[1],beta2.pooled = beta2.pooled[1],tau=tau[2],OR=TRUE,ns=20,splines = TRUE)
end <- Sys.time()
save(S8ORspline2.5,file='S8ORspline2.5')

set.seed('1286')
S8ORspline2.6 <- simpowerSample(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[2],tau=tau[2],OR=TRUE,ns=20,splines = TRUE)
save(S8ORspline2.6,file='S8ORspline2.6')

set.seed('1287')
S8ORspline2.7 <- simpowerSample(nsim=nsim,beta1.pooled = beta1.pooled[3],beta2.pooled = beta2.pooled[3],tau=tau[2],OR=TRUE,ns=20,splines = TRUE)
save(S8ORspline2.7,file='S8ORspline2.7')

set.seed('1288')
S8ORspline2.8 <- simpowerSample(nsim=nsim,beta1.pooled = beta1.pooled[4],beta2.pooled = beta2.pooled[4],tau=tau[2],OR=TRUE,ns=20,splines = TRUE)
save(S8ORspline2.8,file='S8ORspline2.8')

# Scenario 9: RCS ---> log.function
nsim <- 20
ns <- 20
start <- Sys.time()
#set.seed('125')
S9ORspline2 <- replicate(nsim,OneSimulationSigm(ns=ns,doserange = c(1,10),samplesize = 200,OR=T,splines = T),simplify = T)

#S9ORspline2 <- simpowerLog(nsim=nsim,beta1.pooled = 1,beta2.pooled = 1,tau=0.001,OR=TRUE,ns=ns,splines = TRUE)
end <- Sys.time()
save(S9ORspline2,file='S9ORspline2')
dd <-simulateDRmeta.funLog(splines = T)
OneSimulationLog(OR=T,splines = T)
#

#




