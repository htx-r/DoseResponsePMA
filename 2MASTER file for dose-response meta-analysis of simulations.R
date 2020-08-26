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


# Scenario 5: discrete doses
# see the doses distribution in antidep dataset
#
mydata <-  read.csv('DOSEmainanalysis.csv')
antidep=mydata[mydata$exc==F,]
library(MASS)
truehist(antidep$hayasaka_ddd,h=1,prob = F,las=1,xlab='hayasaka dose')

# !!!! change the dose generation from dose~U(0, 10) to dose~(0, Unif(1,2), Unif(4,5), Unif(7,8), Unif(9,10))

start <- Sys.time()
set.seed('125')
S5ORspline2 <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[3],tau=tau[1],OR=TRUE,ns=ns[2],splines = TRUE)
end <- Sys.time()
save(S5ORspline2,file='S5ORspline2')









