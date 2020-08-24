#load libraries and functions needed
source('Functions needed for dosres MA simulations.R')
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


# Scenario 2: ns=8
start <- Sys.time()
set.seed('123')
S2ORspline2 <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[3],tau=tau[1],OR=TRUE,ns=ns[2],splines = TRUE)
end <- Sys.time()

#











