           ##################################################################################
           #         Master analysis for dose-response meta-analysis of Simulated data
           ##################################################################################

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
ns <- 40


### 1. odds ratio (OR): spline
## %% small tau

# Scenario 1
start <- Sys.time()
set.seed('122')
S1ORspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[1],beta2.pooled = beta2.pooled[1],tau=tau[1],OR=TRUE,ns=ns,splines = TRUE)
end <- Sys.time()

# Scenario 2
set.seed('222')

S2ORspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[2],tau=tau[1],OR=TRUE,ns=ns,splines = TRUE)

# Scenario 3
set.seed('322')

S3ORspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[3],beta2.pooled = beta2.pooled[3],tau=tau[1],OR=TRUE,ns=ns,splines = TRUE)

# Scenario 4
set.seed('422')

S4ORspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[4],beta2.pooled = beta2.pooled[4],tau=tau[1],OR=TRUE,ns=ns,splines = TRUE)


## %% large tau

# Scenario 5
set.seed('522')
S5ORspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[1],beta2.pooled = beta2.pooled[1],tau=tau[2],OR=TRUE,ns=ns,splines = TRUE)

# Scenario 6
set.seed('622')

S6ORspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[2],tau=tau[2],OR=TRUE,ns=ns,splines = TRUE)

# Scenario 7
set.seed('722')

S7ORspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[3],beta2.pooled = beta2.pooled[3],tau=tau[2],OR=TRUE,ns=ns,splines = TRUE)

# Scenario 8:
set.seed('822')

S8ORspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[4],beta2.pooled = beta2.pooled[4],tau=tau[2],OR=TRUE,ns=ns,splines = TRUE)

# combine all results ...
resORspline40sim1000 <- rbind(S1ORspline$res1,S2ORspline$res1,S3ORspline$res1,S4ORspline$res1,
                              S5ORspline$res1,S6ORspline$res1,S7ORspline$res1,S8ORspline$res1)
resORspline40sim1000ALL <- rbind(S1ORspline$res2,S2ORspline$res2,S3ORspline$res2,S4ORspline$res2,
                                 S5ORspline$res2,S6ORspline$res2,S7ORspline$res2,S8ORspline$res2)
# ... and save them
write.csv(resORspline40sim1000,file=paste0(Sys.Date(),"resORspline40sim1000.csv"),sep=" ")
write.csv(resORspline40sim1000ALL,file=paste0(Sys.Date(),"resORspline40sim1000ALL.csv"))
save(resORspline40sim1000ALL,file='resORspline40sim1000ALL')
save(resORspline40sim1000,file='resORspline40sim1000')
save(S1ORspline,S2ORspline,S3ORspline,S4ORspline,
     S5ORspline,S6ORspline,S7ORspline,S8ORspline,file = 'ORspline')
# end of OR spline model

#### 2. risk ratio (RR): spline

## %% small tau
# Scenario 1
set.seed('197')
S1RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[1],beta2.pooled = beta2.pooled[1],tau=tau[1],ns=ns,OR=FALSE,splines = TRUE)

# Scenario 2
set.seed('297')

S2RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[2],tau=tau[1],ns=ns,OR=FALSE,splines = TRUE)

# Scenario 3
set.seed('397')

S3RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[3],beta2.pooled = beta2.pooled[3],tau=tau[1],ns=ns,OR=FALSE,splines = TRUE)

# Scenario 4
set.seed('497')

S4RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[4],beta2.pooled = beta2.pooled[4],tau=tau[1],ns=ns,OR=FALSE,splines = TRUE)

## %% large tau

# Scenario 5
set.seed('597')
S5RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[1],beta2.pooled = beta2.pooled[1],tau=tau[2],ns=ns,OR=FALSE,splines = TRUE)

# Scenario 6
set.seed('697')

S6RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[2],tau=tau[2],ns=ns,OR=FALSE,splines = TRUE)

# Scenario 7
set.seed('797')

S7RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[3],beta2.pooled = beta2.pooled[3],tau=tau[2],ns=ns,OR=FALSE,splines = TRUE)


# Scenario 8:
set.seed('897')

S8RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[4],beta2.pooled = beta2.pooled[4],tau=tau[2],ns=ns,OR=FALSE,splines = TRUE)

# combine the results ...

# resRRspline <- rbind(S1RRspline$res1,S2RRspline$res1,S3RRspline$res1,S4RRspline$res1,S5RRspline$res1,S8RRspline$res1)
# resRRsplineALL <- rbind(S1RRspline$res2,S2RRspline$res2,S3RRspline$res2,S4RRspline$res2,S5RRspline$res2,S8RRspline$res2)

# ... and save them
# write.csv(resRRspline,file=paste0(Sys.Date(),"resRRspline40sim1000.csv"))
# save(resRRsplineALL,file='resRRspline40sim1000ALL')
# save(resRRspline,file='resRRspline40sim1000')
resRRspline40sim1000 <- rbind(S1RRspline$res1,S2RRspline$res1,S3RRspline$res1,S4RRspline$res1,
                              S5RRspline$res1,S6RRspline$res1,S7RRspline$res1,S8RRspline$res1)
resRRspline40sim1000ALL <- rbind(S1RRspline$res2,S2RRspline$res2,S3RRspline$res2,S4RRspline$res2,
                                 S5RRspline$res2,S6RRspline$res2,S7RRspline$res2,S8RRspline$res2)


load('SRR5678')

write.csv(resRRspline40sim1000,file=paste0(Sys.Date(),"resRRspline40sim1000.csv"))
write.csv(resRRspline40sim1000ALL,file=paste0(Sys.Date(),"resRRspline40sim1000ALL.csv"))
save(resRRspline40sim1000ALL,file='resRRspline40sim1000ALL')
save(resRRspline40sim1000,file='resRRspline40sim1000')
save(S1RRspline,S2RRspline,S3RRspline,S4RRspline,
     S5RRspline,S6RRspline,S7RRspline,S8RRspline,file = 'RRspline')

#

########################################################################
#*** Additional simulations as SMMR reviewers asked

# spline + OR
nsim <- 1000
beta1.pooled <- c(0,0.04,0.1,0.2)
beta2.pooled <- c(0,0,0.03,-0.2 )
tau <- c(0.001,0.01)

#** 1. Few trials
# Scenario 2: ns=8
start <- Sys.time()
set.seed('123')
S2ORspline2 <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[3],tau=tau[1],OR=TRUE,ns=8,splines = TRUE)
end <- Sys.time()
save(S2ORspline2,file='S2ORspline2')

# Scenario 7: ns=16
start <- Sys.time()
set.seed('1234')
S7ORspline2 <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[3],tau=tau[1],OR=TRUE,ns=16,splines = TRUE)
end <- Sys.time()
save(S7ORspline2,file='S7ORspline2')

#** 2. Overlapped doses
# Scenario 3: d~ unif(1,6) U unif(4,10)
start <- Sys.time()
set.seed('123')
S3ORspline2 <- simpowerDose(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[3],tau=tau[1],OR=TRUE,ns=20,splines = TRUE)
end <- Sys.time()
save(S3ORspline2,file='S3ORspline2')


#** 3. Discrete doses
# Scenario 6: d~sample(1,10)
start <- Sys.time()
set.seed('126')
S6ORspline2 <- simpowerDiscDose(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[3],tau=tau[1],OR=TRUE,ns=20,splines = TRUE)
end <- Sys.time()
save(S6ORspline2,file='S6ORspline2')


#** 4. small sample size
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

# 5. different dose shapes

#** Scenario 9: half sigmoid & uniform & harrell's knots
nsim <- 1000
set.seed('125')
S9ORspline2 <- replicate(nsim,OneSimulationSigm(ns=ns,doserange = c(1,10),samplesize = 200,OR=TRUE,splines = TRUE),simplify = T)
save(S9ORspline2,file='S9ORspline2')


#** Scenario 10: half sigmoid & chi2 & harrell's knots
nsim <- 1000
set.seed('129')
S10ORspline2 <- replicate(nsim,OneSimulationSigmChi(ns=ns,doserange = c(1,10),samplesize = 200,OR=TRUE,splines = TRUE),simplify = T)
save(S10ORspline2,file='S10ORspline2')

#** Scenario 11: half sigmoid & uniform & my knots
nsim <- 1000
set.seed('1252')
S11ORspline2 <- replicate(nsim,OneSimulationSigm(ns=ns,doserange = c(1,10),samplesize = 200,OR=TRUE,splines = TRUE),simplify = T)
save(S11ORspline2,file='S11ORspline2')


#** Scenario 12: half sigmoid & chi2 & my knots
nsim <- 1000
set.seed('1292')
S12ORspline2 <- replicate(nsim,OneSimulationSigmChi(ns=ns,doserange = c(1,10),samplesize = 200,OR=TRUE,splines = TRUE),simplify = T)
save(S12ORspline2,file='S12ORspline2')
#

#** Scenario 13: Log & uniform & harrell's knots
nsim <- 1000
set.seed('1253')
S13ORspline2 <- replicate(nsim,OneSimulationLog(ns=20,doserange = c(1,10),samplesize = 200,OR=TRUE,splines = TRUE),simplify = T)
save(S13ORspline2,file='S13ORspline2')


#** Scenario 14:Log & chi2 & harrell's knots
nsim <- 1000
set.seed('12934')
S14ORspline2 <- replicate(nsim,OneSimulationLogChi(ns=20,doserange = c(1,10),samplesize = 200,OR=TRUE,splines = TRUE),simplify = T)
save(S14ORspline2,file='S14ORspline2')

#** Scenario 15: Log& uniform & my knots (0,2,4)
nsim <- 1000
set.seed('12515')
S15ORspline2 <- replicate(nsim,OneSimulationLog(ns=20,doserange = c(1,10),samplesize = 200,OR=TRUE,splines = TRUE),simplify = T)
save(S15ORspline2,file='S15ORspline2')


#** Scenario 16: Log & chi2 & my knots
nsim <- 1000
set.seed('1692')
S16ORspline2 <- replicate(nsim,OneSimulationLogChi(ns=20,doserange = c(1,10),samplesize = 200,OR=TRUE,splines = TRUE),simplify = T)
save(S16ORspline2,file='S16ORspline2')
#








