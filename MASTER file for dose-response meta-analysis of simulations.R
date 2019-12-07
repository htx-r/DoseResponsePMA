           ##################################################################################
           #         Master analysis for dose-response meta-analysis of Simulated data
           ##################################################################################

#load libraries and functions needed
source('Functions needed for dose-response meta-analysis of simulations.R')
library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force = T)
library(DoseResponseNMA)
library('rsimsum')
library(tidyr)

###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Restricted Cubic Splines (RCS): 1. odds ratio (OR) 2. risk ratio (RR)
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# spline simulation settings for OR and RR
nsim <- 1
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
end-start
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
save(resORspline40sim1000ALL,file='resORspline40sim1000ALL')
save(resORspline40sim1000,file='resORspline40sim1000')

# end of OR spline model

#### 2. risk ratio (RR): spline

## %% small tau
# Scenario 1
set.seed('197')
S1RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[1],beta2.pooled = beta2.pooled[1],tau=tau[1],OR=FALSE,splines = TRUE)

# Scenario 2
set.seed('297')

S2RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[2],tau=tau[1],OR=FALSE,splines = TRUE)

# Scenario 3
set.seed('397')

S3RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[3],beta2.pooled = beta2.pooled[3],tau=tau[1],OR=FALSE,splines = TRUE)

# Scenario 4
set.seed('497')

S4RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[4],beta2.pooled = beta2.pooled[4],tau=tau[1],OR=FALSE,splines = TRUE)

## %% large tau

# Scenario 5
set.seed('597')
S5RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[1],beta2.pooled = beta2.pooled[1],tau=tau[2],OR=FALSE,splines = TRUE)

# Scenario 6
set.seed('697')

S6RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[2],tau=tau[2],OR=FALSE,splines = TRUE)

# Scenario 7
set.seed('797')

S7RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[3],beta2.pooled = beta2.pooled[3],tau=tau[2],OR=FALSE,splines = TRUE)


# Scenario 8:
set.seed('897')

S8RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[4],beta2.pooled = beta2.pooled[4],tau=tau[2],OR=FALSE,splines = TRUE)

# combine the results ...

resRRspline <- rbind(S1RRspline$res1,S2RRspline$res1,S3RRspline$res1,S4RRspline$res1,S5RRspline$res1,S8RRspline$res1)
resRRsplineALL <- rbind(S1RRspline$res2,S2RRspline$res2,S3RRspline$res2,S4RRspline$res2,S5RRspline$res2,S8RRspline$res2)

# ... and save them
write.csv(resRRspline,file=paste0(Sys.Date(),"resRRspline40sim1000.csv"))
save(resRRsplineALL,file='resRRspline40sim1000ALL')
save(resRRspline,file='resRRspline40sim1000')

# end of risk ratio spline model












###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           # Linear: 1. odds ratio (OR) 2. risk ratio (RR)
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# linear simulation settings
nsim <- 3
beta.pooled <- c(0,0.03,0.1,0.3)
tau <- c(0.001,0.01)
ns <- 40

####### 1. odds ratio (OR): linear

## %% small tau
# Scenario 1
#set.seed('145')
S1ORlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[1],ns=ns,OR=TRUE,splines=FALSE)


# Scenario 2
#set.seed('245')
S2ORlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[2],tau=tau[1],ns=ns,OR=TRUE,splines=FALSE)

# Scenario 3
#set.seed('345')

S3ORlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[3],tau=tau[1],ns=ns,OR=TRUE,splines=FALSE)


# Scenario 4
#set.seed('445')

S4ORlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[4],tau=tau[1],ns=ns,OR=TRUE,splines=FALSE)

## %% large tau

# Scenario 5
#set.seed('545')
S5ORlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[2],ns=ns,OR=TRUE,splines=FALSE)


# Scenario 6:
#set.seed('645')

S6ORlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[2],tau=tau[2],ns=ns,OR=TRUE,splines=FALSE)

# Scenario 7:
#set.seed('745')

S7ORlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[3],tau=tau[2],ns=ns,OR=TRUE,splines=FALSE)

# Scenario 8:
#set.seed('845')
S8ORlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[4],tau=tau[2],ns=ns,OR=TRUE,splines=FALSE)


# combine the results ...
resORlinear <- rbind(S1ORlinear$res1,S2ORlinear$res1,S3ORlinear$res1,S4ORlinear$res1,S5ORlinear$res1,S6ORlinear$res1,S7ORlinear$res1,S8ORlinear$res1)
resORlinearALL <- rbind(S1ORlinear$res2,S2ORlinear$res2,S3ORlinear$res2,S4ORlinear$res2,S5ORlinear$res2,S6ORlinear$res2,S7ORlinear$res2,S8ORlinear$res2)

# ... and save them
save(resORlinearALL,file='resORlinear40sim1000ALL')
save(resORlinear,file='resORlinear40sim1000')
write.csv(resORlinear,file=paste0(Sys.Date(),'resORlinear40sim1000.csv')) # keeps the rownames
# end for linear OR model

#### 2. risk ratio (RR)
## %% small tau
# Scenario 1
#set.seed('123')
S1RRlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[1],OR=FALSE,splines = FALSE)


# Scenario 2
#set.seed('223')
S2RRlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[2],tau=tau[1],OR=FALSE,splines = FALSE)

# Scenario 3
#set.seed('323')

S3RRlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[3],tau=tau[1],OR=FALSE,splines = FALSE)


# Scenario 4
#set.seed('423')
S4RRlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[4],tau=tau[1],OR=FALSE,splines = FALSE)

## %% large tau

# Scenario 5
#set.seed('523')
S5RRlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[2],OR=FALSE,splines = FALSE)


# Scenario 6:
#set.seed('623')

S6RRlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[2],tau=tau[2],OR=FALSE,splines = FALSE)

# Scenario 7:
# set.seed('723')
S7RRlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[3],tau=tau[2],OR=FALSE,splines = FALSE)

# Scenario 8:
#set.seed('823')

S8RRlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[4],tau=tau[2],OR=FALSE,splines = FALSE)


# combine the results ...
resRRlinear <- rbind(S1RRlinear$res1,S2RRlinear$res1,S3RRlinear$res1,S4RRlinear$res1,S5RRlinear$res1,S6RRlinear$res1,S8RRlinear$res1)
resRRlinearALL <- rbind(S1RRlinear$res2,S2RRlinear$res2,S3RRlinear$res2,S4RRlinear$res2,S5RRlinear$res2,S6RRlinear$res2,S8RRlinear$res2)

# ... and save them
save(resRRlinearALL,file='resRRlinear40sim1000ALL')
save(resRRlinear,file='resRRlinear40sim1000')
write.csv(resRRlinear,file=paste0(Sys.Date(),"resRRlinear40sim1000.csv")) # keeps the rownames

# end for RR linear model


