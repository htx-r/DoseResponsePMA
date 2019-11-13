source('FunctionsForSimulations.R')
library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force = T)
library(DoseResponseNMA)
library('rsimsum')
library(tidyr)

nsim <- 100
beta1.pooled <- c(0,0.04,0.1,0.2)
beta2.pooled <- c(0,0,0.03,-0.2 )
tau <- c(0.001,0.01)
ns <- 40


### 1. odds ratio (OR)
## %% smaller tau
# Scenario 1
set.seed('122')
S1ORspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[1],beta2.pooled = beta2.pooled[1],tau=tau[1],OR=TRUE,ns=ns,splines = TRUE)

# Scenario 2
set.seed('222')

S2ORspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[2],tau=tau[1],OR=TRUE,ns=ns,splines = TRUE)

# Scenario 3
set.seed('322')

S3ORspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[3],beta2.pooled = beta2.pooled[3],tau=tau[1],OR=TRUE,ns=ns,splines = TRUE)

# Scenario 4
set.seed('422')

S4ORspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[4],beta2.pooled = beta2.pooled[4],tau=tau[1],OR=TRUE,ns=ns,splines = TRUE)

## %% Larger tau

# Scenario 5
set.seed('522')
S5ORspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[1],beta2.pooled = beta2.pooled[1],tau=tau[2],OR=TRUE,ns=ns,splines = TRUE)

# Scenario 6
set.seed('622')

S6ORspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[2],tau=tau[2],OR=TRUE,ns=ns,splines = TRUE)

# Scenario 7
set.seed('722')

S7ORspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[3],beta2.pooled = beta2.pooled[3],tau=tau[3],OR=TRUE,ns=ns,splines = TRUE)

# Scenario 8:
set.seed('822')

S8ORspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[4],beta2.pooled = beta2.pooled[4],tau=tau[2],OR=TRUE,ns=ns,splines = TRUE)






resORspline40sim1000 <- rbind(S1ORspline$res1,S2ORspline$res1,S3ORspline$res1,S4ORspline$res1,S5ORspline$res1,S6ORspline$res1,S7ORspline$res1,S8ORspline$res1,S9ORspline$res1,S10ORspline$res1)
resORspline40sim1000ALL <- rbind(S1ORspline$res2,S2ORspline$res2,S3ORspline$res2,S4ORspline$res2,S5ORspline$res2,S6ORspline$res2,S7ORspline$res2,S8ORspline$res2,S9ORspline$res2,S10ORspline$res2)

write.csv(resORspline40sim1000,file=paste0(Sys.Date(),"resORspline40sim1000.csv")) # keeps the rownames
save(resORspline40sim1000ALL,file='resORspline40sim1000ALL')
save(resORspline40sim1000,file='resORspline40sim1000')

# end of OR spline model

