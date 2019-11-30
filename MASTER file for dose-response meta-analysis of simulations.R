           ##################################################################################
           #         Master analysis for dose-response meta. analysis of Simulated data
           ##################################################################################

#load libraries
source('Functions needed for dose-response meta-analysis of simulations.R')
library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force = T)
library(DoseResponseNMA)
library('rsimsum')
library(tidyr)

###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           # Linear: 1. odds ratio (OR) 2. risk ratio (RR)
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsim <- 3
beta.pooled <- c(0,0.02,0.04,0.06,0.1,0.2,0.3)
tau <- c(0.001,0.01)

####### 1. odds ratio (OR)

## %% smaller tau
# Scenario 1
#set.seed('145')
S1ORlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[1],OR=TRUE,splines=FALSE)


# Scenario 2
#set.seed('245')
S2ORlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[2],tau=tau[1],OR=TRUE,splines=FALSE)

# Scenario 3
#set.seed('345')

S3ORlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[3],tau=tau[1],OR=TRUE,splines=FALSE)


# Scenario 4
#set.seed('445')

S4ORlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[4],tau=tau[1],OR=TRUE,splines=FALSE)


# Scenario 5
#set.seed('545')
S5ORlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[5],tau=tau[1],OR=TRUE,splines=FALSE)


# Scenario 6:
#set.seed('645')

S6ORlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[6],tau=tau[1],OR=TRUE,splines=FALSE)

# Scenario 7:
#set.seed('745')

S7ORlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[7],tau=tau[1],OR=TRUE,splines=FALSE)

## %% larger tau
# Scenario 8:
#set.seed('845')

S8ORlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[2],OR=TRUE,splines=FALSE)

# Scenario 9:
#set.seed('945')
S9ORlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[2],tau=tau[2],OR=TRUE,splines=FALSE)

# Scenario 10:
#set.seed('1045')

S10ORlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[3],tau=tau[2],OR=TRUE,splines=FALSE)

# Scenario 11:
#set.seed('1145')

S11ORlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[4],tau=tau[2],OR=TRUE,splines=FALSE)

# Scenario 12:
#set.seed('1245')

S12ORlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[5],tau=tau[2],OR=TRUE,splines=FALSE)

# Scenario 13:
#set.seed('1345')
S13ORlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[6],tau=tau[2],OR=TRUE,splines=FALSE)

# Scenario 14:
#set.seed('1445')
S14ORlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[7],tau=tau[2],OR=TRUE,splines=FALSE)

# Save the results in a file
resORlinear <- rbind(S1ORlinear,S2ORlinear,S3ORlinear,S4ORlinear,S5ORlinear,S6ORlinear,S7ORlinear,
                     S8ORlinear,S9ORlinear,S10ORlinear,S11ORlinear,S12ORlinear,S13ORlinear,S14ORlinear)
write.csv(resORlinear,file=paste0(Sys.Date(),'resORlinear.csv')) # keeps the rownames
# end for linear OR model

#### 2. risk ratio (RR)
## %% smaller tau
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


# Scenario 5
#set.seed('523')
S5RRlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[5],tau=tau[1],OR=FALSE,splines = FALSE)


# Scenario 6:
#set.seed('623')

S6RRlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[6],tau=tau[1],OR=FALSE,splines = FALSE)

# Scenario 7:
# set.seed('723')
S7RRlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[7],tau=tau[1],OR=FALSE,splines = FALSE)

## %% larger tau
# Scenario 8:
#set.seed('823')

S8RRlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[1],tau=tau[2],OR=FALSE,splines = FALSE)

# Scenario 9:
#set.seed('923')

S9RRlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[2],tau=tau[2],OR=FALSE,splines = FALSE)

# Scenario 10:
#set.seed('1023')

S10RRlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[3],tau=tau[2],OR=FALSE,splines = FALSE)

# Scenario 11:
#set.seed('1123')

S11RRlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[4],tau=tau[2],OR=FALSE,splines = FALSE)

# Scenario 12:
#set.seed('1223')

S12RRlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[5],tau=tau[2],OR=FALSE,splines = FALSE)

# Scenario 13:
#set.seed('131323')
S13RRlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[6],tau=tau[2],OR=FALSE,splines = FALSE)

# Scenario 14:
# set.seed('1423')
 S14RRlinear <- simpower(nsim=nsim,beta1.pooled=beta.pooled[7],tau=tau[2],OR=FALSE,splines = FALSE)

# Save the results in a file
#eval(parse(text=paste0('S',1:2,'RRlinear')))
resRRlinear <- rbind(S1RRlinear,S2RRlinear,S3RRlinear,S4RRlinear,S5RRlinear,S6RRlinear,S8RRlinear,S9RRlinear,S10RRlinear,S11RRlinear,S12RRlinear,S14RRlinear)#,S7RRlinear,S14RRlinear)
write.csv(resRRlinear,file=paste0(Sys.Date(),"RRlinear.csv")) # keeps the rownames

# end for RR linear model

#simulateDRmeta.fun(beta1.pooled=beta.pooled[6],tau=tau[2],OR=FALSE,splines = FALSE)

###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         # Restricted Cubic Splines (RCS): 1. odds ratio (OR) 2. risk ratio (RR)
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

S7ORspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[3],beta2.pooled = beta2.pooled[3],tau=tau[2],OR=TRUE,ns=ns,splines = TRUE)

# Scenario 8:
set.seed('822')

S8ORspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[4],beta2.pooled = beta2.pooled[4],tau=tau[2],OR=TRUE,ns=ns,splines = TRUE)






resORspline40sim1000 <- rbind(S1ORspline$res1,S2ORspline$res1,S3ORspline$res1,S4ORspline$res1,
                              S5ORspline$res1,S6ORspline$res1,S7ORspline$res1,S8ORspline$res1)
resORspline40sim1000ALL <- rbind(S1ORspline$res2,S2ORspline$res2,S3ORspline$res2,S4ORspline$res2,
                                 S5ORspline$res2,S6ORspline$res2,S7ORspline$res2,S8ORspline$res2)

write.csv2(resORspline40sim1000,file=paste0(Sys.Date(),"resORspline40sim1000.csv"),sep=" ") # keeps the rownames
write.ftable(resORspline40sim1000,file=paste0(Sys.Date(),"resORspline40sim1000.txt"))
save(resORspline40sim1000ALL,file='resORspline40sim1000ALL')
save(resORspline40sim1000,file='resORspline40sim1000')

# end of OR spline model

#### 2. risk ratio (RR)

## %% smaller tau
# Scenario 1
#set.seed('197')
S1RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[1],beta2.pooled = beta2.pooled[1],tau=tau[1],OR=FALSE,splines = TRUE)

# Scenario 2
#set.seed('297')

S2RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[2],tau=tau[1],OR=FALSE,splines = TRUE)

# Scenario 3
#set.seed('397')

S3RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[3],beta2.pooled = beta2.pooled[3],tau=tau[1],OR=FALSE,splines = TRUE)

# Scenario 4
#set.seed('497')

S4RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[4],beta2.pooled = beta2.pooled[4],tau=tau[1],OR=FALSE,splines = TRUE)

# Scenario 5
#set.seed('597')
S5RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[5],beta2.pooled = beta2.pooled[5],tau=tau[1],OR=FALSE,splines = TRUE)

# Scenario 6
#set.seed('697')

S6RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[6],beta2.pooled = beta2.pooled[6],tau=tau[1],OR=FALSE,splines = TRUE)

# Scenario 7
#set.seed('797')

S7RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[7],beta2.pooled = beta2.pooled[7],tau=tau[1],OR=FALSE,splines = TRUE)

## %% Larger tau

# Scenario 8:
#set.seed('897')

S8RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[1],beta2.pooled = beta2.pooled[1],tau=tau[2],OR=FALSE,splines = TRUE)

# Scenario 9
#set.seed('997')

S9RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[2],beta2.pooled = beta2.pooled[2],tau=tau[2],OR=FALSE,splines = TRUE)


# Scenario 10
#set.seed('1097')

S10RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[3],beta2.pooled = beta2.pooled[3],tau=tau[2],OR=FALSE,splines = TRUE)


#Scenario 11
#set.seed('1197')

S11RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[4],beta2.pooled = beta2.pooled[4],tau=tau[2],OR=FALSE,splines = TRUE)


# Scenario 12
#set.seed('1279')

S12RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[5],beta2.pooled = beta2.pooled[5],tau=tau[2],OR=FALSE,splines = TRUE)


# Scenario 13
#set.seed('1397')

S13RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[6],beta2.pooled = beta2.pooled[6],tau=tau[2],OR=FALSE,splines = TRUE)

# Scenario 14
#set.seed('1497')
S14RRspline <- simpower(nsim=nsim,beta1.pooled = beta1.pooled[7],beta2.pooled = beta2.pooled[7],tau=tau[2],OR=FALSE,splines = TRUE)

data("cc_ex")

## Fitting the model
mod1 <- dosresmeta(formula = logrr ~ dose, type = "cc", cases = case,
                   n = n, lb = lb, ub = ub, data= cc_ex)

# Save the results in a file

resRRspline <- rbind(S1RRspline,S2RRspline,S3RRspline,S4RRspline,S5RRspline,S8RRspline,S9RRspline,S10RRspline, S11RRspline)#, S12RRspline)S6RRspline,S7RRspline,S13RRspline,S14RRspline)
write.csv(resRRspline,file=paste0(Sys.Date(),"resRRspline.csv")) # keeps the rownames

# end of risk ratio spline model

