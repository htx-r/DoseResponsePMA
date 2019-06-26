library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force = T)

library(DoseResponseNMA)

## Cubic Splines

# compare dosresmeta vs true values
m1 <- m2 <- f1<-f2 <- vector()
n.sim.data <- 1000
ns <- 20
start <- Sys.time()


for (j in 1:n.sim.data) {

  sim <- simulateDRsplinedata.fun(beta1.pooled=0.01,beta2.pooled=0.02,tau=0.001,doserange = c(1,10))
  sim.data <- sim$simulatedDRdata
  knots <- sim$knots
  # Freq
  rcsplineDRmetaFreq <- dosresmeta(formula = logRR~rcs(sim.data$dose,knots), id = Study_No,type=type,
                                   se = selogRR, cases = cases, n = cases+noncases, data = sim.data, proc='1stage',method = 'fixed')
  m <- lm(formula = logRR ~  rcs(sim.data$dose,knots)-1, data = sim.data)

  f1[j] <-coef(rcsplineDRmetaFreq)[1]
  f2[j] <-coef(rcsplineDRmetaFreq)[2]
  m1[j] <-coef(m)[1]
  m2[j] <-coef(m)[2]

}


fres1 <-mean(f1)
fres2 <- mean(f2)
mres1 <- mean(m1)
mres2 <- mean(m2)
## Almost as the true value! 0.01 and 0.02


## Linear

# compare dosresmeta vs true values
f1 <- vector()
n.sim.data <- 1000
for (j in 1:n.sim.data) {

  sim.data <- simulateDRlineardata.fun(beta.pooled=0.02,tau=0.001)

  # Freq
  linearDRmetaFreq=dosresmeta(formula = logRR~dose, id = Study_No,type=type,
                              se = selogRR, cases = cases, n = cases+noncases, data = sim.data, proc='2stage',method = 'fixed')
  m <- lm(formula = logRR ~ dose - 1, data = sim.data)
  f1[j] <-coef(linearDRmetaFreq)[1]

}

resf <- mean(f1) ##
bias<- 0.02-resf ## I changed beta.pooled several times and normally I get 0.03 bias!


## Quadratic
m1 <- m2 <-f1<-f2 <- vector()
n.sim.data <- 1000
for (i in 1:n.sim.data) {

  sim.data <- simulateDRquadraticdata.fun(beta1.pooled=0.01,beta2.pooled=0.02,tau=0.001,doserange = c(1,10))

   # Freq
  #quadraticDRmetaFreq=dosresmeta(formula = logRR~dose+I(dose^2), id = Study_No,type=type,
   #                              se = selogRR, cases = cases, n = cases+noncases, data = sim.data, proc='1stage',method = 'fixed')#!!!!!!!!!!!!!!
  m <- lm(formula = logRR ~ dose+I(dose^2) - 1, data = sim.data)

  f1[j] <-coef(quadraticDRmetaFreq)[1]
  f2[j] <-coef(quadraticDRmetaFreq)[2]

  m1[j] <-coef(m)[1]


  m2[j] <-coef(m)[2]


}
mean(f1,na.rm = T)
mean(f2,na.rm = T)
mean(m1,na.rm = T)
mean(m2,na.rm = T)

#

























# # Data: exclude the 9 studies with NA log RR.
# #antiDep <-  read.csv('~/Desktop/TasnimPhD/DoseResponseNMA1/DOSEmainanalysis.csv')
# antiDep <-  read.csv('C:/Users/th19g983/Desktop/DoseResponseNMA/DOSEmainanalysis.csv')
# #antiDep <-  read.csv("~/Google Drive/_mydrive/HTx/HTx-R/DoseResponseNMA/DOSEmainanalysis.csv")
# NAstudyid <- antiDep$Study_No[is.na(antiDep$logRR)]
# antiDep1 <-antiDep[!antiDep$Study_No %in% NAstudyid,]
