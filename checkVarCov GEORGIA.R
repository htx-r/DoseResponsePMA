library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force = T)

library(DoseResponseNMA)



diffincov<-c()
diffincorr<-c()
for(i in 1:100){
# Data
sim.data <- simulateDRsplinedata.fun(beta1.pooled=0.02,beta2.pooled=0.05,tau=0.0001,doserange = c(1,10))
data <- sim.data$simulatedDRdata
data<-data[data$Study_No==1,]# WE DO NOT NEED MANY STUDIES, ONE IS ENOUGH TO CHECK!
data$studyid <- data$Study_No
data$n<-data$cases+data$noncases
data$pr<-data$cases/(data$n)
data$mylogRR<-log(data$pr/data$pr[1])
data$myselogRR<-sqrt(1/data$cases-1/data$n+1/data$cases[1]-1/data$n[1])
data$myselogRR[1]<-NA

data#####!!!!! NOTICE THAT YOUR LOGRR AND MYLOGRR ARE DIFFERENT. I USE THE SAMPLE LOG RR AND YOU USE THE UNDERLYING - THIS IS NOT CORRECT.
#CONSEQUENTLY, THE GL METHOD TO ESTIMATE CORRELATIONS THINKS THAT YOUR LOGRR ARE ADJUSTED FOR SOME FACTORS AND GIVES OF COURSE DIFFERENT RESULTS


#  Compute the var-cov "by hand"
covariance<-1/data$cases[1]-1/data$n[1]
varcovar1<-matrix(covariance,nrow=2,ncol=2)
diag(varcovar1)<-(data$selogRR^2)[-1]


# Compute the var-cov from doseresmeta
rcsplineDRmetaFreq <- dosresmeta(formula = logRR~dose1+dose2, id = Study_No,type=type,
                                 se = selogRR, cases = cases, n = n, data = data, covariance = 'gl')
varcovar2<-rcsplineDRmetaFreq$Slist[[1]]
varcovar1-varcovar2

# Compute the var-cov from covar.logrr() function
varcovar3<-covar.logrr(data$cases,data$n,data$logRR,data$selogRR^2,data$type)
varcovar1-varcovar2
varcovar1-varcovar3
varcovar3-varcovar2

diffincov<-c(diffincov,(varcovar1[1,2]-varcovar2[1,2]))
diffincorr<-c(diffincorr,(varcovar1[1,2]-varcovar2[1,2])/sqrt(varcovar1[1,1]*varcovar2[2,2])  )
             }

mean(diffincov)
mean(diffincorr)


