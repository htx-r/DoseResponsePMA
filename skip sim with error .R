
freqfun <- function(tau,OR=TRUE){

  sim.data <- try(simulateDRmeta.fun(beta1.pooled=0,beta2.pooled=0,tau=tau,ns=20,doserange=c(1, 5),samplesize=200,OR=TRUE,splines = TRUE),silent = TRUE)
if(class(sim.data)=='try-error'){
  coefFreq <- c(NA,NA)
}else{
  coefFreq<-dosresmeta(formula = logrr~dose1+dose2, id = Study_No,type=type,
                       se = selogrr, cases = cases, n = cases+noncases, data = sim.data, proc='1stage',method = 'reml',covariance = 'gl')$coefficients
  }
#   m <- try(simulateDRmeta.fun(beta1.pooled=0,beta2.pooled=0,tau=tau,ns=20,doserange=c(1, 5),samplesize=200,OR=TRUE,splines = TRUE),silent = TRUE)
#
#   if(class(m)=='try-error') {
#     m <- try(simulateDRmeta.fun(beta1.pooled=0,beta2.pooled=0,tau=tau,ns=20,doserange=c(1, 5),samplesize=200,OR=TRUE,splines = TRUE),silent = TRUE)
# print('na')
#   }


  return(coefFreq)
}


set.seed('122')
nrep <- 10
tau <- 0.001
rval <-t(replicate(nrep,freqfun(tau = tau),simplify = TRUE))
colMeans(rval,na.rm = T)

#
set.seed('124')
nrep=3
replicate(n=2,OneSimulation(beta1.pooled=0,beta2.pooled=0,tau=tau,ns=20,doserange=c(1, 5),samplesize=200,OR=TRUE,splines = TRUE))
res <- simpower(nsim=2,beta1.pooled=0,beta2.pooled=0,tau=tau,ns=20,doserange=c(1, 4.5),samplesize=200,OR=TRUE,splines = TRUE)
res

v <- 'try-error'
c <- 0
while (v=='try-error') {
  zz <- try(simulateDRmeta.fun(beta1.pooled=0,beta2.pooled=0,tau=tau,ns=40,doserange=c(1, 4),samplesize=200,OR=TRUE,splines = TRUE),silent = TRUE)
v<- class(zz)
c <- c+1
  }
c


inputs = list(1, 2, 4, -5, 'oops', 0, 10)

for(input in inputs) {
     print(paste("log of", input, "=", log(input)))
}

# skip the error and execute the rest
for(input in inputs) {
     try(print(paste("log of", input, "=", log(input))))
   }

# you might want substitute your own return value when errors (or warnings) are returned.
for(input in inputs) {
       tryCatch(print(paste("log of", input, "=", log(input))),
                               warning = function(w) {print(paste("negative argument", input));
                   log(-input)},
                             error = function(e) {print(paste("non-numeric argument", input));
                   NaN})
   }




m <- list(list(name=1, data=data.frame()),list(name=2, data=data.frame()))

sim.data <- simulateDRmeta.fun(beta1.pooled=0.01,beta2.pooled=0.02,tau=tau,ns=20,doserange=c(1, 10),samplesize=200,OR=TRUE,splines = TRUE)
jagsdata<- makejagsDRmeta(Study_No,logrr,dose1,dose2,cases,noncases,se=selogrr,type=type,data=sim.data,Splines=T,new.dose.range = c(1,10))
jagsdata$rr <- jagsdata$r
knots<-unlist(round(quantile(jagsdata$new.dose,c(0.25,0.5,0.75))))
jagsdata$f.new.dose <- rcspline.eval(jagsdata$new.dose,knots,inclx = T)[,2]
jagsdata$nd.new <- length(jagsdata$new.dose)
splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau','Z','p.drug','OR'),model.file = modelBinSplineDRmetaOR,
                                          n.chains=2,n.iter = 1000,n.burnin =200,DIC=F,n.thin = 1)

splineDRmetaJAGSmodelBin$BUGSoutput$mean$Z
exp(Z)
p.drug <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$p.drug
c(jagsdata$X1[,2:3])
c(odds)
par(mfrow=c(1,1))
myd <- data.frame(dose=c(jagsdata$X1[,2:3]),odds=c(odds))

myd <- myd[order(myd$dose),]
plot(sort(c(odds)), ylim = c(0.9,1.4))
abline(h=1)
np <- nrow(sim.data[sim.data$dose1==0,])


## antidepressant
knots = c(10,20,50)
new.dose <- seq(0,80,1)
new.dose1 <- c(rcs(new.dose,knots)[,1])
new.dose2 <- c(rcs(new.dose,knots)[,2])

jagsdataORspline<- makejagsDRmeta(studyid=studyid,logOR,dose1=dose1,dose2=dose2,cases=Responders,noncases=nonResponders,se=selogOR,type=type,data=antidep,Splines=T,new.dose.range = c(5,10))
jagsdataORspline$np <- 58# sum(jagsdataORspline$X1[,1]==0)
jagsdataORspline$nn <- jagsdataORspline$n[-56,]
jagsdataORspline$rr <- jagsdataORspline$r[-56,]
jagsdataORspline$new.dose <-  seq(1,80,1)
jagsdataORspline$f.new.dose <- rcspline.eval(jagsdataORspline$new.dose,knots,inclx = T)[,2]
jagsdataORspline$nd.new <- length(jagsdataORspline$new.dose)

doseresORsplineBin <- jags.parallel(data = jagsdataORspline,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau','Z','p.drug','OR'),model.file = modelBinSplineDRmetaOR,
                                    n.chains=3,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 1)

p.placebo<- exp(doseresORsplineBin$BUGSoutput$mean$Z)/(1+exp(doseresORsplineBin$BUGSoutput$mean$Z))
exp(doseresORsplineBin$BUGSoutput$summary['Z',])/(1+exp(doseresORsplineBin$BUGSoutput$summary['Z',]))

p.drug <- doseresORsplineBin$BUGSoutput$mean$p.drug
plot(jagsdataORspline$new.dose,p.drug,type='l',ylim = c(0.3,0.5),xlab = 'dose')
abline(h=p.placebo,lty=2)
mean(jagsdataORspline$rr[,1]/jagsdataORspline$nn[,1])







dose <- jagsdataORspline$new.dose
p.placebo<- exp(doseresORsplineBin$BUGSoutput$mean$Z)/(1+exp(doseresORsplineBin$BUGSoutput$mean$Z))
p.drug <- doseresORsplineBin$BUGSoutput$mean$p.drug
l.ci <-  doseresORsplineBin$BUGSoutput$summary[-c(1,2,3,84),'2.5%']
u.ci <- doseresORsplineBin$BUGSoutput$summary[-c(1,2,3,84),'97.5%']

plot(dose,p.drug,type='l',ylim = c(0.3,2),lwd=2,las=1,ylab='response')
lines(dose,l.ci,lty=2,lwd=2)
lines(dose,u.ci,lty=2,lwd=2)

abline(h=p.placebo,col=2,lwd=2)


doseresORsplineNor$BUGSoutput$summary['beta1.pooled','sd']
doseresORsplineNor$BUGSoutput$summary['beta2.pooled','sd']

doseresORsplineNor00 <- jags.parallel(data = jagsdataORspline,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau','mean'),model.file = modelNorSplineDRmeta,
                                    n.chains=3,n.iter = 1000,n.burnin = 20,DIC=F,n.thin = 1)







