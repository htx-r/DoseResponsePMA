# load libraries
library(rms) # for rcs()
library(MASS) # for truehist()
library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponsePMA",force=TRUE)
library(DoseResponseNMA)
library(meta)
devAskNewPage(ask=F)
source('Functions needed for dosres MA antidep.R')

########################################
#     load data an prepare

# load and exclude single arm studies
mydata <-  read.csv('~/Google Drive/DoseResponseNMA/DoseResponsePMA/DOSEmainanalysis.csv')
antidep=mydata[mydata$exc==F,]
sum(antidep$No_randomised)
#
antidep$studyid <- as.numeric(as.factor(antidep$Study_No))
antidep$nonResponders <- antidep$No_randomised- antidep$Responders

# apply the function above to all studies
logORmat <- sapply(unique(antidep$studyid),function(i) createORreference.fun(antidep$Responders[antidep$studyid==i],antidep$No_randomised[antidep$studyid==i]),simplify = FALSE)
logORmat <- do.call(rbind,logORmat)
antidep$logOR <- c(logORmat[,1])
antidep$selogOR <- c(logORmat[,2])

# restricted cubic spline transformation doses
knots = c(10,20,50)
antidep$dose1 <- as.matrix(rcs(antidep$hayasaka_ddd,knots))[,1]
antidep$dose2 <- as.matrix(rcs(antidep$hayasaka_ddd,knots))[,2]

# transform data into jags format
jagsdataRRspline<- makejagsDRmeta(studyid=studyid,logRR,dose1=dose1,dose2=dose2,cases=Responders,noncases=nonResponders,se=selogRR,type=type,data=antidep,splines=T)
jagsdataORspline<- makejagsDRmeta(studyid=studyid,logOR,dose1=dose1,dose2=dose2,cases=Responders,noncases=nonResponders,se=selogOR,type=type,data=antidep,splines=T)

# compute sigma matrix per study
nd    <- as.numeric(table(antidep$studyid)) ## number of all doses (with zero dose)
max.nd <- max(nd)                        ## maximum number of doses
ns <- length(unique(antidep$studyid))       ## number of studies
tncomp <- sum(as.numeric(table(antidep$studyid))-1) ## total number of non-zero comparisons
sigma <- sapply(1:ns, sigmaMat)

# and residual het matrix per study
resMat <- sapply(1:ns, function(i) diag(1,nrow = nd[i]-1)+(1-diag(1,nrow = nd[i]-1))*0.5)

#
sigmamat <- matrix(NA,tncomp,max.nd-1)
s <- matrix(NA, ns,max.nd-1)
b <- no.d <-vector()
index <- 1:tncomp
for (i in 1:ns) {
  b[1] <- 0
  no.d[i] <- as.numeric(table(antidep$studyid)[i])-1
  sigmamat[(b[i]+1):(b[i]+no.d[i]),1:(no.d[i])] <- sigma[[i]]
  s[i,1:no.d[i]] <- index[(b[i]+1):(b[i]+no.d[i])]
  b[i+1] <- b[i]+ no.d[i]
  sigmamat
}

resmat <- matrix(NA,tncomp,max.nd-1)
s <- matrix(NA, ns,max.nd-1)
b <- no.d <-vector()
index <- 1:tncomp
for (i in 1:ns) {
  b[1] <- 0
  no.d[i] <- as.numeric(table(antidep$studyid)[i])-1
  resmat[(b[i]+1):(b[i]+no.d[i]),1:(no.d[i])] <- resMat[[i]]
  s[i,1:no.d[i]] <- index[(b[i]+1):(b[i]+no.d[i])]
  b[i+1] <- b[i]+ no.d[i]
  resmat
}

jagsdataORspline$sigmamat <- sigmamat
jagsdataORspline$resmat <- resmat


doseresORsplineBin2 <- jags.parallel(data = jagsdataORspline,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau','tau.res'),model.file = modelBinSplineDRmetaOR2,
                                    n.chains=2,n.iter = 100,n.burnin = 20,DIC=F,n.thin = 1)

doseresORsplineBin2$BUGSoutput$mean$tau





