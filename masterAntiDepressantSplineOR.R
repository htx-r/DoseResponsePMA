# Libraries
library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force=TRUE)
library(DoseResponseNMA)
# Data
antiDep <-  read.csv('~/Desktop/TasnimPhD/DoseResponseNMA1/DOSEmainanalysis.csv')
#antiDep <-  read.csv("~/Google Drive/_mydrive/HTx/HTx-R/DoseResponseNMA/DOSEmainanalysis.csv")
NAstudyid <- antiDep$Study_No[is.na(antiDep$logRR)]
antiDep1 <-antiDep[!antiDep$Study_No %in% NAstudyid,]

# add noncases to the data and studyid
antiDep1$nonResponders <- antiDep1$No_randomised- antiDep1$Responders
antiDep1$studyid <- as.numeric(as.factor(antiDep1$Study_No))

# Compute the OR and selogOR
ORfun <- function(p0,p1) p1*(1-p0)/(p0*(1-p1))
selogORfun <- function(c0,n0,c1,n1) sqrt(1/c0+1/n0+1/c1+1/n1)
OR <- sapply(unique(antiDep1$studyid), function(i) c(1,ORfun(p0=antiDep1$Responders[antiDep1$studyid==i&antiDep1$hayasaka_ddd==0]/antiDep1$No_randomised[antiDep1$studyid==i&antiDep1$hayasaka_ddd==0],p1=antiDep1$Responders[antiDep1$studyid==i&antiDep1$hayasaka_ddd!=0]/antiDep1$No_randomised[antiDep1$studyid==i&antiDep1$hayasaka_ddd!=0])),simplify = F)
selogOR <- sapply(unique(antiDep1$studyid), function(i) c(NA,selogORfun(c0=antiDep1$Responders[antiDep1$studyid==i&antiDep1$hayasaka_ddd==0],n0=antiDep1$No_randomised[antiDep1$studyid==i&antiDep1$hayasaka_ddd==0],c1=antiDep1$Responders[antiDep1$studyid==i&antiDep1$hayasaka_ddd!=0],n1=antiDep1$No_randomised[antiDep1$studyid==i&antiDep1$hayasaka_ddd!=0])),simplify = F)

# compute OR and SE for study '56' by hand, since it does not have 0 dose, I will do it later more generic
p0 <- 47/111
p1 <- 90/110
OR56 <- ORfun(p0,p1)
OR[[56]] <- c(1, OR56)
selogOR[[56]] <- c(NA, selogORfun(47,11,90,110))

antiDep1$logOR <- log(unlist(OR))
antiDep1$selogOR <- unlist(selogOR)

# dose rcs transformation
dose  <- as.matrix(rcs(antiDep1$hayasaka_ddd))
antiDep1$dose1 <- dose[,1]
antiDep1$dose2 <- dose[,2]
antiDep1$dose3 <-dose[,3]
attr(dose,'ATT') <- NULL
# Jags data
jagsdataLinear <- makeJAGSDRmeta(studyid=studyid,logrr =  logOR,dose1 =  dose1,dose2=dose2,cases=Responders,noncases=nonResponders,data=antiDep1,Splines=T)

# Anaylsis:
# a.Frequentist: dosresmeta

splineDRmetaFreq <-dosresmeta(formula = logRR ~ rcs(hayasaka_ddd,knots=c(10,20,30)), type = type, id = Study_No,
                              se = selogRR, cases = Responders, n = No_randomised  , data = antiDep1,covariance = 'gl',proc = '1stage')#!!!!!!!!!!!!!!

# b. Bayes: JAGS model

prec <-  sapply(linearDRmetaFreq$Slist,solve,simplify = F)

tncomp <- sum(as.numeric(table(antiDep1$Study_No))-1) ## total number of non-zero comparisons
max.nd <- max(as.numeric(table(antiDep1$Study_No)))

precmat <- matrix(NA,tncomp,max.nd-1)
b <- vector()

for (i in 1:ns) {
  b[1] <- 0
  nd[i] <- as.numeric(table(antiDep1$Study_No)[i])-1
  precmat[(b[i]+1):(b[i]+nd[i]),1:(nd[i])] <- prec[[i]]
  b[i+1] <- b[i]+ nd[i]
  precmat
}
jagsdataLinear$prec <- precmat
jagsdataLinear$nd <- jagsdataLinear$nd-1
linearDRmetaJAGSmodel <- jags(data = jagsdataLinear,inits=NULL,parameters.to.save = c('beta.pooled','tau'),model.file = modelLinearDRmeta,
                              n.chains=2,n.iter = 100000,n.burnin = 20000,DIC=F,n.thin = 10)
traceplot(linearDRmetaJAGSmodel,varname='beta.pooled') ## looks good
traceplot(linearDRmetaJAGSmodel,varname='tau') ## ALWAYS check the tau for convergence

# c.Compare Bayesian and Frequentist
linearDRmetaJAGSmodel$BUGSoutput$mean$beta.pooled - coef(linearDRmetaFreq)
jagsdataLinearBin <- makeBinomialJAGSDRmeta(studyid=Study_No,dose1 = hayasaka_ddd,dose2=NULL,cases=Responders,noncases=nonResponders,data=antiDep1,Splines=F)

linearDRmetaJAGSmodelBin <- jags.parallel(data = jagsdataLinearBin,inits=NULL,parameters.to.save = c('beta.pooled','beta','tau'),model.file = modelBinomialLinearDRmeta,
                                          n.chains=2,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 1)
linearDRmetaJAGSmodelBin$BUGSoutput$mean$beta.pooled


