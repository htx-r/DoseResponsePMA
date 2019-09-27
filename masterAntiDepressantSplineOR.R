# Libraries
library(rms)
library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force=TRUE)
library(DoseResponseNMA)

#########################################
# Data
#########################################

antiDep <-  read.csv('~/Desktop/TasnimPhD/DoseResponseNMA/DoseResponseNMA/DOSEmainanalysis.csv')
#antiDep <-  read.csv("~/Google Drive/_mydrive/HTx/HTx-R/DoseResponseNMA/DOSEmainanalysis.csv")
NAstudyid <- antiDep$Study_No[is.na(antiDep$logRR)]
antiDep1 <-antiDep[!antiDep$Study_No %in% NAstudyid,]

# add noncases to the data and studyid
antiDep1$nonResponders <- antiDep1$No_randomised- antiDep1$Responders
antiDep1$studyid <- as.numeric(as.factor(antiDep1$Study_No))

max.nd <- max(as.numeric(table(antiDep1$Study_No)))
ns <- length(unique(antiDep1$Study_No))
nd <- as.numeric(table(antiDep1$Study_No))-1

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
dose  <- as.matrix(rcs(antiDep1$hayasaka_ddd,knots=c(10,20,30,50)))
rcspline.eval(antiDep1$hayasaka_ddd,knots=c(10,20,30,50),inclx = T)
antiDep1$dose1 <- dose[,1]
antiDep1$dose2 <- dose[,2]
antiDep1$dose3 <-dose[,3]
study_id <- unique(antiDep1$studyid)

#########################################
# Anaylsis:
#########################################
# a.Frequentist: dosresmeta

splineDRmetaFreq <-dosresmeta(formula = logOR ~ dose1+dose2+dose3, type = type, id = Study_No,
                              se = selogOR, cases = Responders, n = No_randomised  , data = antiDep1,covariance = 'gl',proc = '1stage')#!!!!!!!!!!!!!!
splineDRmetaFreq2 <-dosresmeta(formula = logOR ~ rcs(hayasaka_ddd,knots=c(10,20,30)), type = type, id = Study_No,
                              se = selogOR, cases = Responders, n = No_randomised  , data = antiDep1,covariance = 'gl',proc = '1stage')#!!!!!!!!!!!!!!

# b. Bayes: JAGS model
jagsdataSpline <- makeJAGSDRmeta(studyid=studyid,logrr =  logOR,dose1 =  dose1,dose2=dose2,cases=Responders,noncases=nonResponders,data=antiDep1,Splines=T)

prec <-  sapply(splineDRmetaFreq$Slist,solve,simplify = F)
tncomp <- sum(as.numeric(table(antiDep1$Study_No))-1) ## total number of non-zero comparisons
precmat <- matrix(NA,tncomp,max.nd-1)
b <- vector()
for (i in 1:ns) {
  b[1] <- 0
  nd[i] <- as.numeric(table(antiDep1$Study_No)[i])-1
  precmat[(b[i]+1):(b[i]+nd[i]),1:(nd[i])] <- prec[[i]]
  b[i+1] <- b[i]+ nd[i]
  precmat
}
# Jags data
X3mat <- matrix(NA,ns,max.nd)
for (i in 1:ns) {
  X3mat[i,1:as.numeric(table(antiDep1$studyid)[i])] <- antiDep1$dose3[antiDep1$studyid == study_id[i]]
}
jagsdataSpline$X3 <- X3mat

jagsdataSpline$prec <- precmat
jagsdataSpline$nd <- jagsdataSpline$nd-1
splineDRmetaJAGSmodel <- jags(data = jagsdataSpline,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','beta3.pooled','tau'),model.file = modelRCSplineDRmeta,
                              n.chains=2,n.iter = 100000,n.burnin = 20000,DIC=F,n.thin = 10)

# c. Binomial
jagsdataSplineBin <- makeBinomialJAGSDRmeta(studyid=studyid,dose1 = dose1,dose2=dose2,cases=Responders,noncases=nonResponders,data=antiDep1,Splines=T)
jagsdataSplineBin$X3 <- X3mat
splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdataSplineBin,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','beta3.pooled'),model.file = modelBinomialRCSsplineDRmeta,
                                          n.chains=2,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 1)

#########################################
# Compare the three approaches
#########################################
cbind(Bin=splineDRmetaJAGSmodelBin$BUGSoutput$mean,Nor=splineDRmetaJAGSmodel$BUGSoutput$mean[1:3],Freq=coef(splineDRmetaFreq))

## Plots

# new dose range
dosepred <- rcs(seq(0,80),knots=c(10,20,30))


# Freq
beta1.pooled <- coef(splineDRmetaFreq)[1]
beta2.pooled <- coef(splineDRmetaFreq)[2]
beta3.pooled <- coef(splineDRmetaFreq)[3]

ypred <- coef(splineDRmetaFreq)[1]*c(dosepred[,1])+coef(splineDRmetaFreq)[2]*c(dosepred[,2])+coef(splineDRmetaFreq)[3]*c(dosepred[,3])
ORpredF <- exp(ypred)

# Normal bayes
jagsdataSpline$new.n <- nrow(dosepred)
jagsdataSpline$new.dose1 <- c(dosepred[,1])
jagsdataSpline$new.dose2 <- c(dosepred[,2])
jagsdataSpline$new.dose3 <- c(dosepred[,3])

splineDRmetaJAGSmodelPred <- jags(data = jagsdataSpline,inits=NULL,parameters.to.save = c('newY','newOR'),model.file = modelRCSplineDRmeta,
                              n.chains=2,n.iter = 100000,n.burnin = 20000,DIC=F,n.thin = 10)

ORpredNor <- splineDRmetaJAGSmodelPred$BUGSoutput$mean$newOR

# Binomial bayes
jagsdataSplineBin$new.n <- nrow(dosepred)
jagsdataSplineBin$new.dose1 <- c(dosepred[,1])
jagsdataSplineBin$new.dose2 <- c(dosepred[,2])
jagsdataSplineBin$new.dose3 <- c(dosepred[,3])

splineDRmetaJAGSmodelBinPred <- jags.parallel(data = jagsdataSplineBin,inits=NULL,parameters.to.save = c('newY','newOR'),model.file = modelBinomialRCSsplineDRmeta,
                                          n.chains=2,n.iter = 10000,n.burnin = 2000,DIC=F,n.thin = 1)

ORpredBin <- splineDRmetaJAGSmodelBinPred$BUGSoutput$mean$newOR

plot(dosepred[,1],ORpredF,type = 'l',ylim=c(1,2),xlim = c(0,50),ylab = 'OR',xlab='dose',las=1)
lines(0:80,ORpredNor,col=2)
lines(0:80,ORpredBin,col=3)
#with(antiDep1,points(hayasaka_ddd[logOR != 0], exp(logOR[logOR != 0])))
legend('topleft',legend=c('frequentist','normal','binomial'), col=1:3,lty=1,bty = 'n')
#
















# nd <- as.numeric(table(antiDep1$Study_No))-1
# ypred <- c()
# k <- 1:tncomp
# x <- rcs(seq(0,160),knots=c(10,20,30))
# x1 <- c(x[,1])
# x2 <- c(x[,2])
# ypred <- matrix(NA,ns,4)
# plot(0:160, beta1.pooled*x1+beta2.pooled*x2,type='l')
#
# for (i in 1:ns) {
#   ypred[i,1:nd[i]] <-  beta1.pooled*(jagsdataSpline$X1[i, 2:(nd[i]+1)]-jagsdataSpline$X1[i, 1])+ beta2.pooled*(jagsdataSpline$X2[i, 2:(nd[i]+1)]-jagsdataSpline$X2[i,1])   #+beta3[i]*(X3[i, 2:(nd[i]+1)]-X3[i,1])
# }
# plot(sort(splineDRmetaFreq1$fitted.values),type = 'l')
#
# newdata = data.frame(hayasaka_ddd = seq(0, 160))
# xref = min(antiDep1$hayasaka_ddd)
# predict((splineDRmetaFreq2), newdata,xref, exp = TRUE)
# with(predict(doseresRR, newdata, xref, exp = TRUE), {
#   plot(get("rcs(dose, knots)dose"), pred, log = "y", type = "l",
#        xlim = c(0, 160), ylim = c(0.95, 2.5),
#        xlab = "Dose", ylab = "RR", main = c("Splines", text))
#   matlines(get("rcs(dose, knots)dose"), cbind(ci.ub, ci.lb),
#            col = 1, lty = "dashed")
# })
# with(AntidepressantsDOSE, points(dose[logRR != 0], exp(logRR[logRR != 0])))
# with(AntidepressantsDOSE, rug(dose, quiet = TRUE))
# antiDep1$ypred <- ypred
# rcspline.plot(x=antiDep1$hayasaka_ddd,y=ypred)
