library(R2jags)
library(dosresmeta)
library(devtools)
install_github("htx-r/DoseResponseNMA",force = T)

library(DoseResponseNMA)

# Greenland-Longnecker Variance-coavriance matrix for RR
logRRvarcov <- function(cases,casesRef){
  varcov <- matrix(1/casesRef , nrow = length(cases),ncol = length(cases))
  diag(varcov) <- 1/casesRef+ 1/cases
  return(varcov)
}

# Greenland-Longnecker Variance-coavriance matrix for OR
logORvarcov <- function(cases,casesRef,noncases,noncasesRef){
  varcov <- matrix(1/casesRef -1/noncasesRef, nrow = length(cases),ncol = length(cases))
  diag(varcov) <- 1/casesRef+ 1/cases -1/noncasesRef- 1/noncases
  return(varcov)
}

# Data
sim.data <- simulateDRsplinedata.fun(beta1.pooled=0.1,beta2.pooled=0.05,tau=0.001,doserange = c(1,10))
data <- sim.data$simulatedDRdata
data<-data[data$Study_No==1,]
data$studyid <- data$Study_No
data$pr<-data$cases/(data$cases+data$noncases)
data$mylogRR<-log(data$pr/data$pr[1])
data$myselogRR<-sqrt(1/data$cases-1/(data$cases+data$noncases)+1/data$cases[1]-1/(data$cases[1]+data$noncases[1]))
data$myselogRR[1]<-NA
data#####!!!!! NOTICE THAT YOUR LOGRR AND MYLOGRR ARE DIFFERENT. I USE THE SAMPLE LOG RR AND YOU USE THE UNDERLYING - THIS IS NOT CORRECT.
#CONSEQUENTLY, THE GL METHOD TO ESTIMATE CORRELATIONS THINKS THAT YOUR LOGRR ARE ADJUSTED FOR SOME FACTORS AND GIVES OF COURSE DIFFERENT RESULTS


#  Compute the var-cov using the three approaches

# 1
prRR <- sapply(unique(data$studyid), function(i) logRRvarcov(cases=data$cases[data$studyid==i&data$dose1!=0],
                                                              casesRef =data$cases[data$studyid==i&data$dose1==0]),simplify = F)
# 2
prOR <- sapply(unique(data$studyid), function(i) logORvarcov(cases=data$cases[data$studyid==i&data$dose1!=0], noncases=data$noncases[data$studyid==i&data$dose1!=0],
         casesRef =data$cases[data$studyid==i&data$dose1==0],noncasesRef =data$noncases[data$studyid==i&data$dose1==0]),simplify = F)

# 3
rcsplineDRmetaFreq <- dosresmeta(formula = logRR~dose1+dose2, id = Study_No,type=type,
                                 se = sqrt(selogRR), cases = cases, n = cases+noncases, data = sim.data$simulatedDRdata, proc='1stage',covariance = 'gl')


# 4 I use covar.logrr() function directly from dosresmeta which is basically used in 'dosresmeta' source file
   # to find out Slist (you can check the source file)
ns <- length(unique(sim.data$simulatedDRdata$Study_No))
res <- (sapply(1:ns, function(i) with(sim.data$simulatedDRdata[sim.data$simulatedDRdata$Study_No==i,],covar.logrr(cases,cases+noncases,logRR,selogRR^2,type)),simplify = F))


##
# Compare the results of the four approaches
prRR[[1]]
prOR[[1]]
rcsplineDRmetaFreq$Slist[[1]]
res[[1]]
## I expected at least res and Slist to be exactly the same because in both of them covar.logrr()function is used!
# However, among all, res is closest to prOR.


## Run JAGS model under each one of the 4 approaches and see the results.
jagsdataRCS<- makeJAGSDRmeta(Study_No,logRR,dose1,dose2,cases,noncases,data=sim.data$simulatedDRdata,Splines=T,new.dose.range = c(5,10))

# 1. use prRR
# I computed invS_i as 1.11 formula in page 10
invS_i <- sapply(1:ns, function(i) with(sim.data$simulatedDRdata[sim.data$simulatedDRdata$Study_No==i,],t(cbind(dose1,dose2)[2:3,])%*%solve(prRR[[i]])%*%cbind(dose1,dose2)[2:3,]),simplify = F)

jagsdataRCS$prec <-  matrix(unlist(invS_i),40,2,byrow = T)
rcsplineDRmetaJAGSmodel <- jags.parallel(data = jagsdataRCS,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau1','tau2','newRR'),model.file = modelRCSplineDRmeta,
                                         n.chains=2,n.iter = 10000,n.burnin = 200,DIC=F,n.thin = 1)
# Error! the values in invS_i are large!
# 2. use prOR
invS_i <- sapply(1:ns, function(i) with(sim.data$simulatedDRdata[sim.data$simulatedDRdata$Study_No==i,],t(cbind(dose1,dose2)[2:3,])%*%solve(prOR[[i]])%*%cbind(dose1,dose2)[2:3,]),simplify = F)

jagsdataRCS$prec <-  matrix(unlist(invS_i),40,2,byrow = T)
rcsplineDRmetaJAGSmodel <- jags.parallel(data = jagsdataRCS,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau1','tau2','newRR'),model.file = modelRCSplineDRmeta,
                                         n.chains=2,n.iter = 10000,n.burnin = 200,DIC=F,n.thin = 1)
# Error again! the values in invS_i are large!

# 3. use Slist
jagsdataRCS$prec <-  matrix(unlist(sapply(rcsplineDRmetaFreq$Slist,solve,simplify = F)),40,2,byrow = T)
rcsplineDRmetaJAGSmodel <- jags.parallel(data = jagsdataRCS,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau1','tau2','newRR'),model.file = modelRCSplineDRmeta,
                                         n.chains=2,n.iter = 10000,n.burnin = 200,DIC=F,n.thin = 1)

# 4. use covar.logrr()
invS_i <- sapply(1:ns, function(i) with(sim.data$simulatedDRdata[sim.data$simulatedDRdata$Study_No==i,],t(cbind(dose1,dose2)[2:3,])%*%solve(res[[i]])%*%cbind(dose1,dose2)[2:3,]),simplify = F)

jagsdataRCS$prec <-  matrix(unlist(invS_i),40,2,byrow = T)
rcsplineDRmetaJAGSmodel <- jags.parallel(data = jagsdataRCS,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau1','tau2','newRR'),model.file = modelRCSplineDRmeta,
                                         n.chains=2,n.iter = 10000,n.burnin = 200,DIC=F,n.thin = 1)
# Error again! the values in invS_i are large!



# Until now I coudlnt handle this problem,
  # but discussing with orsini and see how they exactly compute Slist would help!
# I plan to include the prec matrix into makejags() function as you asked me after I have a clear vision how
 # it is computed. I dont want to mess up with functions and create for each appraoch new function!



























dose <- sim.data$all.dose[2:3,]
solve(t(dose)%*%solve(prOR[[1]])%*%dose)

## Greenland-longnecker covaiance-variance function as used in 'covariance' file in dosresmeta (It was difficult for me to follow the code!)

y <- sim.data$simulatedDRdata$logRR
v <- sim.data$simulatedDRdata$selogRR^2
cases <- sim.data$simulatedDRdata$cases
n <- sim.data$simulatedDRdata$cases+sim.data$simulatedDRdata$noncases
type <- sim.data$simulatedDRdata$type
sim.data$simulatedDRdata$n <- sim.data$simulatedDRdata$cases+sim.data$simulatedDRdata$noncases
grl(y=logRR,v=selogRR^2,cases=cases,n=n,type,data=subset(sim.data$simulatedDRdata,Study_No==1))
dd <- subset(sim.data$simulatedDRdata,Study_No==1)
1/dd$cases[1] + 1/dd$cases[2:3]

grl <- function(y, v, cases, n, type, data, tol = 1e-05){
  if (missing(data))
    data <- NULL
  if (is.null(data)) {
    data <- sys.frame(sys.parent())
  }
  else {
    if (!is.data.frame(data)) {
      data <- data.frame(data)
    }
  }
  mf <- match.call(expand.dots = FALSE)
  mf.y <- mf[[match("y", names(mf))]]
  mf.v <- mf[[match("v", names(mf))]]
  mf.cases <- mf[[match("cases", names(mf))]]
  mf.n <- mf[[match("n", names(mf))]]
  mf.type <- mf[[match("type", names(mf))]]
  y <- eval(mf.y, data, enclos = sys.frame(sys.parent()))
  v <- eval(mf.v, data, enclos = sys.frame(sys.parent()))
  v[is.na(v)] <- 0
  cases <- eval(mf.cases, data, enclos = sys.frame(sys.parent()))
  n <- eval(mf.n, data, enclos = sys.frame(sys.parent()))
  type <- eval(mf.type, data, enclos = sys.frame(sys.parent()))
  if (is.null(type))
    type <- as.vector(mf.type)
  Ax <- Axp <- cases
  repeat{
    A0 <- sum(cases) - sum(Ax[v!=0],na.rm = T) ## sum of cases0
    cx <- if (!as.character(type[1]) == "ir")
      1/Ax + 1/(n - Ax)                  #cx=1/Ax + 1/(n - Ax)

    else
      1/Ax
    e <- if (!as.character(type[1]) == "ir")
      y[v!=0] + log(A0) + log(n[v!=0]-Ax[v!=0]) - log(Ax[v!=0]) -
      log(n[v==0] - A0)
    else
      y[v!=0] + log(A0) + log(n[v!=0]) - log(Ax[v!=0]) - log(n[v==0])
    H <- diag(cx[v!=0] + cx[v==0], nrow = sum(v!=0))  ## diagonal matrix with (cases1+case0)
    H[upper.tri(H)] <- H[lower.tri(H)] <- cx[v==0]   ## off diagonal
    Axp[v==0] <- A0
    Axp[v!=0] <- Ax[v!=0] + solve(H) %*% e
    delta <- sum((Axp[v!=0] - Ax[v!=0])^2)
    if (delta < tol)
      break
    Ax <- Axp
  }
  cbind(A = Axp, N = n)
}
