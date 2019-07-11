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
data$studyid <- data$Study_No

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

# Compare the results of the three approaches
prRR[[1]]
prOR[[1]]
rcsplineDRmetaFreq$Slist[[1]]

# I think Slist refers to C_i but S_i = Var(theta_hat) = (X'C_i^-1 X)^-1
# For example, for the first study is computed as follow
 dose <- sim.data$all.dose[2:3,]
solve(t(dose)%*%solve(prOR[[1]])%*%dose)
##













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
