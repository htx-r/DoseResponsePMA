logisticMA <- function(){
  for (i in 1:ns) {

    # placebo effect
rp[i]~dbinom(pp[i], np[i])
    logit(pp[i]) <- BR1[i]

    # drug effect
rd[i] ~dbinom(pd[i], nd[i])
logit(pd[i]) <- BR1[i]+delta[i]
delta[i] <- a+b*(BR2[i]-BR2.bar) ## fixed effect not random
  }

  # prior

  for (i in 1:ns) {
    BR1[i]~dnorm(0,1)
  }
   a ~dnorm(0,1)
   b ~dnorm(0,1)

}

simdata <- function(pp=0.2,pd=0.6,ns=20){
  np <-nd <-  round(runif(ns,180,220))
  rp <- rbinom(ns,np, pp)
  rd <- rbinom(ns,nd,pd)
sim.dat <- data.frame(studyid=1:ns,np=np,rp=rp,nd=nd,rd=rd)
return(sim.dat)
}
simdata <- simdata()
library(meta)

(0.6*0.8)/(0.4*0.2)

m <- metabin(rd,np,rp,np,sm="OR",data=simdata,method = 'GLMM')
metareg(m,BR)

createORreference.fun=function(rp,np,rd,nd)
{


    calculate=metabin(rd,nd,rp,np,sm="OR")
    logOR=c(calculate$TE)
    selogOR=c(calculate$seTE)


  return(cbind(logOR=logOR,selogOR=selogOR))
}
pp <- simdata$rp/simdata$np
BR <- pp/(1-pp)
simdata$BR <- log(BR)
 logORmat <- sapply(unique(simdata$studyid),function(i) createORreference.fun(simdata$rp[simdata$studyid==i],simdata$np[simdata$studyid==i],simdata$rd[simdata$studyid==i],simdata$nd[simdata$studyid==i]),simplify = FALSE)
logORmat <- do.call(rbind,logORmat)
simdata$logOR <- c(logORmat[,1])
simdata$selogOR <- c(logORmat[,2])







