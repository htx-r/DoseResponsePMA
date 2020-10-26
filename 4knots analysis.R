
# 1. simulate
simulateDRmeta.funSigm4knots=function(ns=20,doserange=c(1, 10),samplesize=200,p0=0.1,OR=TRUE,splines=TRUE){ #

  require(Hmisc)

  # 1. generate the doses from uniform distribution within the doserange specified in the arguments above
  d<-cbind(rep(0,ns),matrix(round(c(runif(2*ns,doserange[1],doserange[2])),2),nrow=ns))##
  d<-t(apply(d,1,sort))
  dose <- c(t(d))


  # 2. Create the study-specific logOR and logRR either for splines or linear
    # find the dose cubic transformation
    knots<-unlist(round(quantile(d,c(0.05,0.33,0.41,0.59,0.77,0.95))))
    #knots <- c(0.5,1,3)
    trans.d<-rcspline.eval(c(t(d)),knots,inclx = T)
    dose1 <- c(trans.d[,1])
    dose2 <- c(trans.d[,2])
    dose3 <- c(trans.d[,3])
    dose4 <- c(trans.d[,4])
    dose5 <- c(trans.d[,5])
    #trans.d<-rcspline.eval(c(t(d)),knots,inclx = T)
    # dose2 <- log(dose)#ifelse(dose==0,0,log(dose))
    # implement the dose-response model
    # maxlogRR<- (beta1.pooled+2*tau)*max(trans.d[,1]) +(beta2.pooled+2*tau)*max(trans.d[,2]) # (only for RR) compute the maximum value of the logRR to choose p0 such that RR does not exceed 1
    # beta1.pooled<-c(sapply(rnorm(ns,beta1.pooled,sd=tau),rep,3)) # random effects of the linear coeff: generate the study-specific linear coeff from a normal distribution
    # beta2.pooled<-c(sapply(rnorm(ns,beta2.pooled,sd=tau),rep,3)) # random effects of the rcs-coeff: generate the study-specific rcs-coeff from a normal distribution
    # logrr <- beta1.pooled*trans.d[,1]+beta2.pooled*trans.d[,2]   # derive study-specific logOR by implementing the spline dose-response model
    #logrr <- ifelse(dose==0,0,log(log(dose)+1))
    logrr <- dose/(sqrt(1+dose^2))
  # 3. Generate the dose-specific logOR and logRR
   ## odds ratio (OR)
    # find the probabilities of the event to occur in zero dose;p0 and in the nonzero dose; p1
    odds0 <- p0/(1-p0)  # as p0 is provided as an argument we compute the odds in zero dose
    odds1 <- exp(logrr)*odds0 # from above we get the study-specific OR or RR
    p1<- odds1/(1+odds1) #compute p1 the probability of a case at any dose level

    # for the given probabilities above p0 and p1 and the uniformaly generated sample sizes, generate the dose-specific 'cases'
    uniquess<-round(runif(ns,samplesize-20,samplesize+20)) # sample size of each study arm
    ss<-c(sapply(uniquess,rep,3)) # sample size for all study arms
    cases<- matrix(rbinom(ns*3,ss,p1),nrow=3) # generate the cases for all dose levels
    noncases<-matrix(c(ss-cases),nrow=3) # calculate the noncases = n - cases
    #%%%!!!!! To avoid errors that occur due to having zero or n cases I replace them by 1 or n-1, respectively.
    cases[cases==0] <- 1
    noncases[noncases==0] <- 1

    # compute the estimated odds ratio (hatOR) and its standard error for each dose level
    hatodds<-cases/noncases # hat odds
    hatlogOR<-t(log(apply(hatodds,1,"/",hatodds[1,]))) #  log hat OR
    SEhatlogOR<-sqrt(1/cases[1,]+1/cases[c(2,3),]+1/(noncases[1,]) + 1/(noncases[c(2,3),])) #SE of log hat OR
    SEhatlogOR<-rbind(rep(NA,ns),SEhatlogOR) # add NA for SE in the zero dose1

    # the final results go to the the returned data
    hatlogrr <- hatlogOR
    SEhatlogrr <- SEhatlogOR
    type=rep('cc',3*ns) # type of the data is cc= case-control for OR. This is needed in dosresmeta function to choose how to compute the SE

    p1 <- ifelse(p1>1, 0.97,p1)
    # for the given probabilities above p0 and p1 and the uniformaly generated sample sizes, generate the dose-specific 'cases'
    uniquess<-round(runif(ns,samplesize-20,samplesize+20)) # sample size in study arm of zero dose
    ss<-c(sapply(uniquess,rep,3)) # sample size per study arm
    cases<-matrix(rbinom(ns*3,ss,p1),nrow = 3)     # events per study at zero dose
    noncases<-matrix(c(ss-cases),nrow = 3)     # events per study at zero dose
    cases[cases==0] <- 1
    noncases[noncases==0] <- 1
    # compute the estimated risk ratio (hatRR) and its standard error for each dose level
    hatRR <- cases[2:3,]/cases[1,]
    hatlogRR <- log(rbind(rep(1,ns),hatRR))
    SEhatlogRR<-sqrt(1/cases[2:3,]+1/cases[1,]-2/uniquess)
    selogRR<-c(rbind(NA,SEhatlogRR))

    # the final results go to the the returned data object
    hatlogrr <- hatlogRR
    SEhatlogrr <- selogRR
    type=rep('ci',3*ns) # type of the data is ci= cumulative incidance for RR. This is needed in dosresmeta function to choose how to compute the SE


  Study_No<-rep(1:ns,each=3) # label the studies in numbers

  # a data frame of the simulated data that we need to get from this function
  simulatedDRdata<-cbind.data.frame(Study_No=Study_No,logrr=c(hatlogrr),dose1=dose1,dose2=dose2,dose3=dose3,dose4=dose4,dose5=dose5,cases=c(cases),noncases=c(noncases),
                                    selogrr =c(SEhatlogrr), type=type)

  return(simulatedDRdata=simulatedDRdata)


}
# 2. make jags data
makejagsDRmeta4knots <- function(studyid, y,dose1,dose2,dose3,dose4,dose5,cases,noncases,se,type,data){
  # This function represent the dataset format so it can be input in all dose-response meta-analysis JAGS model
  # different settings are covered in this function: linear, spline, OR or RR

  # The arguments: all of them are vectors that are assumed to be elements in a dataframe ''data''
  # studyid: study id or number
  # y: a dose-specific outomce as log RR = log p_nonzero/p_zero or log OR
  # dose1: the doses in all studies
  # dose2: the spline transformed doses in all studies
  # cases: the dose-specific number of events
  # noncases: the dose-specific number of non events
  # se: standard error of y: log RR or log OR
  # type: (characterstic) either 'cc' when y= logOR or 'ci' when y= logRR
  # data: a dataframe that contains all the above arguments
  # splines: logical (T/F) to indicate whether we want a jags data for splines, or linear

  library(dosresmeta)

  # evaluate the arguments in the data
  data$studyid <-  eval(substitute(studyid), data)
  data$y <- eval(substitute(y), data)
  data$v <- eval(substitute(se), data)^2
  data$type <-eval(substitute(type), data)
  data$dose1 <- eval(substitute(dose1), data)
  data$dose2 <- eval(substitute(dose2), data)
  data$dose3 <- eval(substitute(dose3), data)
  data$dose4 <- eval(substitute(dose4), data)
  data$dose5 <- eval(substitute(dose5), data)

  data$cases <- eval(substitute(cases), data)
  data$noncases <- eval(substitute(noncases), data)
  data$n <- data$cases + data$noncases


  # additional needed quantities
  study_id <- unique(data$studyid)         ## a vector of the study id
  nd    <- as.numeric(table(data$studyid)) ## number of all doses (with zero dose)
  max.nd <- max(nd)                        ## maximum number of doses
  ns <- length(unique(data$studyid))       ## number of studies
  tncomp <- sum(as.numeric(table(data$studyid))-1) ## total number of non-zero comparisons

  ######################################################################  ######################################################################
  # convert the information from vectors format into matrices with study-specific values in columns and dose levels in rows
  ######################################################################  ######################################################################

  ## Matrix for the number of cases/events 'r' where each row refers to study and the columns refers to the dose levels.
  rmat <- matrix(NA,ns,max.nd)
  for (i in 1:ns) {
    rmat[i,1:as.numeric(table(data$studyid)[i])] <- data$cases[data$studyid == study_id[i]]
  }

  ## Matrix for the sample size 'n' where each row refers to study and the columns refers to the dose levels.

  nmat <- matrix(NA,ns,max.nd)
  for (i in 1:ns) {
    nmat[i,1:as.numeric(table(data$studyid)[i])] <- data$n[data$studyid == study_id[i]]
  }

  ## Matrix for the effects log RR or log OR where each row refers to study and the columns refers to the dose levels.
  Ymat <- matrix(NA,ns,max.nd)
  for (i in 1:ns) {
    Ymat[i,1:as.numeric(table(data$studyid)[i])] <- data$y[data$studyid == study_id[i]]
  }

  ## Matrix for the doses  where each row refers to a study and the columns refers to the dose levels.
  Xmat <- matrix(NA,ns,max.nd)
  for (i in 1:ns) {
    Xmat[i,1:as.numeric(table(data$studyid)[i])] <- data$dose1[data$studyid == study_id[i]]
  }

  ## Find the inverse of the variance-covariance matrix of the dose-specific effects y within each study;
  # so within each study we get a matrix (nd-1)x(nd-1) then we bind all the matrices into one big matrix so we can input it in jags model

  Slist <- sapply(unique(data$studyid), function(i) covar.logrr( cases = cases, n = cases+noncases, y=y,v=v,type = type,data = data[data$studyid==i,])
                  ,simplify = F)

  precmat <- matrix(NA,tncomp,max.nd-1)
  s <- matrix(NA, ns,max.nd-1)
  b <- no.d <-vector()
  index <- 1:tncomp
  for (i in 1:ns) {
    b[1] <- 0
    no.d[i] <- as.numeric(table(data$studyid)[i])-1
    precmat[(b[i]+1):(b[i]+no.d[i]),1:(no.d[i])] <- solve(Slist[[i]])
    s[i,1:no.d[i]] <- index[(b[i]+1):(b[i]+no.d[i])]
    b[i+1] <- b[i]+ no.d[i]
    precmat
  }

  # for bivariate model, I need these two matrices
  idmat <- diag(1,2)
  idmati <- matrix(c(0,1,1,0), nrow = 2,ncol = 2)

  #


    ## add the spline dose transformation to the list
    X2mat <- matrix(NA,ns,max.nd)
    for (i in 1:ns) {
      X2mat[i,1:as.numeric(table(data$studyid)[i])] <- data$dose2[data$studyid == study_id[i]]
    }

    X3mat <- matrix(NA,ns,max.nd)
    for (i in 1:ns) {
      X3mat[i,1:as.numeric(table(data$studyid)[i])] <- data$dose3[data$studyid == study_id[i]]
    }
    X4mat <- matrix(NA,ns,max.nd)
    for (i in 1:ns) {
      X4mat[i,1:as.numeric(table(data$studyid)[i])] <- data$dose4[data$studyid == study_id[i]]
    }

    X5mat <- matrix(NA,ns,max.nd)
    for (i in 1:ns) {
      X5mat[i,1:as.numeric(table(data$studyid)[i])] <- data$dose5[data$studyid == study_id[i]]
    }

    JAGSdata <- list(Y=Ymat[,-1],r=rmat,n=nmat,X1=Xmat,X2=X2mat,X3=X3mat,X4=X4mat,X5=X5mat,nd=nd,ns=ns,prec=precmat,s=s,idmat=idmat,idmati=idmati) # X3=X3mat[,-1], X3ref=X3mat[,1],

  return(JAGSdata)

}

# 3.jags model
#******* jags model of spline dose-response model with binomial likelihood for OR

modelBinSplineDRmetaOR4knots <- function(){
  for (i in 1:ns) { ## for each study
    # binomial likelihood of number of events in the *refernce* dose level in a study i
    r[i,1] ~ dbinom(p[i,1],n[i,1])

    # logit parametrization of probabilities at each *refernce* dose level: by that exp(beta)= OR
    logit(p[i,1]) <- u[i]

    for (j in 2:(nd[i])) { ## for each dose
      # binomial likelihood of number of events for the *non-refernce* dose in a study i
      r[i,j] ~ dbinom(p[i,j],n[i,j])

      # logit parametrization of probabilities at each *non-refernce* dose level: by that exp(beta)= OR
      logit(p[i,j]) <- u[i] + delta[i,j]
      delta[i,j] <-   beta1[i]*(X1[i,j]-X1[i,1]) + beta2[i]*(X2[i,j]-X2[i,1])+beta3[i]*(X3[i,j]-X3[i,1])+beta4[i]*(X4[i,j]-X4[i,1])+beta5[i]*(X5[i,j]-X5[i,1])
    }

  }

  # distribution of random effects
  for(i in 1:ns) {
    beta1[i]~dnorm(beta1.pooled,prec.beta)
    beta2[i]~dnorm(beta2.pooled,prec.beta)
    beta3[i]~dnorm(beta3.pooled,prec.beta)
    beta4[i]~dnorm(beta4.pooled,prec.beta)
    beta5[i]~dnorm(beta5.pooled,prec.beta)

    u[i]~dnorm(0,0.001)
  }

  # prior distribution for heterogenity
  prec.beta<-1/variance
  variance<-tau*tau
  tau~ dnorm(0,1)%_%T(0,)

  # log.tau~ dunif(0,100)#%_%T(0,)
  # tau <- exp(log.tau)

  # prior distribution for both regression coeff beta1 and beta2
  beta1.pooled ~ dnorm(0,0.001)
  beta2.pooled ~ dnorm(0,0.001)
  beta3.pooled ~ dnorm(0,0.001)
  beta4.pooled ~ dnorm(0,0.001)
  beta5.pooled ~ dnorm(0,0.001)

  # for(j in 1:nd.new){
  #   OR[j]<- exp(beta1.pooled*new.dose[j]+ beta2.pooled*f.newdose[j])
  # }

  # This part below is to obtain the absolute response over newdose range: 1 to 80, only for antidepressant not simulation

  #   for (i in 1:np) { ## for each study
  #     rr[i,1] ~ dbinom(p0[i],nn[i,1])
  #     logit(p0[i]) <- z[i]
  #     z[i] ~ dnorm(Z, prec.z)
  #   }
  #   # priors
  #   Z ~ dnorm(0, 0.001)
  #   prec.z <- 1/v.z
  #   v.z <- sigma.z * sigma.z
  #   sigma.z ~ dnorm(0,1)%_%T(0,)
  #
  #   for( j in 1:nd.new){
  #    OR[j] <- exp(beta1.pooled*new.dose[j]+ beta2.pooled*f.new.dose[j])
  #    odds.drug[j] <- OR[j]*exp(Z)
  #    p.drug[j] <- odds.drug[j]/(1+odds.drug[j])
  #
  # }
  # p.drug3020 <- step(p.drug[30]-p.drug[20])
  # p.drug4030 <- step(p.drug[40]-p.drug[30])

}

# 4. run jags for one simulation
OneSimulationSigm4knots <- function(ns=20,doserange=c(1, 10),samplesize=200){

  #** 1. simulate the data;
  #!!! I used this while loop to resimulate the data again if we got an error in the simulations.
  #!!! we did it becuase sometimes 'by chance' the function is failed to produce the spline transformation
  #!!! for specific doses espacially for narrow dose range.
  v <- 'try-error'
  while (v=='try-error') {
    sim.data <- try(simulateDRmeta.funSigm4knots(ns=ns,doserange = doserange,samplesize = samplesize),silent = TRUE)
    v<- class(sim.data)
  }

    # Freq: dosresmeta
    rcsplineDRmetaFreq <- dosresmeta(formula = logrr~dose1+dose2, id = Study_No,type=type,
                                     se = selogrr, cases = cases, n = cases+noncases, data = sim.data, proc='1stage',method = 'reml',covariance = 'gl')

    # Bayes Normal: jags
     jagsdata<- makejagsDRmeta4knots(Study_No,logrr,dose1,dose2,dose3,dose4,dose5,cases,noncases,se=selogrr,type=type,data=sim.data)

     rcsplineDRmetaJAGSmodel <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','tau'),model.file = modelNorSplineDRmeta,
                                             n.chains=3,n.iter = 100,n.burnin = 10,DIC=F,n.thin = 1)
    # Bayes Binomial: jags
      splineDRmetaJAGSmodelBin <- jags.parallel(data = jagsdata,inits=NULL,parameters.to.save = c('beta1.pooled','beta2.pooled','beta3.pooled','beta4.pooled','beta5.pooled','tau'),model.file = modelBinSplineDRmetaOR4knots,
                                                n.chains=3,n.iter = 10000,n.burnin =3000,DIC=F,n.thin = 1)


    #** 3s. spline results

    # beta1
    # mean
    f1 <-coef(rcsplineDRmetaFreq)[1]
    b1n <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta1.pooled
    b1b <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$beta1.pooled

    # standard error
    sdF1 <- sqrt(rcsplineDRmetaFreq$vcov[1,1])
    sdNor1 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta1.pooled','sd']
    sdBin1 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta1.pooled','sd']

    # measure to check the convergence: Rhat gelamn statistic
    RhatN1 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta1.pooled','Rhat']
    RhatB1 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta1.pooled','Rhat']

    # heterogenity tau
    # mean
    tf1 <- sqrt(rcsplineDRmetaFreq$Psi[1,1])
    tn <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$tau
    tb <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$tau

    # standard error
    sdF1 <- sqrt(rcsplineDRmetaFreq$vcov[1,1])
    sdNor1 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta1.pooled','sd']
    sdBin1 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta1.pooled','sd']


    # beta2
    # mean
    f2 <-coef(rcsplineDRmetaFreq)[2]
    b2n <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$beta2.pooled
    b2b <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$beta2.pooled

    # standard error
    sdF2 <- sqrt(rcsplineDRmetaFreq$vcov[2,2])
    sdBin2 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta2.pooled','sd']
    sdNor2 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta2.pooled','sd']

    # measure to check the convergence: Rhat gelamn statistic
    RhatN2 <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['beta2.pooled','Rhat']
    RhatB2 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta2.pooled','Rhat']


    # beta3
    b3b <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$beta3.pooled

    # standard error
    sdBin3 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta3.pooled','sd']

    # measure to check the convergence: Rhat gelamn statistic
    RhatB3 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta3.pooled','Rhat']

    # beta4
    b4b <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$beta4.pooled

    # standard error
    sdBin4 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta4.pooled','sd']

    # measure to check the convergence: Rhat gelamn statistic
    RhatB4 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta4.pooled','Rhat']

    # beta5
    b5b <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$beta5.pooled

    # standard error
    sdBin5 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta5.pooled','sd']

    # measure to check the convergence: Rhat gelamn statistic
    RhatB5 <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['beta5.pooled','Rhat']

    # commomn heterogenity tau in beta1, beta2 and beta3
    # mean
    tf <- sqrt(rcsplineDRmetaFreq$Psi[1,1])
    tn <- rcsplineDRmetaJAGSmodel$BUGSoutput$mean$tau
    tb <- splineDRmetaJAGSmodelBin$BUGSoutput$mean$tau

    # standard error
    sdBintau <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['tau','sd']
    sdNortau <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['tau','sd']

    # measure to check the convergence: Rhat gelamn statistic
    RhatNtau <- rcsplineDRmetaJAGSmodel$BUGSoutput$summary['tau','Rhat']
    RhatBtau <- splineDRmetaJAGSmodelBin$BUGSoutput$summary['tau','Rhat']

    # the return object is vector that combine all results: spline
    rval <- c(BayesB1=b1b,BayesN1=b1n,Freq1=unname(f1),sdF1=sdF1,sdNor1=sdNor1,sdBin1=sdBin1,tauN=tn,tauB=tb,tauF1=tf,RhatN1=RhatN1,RhatB1=RhatB1,
              BayesB2=b2b,BayesN2=b2n,Freq2=unname(f2),sdF2=sdF2,sdNor2=sdNor2,sdBin2=sdBin2,tauN=tn,tauB=tb,RhatN2=RhatN2,RhatB2=RhatB2,
              BayesB3=b3b,sdBin3=sdBin3,RhatB3=RhatB3,
              BayesB4=b4b,sdBin4=sdBin4,RhatB4=RhatB4,
              BayesB5=b5b,sdBin5=sdBin5,RhatB5=RhatB5,
              sdBintau=sdBintau,sdNortau=sdNortau,RhatBtau=RhatBtau,RhatNtau=RhatNtau)


  return(rval)

}

# 5. run jags for multiple simulations
nsim <- 100
spline4knots <- replicate(nsim,OneSimulationSigm4knots(ns=20,doserange = c(1,10),samplesize = 200),simplify = T)

rval <- colMeans(t(spline4knots))

#  our x : doses
d <- seq(0,10,l=1000)
knots<-unlist(round(quantile(d,c(0.05,0.33,0.41,0.59,0.77,0.95))))

new.dose<-rcspline.eval(d,knots,inclx = T)
new.dose1 <- c(new.dose[,1])
new.dose2 <- c(new.dose[,2])
new.dose3 <- c(new.dose[,3])
new.dose4 <- c(new.dose[,4])
new.dose5 <- c(new.dose[,5])

n.d <- length(d)




# our beta1
beta1 <- rval[c('BayesB1')]

# our beta2
beta2 <-rval[c('BayesB2')]

# our beta3
beta3 <-rval[c('BayesB3')]

# our beta3
beta4 <-rval[c('BayesB4')]

# our beta3
beta5 <-rval[c('BayesB5')]

y_bin <- exp(beta1*new.dose1+beta2*new.dose2+beta3*new.dose3+beta4*new.dose4+beta5*new.dose5)
y_true <-  exp(new.dose1/(sqrt(1+new.dose1^2)))
df11 <- data.frame(new.dose1=new.dose1,y1=y_bin,y2=y_true)

theme_set(
  theme_minimal() +
    theme(legend.position = "right",axis.text.x = element_text(size=12),axis.text.y = element_text(size=12))
)

# S1
ggplot(data = df11,aes(x=new.dose1)) +
  geom_line(size=1.3,aes(y=y1),col='steelblue')+
  geom_line(size=1.3,aes(y=y2),col='darkred')+
  xlab('')+
  ylab('')+
  # xlim(2.5,10)+
  ylim(0,3.5)+
  # scale_color_manual(values=c('lightblue2','steelblue','darkred','black'))+
  theme(panel.background = element_rect(fill = 'snow1',colour = 'white'),
        legend.position = c(0.72,0.90), legend.key.size = unit(3, "cm"), legend.text = element_text(size=12),
        legend.key.height  = unit(0.5,"cm"),legend.title = element_blank(), legend.background = element_blank(),
        axis.text.x = element_text(face='bold',size=14),axis.text.y = element_text(face='bold',size=14))

# area between curves
y_binF <- function(new.dose1){
  exp(beta1*new.dose1+beta2*new.dose2+beta3*new.dose3+beta4*new.dose4+beta5*new.dose5)
  }

integrate(y_binF,0,10)
area.between.curves(df11$new.dose1,df11$y1, df11$y2, xrange = range(df11$new.dose1))








