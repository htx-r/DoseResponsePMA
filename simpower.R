source('oneSimulation.R')
simpower <- function(nrep=3,beta1.pooled=0.02,beta2.pooled=NULL,tau=0.001,ns=20,doserange=c(1, 10),samplesize=200,OR=FALSE,splines = FALSE){

  # Repeat the simulation nrep times
  res <- replicate(nrep,OneSimulation(beta1.pooled=beta1.pooled,beta2.pooled = beta2.pooled,tau=tau,ns=ns,doserange = doserange,samplesize = samplesize,OR=OR,splines = splines),simplify = T)
if(splines==FALSE){
  # Biases
  biasBbin <- colMeans(t(res))['BayesB']-beta1.pooled
  biasBnor <- colMeans(t(res))['BayesN']-beta1.pooled
  biasF <- colMeans(t(res))['Freq']- beta1.pooled

  # heterogenity
  tauN.hat <- colMeans(t(res))['tauN']
  tauB.hat <- colMeans(t(res))['tauB']

  # MSE: Mean square error
  mseBbin <- mean((t(res)[,'BayesB']-beta1.pooled)^2)
  mseBnor <- mean((t(res)[,'BayesN']-beta1.pooled)^2)
  mseF <- mean((t(res)[,'Freq']-beta1.pooled)^2)

  # Type 1 error
  alphaNor <- ifelse(beta1.pooled==0,mean(beta1.pooled==0&res['sig.testNor',]==1),NA)
  alphaBin <- ifelse(beta1.pooled==0,mean(beta1.pooled==0&res['sig.testBin',]==1),NA)
  alphaF <- ifelse(beta1.pooled==0,mean(beta1.pooled==0&res['sig.testF',]==1),NA)

  # Type 2 error
  betaNor <- ifelse(beta1.pooled!=0,mean(beta1.pooled!=0 & res['sig.testNor',]==0),NA)
  betaBin <- ifelse(beta1.pooled!=0,mean(beta1.pooled!=0 & res['sig.testBin',]==0),NA)
  betaF <- ifelse(beta1.pooled!=0,mean(beta1.pooled!=0 & res['sig.testF',]==0),NA)

  ## Monte Carlo SE of estimate
  MCseBin <- sqrt(sum((t(res)[,'BayesB']-colMeans(t(res))['BayesB'])^2)/(ncol(res)*(ncol(res)-1)))
  MCseNor <- sqrt(sum((t(res)[,'BayesN']-colMeans(t(res))['BayesN'])^2)/(ncol(res)*(ncol(res)-1)))
  MCseF <- sqrt(sum((t(res)[,'Freq']-colMeans(t(res))['Freq'])^2)/(ncol(res)*(ncol(res)-1)))

  rval <- cbind(beta.pooled=beta1.pooled,tau=tau, # true values
                   tauN.hat=tauN.hat,tauB.hat=tauB.hat, # estimation of tau
                   biasBnor=biasBnor,biasBbin=biasBbin,biasF=biasF, # bias of beta
                   mseBnor=mseBnor,mseBbin=mseBbin,mseF=mseF, # mean squared error for beta
                   alphaBin=alphaBin,alphaNor=alphaNor,alphaF=alphaF, # type 1 error (alpha)
                   betaBin=betaBin,betaNor=betaNor,betaF=betaF, # type 2 error (beta)
                   MCseBin=MCseBin,MCseNor=MCseNor,MCseF=MCseF) # monte carlo standard error
  row.names(rval) <- NULL
}else{
  biasBbin1 <- colMeans(t(res))['BayesB1']-beta1.pooled
  biasBnor1 <- colMeans(t(res))['BayesN1']-beta1.pooled
  biasF1 <- colMeans(t(res))['Freq1']- beta1.pooled

  # heterogenity
  tauN.hat1 <- colMeans(t(res))['tauN1']
  tauB.hat1 <- colMeans(t(res))['tauB1']

  # MSE: Mean square error
  mseBbin1 <- mean((t(res)[,'BayesB1']-beta1.pooled)^2)
  mseBnor1 <- mean((t(res)[,'BayesN1']-beta1.pooled)^2)
  mseF1 <- mean((t(res)[,'Freq1']-beta1.pooled)^2)

  # Type 1 error
  alphaNor1 <- ifelse(beta1.pooled==0,mean(beta1.pooled==0&res['sig.testNor1',]==1),NA)
  alphaBin1 <- ifelse(beta1.pooled==0,mean(beta1.pooled==0&res['sig.testBin1',]==1),NA)
  alphaF1 <- ifelse(beta1.pooled==0,mean(beta1.pooled==0&res['sig.testF1',]==1),NA)

  # Type 2 error
  betaNor1 <- ifelse(beta1.pooled!=0,mean(beta1.pooled!=0 & res['sig.testNor1',]==0),NA)
  betaBin1 <- ifelse(beta1.pooled!=0,mean(beta1.pooled!=0 & res['sig.testBin1',]==0),NA)
  betaF1 <- ifelse(beta1.pooled!=0,mean(beta1.pooled!=0 & res['sig.testF1',]==0),NA)

  ## Monte Carlo SE of estimate
  MCseBin1 <- sqrt(sum((t(res)[,'BayesB1']-colMeans(t(res))['BayesB1'])^2)/(ncol(res)*(ncol(res)-1)))
  MCseNor1 <- sqrt(sum((t(res)[,'BayesN1']-colMeans(t(res))['BayesN1'])^2)/(ncol(res)*(ncol(res)-1)))
  MCseF1 <- sqrt(sum((t(res)[,'Freq1']-colMeans(t(res))['Freq1'])^2)/(ncol(res)*(ncol(res)-1)))

  #### #### #### ####
  #### beta2.pooled
  #### #### #### ####
  # Biases
  biasBbin2 <- colMeans(t(res))['BayesB2']-beta2.pooled
  biasBnor2 <- colMeans(t(res))['BayesN2']-beta2.pooled
  biasF2 <- colMeans(t(res))['Freq2']- beta2.pooled

  # heterogenity
  tauN.hat2 <- colMeans(t(res))['tauN2']
  tauB.hat2 <- colMeans(t(res))['tauB2']

  # MSE: Mean square error
  mseBbin2 <- mean((t(res)[,'BayesB2']-beta2.pooled)^2)
  mseBnor2 <- mean((t(res)[,'BayesN2']-beta2.pooled)^2)
  mseF2 <- mean((t(res)[,'Freq2']-beta2.pooled)^2)

  # Type 1 error
  alphaNor2 <- ifelse(beta2.pooled==0,mean(beta2.pooled==0&res['sig.testNor2',]==1),NA)
  alphaBin2 <- ifelse(beta2.pooled==0,mean(beta2.pooled==0&res['sig.testBin2',]==1),NA)
  alphaF2 <- ifelse(beta2.pooled==0,mean(beta2.pooled==0&res['sig.testF2',]==1),NA)

  # Type 2 error
  betaNor2 <- ifelse(beta2.pooled!=0,mean(beta2.pooled!=0 & res['sig.testNor2',]==0),NA)
  betaBin2 <- ifelse(beta2.pooled!=0,mean(beta2.pooled!=0 & res['sig.testBin2',]==0),NA)
  betaF2 <- ifelse(beta2.pooled!=0,mean(beta2.pooled!=0 & res['sig.testF2',]==0),NA)

  ## Monte Carlo SE of estimate
  MCseBin2 <- sqrt(sum((t(res)[,'BayesB2']-colMeans(t(res))['BayesB2'])^2)/(ncol(res)*(ncol(res)-1)))
  MCseNor2 <- sqrt(sum((t(res)[,'BayesN2']-colMeans(t(res))['BayesN2'])^2)/(ncol(res)*(ncol(res)-1)))
  MCseF2 <- sqrt(sum((t(res)[,'Freq2']-colMeans(t(res))['Freq2'])^2)/(ncol(res)*(ncol(res)-1)))


  # Return
  rval <- c( # beta1.pooled
    beta1.pooled=beta1.pooled,tau=tau, # true values
    tauN.hat1=tauN.hat1,tauB.hat1=tauB.hat1, # estimation of tau
    biasBnor1=biasBnor1,biasBbin1=biasBbin1,biasF1=biasF1, # bias of beta
    mseBnor1=mseBnor1,mseBbin1=mseBbin1,mseF1=mseF1, # mean squared error for beta
    alphaBin1=alphaBin1,alphaNor1=alphaNor1,alphaF1=alphaF1, # type 1 error (alpha)
    betaBin1=betaBin1,betaNor1=betaNor1,betaF1=betaF1, # type 2 error (beta)
    MCseBin1=MCseBin1,MCseNor1=MCseNor1,MCseF1=MCseF1,
    # beta2.pooled
    beta2.pooled=beta2.pooled,tau=tau, # true values
    tauN.hat2=tauN.hat2,tauB.hat2=tauB.hat2, # estimation of tau
    biasBnor2=biasBnor2,biasBbin2=biasBbin2,biasF2=biasF2, # bias of beta
    mseBnor2=mseBnor2,mseBbin2=mseBbin2,mseF2=mseF2, # mean squared error for beta
    alphaBin2=alphaBin2,alphaNor2=alphaNor2,alphaF2=alphaF2, # type 1 error (alpha)
    betaBin2=betaBin2,betaNor2=betaNor2,betaF2=betaF2, # type 2 error (beta)
    MCseBin2=MCseBin2,MCseNor2=MCseNor2,MCseF2=MCseF2) # monte carlo standard error
  row.names(rval) <- NULL
}
  return(rval)
}
