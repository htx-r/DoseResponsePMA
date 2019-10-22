#******* Linear dose-response jags model with normal likelihood for logRR=Y

modelNorLinearDRmeta <- function(){
  b[1] <-0
  for (i in 1:ns) { ## for each study

    # Likelihood

    Y[i,1:(nd[i])]  ~ dmnorm(mean[i,1:(nd[i])], prec[(b[i]+1):(b[i]+nd[i]),1:(nd[i])])

    mean[i,1:(nd[i])] <-  beta[i]*(X[i, 2:(nd[i]+1)]-X[i, 1])

    b[i+1] <- b[i]+ nd[i]
  }

  # Random effect

  for(i in 1:ns) {
    # beta[i] <- xi*eta[i]
    # eta[i]~dnorm(beta.pooled,tau.eta)
    beta[i] ~dnorm(beta.pooled,prec.tau)

  }

  # Priors
   prec.tau<-1/variance
   variance<-tau*tau
   tau~dunif(0,1)

  # xi~dnorm(0,0.1)
  # tau.eta~dgamma(5,5)
  # tau <- abs(xi)/sqrt(tau.eta)

  beta.pooled ~ dnorm(0,0.1)

  # log(tau) <- log.tau
  # log.tau~ dunif(-20,20)

  ## Predictions

 #   for (i in 1:new.n) {
 #  newbeta[i]~dnorm(beta.pooled,prec.tau)
 # newY[i] <-newbeta[i]*new.dose[i]
 # newRR[i]<-exp(newY[i])
 #   }
}

