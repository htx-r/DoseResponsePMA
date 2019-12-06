#******* jags model of linear dose-response model with normal likelihood for OR and RR

modelNorLinearDRmeta <- function(){
  b[1] <-0
  for (i in 1:ns) { ## for each study

    # multivariate normal likelihood of the response Y = logOR or logRR for the *reference* and *non-refernce* dose in a study i
    Y[i,1:(nd[i]-1)]  ~ dmnorm(mean[i,1:(nd[i]-1)], prec[(b[i]+1):(b[i]+nd[i]-1),1:(nd[i]-1)])

    # the linear dose-response model
    mean[i,1:(nd[i]-1)] <-  beta[i]*(X[i, 2:(nd[i])]-X[i, 1])

    # window to change the index in precision matrix
    b[i+1] <- b[i]+ nd[i]-1
  }

  # distribution of random effects

  for(i in 1:ns) {
    beta[i] ~dnorm(beta.pooled,prec.tau)
  }

  # prior distribution for heterogenity
   prec.tau<-1/variance
   variance<-tau*tau
   tau~dnorm(0,1)%_%T(0,)

   # prior distribution for the regression coeff beta
   beta.pooled ~ dnorm(0,0.001)
}
#


























# xi~dnorm(0,0.1)
# tau.eta~dgamma(5,5)
# tau <- abs(xi)/sqrt(tau.eta)


# log(tau) <- log.tau
# log.tau~ dunif(-20,20)
