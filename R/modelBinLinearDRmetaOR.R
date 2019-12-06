#******* jags model of linear dose-response model with binomial likelihood for OR

modelBinLinearDRmetaOR <- function(){

  for (i in 1:ns) { ## for each study

    # binomial likelihood of number of events in the *refernce* dose level in a study i
    r[i,1] ~ dbinom(p[i,1],n[i,1])
    # logit parametrization of probabilities at each *refernce* dose level: by that exp(beta)= OR
    logit(p[i,1])<- u[i]

    for (j in 2:(nd[i])) { ## for each dose

      # binomial likelihood of number of events for the *non-refernce* dose in a study i
      r[i,j] ~ dbinom(p[i,j],n[i,j])

      # logit parametrization of probabilities at each *non-refernce* dose level: by that exp(beta)= OR
      logit(p[i,j])<- u[i] + delta[i,j]
      delta[i,j] <-  beta[i]*(X[i,(j)]-X[i,1])
    }
  }

  # distribution of random effects
  for(i in 1:ns) {
    beta[i] ~dnorm(beta.pooled,prec.tau)
    u[i]~dnorm(0,0.001)
  }

  # prior distribution for heterogenity
  prec.tau<-1/variance
  variance<-tau*tau
  tau~dnorm(0,1)%_%T(0,)

  # prior distribution for the regression coeff beta
  beta.pooled ~ dnorm(0,0.001)
}






























# xi~dnorm(0,0.1)
# tau.eta~dgamma(5,5)
# tau <- abs(xi)/sqrt(tau.eta)
# log(tau) <- log.tau
# log.tau~ dunif(-20,20)
