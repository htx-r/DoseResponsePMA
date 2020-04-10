#******* jags model of spline dose-response model with binomial likelihood for OR

modelBinSplineDRmetaORdrugcluster <- function(){
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
      delta[i,j] <-   beta1[i,drug[i]]*(X1[i,j]-X1[i,1]) + beta2[i,drug[i]]*(X2[i,j]-X2[i,1])
    }

  }

  # distribution of random effects
  for(i in 1:ns) {
    beta1[i,drug[i]]~dnorm(beta1D[drug[i]],prec.beta)
    beta2[i,drug[i]]~dnorm(beta2D[drug[i]],prec.beta)
  }

  for (k in c(1:4,6)) {
    beta1D[k] ~dnorm(b1, prec.betaD)
    beta2D[k] ~dnorm(b2, prec.betaD)
  }

for (i in 1:ns) {
  u[i]~dnorm(0,0.001)
}
  # prior distribution for heterogenity
  prec.beta<-1/variance
  variance<-tau*tau
  tau~ dnorm(0,1)%_%T(0,)


  # prior distribution for both regression coeff beta1 and beta2
  b1 ~ dnorm(0,0.001)
  b2 ~ dnorm(0,0.001)

  # prior distribution for heterogenity
  prec.betaD<-1/varianceD
  varianceD<-sigma*sigma
  sigma~ dnorm(0,1)%_%T(0,)

}
