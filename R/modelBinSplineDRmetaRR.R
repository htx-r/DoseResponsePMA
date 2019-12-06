#******* jags model of spline dose-response model with binomial likelihood for RR

modelBinSplineDRmetaRR <- function(){

  for (i in 1:ns) { ## for each study

    # binomial likelihood of number of events in the *reference* dose level in a study i
    r[i,1] ~ dbinom(p[i,1],n[i,1])

    # log parametrization of probabilities at each *reference* dose level: by that exp(beta)= RR
    log(p[i,1])<- u[i]

    for (j in 2:(nd[i])) { ## for each dose

      # binomial likelihood of number of events for the *non-refernce* dose in a study i
      r[i,j] ~ dbinom(p[i,j],n[i,j])

      # log parametrization of probabilities at each *non-refernce* dose level: by that exp(beta)= RR
      log(p[i,j])<- u[i] + delta[i,j]
      delta[i,j] <-  beta1[i]*(X1[i,j]-X1[i,1]) + beta2[i]*(X2[i,j]-X2[i,1])#+beta3[i]*(X3[i,j]-X3[i,1])

    }
  }

  # distribution of random effects
  for(i in 1:ns) {
    beta1[i]~dnorm(beta1.pooled,prec.beta)
    beta2[i]~dnorm(beta2.pooled,prec.beta)
    u[i]~dnorm(0,0.001)%_%T(,0)
  }

  # prior distribution for heterogenity
  prec.beta<-1/variance
  variance<-tau*tau
  tau~ dnorm(0,1)%_%T(0,)

  # prior distribution for both regression coeff beta1 and beta2
  beta1.pooled ~ dnorm(0,0.001)
  beta2.pooled ~ dnorm(0,0.001)
}

#

























# for (i in 1:new.n) {
#   newbeta1[i]~dnorm(beta1.pooled,prec.beta)
#   newbeta2[i]~dnorm(beta2.pooled,prec.beta)
#   newbeta3[i]~dnorm(beta3.pooled,prec.beta)
#
#   newY[i] <-newbeta1[i]*new.dose1[i]+newbeta2[i]*new.dose2[i]+newbeta3[i]*new.dose3[i]
#   newOR[i]<-exp(newY[i])
# }

