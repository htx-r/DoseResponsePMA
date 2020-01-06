#******* jags model of spline dose-response model with binomial likelihood for OR

modelBinSplineDRmetaOR <- function(){
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
      delta[i,j] <-   beta1[i]*(X1[i,j]-X1[i,1]) + beta2[i]*(X2[i,j]-X2[i,1])
    }

  }

  # distribution of random effects
  for(i in 1:ns) {
    beta1[i]~dnorm(beta1.pooled,prec.beta)
    beta2[i]~dnorm(beta2.pooled,prec.beta)
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
#   p.drug3020 <- step(p.drug[30]-p.drug[20])
#   p.drug4030 <- step(p.drug[40]-p.drug[30])

}

#
















#   for (i in 1:np) { ## for each study
#     rr[i,1] ~ dbinom(p0[i],nn[i,1])
#     logit(p0[i]) <- z[i]
#     z[i] ~ dnorm(Z, prec.z)
#   }
# #
#   Z ~ dnorm(0, 0.001)
#   prec.z <- 1/v.z
#   v.z <- sigma.z * sigma.z
#   sigma.z ~ dnorm(0,1)%_%T(0,)
