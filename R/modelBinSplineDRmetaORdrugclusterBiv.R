#******* jags model of spline dose-response model with binomial likelihood for OR

modelBinSplineDRmetaORdrugclusterBiv <- function(){
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
      delta[i,j] <-   beta_c[i,drug[i],1]*(X1[i,j]-X1[i,1]) + beta_c[i,drug[i],2]*(X2[i,j]-X2[i,1])
    }

  }

  # distribution of random effects
   # within clusters
  for(i in 1:ns) {
    beta_c[i,drug[i],1:2]~dmnorm(b_c[drug[i],1:2],inv.det*(tau.sq.with*idmat - rho.with*idmati))
  }

  # between clusters
  for (c in c(1:4,6)) {
    b_c[c,1:2] ~dmnorm(b[1:2], inv.det_c*(tau.sq.betw*idmat - rho.betw*idmati))
  }

  # minimal dose effect
  for (i in 1:ns) {
    u[i]~dnorm(0,0.001)
  }

  # prior distribution to heterogenity within clusters
  prec.beta.with<-1/tau.sq.with
  tau.sq.with<-tau.with*tau.with
  tau.with~ dnorm(0,1)%_%T(0,)
  inv.det <- 1/(tau.sq.with^2 + rho.with^2) # inverse of the determinant of the matrix

  # prior to covariances within clusters
  rho.with <- cov.with/tau.sq.with
  cov.with ~ dunif(-10,10)

  # prior distributions to b1 and b2
  b[1] <- b1
  b[2] <- b2
  b1 ~ dnorm(0,0.001)
  b2 ~ dnorm(0,0.001)

  # prior distribution to heterogenity between clusters
  prec.beta.betw<-1/tau.sq.betw
  tau.sq.betw<-tau.betw*tau.betw
  tau.betw~ dnorm(0,1)%_%T(0,)
  inv.det_c <- 1/(tau.sq.betw^2 + rho.betw^2) # inverse of the determinant of the matrix

  # prior to covariances between clusters
  rho.betw <- cov.betw/tau.sq.betw
  cov.betw ~ dunif(-10,10)

}
