  model {

    for (i in 1:n.obs){
      # We assume our ecological response(s) follow a normal distribution
      y[i] ~ dnorm(mu[i], w)
      
      # Here, year is treated separately from the climate covariates because
      # we always want to include year in the model because we are working
      # with time series data. Similarly, if sampling data were available for
      # our case studies, we would include an effect for sampling effort with
      # its own coefficient to ensure it was always included in the model
      mu[i] <- beta0 + inprod(gamma[]*beta[], covar[i,]) + beta.year*years[i]
    } #betas come from the horseshoe prior

    for (k in 1:n.covars){
      # Priors for our climate covariates
      beta[k] ~ dnorm(0, lambda * v.fix[k] * vs.fix)

      # We use a half-Cauchy prior for v[k]
      xi[k] ~ dnorm(0, 1)
      tau.eta[k] ~ dgamma(0.5, 0.5)
      v[k] <- abs(xi[k]) / sqrt(tau.eta[k])
      v.fix[k] <- ifelse(v[k] < 0.0001, 0.0001, v[k])
      gamma[k] ~ dbern(0.5)
    }

    # Half-Cauchy prior for vs
    xis ~ dnorm(0, 1)
    tau.etas ~ dgamma(0.5, 0.5) # chi^2 with 1 d.f.
    vs <- abs(xis) / sqrt(tau.etas)
    vs.fix <- ifelse(vs < 0.001, 0.001, vs)
    
    # Global priors
    w ~ dgamma(0.01, 0.01) # gamma prior for error term
    beta0 ~ dnorm(0, 0.001) # normal prior for intercept term
    beta.year ~ dnorm(0, 0.001) # normal prior for the effect of year

  }

  
  
  