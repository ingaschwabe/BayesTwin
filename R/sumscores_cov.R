#==========================================================
# Analysis of twin data (ACE model), using sum scores
# Inclusion of covariates is possible.  
# function sumscores_cov.R for BayesTwin package
# 
# This function writes a JAGS script for an ACE model 
# with covaristes. A prior probability distribution is 
# used for missing data in the covariates (Schwabe & van 
# den Berg, in preparation)
#==========================================================

sumscores_cov <- function(data_mz, data_dz, 
                          X_mz_twin1, X_mz_twin2, 
                          X_dz_twin2, X_dz_twin2,
                          data_mz_cov, data_dz_cov, 
                          n_burnin, n_iter){
  
  #Install packages when necessary (later niet meer nodig!!)
  if(!require(MCMCpack)){ install.packages('MCMCpack'); require(MCMCpack)}
  if(!require(R.utils)){ install.packages('R.utils'); require(R.utils)}
  if(!require(MCMCpack)){ install.packages('MASS'); require(MASS)}    
  if(!require(rjags)){ install.packages('rjags'); require(rjags)}    
  
  # Determine number of twin pairs
  n_mz <- nrow(data_mz); n_dz <- nrow(data_dz)
  
  # Determine number of environmental covariates: 
  N = ncol(X_mz_twin1)
  
  #==========================================================
  # I. Write JAGS model file
  #==========================================================
  jags_model_sumscores_cov <- "model{
    ##MZ twins
    for (fam in 1:n_mz){ #for each MZ family: 
        c_mz[fam] ~ dnorm(mu, tau_c)
        f_mz[fam] ~ dnorm(c_mz[fam], tau_a) 
  
        #Response variables:
        y_mz[fam,1] ~ dnorm(f_mz[fam] + inprod(X_mz_twin1[fam,], b), tau_e) 
        y_mz[fam,2] ~ dnorm(f_mz[fam] + inprod(X_mz_twin2[fam,], b), tau_e) 
    }
  
    ##DZ twins
    for (fam in 1:n_dz){ #for each DZ family: 
        c_dz[fam] ~ dnorm(mu, tau_c) 
        f1_dz[fam] ~ dnorm(c_dz[fam], doubletau_a) 
  
        for (twin in 1:2){
            f2_dz[fam,twin] ~ dnorm(f1_dz[fam], doubletau_a)
        }
  
        #Response variables: 
        y_dz[fam,1] ~ dnorm(f2_dz[fam,1] + inprod(X_dz_twin1[fam,], b), tau_e) 
        y_dz[fam,2] ~ dnorm(f2_dz[fam,2] + inprod(X_dz_twin2[fam,], b), tau_e)
    }
  
  
    #For DZ twins:
    doubletau_a <- 2*tau_a    
  
    #Priors:
    tau_a ~ dgamma(1,1) 
    tau_e ~ dgamma(1,1)
    tau_c ~ dgamma(1,.5)
     mu ~ dnorm(0, .1)
    
    #Prior for regression coefficients: 
    #Multivariate normal 
    b[1:N] ~ dmnorm(mu_b[1:N], tau_b[,]) 
    
    ##Latent variable psi for MZ and DZ twins
    for(fam in 1:n_mz){
        psi_mz[fam, 1:N] ~ dmnorm(gamma_cov[1:N], tau_b_mz) 
    }

    #for(fam in 1:n_dz){
    #    psi_dz[fam, 1:N] ~ dmnorm(gamma_cov[1:N], tau_b_dz) 
    #}

    #We use independent standard normal distributions as priors for 
    #the expected values of the covariates
    for (i in 1:N){
        gamma_cov[i] ~ dnorm(0, .1) 
    }

    #As a prior distribution for tau b, we use a 
    #wishart distribution. 
    tau_b_mz[1:N, 1:N] ~ dwish(omgea_tau_b[,], N)  
    tau_b_dz[1:N, 1:N] ~ dwish(omega_tau_b[,], N)

    #Final prior distribution for the covariates: 
    for (fam in 1:n_mz){
        exp_mz_twin1[fam,1:N] ~ dmnorm(psi_mz[fam,], tau_w_mz)
        exp_mz_twin2[fam,1:N] ~ dmnorm(psi_dz[fam,], tau_w_mz)
    }
  
    for (fam in 1:n_dz){
        exp_dz_twin1[fam,1:N] ~ dmnorm(psi_dz[fam,], tau_w_dz)
        exp_dz_twin2[fam,1:N] ~ dmnorm(psi_dz[fam,], tau_w_dz)
    }
  
    #As a prior distribution for tau w, we use a 
    #wishart distribution. 
    tau_w_mz[1:N, 1:N] ~ dwish(omgea_tau_w[,], N)  
    tau_w_dz[1:N, 1:N] ~ dwish(omega_tau_w[,], N)
  
    #As JAGS does not allow to define anything other 
    #than probability distributions for observed data,
    #we place a normally distributed prior distribution 
    #on the covariates with a very samll residual 
    #Variance of 0.01
  
    for (fam in 1:n_mz){
        for (cov in 1:N){
            X_mz_twin1[fam,cov] ~ dnorm(exp_mz_twin1[fam,cov], 100)
            X_mz_twin2[fam,cov] ~ dnorm(exp_mz_twin2[fam,cov], 100)
        }
    }
    
    for (fam in 1:n_dz){
        for (cov in 1:N){
            X_dz_twin1[fam,cov] ~ dnorm(exp_dz_twin1[fam,cov], 100)
            X_dz_twin2[fam,cov] ~ dnorm(exp_dz_twin2[fam,cov], 100)
        }
    }
  }"

  jags_file_sumscores_cov <- tempfile(fileext=".txt")
  write(jags_model_sumscores_cov, jags_file_sumscores_cov)
  
  #==========================================================
  # II. Run JAGS analysis
  #==========================================================
  inits = list(tau_a = 2, tau_c = 5)
  jags_data <- list("y_mz" = data_mz, "y_dz" = data_dz, 
                    "X_mz_twin1" = X_mz_twin1, "X_mz_twin2" = X_mz_twin2, 
                    "X_dz_twin1" = X_dz_twin1, "X_dz_twin2" = X_dz_twin2,
                    "n_mz" = n_mz, "n_dz" = n_dz, "N" = N,
                    "mu_b" = rep(0,N), "tau_b" = diag(1,N),
                    "omega_tau_b" = diag(1,N), "omega_tau_w" = diag(1,N))
  jags <- jags.model(jags_file_sumscores_cov, jags_data, inits, n.chains = 1, quiet=FALSE)
  update(jags, n_burnin)
  out <- jags.samples(jags, c("tau_a", "tau_c", "tau_e"), n_iter)
  
  
  #==========================================================
  # III. Organize results  
  #==========================================================
  samples_var_a = 1/out$tau_a[,,1]; samples_var_c = 1/out$tau_c[,,1]
  samples_var_e = 1/out$tau_e[,,1]
  
  #Calculate mean values: 
  var_a = mean(samples_var_a); var_c = mean(samples_var_c); var_e = mean(samples_var_e)
  
  #Calculate SDs: 
  sd_var_a = sd(samples_var_a); sd_var_c = sd(samples_var_c); sd_var_e = sd(samples_var_e)
  
  #Calculate HPD: 
  hpd_var_a = HPD(samples_var_a, 0.95); hpd_var_c = HPD(samples_var_c, 0.95)
  hpd_var_e = HPD(samples_var_e, 0.95)
  
  #Put results in a table 
  results = matrix(c(var_a, sd_var_a, var_c, sd_var_c, var_e, sd_var_e), 2, 3)
  colnames(results) <- c("varA","varC","varE")
  rownames(results) <- c("Posterior mean","Posterior standard deviation")
  results = as.table(results) 
  
  #And save output in a list: 
  output = list(results = results, samples_var_a = samples_var_a, samples_var_c = samples_var_c, 
                samples_var_e = samples_var_e,
                hpd_var_a = hpd_var_a, hpd_var_c = hpd_var_c, hpd_var_e = hpd_var_e)
  
  #Change class of objects in order to use right plot method 
  class(output) <- c("sumscores","list")    
  class(output$samples_var_a) <- "samples"
  class(output$samples_var_c) <- "samples"
  class(output$samples_var_e) <- "samples"
  
  return(output)
}


# #Test function: 
# library(rjags)
source("simulate_twin_data.R")
source("HPD.R")
source("plot.samples.R")
data = simulate_twin_data(50, 20, n_items = 0, n_var = 4)
data_mz = data$y_mz
data_dz = data$y_dz
X_mz_twin1 = data$cov_mz[,1:3]
X_mz_twin2 = data$cov_mz[,4:6]

X_dz_twin1 = data$cov_dz[,1:3]
X_dz_twin2 = data$cov_dz[,4:6]

X_dz_twin1[]

xx = sumscores_cov(data_mz = data_mz, data_dz = data_dz, 
                   X_mz_twin1 = X_mz_twin1, X_mz_twin2 = X_mz_twin2, 
                   X_dz_twin1 = X_dz_twin1, X_dz_twin2 = X_dz_twin2,
                   n_burnin = 200, n_iter = 200)
plot(xx$samples_var_a)