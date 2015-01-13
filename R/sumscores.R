#==========================================================
# Analysis of twin data (ACE model), using sum scores
# function irt_ge.R for BayesTwin package
# 
# This function writes a JAGS script for the model with 
# GxE interaction, draws from the posterior distribution
# and then organizes the output
#==========================================================

sumscores <- function(data_mz, data_dz, n_burnin, n_iter){
  
  #Install packages when necessary (later niet meer nodig!!)
  if(!require(MCMCpack)){ install.packages('MCMCpack'); require(MCMCpack)}
  if(!require(R.utils)){ install.packages('R.utils'); require(R.utils)}
  if(!require(MCMCpack)){ install.packages('MASS'); require(MASS)}    
  if(!require(MCMCpack)){ install.packages('rjags'); require(rjags)}    
  
  # Determine number of twin pairs
  n_mz <- nrow(data_mz); n_dz <- nrow(data_dz)
  
  #==========================================================
  # I. Write JAGS model file
  #==========================================================
  jags_model_sumscores <- "model{
    ##MZ twins
    for (fam in 1:n_mz){
        c_mz[fam] ~ dnorm(mu, tau_c)
        f_mz[fam] ~ dnorm(c_mz[fam], tau_a) 

        for (twin in 1:2){
            data_mz[fam,twin]~ dnorm(f_mz[fam],tau_e)   
        }
    }
  
    ##DZ twins
    for (fam in 1:n_dz){
        c_dz[fam] ~ dnorm(mu, tau_c)
        f0_dz[fam] ~ dnorm(c_dz[fam], doubletau_a)
  
        for (twin in 1:2){										
            f_dz[fam,twin] ~ dnorm(f0_dz[fam], doubletau_a)
            data_dz[fam,twin] ~ dnorm(f_dz[fam,twin], tau_e)
        }
    }
  
    doubletau_a <- 2*tau_a

    #Priors
    mu ~ dnorm(0,.1)
    beta0 ~ dnorm(-1,.5)
    beta1 ~ dnorm(0,.1)
    tau_a ~ dgamma(1,1)   
    tau_c ~ dgamma(1,.5)
    tau_e ~ dgamma(1,1)
    }"

  jags_file_sumscores <- tempfile(fileext=".txt")
  write(jags_model_sumscores,jags_file_sumscores)
  
  #==========================================================
  # II. Run JAGS analysis
  #==========================================================
  inits = list(tau_a = 2, tau_c = 5)
  jags_data <- list(data_mz, data_dz, n_mz, n_dz)
  names(jags_data)<- c("data_mz", "data_dz", "n_mz", "n_dz") 
  jags <- jags.model(jags_file_sumscores, jags_data, inits, n.chains = 1, quiet=FALSE)
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
'source("simulate_twin_data.R")
source("HPD.R")
source("plot.samples.R")
data = simulate_twin_data(50, 20, n_items = 0)
data_mz = data$y_mz
data_dz = data$y_dz
xx = sumscores(data_mz = data_mz, data_dz = data_dz, n_burnin = 200, n_iter = 200)
xx
plot(xx$samples_var_a)'
