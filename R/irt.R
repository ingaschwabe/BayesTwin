#==========================================================
# Analysis of twin data. 
# With IRT measurement model, without GE
# function irt.R for BayesTwin package
#==========================================================

irt <- function(data_mz, data_dz, n_burnin, n_iter){
    
    #Install packages when necessary (later niet meer nodig!!)
    #if(!require(MCMCpack)){ install.packages('MCMCpack'); require(MCMCpack)}
    #if(!require(R.utils)){ install.packages('R.utils'); require(R.utils)}
    #if(!require(MCMCpack)){ install.packages('MASS'); require(MASS)}    
    #if(!require(MCMCpack)){ install.packages('rjags'); require(rjags)}    
    
    # determine number of twin pairs
    n_mz <- nrow(data_mz) ; n_dz<- nrow(data_dz)
    
    # determine number of phenotypic variables
    n_var <- ncol(data_mz)/2
    
    
    #==========================================================
    # I. Write JAGS model file
    #==========================================================
    jags_model_irt <- "model{
      ##MZ twins
      for (fam in 1:n_mz){
        c_mz[fam] ~ dnorm(mu, tau_c)
        f_mz[fam] ~ dnorm(c_mz[fam], tau_a) 
          
          for (twin in 1:2){
            pheno_mz[fam,twin]~ dnorm(f_mz[fam],tau_e)   
          }
    
            #1pl model twin1
            for (k in 1:n_var){
              logit(p[fam,k]) <- pheno_mz[fam,1] - b[k]
              data_mz[fam,k] ~ dbern(p[fam,k])
            }      	
    
            #1pl model twin2
            for (k in (n_var+1):(2*n_var)){
              logit(p[fam,k]) <- pheno_mz[fam,2] - b[k-n_var]
              data_mz[fam,k] ~ dbern(p[fam,k])
            }
      }
    
      ##DZ twins
      for (fam in 1:n_dz){
        c_dz[fam] ~ dnorm(mu, tau_c)
        f0_dz[fam] ~ dnorm(c_dz[fam], doubletau_a)
    
          for (twin in 1:2){										
            f_dz[fam,twin] ~ dnorm(f0_dz[fam], doubletau_a)
            pheno_dz[fam,twin] ~ dnorm(f_dz[fam,twin], tau_e)
          }
    
            #1pl model twin1 (DZ)
            for (k in 1:n_var){
              logit(p2[fam,k]) <- pheno_dz[fam,1] - b[k]
              data_dz[fam,k] ~ dbern(p2[fam,k])
            }
    
            #1pl model twin2 (DZ)
            for (k in (n_var+1):(2*n_var)){
              logit(p2[fam,k]) <- pheno_dz[fam,2] - b[k-n_var]
              data_dz[fam,k] ~ dbern(p2[fam,k])
            }
    
        }
    
    
    #Priors
    mu <- 0 #to identify scale 
    doubletau_a <- 2*tau_a
    tau_a ~ dgamma(1,1)   
    tau_c ~ dgamma(1,.5)
    tau_e ~ dgamma(1,1)
    
    for (i in 1:n_var){
      b[i] ~ dnorm(0,.1)
    }
    }"

    jags_file_irt <- tempfile(fileext=".txt")
    write(jags_model_irt,jags_file_irt)
    
    #==========================================================
    # II. Run JAGS analysis
    #==========================================================
    inits = list(tau_a = 2, tau_c = 5, tau_e = 2)
    jags_data <- list(data_mz, data_dz, n_mz, n_dz, n_var)
    names(jags_data)<- c("data_mz", "data_dz", "n_mz", "n_dz", "n_var") 
    jags <- jags.model(jags_file_irt, jags_data, inits, n.chains = 1, quiet=FALSE)
    update(jags, n_burnin)
    out <- jags.samples(jags, c("tau_a", "tau_c", "tau_e", "b"), n_iter)
  
    
    #==========================================================
    # III. Organize results  
    #==========================================================
    
    samples_var_a = 1/out$tau_a[,,1]
    samples_var_c = 1/out$tau_c[,,1]
    samples_var_e = 1/out$tau_e[,,1]
    samples_b = out$b[,,1]
    
    #Calculate mean values: 
    var_a = mean(samples_var_a); var_c = mean(samples_var_c); var_e = mean(samples_var_e)
    b = apply(samples_b,1,mean)
    
    #Calculate SDs: 
    sd_var_a = sd(samples_var_a); sd_var_c = sd(samples_var_c); sd_var_e = sd(samples_var_e)
    sd_b = apply(samples_b,1, sd)
    
    #Calculate HPD: 
    hpd_var_a = HPD(samples_var_a, 0.95); hpd_var_c = HPD(samples_var_c, 0.95)
    hpd_var_e = HPD(samples_var_e, 0.95); hpd_b = apply(samples_b, 1,function(x) HPD(x,0.95))
    
    #Put results in a table 
    results = matrix(c(var_a, sd_var_a, var_c, sd_var_c, var_e, sd_var_e), 2, 3)
    colnames(results) <- c("varA","varC","varE")
    rownames(results) <- c("Posterior mean","Posterior standard deviation")
    results = as.table(results)  
    
    #And save output in a list: 
    output = list(results = results, samples_var_a = samples_var_a, samples_var_c = samples_var_c,
                  samples_var_e = samples_var_e, samples_b = samples_b,
                  hpd_var_a = hpd_var_a, hpd_var_c = hpd_var_c, hpd_var_e = hpd_var_e, hpd_b = hpd_b)
    
    #Change class of objects in order to use right plot method 
    class(output) <- c("ge_irt","list")
    class(output$samples_var_a) <- "samples"
    class(output$samples_var_c) <- "samples"
    class(output$samples_var_e) <- "samples"
    class(output$samples_b) <- "samples"  
    
    return(output)    
}


#Test function: 
'library(rjags)
source("simulate_twin_data.R")
data = simulate_twin_data(50, 20, n_items = 3)
data_mz = data$y_mz
data_dz = data$y_dz
n_mz = 50; n_dz = 20; n_items = 3; burnin = 200; n_iterations = 200
xx = irt(data_mz=data_mz, data_dz=data_dz, n_burnin=200, n_iter=200)
x = xx$samples_var_c
plot(x, type = "Sampling plot")'