#==========================================================
# Master function that calls subroutines
# BayesTwin package - Bayesian Analysis of Twin Data 
#==========================================================
#setwd("C:/Users/schwabei/Dropbox/International student performance_IngaStephanie/R/BayesTwin3/R")
#setwd("/Users/stephanievandenberg/Dropbox/International student performance_IngaStephanie/R/BayesTwin3/R")
#setwd("C:/Users/inga/Dropbox/International student performance_IngaStephanie/R/BayesTwin3/R")
#setwd("/Users/inga/Dropbox/International student performance_IngaStephanie/R/BayesTwin3/R")

twin_mcmc = function(data_mz, data_dz, twin1_datacols_p, twin2_datacols_p,
                         twin1_datacols_cov = NA, twin2_datacols_cov = NA,
                         #multivariate = F, ordinal = F, model = "ACE", common = T,
                         n_iter = 10000, n_burnin = 5000, ge = FALSE){
    
    ### At this moment: 
    ### ACE model, with possible inclusion of environmental covariates (missing data possible) and 
    ### estimation of genotype by environment interaction. Sumscores as well as item scores, measurement
    ### model is item response theory (IRT) model, Rasch + GPCM for ordinal data. 
    ###
    ### Coming in the future: 
    ### Multivariate twin analysis (more than 1 phenotype) - independent and common pathway models. 
    ### Estimation of ADE model: modelling dominance effects instead of common-shared environmental variance.

    #Install packages (for debugging, not necessary anymore when R package is fully working)
    if(!require(MCMCpack)){ install.packages('MCMCpack'); require(MCMCpack)}
    if(!require(R.utils)){ install.packages('R.utils'); require(R.utils)}
    if(!require(MCMCpack)){ install.packages('MASS'); require(MASS)}    
    if(!require(MCMCpack)){ install.packages('rjags'); require(rjags)}      
                                    
    #Load HPD function: 
    source("HPD.R")                                                         
                                                                    
    #Calculate number of MZ and DZ twins: 
    n_mz = nrow(data_mz)
    n_dz = nrow(data_dz)
    
    #Make sure that first columns are only twin 1, next the twin 2 columns
    data_mz = data_mz[,c(twin1_datacols_p, twin2_datacols_p)]
    data_dz = data_dz[,c(twin1_datacols_p, twin2_datacols_p)]
    
    #Test if its item data: 
    itemdata = NA
    if(length(twin1_datacols) == 1){
        itemdata = FALSE
    } else {itemdata = TRUE}
    
    #Test if covariates are used
    covariates = NA
    if(is.na(twin1_datacols_cov) == FALSE){
        X_mz_twin1 = data_mz[,twin1_datacols_cov]
        X_mz_twin2 = data_mz[,twin2_datacols_cov]
        X_dz_twin1 = data_dz[,twin1_datacols_cov]
        X_dz_twin2 = data_dz[,twin2_datacols_cov]
    } else {covariates = FALSE}
    
    #==========================================================
    #==========================================================
    # Subroutines for ACE model
    #==========================================================

    #==========================================================
    # Call subroutines for analyses based on sumscores
    #==========================================================
    #No covariates, no GxE
    if(itemdata == FALSE && ge == FALSE){
        source("sumscores.R")
        output = sumscores(data_mz = data_mz, data_dz = data_dz,
                           n_burnin = n_burnin, n_iter = n_iter)            
    } 
    
    #No covariates, GxE
    if(itemdata == FALSE && ge == TRUE){
        source("sumscores_ge.R")
        output = sumscores_ge(data_mz = data_mz, data_dz = data_dz,
                              n_burnin = n_burnin, n_iter = n_iter)            
    } 
    
    #Covariates, no GxE
    if(itemdata == FALSE && ge == FALSE && covariates == TRUE ){
        source("sumscores_cov.R")
        output = sumscores_cov(data_mz = data_mz, data_dz = data_dz,
                              X_mz_twin1 = X_mz_twin1, X_mz_twin2 = X_mz_twin2,
                              X_dz_twin1 = X_dz_twin1, X_dz_twin2 = X_dz_twin2,
                              n_burnin = n_burnin, n_iter = n_iter) 
    }
    
    #Covariates, GxE
    #if(itemdata == FALSE && ge == TRUE && covariates == TRUE ){
    #    source("sumscores_cov_ge.R")
    #    output = sumscores_cov_ge(data_mz = data_mz, data_dz = data_dz,
    #                              X_mz_twin1 = X_mz_twin1, X_mz_twin2 = X_mz_twin2,
    #                              X_dz_twin1 = X_dz_twin1, X_dz_twin2 = X_dz_twin2,
    #                              n_burnin = n_burnin, n_iter = n_iter) 
    #}
    

    #==========================================================
    # Call subroutines for item data analyses
    #==========================================================
    #No covariates, no GxE
    if(itemdata == TRUE && ge == FALSE){
        source("irt.R")
        output = irt(data_mz = data_mz, data_dz = data_dz, 
                     n_burnin = n_burnin, n_iter = n_iter)              
    } 
    
    #No covariates, GxE
    if(itemdata == TRUE && ge == TRUE){
        source("irt_ge.R")
        output = irt_ge(data_mz = data_mz, data_dz = data_dz,
                        n_burnin = n_burnin, n_iter = n_iter)            
    } 
    
    #Covariates, no GxE
    
    #Covariates, GxE
    
    
    #==========================================================
    #==========================================================
    # Subroutines for concordance rate analysis
    #==========================================================
    
    #==========================================================
    # Output
    #==========================================================
    source("plot.samples.R")
    
    #Remind user of convergence-issue
    cat("\n Here give some advice over convergence")
    cat("\n Here explain how to get results etc")
    
    #Print results on the fly: 
    cat("\n") 
    print(output$results) #table that is made in subroutine (always give same name!)
    return(output) 
}

#Run function: 
#Simulate data: 
#source("simulate_twin_data.R")
#data <- simulate_twin_data(nmz = 10, ndz = 20, var_a = 0.5, var_c = 0.3,  var_e = 0.2,
#                           model = "ACE", n_items = 8, n_var = 0)
#IRT + GE
#twin_analysis2 = twin_analysis(data_mz = data_mz, data_dz = data_dz, 
#                               twin1_datacols = 1:8, twin2_datacols = 9:16, 
#                               ordinal = T, model = "ACE", common = T,
#                               n_iter = 100, n_burnin = 100, ge = TRUE)
#twin_analysis2
#plot(twin_analysis2$samples_var_a, type = "Sampling plot")