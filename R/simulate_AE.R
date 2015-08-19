#==========================================================
# simulate_AE.R
# Simulation of twin data under the AE model 
# Subroutine for master function simulate_twin_data()
# BayesTwin package
#==========================================================

simulate_AE = function(n_mz, n_dz, var_a, var_e,
                       n_items, ge, ge_beta0, ge_beta1,
                       irt_model){
    ## MZ twins:
    #Simulate genetic additive effects 
    a_mz <- rnorm(n_mz, 0, sqrt(var_a)) 
    
    ## Do the same for DZ twins
    a1_dz <- rnorm(n_dz, 0, sqrt(0.5 * var_a)) #genetic correlation of 0.5
    a2_dz <- cbind(rnorm(n_dz, a1_dz, sqrt(0.5 * var_a)),
                   rnorm(n_dz, a1_dz, sqrt(0.5 * var_a)))
    
    #Change residual variance in case of GE: 
    if(ge == TRUE){
        var_e_mz = exp(ge_beta0 + (ge_beta1 * a_mz))
        var_e_dz_twin1 = exp(ge_beta0 + (ge_beta1 * a2_dz[,1]))
        var_e_dz_twin2 = exp(ge_beta0 + (ge_beta1 * a2_dz[,2]))     
        print("Simulating model with GxE interaction...")
    } else {
        #In case of no GE:
        var_e_mz = var_e
        var_e_dz_twin1 = var_e
        var_e_dz_twin2 = var_e
    }
    
    #Phenotype data: 
    pheno_mz = cbind(rnorm(n_mz, a_mz, sqrt(var_e_mz)), 
                     rnorm(n_mz, a_mz, sqrt(var_e_mz)))
    
    pheno_dz = cbind(rnorm(n_dz, a2_dz[,1], sqrt(var_e_dz_twin1)),
                     rnorm(n_dz, a2_dz[,2], sqrt(var_e_dz_twin2)))
    #cor(pheno_mz[,1], pheno_mz[,2]) #to check, must be ~ var_a 
    #cor(pheno_dz[,1], pheno_dz[,2]) #to check, must be ~ 1/2 var_a 
    
    #Trait values
    traits_mz_twin1 <- matrix(pheno_mz[,1], n_mz, n_items)
    traits_mz_twin2 <- matrix(pheno_mz[,2], n_mz, n_items)
    
    traits_dz_twin1 <- matrix(pheno_dz[,1], n_dz, n_items)#DZ twins
    traits_dz_twin2 <- matrix(pheno_dz[,2], n_dz, n_items)
    
    #Generate answers patterns under IRT model
    
    #Simulate data for the betas 
    bp <- as.matrix(rnorm(n_items, 0,1))
    bp_mz <- t(matrix(bp, n_items, n_mz))
    bp_dz <- t(matrix(bp, n_items, n_dz))
    
    if(irt_model == "1PL"){
        print("Using 1 PL model to generate item data...")
        #Calculate p (simple Rasch model) for item data
        p_mz_twin1 <- (exp(traits_mz_twin1-bp_mz))/(1+(exp(traits_mz_twin1-bp_mz)))#MZ twins
        p_mz_twin2 <- (exp(traits_mz_twin2-bp_mz))/(1+(exp(traits_mz_twin2-bp_mz)))
        mz_twin1_itemdata <- t(apply(p_mz_twin1, 1, function(x) rbinom (n_items, 1, x)))
        mz_twin2_itemdata <- t(apply(p_mz_twin2, 1, function(x) rbinom (n_items, 1, x)))
        
        #Organize the data, suitable for MCMC analysis: 
        y_mz <- matrix(0,n_mz,(n_items*2))
        y_mz[,1:(n_items)] <- mz_twin1_itemdata
        y_mz[,(n_items+1):(n_items*2)] <- mz_twin2_itemdata
            
        #Calculate p (simple Rasch model) for item data
        p_dz_twin1 <- (exp(traits_dz_twin1-bp_dz))/(1+(exp(traits_dz_twin1-bp_dz)))
        p_dz_twin2 <- (exp(traits_dz_twin2-bp_dz))/(1+(exp(traits_dz_twin2-bp_dz)))
        dz_twin1_itemdata <- t(apply(p_dz_twin1, 1, function(x) rbinom (n_items, 1, x))) 
        dz_twin2_itemdata <- t(apply(p_dz_twin2, 1, function(x) rbinom (n_items, 1, x)))
        
        #Organize the data, suitable for MCMC analysis: 
        y_dz <- matrix(0,n_dz,(n_items*2))
        y_dz[,1:(n_items)] <- dz_twin1_itemdata
        y_dz[,(n_items+1):(n_items*2)] <- dz_twin2_itemdata
    } else {
        print("Using 2 PL model to generate item data...")
        alpha <- as.matrix(runif(n_items, .75,1.25))
        alpha_mz <- t(matrix(alpha, n_items, n_mz))   
        alpha_dz <- t(matrix(alpha, n_items, n_dz))
        
        #Calculate p (2PL) for item data
        p_mz_twin1 <- (exp(traits_mz_twin1-bp_mz))/(1+(exp(traits_mz_twin1-bp_mz)))
        p_mz_twin2 <- (exp(traits_mz_twin2-bp_mz))/(1+(exp(traits_mz_twin2-bp_mz)))
        mz_twin1_itemdata <- t(apply(p_mz_twin1, 1, function(x) rbinom (n_items, 1, x)))
        mz_twin2_itemdata <- t(apply(p_mz_twin2, 1, function(x) rbinom (n_items, 1, x)))
        
        #Organize the data, suitable for MCMC analysis: 
        y_mz <- matrix(0,n_mz,(n_items*2))
        y_mz[,1:(n_items)] <- mz_twin1_itemdata
        y_mz[,(n_items+1):(n_items*2)] <- mz_twin2_itemdata
        
        #Calculate p (simple Rasch model) for item data
        p_dz_twin1 <- (exp(alpha_dz*(traits_dz_twin1-bp_dz)))/(1+(exp(alpha_dz*(traits_dz_twin1-bp_dz))))#DZ twins
        p_dz_twin2 <- (exp(alpha_dz*(traits_dz_twin2-bp_dz)))/(1+(exp(alpha_dz*(traits_dz_twin2-bp_dz))))#DZ twins
        dz_twin1_itemdata <- t(apply(p_dz_twin1, 1, function(x) rbinom (n_items, 1, x))) 
        dz_twin2_itemdata <- t(apply(p_dz_twin2, 1, function(x) rbinom (n_items, 1, x)))
        
        #Organize the data, suitable for MCMC analysis: 
        y_dz <- matrix(0,n_dz,(n_items*2))
        y_dz[,1:(n_items)] <- dz_twin1_itemdata
        y_dz[,(n_items+1):(n_items*2)] <- dz_twin2_itemdata
    }
    
    return_list = list(y_mz = y_mz, y_dz = y_dz)
    return(return_list)
}