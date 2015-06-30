#==========================================================
# simulate_AE.R
# Simulation of twin data under the AE model 
# Subroutine for simulate_twin_data.R
#==========================================================
simulate_AE = function(n_mz, n_dz, var_a, var_c, var_e,
                        n_items, n_var, ge, ge_beta0, ge_beta1){

    #I. Simulate genetic additive effects for MZ twins: 
    a_mz <- rnorm(n_mz, 0, sqrt(var_a)) 
    
    #II. Do the same for DZ twins: 
    a1_dz <- rnorm(n_dz, 0, sqrt(0.5 * var_a)) 
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

    #Generate answers for each twin, ACE model
    #When n_items > 0, generate item patterns according to an 1 PL Rasch model. 
    if(n_items > 0 && n_var > 0){
        ### I. MZ twins: 
        #Simulate data for the betas
        bp <- as.matrix(rnorm(n_items, 0,1))
        bp_mz <- t(matrix(bp, n_items, n_mz))
    
        #Simulate data for the coefficients: 
        cov_mz = matrix(0, n_mz, 2*n_var) #first 1:n_var col for twin 1, n_var+1 : 2*n_var col for twin 2
        coefficients = runif(n_var, -2, 2)
        beta = as.matrix(coefficients)
    
        #Generate answers to covariates: 
        for (i in 1:(2*n_var)){
            if(i == 1 || i == (n_var+1)){ #Intercept! 
                cov_mz[,i] = rep(1, n_mz)
            }else {cov_mz[,i]=rnorm(n_mz,0,1)} #rest of the colomns are answers to covariates. 
        }
    
        #Generate answers for each twin: (phenotype data)
        pheno_mz = matrix(NA, n_mz, 2)
        for (i in 1:n_mz){
            pheno_mz[i, 1] = rnorm(1, a_mz[i] + (cov_mz[i,1:n_var] %*% beta), sqrt(var_e_mz))
            pheno_mz[i, 2] = rnorm(1, a_mz[i] + (cov_mz[i,(n_var+1):(2*n_var)] %*% beta), sqrt(var_e_mz))
        }
    
        #cor(pheno_mz[,1], pheno_mz[,2]) #to check. must be ~ var_a + var_c
    
        #Trait values
        traits_mz_twin1 <- matrix(pheno_mz[,1], n_mz, n_items)
        traits_mz_twin2 <- matrix(pheno_mz[,2], n_mz, n_items)
    
        #Calculate p (simple Rasch model) for item data
        p_mz_twin1 <- (exp(traits_mz_twin1-bp_mz))/(1+(exp(traits_mz_twin1-bp_mz)))#MZ twins
        p_mz_twin2 <- (exp(traits_mz_twin2-bp_mz))/(1+(exp(traits_mz_twin2-bp_mz)))
        mz_twin1_itemdata <- t(apply(p_mz_twin1, 1, function(x) rbinom (n_items, 1, x)))
        mz_twin2_itemdata <- t(apply(p_mz_twin2, 1, function(x) rbinom (n_items, 1, x)))
    
        #Organize the data, suitable for MCMC analysis: 
        y_mz <- matrix(0,n_mz,(n_items*2))
        y_mz[,1:(n_items)] <- mz_twin1_itemdata
        y_mz[,(n_items+1):(n_items*2)] <- mz_twin2_itemdata
    
        ### II. DZ twins: 
        cov_dz = matrix(0, n_dz, 2*n_var)
    
        #Generate answers to covariates: 
        for (i in 1:(2*n_var)){
            if(i == 1 || i == (n_var+1)){
                cov_dz[,i] = rep(1, n_dz)
            }else {cov_dz[,i]=rnorm(n_dz,0,1)} #rest of the colomns are answers to covariates. 
        }
    
        #Generate answers for each twin: 
        pheno_dz = matrix(NA, n_dz, 2)
        for (i in 1:n_dz){
            pheno_dz[i,1] = rnorm(1, a2_dz[i,1] + (cov_dz[i,1:n_var] %*% beta), sqrt(var_e_dz_twin1))
            pheno_dz[i,2] = rnorm(1, a2_dz[i,2] + (cov_dz[i,(n_var+1):(2*n_var)] %*% beta) , sqrt(var_e_dz_twin2))
        }
    
        #cor(pheno_dz[,1], pheno_dz[,2]) #to check, must be ~0.5*var_a + var_c
        bp_dz <- t(matrix(bp, n_items, n_dz))
    
        #Trait values
        traits_dz_twin1 <- matrix(pheno_dz[,1], n_dz, n_items)#DZ twins
        traits_dz_twin2 <- matrix(pheno_dz[,2], n_dz, n_items)
    
        #Calculate p (simple Rasch model) for item data
        p_dz_twin1 <- (exp(traits_dz_twin1-bp_dz))/(1+(exp(traits_dz_twin1-bp_dz)))
        p_dz_twin2 <- (exp(traits_dz_twin2-bp_dz))/(1+(exp(traits_dz_twin2-bp_dz)))
        dz_twin1_itemdata <- t(apply(p_dz_twin1, 1, function(x) rbinom (n_items, 1, x))) 
        dz_twin2_itemdata <- t(apply(p_dz_twin2, 1, function(x) rbinom (n_items, 1, x)))
    
        #Organize the data, suitable for MCMC analysis: 
        y_dz <- matrix(0,n_dz,(n_items*2))
        y_dz[,1:(n_items)] <- dz_twin1_itemdata
        y_dz[,(n_items+1):(n_items*2)] <- dz_twin2_itemdata
    
    } else if (n_items > 0 && n_var == 0) {
        ### I. MZ twins: 
        #Simulate data for the betas
        bp <- as.matrix(rnorm(n_items, 0,1))
        bp_mz <- t(matrix(bp, n_items, n_mz))    
    
        #Phenotype data: 
        pheno_mz = cbind(rnorm(n_mz, a_mz, sqrt(var_e_mz)), 
                         rnorm(n_mz, a_mz, sqrt(var_e_mz)))
        
        #cor(pheno_mz[,1], pheno_mz[,2]) #to check. must be ~ var_a + var_c
    
        #Trait values
        traits_mz_twin1 <- matrix(pheno_mz[,1], n_mz, n_items)
        traits_mz_twin2 <- matrix(pheno_mz[,2], n_mz, n_items)
    
        #Calculate p (simple Rasch model) for item data
        p_mz_twin1 <- (exp(traits_mz_twin1-bp_mz))/(1+(exp(traits_mz_twin1-bp_mz)))#MZ twins
        p_mz_twin2 <- (exp(traits_mz_twin2-bp_mz))/(1+(exp(traits_mz_twin2-bp_mz)))
        mz_twin1_itemdata <- t(apply(p_mz_twin1, 1, function(x) rbinom (n_items, 1, x)))
        mz_twin2_itemdata <- t(apply(p_mz_twin2, 1, function(x) rbinom (n_items, 1, x)))
    
        #Organize the data, suitable for MCMC analysis: 
        y_mz <- matrix(0,n_mz,(n_items*2))
        y_mz[,1:(n_items)] <- mz_twin1_itemdata
        y_mz[,(n_items+1):(n_items*2)] <- mz_twin2_itemdata
    
        ### II. DZ twins: 
        pheno_dz = cbind(rnorm(n_dz, a2_dz[,1], sqrt(var_e_dz_twin1)),
                         rnorm(n_dz, a2_dz[,2], sqrt(var_e_dz_twin2)))
    
        #cor(pheno_dz[,1], pheno_dz[,2]) #to check, must be ~0.5*var_a + var_c
        bp_dz <- t(matrix(bp, n_items, n_dz))
    
        #Trait values
        traits_dz_twin1 <- matrix(pheno_dz[,1], n_dz, n_items)#DZ twins
        traits_dz_twin2 <- matrix(pheno_dz[,2], n_dz, n_items)
    
        #Calculate p (simple Rasch model) for item data
        p_dz_twin1 <- (exp(traits_dz_twin1-bp_dz))/(1+(exp(traits_dz_twin1-bp_dz)))
        p_dz_twin2 <- (exp(traits_dz_twin2-bp_dz))/(1+(exp(traits_dz_twin2-bp_dz)))
        dz_twin1_itemdata <- t(apply(p_dz_twin1, 1, function(x) rbinom (n_items, 1, x))) 
        dz_twin2_itemdata <- t(apply(p_dz_twin2, 1, function(x) rbinom (n_items, 1, x)))
    
        #Organize the data, suitable for MCMC analysis: 
        y_dz <- matrix(0,n_dz,(n_items*2))
        y_dz[,1:(n_items)] <- dz_twin1_itemdata
        y_dz[,(n_items+1):(n_items*2)] <- dz_twin2_itemdata
    
    } else if (n_items == 0 && n_var > 0){
        ### I.MZ twins
        #Simulate data for the coefficients: 
        cov_mz = matrix(0, n_mz, 2*n_var) #first 1:n_var col for twin 1, n_var+1 : 2*n_var col for twin 2
        coefficients = runif(n_var, -2, 2)
        beta = as.matrix(coefficients)
    
        #Generate answers to covariates: 
        for (i in 1:(2*n_var)){
            if(i == 1 || i == (n_var+1)){ #Intercept! 
                cov_mz[,i] = rep(1, n_mz)
            }else {cov_mz[,i]=rnorm(n_mz,0,1)} #rest of the colomns are answers to covariates. 
        }
    
        #Generate answers for each twin: (phenotype data)
        y_mz = matrix(NA, n_mz, 2)
        for (i in 1:n_mz){
            y_mz[i, 1] = rnorm(1, a_mz[i] + (cov_mz[i,1:n_var] %*% beta), sqrt(var_e_mz))
            y_mz[i, 2] = rnorm(1, a_mz[i] + (cov_mz[i,(n_var+1):(2*n_var)] %*% beta), sqrt(var_e_mz))
        }
    
        #cor(pheno_mz[,1], pheno_mz[,2]) #to check. must be ~ var_a + var_c
    
        ### II.DZ twins
        cov_dz = matrix(0, n_dz, 2*n_var)
    
        #Generate answers to covariates: 
        for (i in 1:(2*n_var)){
            if(i == 1 || i == (n_var+1)){
                cov_dz[,i] = rep(1, n_dz)
            }else {cov_dz[,i]=rnorm(n_dz,0,1)} #rest of the colomns are answers to covariates. 
        }
    
        #Generate answers for each twin: 
        y_dz = matrix(NA, n_dz, 2)
        for (i in 1:n_dz){
            y_dz[i,1] = rnorm(1, a2_dz[i,1] + (cov_dz[i,1:n_var] %*% beta), sqrt(var_e_dz_twin1))
            y_dz[i,2] = rnorm(1, a2_dz[i,2] + (cov_dz[i,(n_var+1):(2*n_var)] %*% beta) , sqrt(var_e_dz_twin2))
        }
    
    }   else { 
        ### I. MZ twins: 
        y_mz = cbind(rnorm(n_mz, a_mz,  sqrt(var_e_mz)), 
                     rnorm(n_mz, a_mz, sqrt(var_e_mz)))
    
        ### II. DZ twins: 
        y_dz = cbind(rnorm(n_dz, a2_dz[,1], sqrt(var_e_dz_twin1)),
                     rnorm(n_dz, a2_dz[,2], sqrt(var_e_dz_twin2)))
    
        #cor(y_mz[,1], y_mz[,2]) #to check. must be ~ var_a + var_c
        #cor(y_dz[,1], y_dz[,2]) #to check, must be ~ 0.5*var_a + var_c
        
        if(n_var > 0){return_list = list(y_mz = y_mz, y_dz = y_dz, cov_mz = cov_mz[,c(2:n_var, (n_var + 2):(2*n_var))], 
                                         cov_dz = cov_dz[,c(2:n_var, (n_var + 2):(2*n_var))])
        } else {return_list = list(y_mz = y_mz, y_dz = y_dz)}
        
        return(return_list)
        }
}