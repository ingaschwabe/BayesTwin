#==========================================================
# simulate_ADE.R
# Simulation of twin data under the ADE model 
# Subroutine for master function simulate_twin_data()
# BayesTwin package
#==========================================================

simulate_ADE = function(n_mz, n_dz, var_a, var_d, var_e, n_items, 
                        ge, ge_beta0, ge_beta1, irt_model,
                        n_cat){
    
    ## MZ twins: 
    # Simulate gentic additive and dominance effects:
    a_mz = rnorm(n_mz, 0, sqrt(var_a)) #Phenotypic population mean set to zero to
    g_mz <- rnorm(n_mz, a_mz, sqrt(var_d)) 
    #var(g_mz) #to check, should be ~ var_a + var_D
    
    ##Do the same for DZ twins 
    #Additive genetic correlation is 0.5, nonadditive correlation is .25
    #Between VAR
    a_dz = rnorm(n_dz, 0, sqrt(0.5 * var_a))
    g_dz = rnorm(n_dz, a_dz, sqrt(0.25 * var_d))
    
    #Within VAR
    g1_dz =  cbind(rnorm(n_dz, g_dz, sqrt(0.5 * var_a)),
                   rnorm(n_dz, g_dz, sqrt(0.5 * var_a)))
    g2_dz = cbind(rnorm(n_dz, g1_dz[,1], sqrt(0.75 * var_d)),
                  rnorm(n_dz, g1_dz[,2], sqrt(0.75 * var_d)))
    #cor(g2_dz[,1], g2_dz[,2]) #to check, should be ~ 1/2 var_a + 1/4 var_d
    
    if(ge == TRUE){
        var_e_mz = exp(ge_beta0 + (ge_beta1 * g_mz))
        var_e_dz_twin1 = exp(ge_beta0 + (ge_beta1 * g2_dz[,1]))
        var_e_dz_twin2 = exp(ge_beta0 + (ge_beta1 * g2_dz[,2]))   
        print("Simulating model with GxE interaction...")
    } else {
        #In case of no GE:
        var_e_mz = var_e
        var_e_dz_twin1 = var_e
        var_e_dz_twin2 = var_e
    }  
    
    #Phenotype data: 
    pheno_mz = cbind(rnorm(n_mz, g_mz, sqrt(var_e_mz)), 
                     rnorm(n_mz, g_mz, sqrt(var_e_mz)))
    pheno_dz = cbind(rnorm(n_dz, g2_dz[,1], sqrt(var_e_dz_twin1)),
                     rnorm(n_dz, g2_dz[,2], sqrt(var_e_dz_twin2)))
    #cor(pheno_mz[,1], pheno_mz[,2]) #to check, must be ~ var_a + var_d
    #cor(pheno_dz[,1], pheno_dz[,2]) #to check, must be ~ 1/2 var_a + 1/4 var_d
    
    if(irt_model == "1PL"){
        print("Using 1 PL model to generate item data...")
        
        #Trait values
        traits_mz_twin1 <- matrix(pheno_mz[,1], n_mz, n_items)
        traits_mz_twin2 <- matrix(pheno_mz[,2], n_mz, n_items)
        
        traits_dz_twin1 <- matrix(pheno_dz[,1], n_dz, n_items)#DZ twins
        traits_dz_twin2 <- matrix(pheno_dz[,2], n_dz, n_items)
        
        #Generate item patterns according to IRT model
        #Simulate data for the betas
        bp <- as.matrix(rnorm(n_items, 0,1))
        bp_mz <- t(matrix(bp, n_items, n_mz))  
        bp_dz <- t(matrix(bp, n_items, n_dz))  
        
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
        
    }else if (irt_model == "PCM"){
        print("Using Partial Credit model to generate item data...")        
        
        #The first threshold is always equal to zero 
        #The next (second) threshold is drawn from a standard normal distribution. 
        #Rest of the threshold is placed with equal distance to prevent unreasonable frequencies. 
        bp = matrix(NA, n_items, n_cat)
        bp[,1] = 0 #moet altijd 0 zijn 
        bp[,2] = rnorm(n_items, 0,1)#1 itemlocatie parameter, rest met gelijke afstanden hieromheen 
        
        for(i in 3:n_cat){
            bp[,i] = bp[,(i-1)] + 0.5 
        }
        
        alpha = 1
        
        #MZ twins first: 
        mz_twin1_itemdata = matrix(NA, n_mz, n_items)
        mz_twin2_itemdata = matrix(NA, n_mz, n_items)
        
        for (twin in 1:n_mz){
            for (item in 1:n_items){   
                #print(paste("twin",twin))
                #print(paste("item",item))
                mz_twin1_itemdata[twin,item] = (sum(cprm(alpha,bp[item,],pheno_mz[twin,1])))
                mz_twin2_itemdata[twin,item] = (sum(cprm(alpha,bp[item,],pheno_mz[twin,2])))
            }
        }
        
        #DZ twins: 
        dz_twin1_itemdata = matrix(NA, n_dz, n_items)
        dz_twin2_itemdata = matrix(NA, n_dz, n_items)
        
        for (twin in 1:n_dz){
            for (item in 1:n_items){   
                #print(paste("twin",twin))
                #print(paste("item",item))
                dz_twin1_itemdata[twin,item] = (sum(cprm(alpha,bp[item,],pheno_dz[twin,1])))
                dz_twin2_itemdata[twin,item] = (sum(cprm(alpha,bp[item,],pheno_dz[twin,2])))
            }
        }
        
        #Organize the data, suitable for MCMC analysis: 
        y_mz <- matrix(0,n_mz,(n_items*2))
        y_mz[,1:(n_items)] <- mz_twin1_itemdata
        y_mz[,(n_items+1):(n_items*2)] <- mz_twin2_itemdata
        
        y_dz <- matrix(0,n_dz,(n_items*2))
        y_dz[,1:(n_items)] <- dz_twin1_itemdata
        y_dz[,(n_items+1):(n_items*2)] <- dz_twin2_itemdata
        
        y_mz = round(y_mz, 0)
        y_dz = round(y_dz, 0)
        
    } else if (irt_model == "GPCM"){
        print("Using Generalized Partial Credit model to generate item data...")
        
        #The first threshold is always equal to zero 
        #The next (second) threshold is drawn from a standard normal distribution. 
        #Rest of the threshold is placed with equal distance to prevent unreasonable frequencies. 
        bp = matrix(NA, n_items, n_cat)
        bp[,1] = 0 #moet altijd 0 zijn 
        bp[,2] = rnorm(n_items, 0,1)#1 itemlocatie parameter, rest met gelijke afstanden hieromheen 
        
        for(i in 3:n_cat){
            bp[,i] = bp[,(i-1)] + 0.5 
        }
        
        alpha <- as.matrix(runif(n_items, .75,1.25))
        
        #MZ twins first: 
        mz_twin1_itemdata = matrix(NA, n_mz, n_items)
        mz_twin2_itemdata = matrix(NA, n_mz, n_items)
        
        for (twin in 1:n_mz){
            for (item in 1:n_items){   
                #print(paste("twin",twin))
                #print(paste("item",item))
                mz_twin1_itemdata[twin,item] = (sum(cprm(alpha[item],bp[item,],pheno_mz[twin,1])))
                mz_twin2_itemdata[twin,item] = (sum(cprm(alpha[item],bp[item,],pheno_mz[twin,2])))
            }
        }
        
        #DZ twins: 
        dz_twin1_itemdata = matrix(NA, n_dz, n_items)
        dz_twin2_itemdata = matrix(NA, n_dz, n_items)
        
        for (twin in 1:n_dz){
            for (item in 1:n_items){   
                #print(paste("twin",twin))
                #print(paste("item",item))
                dz_twin1_itemdata[twin,item] = (sum(cprm(alpha[item],bp[item,],pheno_dz[twin,1])))
                dz_twin2_itemdata[twin,item] = (sum(cprm(alpha[item],bp[item,],pheno_dz[twin,2])))
            }
        }
        
        #Organize the data, suitable for MCMC analysis: 
        y_mz <- matrix(0,n_mz,(n_items*2))
        y_mz[,1:(n_items)] <- mz_twin1_itemdata
        y_mz[,(n_items+1):(n_items*2)] <- mz_twin2_itemdata
        
        y_dz <- matrix(0,n_dz,(n_items*2))
        y_dz[,1:(n_items)] <- dz_twin1_itemdata
        y_dz[,(n_items+1):(n_items*2)] <- dz_twin2_itemdata
        
        y_mz = round(y_mz, 0)
        y_dz = round(y_dz, 0)
        
    } else {
        print("Using 2 PL model to generate item data...")
        
        #Trait values
        traits_mz_twin1 <- matrix(pheno_mz[,1], n_mz, n_items)
        traits_mz_twin2 <- matrix(pheno_mz[,2], n_mz, n_items)
        
        traits_dz_twin1 <- matrix(pheno_dz[,1], n_dz, n_items)#DZ twins
        traits_dz_twin2 <- matrix(pheno_dz[,2], n_dz, n_items)
        
        #Generate item patterns according to IRT model
        #Simulate data for the betas
        bp <- as.matrix(rnorm(n_items, 0,1))
        bp_mz <- t(matrix(bp, n_items, n_mz))  
        bp_dz <- t(matrix(bp, n_items, n_dz))  
        
        alpha <- as.matrix(runif(n_items, .75,1.25))
        alpha_mz <- t(matrix(alpha, n_items, n_mz))   
        alpha_dz <- t(matrix(alpha, n_items, n_dz))
        
        #Calculate p (2PL) for item data
        p_mz_twin1 <- (exp(alpha_mz*(traits_mz_twin1-bp_mz)))/(1+(exp(alpha_mz*(traits_mz_twin1-bp_mz))))#MZ twins
        p_mz_twin2 <- (exp(alpha_mz*(traits_mz_twin2-bp_mz)))/(1+(exp(alpha_mz*(traits_mz_twin2-bp_mz))))#MZ twins
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