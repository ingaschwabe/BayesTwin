#==========================================================
# simulate_ACE.R
# Simulation of twin data under the ACE model 
# Subroutine for master function simulate_twin_data()
# 
# BayesTwin - An R Package for Bayesian Inference of Item-Level Twin Data
# Copyright (C) 2014-2017 Inga Schwabe
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details
#==========================================================

simulate_ACE = function(n_mz, n_dz, var_a, var_c, var_e, n_items, 
                        ge, ge_beta0, ge_beta1, irt_model,
                        n_cat){
    
    ## MZ twins:
    #Simulate common environmental effects, familial effects and retrieve genetic values for
    c_mz = rnorm(n_mz, 0, sqrt(var_c)) #phenotypic population mean set to zero 
    f_mz <- rnorm(n_mz, c_mz, sqrt(var_a)) #familial influences = common environment + genetic influences
    a_mz <- f_mz - c_mz #genetic influences
    #var(f_mz) #to check, should be ~ var_c + var_a
    
    ## Do the same for DZ twins  
    c_dz = rnorm(n_dz, 0, sqrt(var_c)) #as for MZ twins: phenotypic population mean set to zero 
    f1_dz <- rnorm(n_dz, c_dz, sqrt(0.5 * var_a)) #familial influences = common environment + genetic influences
    f2_dz <- cbind(rnorm(n_dz, f1_dz, sqrt(0.5 * var_a)), #genetic correlation of 1/2 in DZ twins
                   rnorm(n_dz, f1_dz, sqrt(0.5 * var_a)))
    
    a_dz_twin1 <- f2_dz[,1] - c_dz #genetic influences DZ twin 1
    a_dz_twin2 <- f2_dz[,2] - c_dz #genetic influences DZ twin 2
    #cor(a_dz_twin1, a_dz_twin2) #to check, should be ~ 1/2 var_a + var_c
    
    #Change residual variance in case of GE: 
    if(ge == TRUE){
        var_e_mz = exp(ge_beta0 + (ge_beta1 * a_mz))
        var_e_dz_twin1 = exp(ge_beta0 + (ge_beta1 * a_dz_twin1))
        var_e_dz_twin2 = exp(ge_beta0 + (ge_beta1 * a_dz_twin2))  
        print("Simulating model with GxE interaction...")
    } else {
        #In case of no GE:
        var_e_mz = var_e
        var_e_dz_twin1 = var_e
        var_e_dz_twin2 = var_e
    }
    
    #Phenotype data: 
    pheno_mz = cbind(rnorm(n_mz, f_mz, sqrt(var_e_mz)), 
                     rnorm(n_mz, f_mz, sqrt(var_e_mz)))
    
    pheno_dz = cbind(rnorm(n_dz, f2_dz[,1], sqrt(var_e_dz_twin1)),
                     rnorm(n_dz, f2_dz[,2], sqrt(var_e_dz_twin2)))
    
    #cor(pheno_mz[,1], pheno_mz[,2]) #to check, must be ~ var_a + var_c
    #cor(pheno_dz[,1], pheno_dz[,2]) #to check, must be ~ 1/2 var_a + var_c
    
    if(irt_model == "1PL"){
        print("Using 1 PL model to generate item data...")
        
        #Trait values
        traits_mz_twin1 <- matrix(pheno_mz[,1], n_mz, n_items)
        traits_mz_twin2 <- matrix(pheno_mz[,2], n_mz, n_items)
        
        traits_dz_twin1 <- matrix(pheno_dz[,1], n_dz, n_items)#DZ twins
        traits_dz_twin2 <- matrix(pheno_dz[,2], n_dz, n_items)
        
        #Generate item patterns 
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
    } else if (irt_model == "PCM"){
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
        
        #Generate item patterns 
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