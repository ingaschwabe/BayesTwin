#' Simulation of twin data
#'
#' This function simulates twin data. The data can be generated as sumscores but also on item level. Item data is simulated under 
#' a 1 PL Rasch model. 
#' @param nmz refers tot the total number of MZ twins. Defaults to 300
#' @param ndz refers tot the total number of DZ twins. Defaults to 500
#' @param var_a refers to variance due to additive genetic effects. Defaults to 0.5
#' @param var_c refers to variance due to shared-environmental effects. Defaults to 0.3
#' @param var_e refers to variance due to unique-environmental effects. Defautls to 0.2
#' @param model refers to used model ("ACE", "ADE"). Defaults to "ACE"
#' @param n_items refers to the number of items. Defaults to 0 (sumscores are simulated)
#' @keywords twin data, ACE, ADE
#' @export
#' @examples
#' simulate_twin_data(n_items = 60)
#' Simualtes twin data on 60 items with the other parameters set to their default values.  
simulate_twin_data <- function(nmz = 300, ndz = 500, var_a = 0.5, var_c = 0.3,  var_e = 0.2,
                               model = "ACE", n_items = 0){
    
    #I. Simulate genetic additive effects and common shared env. effects for MZ twins: 
    c_mz = rnorm(nmz, 0, sqrt(var_c)) 
    a_mz <- rnorm(nmz, c_mz, sqrt(var_a)) 
    
    #II. Do the same for DZ twins: 
    c_dz = rnorm(ndz, 0, sqrt(var_c)) 
    a1_dz <- rnorm(ndz, c_dz, sqrt(0.5 * var_a)) 
    a2_dz <- cbind(rnorm(ndz, a1_dz, sqrt(0.5 * var_a)),
                   rnorm(ndz, a1_dz, sqrt(0.5 * var_a)))
    
    #Generate answers for each twin, ACE model
    #When n_items > 0, generate item patterns according to an 1 PL Rasch model. 
    if(n_items > 0){
        ### I. MZ twins: 
        #Simulate data for the betas
        bp <- as.matrix(rnorm(n_items, 0,1))
        bp_mz <- t(matrix(bp, n_items, nmz))    
        
        #Phenotype data: 
        pheno_mz = cbind(rnorm(nmz, a_mz, sqrt(var_e)), 
                         rnorm(nmz, a_mz, sqrt(var_e)))
        
        #cor(pheno_mz[,1], pheno_mz[,2]) #to check. must be ~ var_a + var_c
        
        #Trait values
        traits_mz_twin1 <- matrix(pheno_mz[,1], nmz, n_items)
        traits_mz_twin2 <- matrix(pheno_mz[,2], nmz, n_items)
        
        #Calculate p (simple Rasch model) for item data
        p_mz_twin1 <- (exp(traits_mz_twin1-bp_mz))/(1+(exp(traits_mz_twin1-bp_mz)))#MZ twins
        p_mz_twin2 <- (exp(traits_mz_twin2-bp_mz))/(1+(exp(traits_mz_twin2-bp_mz)))
        mz_twin1_itemdata <- t(apply(p_mz_twin1, 1, function(x) rbinom (n_items, 1, x)))
        mz_twin2_itemdata <- t(apply(p_mz_twin2, 1, function(x) rbinom (n_items, 1, x)))
        
        #Organize the data, suitable for MCMC analysis: 
        y_mz <- matrix(0,nmz,(n_items*2))
        y_mz[,1:(n_items)] <- mz_twin1_itemdata
        y_mz[,(n_items+1):(n_items*2)] <- mz_twin2_itemdata
        
        ### II. DZ twins: 
        pheno_dz = cbind(rnorm(ndz, a2_dz[,1], sqrt(var_e)),
                         rnorm(ndz, a2_dz[,2], sqrt(var_e)))
        
        #cor(pheno_dz[,1], pheno_dz[,2]) #to check, must be ~0.5*var_a + var_c
        
        bp_dz <- t(matrix(bp, n_items, ndz))
        
        #Trait values
        traits_dz_twin1 <- matrix(pheno_dz[,1], ndz, n_items)#DZ twins
        traits_dz_twin2 <- matrix(pheno_dz[,2], ndz, n_items)
        
        #Calculate p (simple Rasch model) for item data
        p_dz_twin1 <- (exp(traits_dz_twin1-bp_dz))/(1+(exp(traits_dz_twin1-bp_dz)))
        p_dz_twin2 <- (exp(traits_dz_twin2-bp_dz))/(1+(exp(traits_dz_twin2-bp_dz)))
        dz_twin1_itemdata <- t(apply(p_dz_twin1, 1, function(x) rbinom (n_items, 1, x))) 
        dz_twin2_itemdata <- t(apply(p_dz_twin2, 1, function(x) rbinom (n_items, 1, x)))
        
        #Organize the data, suitable for MCMC analysis: 
        y_dz <- matrix(0,ndz,(n_items*2))
        y_dz[,1:(n_items)] <- dz_twin1_itemdata
        y_dz[,(n_items+1):(n_items*2)] <- dz_twin2_itemdata
        
    } else {
        y_mz = cbind(rnorm(nmz, a_mz,  sqrt(var_e)), 
                     rnorm(nmz, a_mz, sqrt(var_e)))
        
        y_dz = cbind(rnorm(ndz, a2_dz[,1], sqrt(var_e)),
                     rnorm(ndz, a2_dz[,2], sqrt(var_e)))
        
        #cor(y_mz[,1], y_mz[,2]) #to check. must be ~ var_a + var_c
        #cor(y_dz[,1], y_dz[,2]) #to check, must be ~ 0.5*var_a + var_c
    }
    
    
    return(list(y_mz = y_mz, y_dz = y_dz))
    
}#end function simulate_twin_data