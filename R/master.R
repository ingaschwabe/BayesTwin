#==========================================================
# Master function that calls subroutines
# for BayesTwin package
#==========================================================


setwd("C:/Users/schwabei/Dropbox/International student performance_IngaStephanie/R/BayesTwin3/R")
setwd("/Users/stephanievandenberg/Dropbox/International student performance_IngaStephanie/R/BayesTwin3/R")

library(rjags)

#Simulate data: 
source("simulate_twin_data.R")

data <- simulate_twin_data(nmz = 300, ndz = 500, var_a = 0.5, var_c = 0.3,  var_e = 0.2,
                           model = "ACE", n_items = 10, n_var = 0)

    
data_mz = data$y_mz    
data_dz = data$y_dz    


twin_analysis = function(data_mz, data_dz, twin1_datacols, twin2_datacols, 
                         multivariate = F, ordinal = F, model = "ACE", common = T,
                         n_iter = 10000, n_burnin = 5000, ge = FALSE, n_items){
    
    #Install packages when necessary (later niet meer nodig!!)
    if(!require(MCMCpack)){ install.packages('MCMCpack'); require(MCMCpack)}
    if(!require(R.utils)){ install.packages('R.utils'); require(R.utils)}
    if(!require(MCMCpack)){ install.packages('MASS'); require(MASS)}    
    if(!require(MCMCpack)){ install.packages('rjags'); require(rjags)}      
                                    
    #Load HPD function: 
     source("HPD.R")                           
                                    
                                    
    #Change data format:                                
    
    #from .dat to R ->
    
    #from .sav to R -> 
                                    
    #Calculate number of MZ and DZ twins: 
    n_mz = nrow(data_mz)
    n_dz = nrow(data_dz)
    
    #Select item data: 
    itemdata_mz = data_mz[,c(twin1_datacols, twin2_datacols]
                        
    
    itemdata_dz = data_dz[,c(twin1_datacols, twin2_datacols]
        
    #GE + item data: 
    if(ge == TRUE){
        source("ge_irt.R")
        out = ge_irt(y_mz = itemdata_mz, y_dz = itemdata_dz, 
                     n_dz = n_dz, n_mz = n_mz, n_items = n_items, 
                     burnin = n_burnin, n_iterations = n_iter)    
    
        
        #==========================================================
        # III. Organize results  
        #==========================================================
        
        #Calculate mean values: 
        var_a = mean(1/out$tau_a[,,1]); var_c = mean(1/out$tau_c[,,1])
        beta0 = mean(exp(out$beta0[,,1])); beta1 = mean(out$beta1[,,1])
        
        #Calculate SDs: 
        sd_var_a = sd(1/out$tau_a[,,1]); sd_var_c = sd(1/out$tau_c[,,1])
        sd_beta0 = sd(exp(out$beta0[,,1])); sd_beta1 = sd(out$beta1[,,1])
        
        #Calculate HPD: 
        hpd_var_a = HPD(1/out$tau_a[,,1], 0.95); hpd_var_c = HPD(1/out$tau_c[,,1], 0.95)
        hpd_beta0 = HPD(exp(out$beta0[,,1]), 0.95); hpd_beta1 = HPD(out$beta1[,,1], 0.95)
        
        #Put results in a table 
        results = matrix(c(var_a, sd_var_a, var_c, sd_var_c, beta0, sd_beta0, beta1, sd_beta1), 2, 4)
        colnames(results) <- c("varA","varC","beta0", "beta1")
        rownames(results) <- c("Posterior point estimate","Standard deviation")
        results = as.table(results)
        
        #Print results on the fly: 
        cat("\n") 
        print(results)
        cat("\n newline")
        
        #Remind the user of the issue of convergence: 
        
        
        #And save output in a list: 
        output = list(results = results, samples_var_a = 1/out$tau_a[,,1], samples_var_c = 1/out$tau_c[,,1], 
                      samples_beta0 = exp(out$beta0[,,1]), samples_beta1 = out$beta1[,,1],
                      hpd_var_a = hpd_var_a, hpd_var_c = hpd_var_c, hpd_beta0 = hpd_beta0, hpd_beta1 = hpd_beta1)
        
        #Change class of objects in order to use right plot method 
        class(output) <- c("ge_irt","list")
        class(output$samples_var_a) <- "samples"
        class(output$samples_var_c) <- "samples"
        class(output$samples_beta0) <- "samples"
        class(output$samples_beta1) <- "samples"        
    } 
    #item data, no GE: 
    if(ge == FALSE){
            source("irt.R")
            out = ge_irt(y_mz = itemdata_mz, y_dz = itemdata_dz, 
                         n_dz = n_dz, n_mz = n_mz, n_items = n_items, 
                         burnin = n_burnin, n_iterations = n_iter)    
            
            
            #==========================================================
            # III. Organize results  
            #==========================================================
            
            #Calculate mean values: 
            var_a = mean(1/out$tau_a[,,1]); var_c = mean(1/out$tau_c[,,1])
            var_e = mean(1/out$tau_e[,,1]); b = apply(out$b,1, mean)
            
            #Calculate SDs: 
            sd_var_a = sd(1/out$tau_a[,,1]); sd_var_c = sd(1/out$tau_c[,,1])
            sd_var_e = sd(1/out$tau_e[,,1]); sd_b = apply(out$b,1,sd)
            
            #Calculate HPD: 
            hpd_var_a = HPD(1/out$tau_a[,,1], 0.95); hpd_var_c = HPD(1/out$tau_c[,,1], 0.95)
            hpd_var_e = HPD(1/out$tau_e[,,1], 0.95)
            
            #Put results in a table 
            results = matrix(c(var_a, sd_var_a, var_c, sd_var_c, var_e, sd_var_e), 2, 3)
            colnames(results) <- c("varA","varC","varE")
            rownames(results) <- c("Posterior point estimate","Standard deviation")
            results = as.table(results)
            
            #Print results on the fly: 
            cat("\n") 
            print(results)
            cat("\n newline")
            
            #Remind the user of the issue of convergence: 
            
            
            #And save output in a list: 
            output = list(results = results, samples_var_a = 1/out$tau_a[,,1], samples_var_c = 1/out$tau_c[,,1], 
                          samples_var_e = 1/out$tau_e[,,1], samples_b = out$b[,,1],
                          hpd_var_a = hpd_var_a, hpd_var_c = hpd_var_c, hpd_var_e=hpd_var_e)
            
            #Change class of objects in order to use right plot method 
            class(output) <- c("irt","list")
            class(output$samples_var_a) <- "samples"
            class(output$samples_var_c) <- "samples"
            class(output$samples_var_e) <- "samples"
                 
    } 
    
    
    
    
    source("plot.samples.R")
    #Output: 
    return(output) 
}

plot(x)

# IRT + GE
twin_analysis2 = twin_analysis(data_mz = data_mz, data_dz = data_dz, 
                               twin1_datacols = 1:10, twin2_datacols = 11:20, 
                               multivariate = F, ordinal = F, model = "ACE", common = T,
                               n_iter = 100, n_burnin = 100, ge = TRUE, n_items = 10)
twin_analysis2
plot(twin_analysis2$samples_var_a, type = "Sampling plot")

# IRT no GE
twin_analysis1 = twin_analysis(data_mz = data_mz, data_dz = data_dz, 
                               twin1_datacols = 1:10, twin2_datacols = 11:20, 
                               multivariate = F, ordinal = F, model = "ACE", common = T,
                               n_iter = 100, n_burnin = 100, ge = FALSE, n_items = 10)
twin_analysis1
plot(twin_analysis1$samples_var_a, type = "Sampling plot")





