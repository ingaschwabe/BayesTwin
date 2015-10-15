###############################################################
### Simulation to validate BayesTwin R packae 
### I Validate all models that use the ACE variance 
### decomposition model 
##############################################################

#Libraries
library(rjags)

#Set working directory
setwd("/Users/Inga/Dropbox/International student performance_IngaStephanie/R/BayesTwin3/R")

## Source functions 
source("simulate_twin_data.R")
source("twinUniv.R")

###############################################################
### Choose simulation values 
##############################################################

### Number of MZ and DZ twin pairs was fixed: 
### as was A,C and E and number of items (60):
n_pairs = 100; n_mz = round(n_pairs * 0.28); n_dz = n_pairs - n_mz
var_a = 0.5; var_c = 0.3; var_e = 0.2
n_burnin = 2; n_iter = 2

### Different conditions, varying values for: 
### 1) IRT model ("1PL/"2PL"/"GPCM"/"PCM")
### 2) Prior for variance components:  "INV GAMMA"/"UNIF"
### 3) Genotype by environment interaction (TRUE/FALSE)
irt_model = "1PL"; 
var_prior = "INV_GAMMA"
ge = TRUE; if(ge == TRUE){ge_beta0 = log(0.5); ge_beta1 = 0}

###########################################################
## Following function runs the simulation,
## taking the simulation values defined as above
## Function is called by for loop at the end of this script
##
## If you want to do an analysis without GxE (with GxE), 
## choose dummy values (e.g. = NA) for ge_beta0 and 
## ge_beta1 (for var_e). 
###########################################################

test_bayestwin_ACE <- function(var_a, var_c, var_e, irt_model, ge, n_items,
                               ge_beta0, ge_beta1){
    
    ### choose simulation input values and analysis based on GE value (TRUE/FALSE)
    if(ge == TRUE){
        ### Simulate data using simulate_twin_data.R
        print("Simulating data for ACE model with GxE...")
        data = simulate_twin_data(n_mz = n_mz, n_dz = n_dz, var_a = var_a, 
                                  var_c = var_c, ge = TRUE, ge_beta0 = ge_beta0,
                                  ge_beta1 = ge_beta1, irt_model = irt_model,
                                  model = "ACE", n_items = n_items)
        data_mz = data$y_mz
        data_dz = data$y_dz
        
        ### Analyse data using masterfunction twinUniv
        print("Analysing data using an ACE model with GxE")
        output = twinUniv(data_mz = data_mz, data_dz = data_dz, 
                          twin1_datacols_p = c(1:n_items), twin2_datacols_p = c((n_items+1):(2*n_items)),
                          decomp_model = "ACE", irt_model = irt_model, ge = TRUE, n_iter = n_iter, 
                          n_burnin = n_burnin)
    } else {
        print("Simulating data for ACE model without GxE...")
        data = simulate_twin_data(n_mz = n_mz, n_dz = n_dz, var_a = var_a, 
                                 var_c = var_c, var_e = var_e, ge = FALSE,
                                 model = "ACE", n_items = n_items, irt_model = irt_model)
        data_mz = data$y_mz
        data_dz = data$y_dz
        
        ### Analyse data using masterfunction twinUniv
        print("Analysing data using an ACE model without GxE")
        output = twinUniv(data_mz = data_mz, data_dz = data_dz, 
                          twin1_datacols_p = c(1:n_items), twin2_datacols_p = c((n_items+1):(2*n_items)),
                          decomp_model = "ACE", irt_model = irt_model, ge = FALSE, n_iter = n_iter, 
                          n_burnin = n_burnin)
    }

    
    return(output)
}

sim_results  =  test_bayestwin_ACE(var_a, var_c, var_e, irt_model, ge, n_items,
                           ge_beta0, ge_beta1)

names(tryx$results)

###########################################################
## For loop to do no_sim - simulations
###########################################################
no_sim = 2
    
## Do NoSim- simulations: 
for (i in 1:no_sim){
    
    #Print number of simulation:    
    print(paste("################## SIMULATION NR", i, "####################"))
    
    #Create empty list / delete results of earlier analysis
    sim_results_list <- list()
        
    #Run the analysis and save results in sim_results list:  
    sim_results <- test_bayestwin_ACE(var_a = var_a, var_c = var_c, var_e = var_e, 
                                      , var_c, var_e = NA, irt_model, ge, n_items,
    ge_beta0 = NA, ge_beta1 = NA)
    
    ############### Organize output ###################################
    if (ge == TRUE && irt_model == "1PL"){
        
        #Mean values --------------------------------------------------
        #-variance components
        var_a = mean(1/sim_results$samples_var_a)
        var_c = mean(1/sim_results$samples_var_c)
        beta0 = mean(exp(sim_results$samples_beta0))
        beta1 = mean(sim_results$samples_beta1)
        
        #-item dificulties
        dim(sim_results$samples_item_b)
        n_items
        sim_results$samples_item_b
        ## SDs
        var_a_sd = sd(1/sim_results$samples_var_a)
        var_c_sd = sd(1/sim_results$samples_var_c)
        beta0_sd = sd(exp(sim_results$samples_beta0))
        beta1_sd = sd(sim_results$samples_beta1)
        
        names(sim_results)
    }
        
        ### Calculate mean values: 
        
        tryx$output
    
        
        #Save results in text file: 
        xy = sim_results
        
        fn<-"resultsMCAR_zesprocent.txt"
        if(file.exists(fn)) {
            write.table(xy, file=fn, row.names=F, col.names=F, append=T)
        } else {
            write.table(xy, file=fn, row.names=F, col.names=T)
        }
        
    }
    