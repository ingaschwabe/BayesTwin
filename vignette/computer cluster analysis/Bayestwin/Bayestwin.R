### Bayestwin.R
### File that includes the analysis 
library(Bayestwin)
simulated_data = simulatetwin(n_mz = 2000, n_dz = 5000, var_a = 0.5, var_c = 0.3, 
                              model = "ACE", n_items = 40, ge = TRUE,
                              ge_beta0 = log(0.2), ge_beta1 = 0,
                              irt_model = "1PL")

itemdata_mz = simulated_data$y_mz
itemdata_dz = simulated_data$y_dz

results = IRTtwin(data_mz = itemdata_mz, data_dz = itemdata_dz, twin1_datacols_p = 1:40, 
                  twin2_datacols_p = 41:80, decomp_model = "ACE", irt_model = "1PL", ge = TRUE, 
                  n_iter = 7000, n_burnin = 5000, n_chains = 1)
