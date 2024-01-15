READ ME

# Bayesian Kriging Approaches for Spatial Functional Data


## 1. Simulation (Section 4)
## 1-1 simulation_SFoFR.R: conduct a SFoFR with simulation data
## 1-2 simulation_PSFoFR.R: conduct a PSFoFR with simulation data
## 1-3 Fdagstat_Simul.R: conduct a fdagstat package with simulation data
## 1-4 Estimated Regression function_Simul.R: draw estimated regression functions (psi) and 95% significant area
## 1-5 Draw Kriging_Simul.R: draw kriging results (SFoFR, PSFoFR, and fdagstat)
## 1-6 Draw Uncertainty_Simul.R: draw an uncertainty map with test sets 

After conducting 1-1),1-2),1-3), then conduct 1-4),1-5),1-6) with 1-1),1-2),1-3) results


## 2. PM2.5 (Section 5.1)
## 2-0 data
 - japan_coord_dat_align.RData: coordinate of locations
 - japan_no2_dat_align.RData: functional covariates (NOx)
 - japan_pm25_dat_align.RData: functional response (PM2.5)
 - japan_train_index.RData: train index
## 2-1 PM2.5_SFoFR.R: conduct a SFoFR with PM2.5 data
## 2-2 PM2.5_PSFoFR_rank5per.R: conduct a PSFoFR with PM2.5 data
## 2-3 Fdagstat_Japan.R: conduct a fdagstat package with PM2.5 data
## 2-4 Estimated Regression function_Japan.R: draw estimated regression functions (psi) and 95% significant area
## 2-5 Draw Kriging_Japan.R: draw kriging results (SFoFR, PSFoFR, and fdagstat)
## 2-6 Draw Uncertainty_Japan.R: draw an uncertainty map with test sets 

After conducting 2-1),2-2),2-3), then conduct 2-4),2-5),2-6) with 2-1),2-2),2-3) results

## 3. mobility (Section 5.2)
## 3-1 mobility_SFoFR.R: conduct a SFoFR with mobility data
## 3-2 mobility_PSFoFR.R: conduct a PSFoFR with mobility data
## 3-3 Estimated Regression function_Mobility.R: draw estimated regression functions (psi) and 95% significant area

After conducting 3-1),3-2) then conduct 3-3) with 3-1),3-2) results



