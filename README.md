# flu_model_LG

flu_model, adapted from https://github.com/mbkoltai/flu_model/tree/main.

Use cluster_code_sep_COUNTRY_num.R to run the parameter inference for each epidemic one-by-one in one exemplar country (on HPC or own computer). This calls /fcns/inference_function.R using custom_inference_sep_bt(), and 2_1b_model_epidemic_yearcross.R (from https://github.com/NaomiWaterlow/NextGen_Flu_Thai/blob/main/Vacc_epi_model/2_1b_model_epidemic_yearcross.R). 

cluster_code_sep_COUNTRY_num.R uses the datasets in /data_for_cluster/. Epidemic identification, contact matrix contruction etc. uses datasets in /data/.

Vaccination data (efficacy, matching, exemplar country coverage) was produced in vaccine_inputs.R, and used in llikelihood (in /fcns/inference_function.R).

Files for running projections in each ITZ are in 'BS/' (bootstrapping). To carry out simulations:
  1. Download all data from /data_for_BS/ folder (100 simulations of the sampled 30-year period in each cluster, demographic parameters)
  2. Download /BS/ folder (in particular BS_data_fcns.R, BS_demography.R, BS_epid_simulations.R, BS_vaccine_doses.R, BS_vaccine_programs.R)
  3. Download /fcns/ files inference_function.R and 2_1b_model_epidemic_yearcross.R 
  4. Change c_number to the right ITZ number, scenario_name to the right sensitivity analysis or base case, and k to match intended coverage targets
  5. Run BS_epid_simulations.R for epidemic projections of weekly cases, stratified by strain, age group, and vaccination-status
  6. Run BS_vaccine_doses.R (changing c_number and scenario_name) to calculate annual vaccine doses given to each age group (unvaccinated and already-vaccinated individuals)

