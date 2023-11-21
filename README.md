# flu_model_LG

flu_model, adapted from https://github.com/mbkoltai/flu_model/tree/main.

Use cluster_code_sep_COUNTRY_num.R to run the parameter inference for each epidemic one-by-one in one exemplar country (on HPC or own computer). This calls /fcns/inference_function.R using custom_inference_sep(), and 2_1b_model_epidemic_yearcross.R (from https://github.com/NaomiWaterlow/NextGen_Flu_Thai/blob/main/Vacc_epi_model/2_1b_model_epidemic_yearcross.R). Use cluster_code.R to run the joint parameter inference for all epidemics in one exemplar country (slower) using custom_inference().

cluster_code_sep_COUNTRY_num.R uses the datasets in /data_for_cluster/. Epidemic identification, contact matrix contruction etc. uses datasets in /data/.

Vaccination data (efficacy, matching, exemplar country coverage) was produced in vaccine_inputs.R, and used in llikelihood (in /fcns/inference_function.R).

Files for running projections in each ITZ are in 'BS/' (bootstrapping). To carry out simulations:
  1. Download all data from /data_for_BS/ folder, including the 100 simulations of the sampled 30-year period 
  2. Run BS_post_samples.R to merge the posterior samples
  3. Run BS_epid_simulations.R, changing 'c_number' to be the chosen ITZ (sources BS_data_fcns.R, BS_demography.R). Can run *some* of the 100 simulations on multiple cores for speed

