# flu_model_LG

Edits to flu_model.

Use cluster_code_sep_COUNTRY_num.R to run the parameter inference for each epidemic one-by-one in one exemplar country (on HPC or own computer). This calls /fcns/inference_function.R using custom_inference_sep(), and 2_1b_model_epidemic_yearcross.R (from https://github.com/NaomiWaterlow/NextGen_Flu_Thai/blob/main/Vacc_epi_model/2_1b_model_epidemic_yearcross.R). Use cluster_code.R to run the joint parameter inference for all epidemics in one exemplar country (slower) using custom_inference().

Vaccination program data is downloaded from vaccine_inputs.R, and used in llikelihood (in /fcns/inference_function.R).



