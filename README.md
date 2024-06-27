# The potential global impact and cost-effectiveness of next-generation influenza vaccines: A modelling analysis.

flu_model, adapted from https://github.com/mbkoltai/flu_model/tree/main.

Data are mainly in ``` /data/ ```, ```data_for_cluster```, and ``` /data_for_BS/ ```.

## Files for FluNet data cleaning, epidemic identification, ITZ expansion:

In ```/prelim/``` and ```/fcns/```. Main files include:

- ``` itz_expansion.R ``` for geographical expansion of the ITZs by k-means clustering
- ```epid_identif_cont_matrs.R``` runs the epidemic identificatino algorithm and calculates the exemplar countries' contact matrices
- ```inference_function.R``` for the MCMC code

## Files for epidemic inference and viewing posterior parameters:

In ``` /epidemic_inference/ ```.

- ``` cluster_code_sep_CHN_1.R ``` runs inference in some epidemics in China. Can change parameters to run inference on different epidemcs and countries
- ``` viewing_sep_bt.R ``` runs plots of inference outputs

## Files for simulating future epidemics, including demography and vaccination models, and viewing outcomes:

In ``` /BS/ ```. Main files include:

- ``` BS_data_fcns.R ``` for the functions that run the simulations
- ``` BS_demography.R ``` for the demographic model
- ``` BS_vaccine_programs.R ``` for the vaccine scenarios
- ``` runARG.R ``` etc. for running the simulations in a given ITZ (here Southern America)
- ``` BS_projection_plots.R ``` for viewing simulation outputs

## Files for running and viewing the economic analysis:

In ```econ```. Main files:

- ```global_econ.R``` runs the economic analysis, including sensitivity analyses
- ```econ_plots.R``` produces figures of economic outcomes
- ```/outcome_calculations/``` has all input data for the economic model, and the files which calculate them when relevant

