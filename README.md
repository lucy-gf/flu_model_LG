# The potential global impact and cost-effectiveness of next-generation influenza vaccines: A modelling analysis.

flu_model, adapted from https://github.com/mbkoltai/flu_model/tree/main.

Data are mainly in ``` /data/ ``` and ``` /data_for_BS/ ```.

## Files for FluNet data cleaning, epidemic identification, ITZ expansion:

In ```/prelim/```. Main files include:

- ``` itz_expansion.R ``` for geographical expansion of the ITZs by k-means clustering
- 

## Files for epidemic inference and viewing posterior parameters:

## Files for simulating future epidemics, including demography and vaccination models, and viewing outcomes:

In ``` /BS/ ```. Main files include 

- ``` BS_data_fcns.R ``` for the functions that run the simulations
- ``` BS_demography.R ``` for the demographic model
- ``` BS_vaccine_programs.R ``` for the vaccine scenarios
- ``` runARG.R ``` etc. for running the simulations in a given ITZ (here Southern America)
- ``` BS_projection_plots.R ``` for viewing simulation outputs


## Files for running and viewing the economic analysis:


