# combining_biologging_with_captures
Data and code for "A capture-recapture framework for combining biologging data with physical captures to decompose and estimate demographic rates: simulations across life cycles and application to polar bears" by Marwan Naciri, Jon Aars, Magnus Andersen, Andrew E. Derocher, Øystein Wiig, and Sarah Cubaynes 

[![DOI](https://zenodo.org/badge/1098295196.svg)](https://doi.org/10.5281/zenodo.17632764)

The files are organized into three main folders:

1. 01_inputs which contains all files required to run analyses on Svalbard polar bears. 

    - CR_biologging_events.csv: CR dataset for the model that makes use of biologging data 

    - CR_events.csv: CR dataset for the model that ignores biologging data 

    - sea_ice_metrics.csv: values of the sea-ice covariate 

    - two_year_old_cubs.csv: status (independent VS still with mother) of two-year-old cubs captured in spring. Used to compute the probability for a female who successfully raised two-year-old cubs to already be alone at capture.


2. 02_scripts which contains all the scripts to perform the analyses and plot the figures presented in the article and supporting information. This folder contains two files, and two folders

    - 01_simulations: a folder which contains the scripts to perfom the simulations study
    
    - 02_polar_bear_case_study: a folder which contains the script to apply the model to Svalbard polar bear data. Also contains a script to run a CR model that ignores biologging data on Svalbar polar bear data.
    
    - 03_plot_main_figures.R: to plot the figures presented in the main text. Requires the outputs from the simulations from the main example (see below), as well as the output of the model that makes use of biologging data applied to the Svalbar polar bear dataset (see below)
    
    - 04_plot_supplementary_figures.R: to plot the figures presented in the supplementary materials. Requires the outputs from all simulations (see below), the output of the model that makes use of biologging data applied to the Svalbard polar bear dataset, and the output of the model that does not make use of biologging data applied to the Svalbard polar bear dataset (see below).

3. 03_outputs which contains the outputs of the case study, as the models take >15h to run. Also contains a folder for the outptus of the simulations, but the outputs themeselves were not included to limit the size of the folder.


Software and package version:
R version 4.5.2
nimble version 1.3.0
tidyverse version 2.0.0
patchwork version 1.3.2
Rstudio version 2025.9.2.418
