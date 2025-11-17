#==============================================================================#
#                                                                              #
#                 Run CR-biologging model on simulated datasets                #
#                                                                              #
#==============================================================================#

library(nimble)

# The purpose of this script is to be launched on the cluster as many times as there 
# are simulated datasets.

path_input <- "03_outputs/simulations/simulated_datasets/polar_bear/no_biologging/"
path_output <- "03_outputs/simulations/model_outputs/polar_bear/"
# path_input <- "data_simulations/polar_bear_no_biologging/"
# path_output <- "res_simulations/polar_bear/"


# ~ 1. Select dataset file -----------------------------------------------------

# Make sure the model wait a randomly drawn period of time to make sure two model 
# runs don't happen on the same dataset
pause_secs <- sample(seq(from = 10, to = 60, by = 0.5), 1)
Sys.sleep(pause_secs)

# Retrieve input file names
file_list <- list.files(path_input)

# Pick a file at random 
set <- sample(file_list, 1)

# Load the file
load(paste0(path_input, set))

# Delete the file from the directory
file.remove(paste0(path_input, set))                                    


# ~ 2. Run the model -----------------------------------------------------------


# Provide initial values to the model and the parameters monitored
inits <- function() list(beta_phi = runif(n = 3, min = 0.7, max = 1),
                         
                         s0 = runif(n = 1, min = 0.25, max = 0.75), 
                         s1 = runif(n = 1, min = 0.5, max = 1), 
                         
                         beta_beta = runif(n = 4, min = 0.4, max = 1),
                         beta_gamma = runif(n = 3, min = 0.4, max = 1),
                         
                         p = runif(n = 1, min = 0, max = 0.5))
                         # zeta = zeta)

inits_values <- list(inits(), inits())

# Parameters monitored 
params <- c("beta_phi", "s0", "s1", "beta_beta", "beta_gamma", "p")

timestamps <- NULL
start <- Sys.time() 
# Create R model
inits_values <- inits()
CR_model_biologging_sim_R <- nimbleModel(code = CR_model_biologging_sim,
                                         constants = my.constants,
                                         data = dat,
                                         inits = inits_values,
                                         calculate = FALSE)

# Compile model (in C++)
CR_model_biologging_sim_C <- compileNimble(CR_model_biologging_sim_R,
                                           showCompilerOutput = FALSE)

# Configure MCMC
MCMC_conf <- configureMCMC(CR_model_biologging_sim_R,
                           monitors = params)

# Compile MCMC
MCMC_R <- buildMCMC(MCMC_conf)

MCMC_C <- compileNimble(MCMC_R, project = CR_model_biologging_sim_R)
timestamps[1] <- Sys.time() - start

start <- Sys.time() 
output <- runMCMC(mcmc = MCMC_C,
                  nchains = 2,
                  niter = 17000,
                  nburnin = 5000,
                  thin = 10,
                  inits = list(inits(), inits()),
                  setSeed = c(1, 2))

timestamps[2] <- Sys.time() - start
timestamps[3] <- timestamps[1] + timestamps[2]
names(timestamps) <- c("compilation-configuration", "MCMC", "total")

output[[3]] <- timestamps

assign(x = paste0("fit_", substr(set, 1, nchar(set)-6)),
       value = output)

save(list = paste0("fit_", substr(set, 1, nchar(set)-6)),
     file = paste0(path_output, "fit_", set))


# load("07_results/model_outputs/CR_GPS_GLS/simulations_revisions/polar_bear_no_GPS/fit_p_0.25_pGPS_1_pError_0_dataset_47.RData")
# output <- fit_p_0.25_pGPS_1_pError_0_dataset_47[1:2]
# 
# source("05_script/functions_for_models.R")
# 
# nrows <- length(params.plot)
# 
# diagnostic_plot <- check_convergence_detailed(params.plot = params.plot,
#                                               nimble_output = output)


