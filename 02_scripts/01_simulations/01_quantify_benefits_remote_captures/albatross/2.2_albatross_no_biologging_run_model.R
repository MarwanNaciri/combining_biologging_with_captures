#==============================================================================#
#                                                                              #
#                 Run CR-biologging model on simulated datasets                #
#                                                                              #
#==============================================================================#

library(nimble)

# The purpose of this script is to be launched on the cluster as many times as there 
# are simulated datasets.

path_input <- "03_outputs/simulations/simulated_datasets/albatross/no_biologging/"
path_output <- "03_outputs/simulations/model_outputs/albatross/"
# path_input <- "data_simulations/albatross_no_biologging/"
# path_output <- "res_simulations/albatross/"



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
inits <- function() list(b_phi = runif(n = 3, min = 0, max = 1),
                         b_psi1 = runif(n = 4, min = 0, max = 1),
                         b_psi2 = runif(n = 2, min = 0, max = 1),
                         b_eta = runif(n = 5, min = 0, max = 1),
                         a1_beta = runif(n = 1, min = 0, max = 1),
                         p = runif(n = 1, min = 0, max = 1),
                         theta = runif(n = 1, min = 0, max = 1),
                         # z = Z) # 
                         zeta = zeta)

inits_values <- list(inits(), inits())

# Parameters monitored 
params <- c("b_phi", "b_psi1", "b_psi2", "b_eta", "a1_beta", "p", "theta") 

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

# Configure of MCMC
MCMC_conf <- configureMCMC(CR_model_biologging_sim_R,
                           monitors = params)

# Compile MCMC
MCMC_R <- buildMCMC(MCMC_conf)

MCMC_C <- compileNimble(MCMC_R, project = CR_model_biologging_sim_R)
timestamps[1] <- Sys.time() - start

start <- Sys.time() 
output <- runMCMC(mcmc = MCMC_C,
                    nchains = 2,
                    niter = 15000,
                    nburnin = 5000,
                    thin = 5,
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



