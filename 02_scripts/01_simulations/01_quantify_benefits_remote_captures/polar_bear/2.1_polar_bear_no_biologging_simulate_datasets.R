#==============================================================================#
#                                                                              #
#       Simulation of datasets for the CR-GPS-GLS model: no GPS & no GLS       #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(nimble)

# How many datasets?
n_sim <- 50
n_sim_start <- 1


start_global <- Sys.time() ; print(paste0("global start: ", start_global))

# A. Build the model ===========================================================

CR_model_biologging_sim <- nimbleCode({
  # -------------------------------------------------+
  # Parameters:
  # phi: juvenile, subadult and adult survival probability
  # eta: denning probability
  # beta: early litter survival
  # gamma: probability of a twin litter
  # s0: first-year cub survival 
  # s1: second-year cub survival probability
  # p: detection probability
  
  # -------------------------------------------------+
  # States (z):
  # 1 = J2: 2-yr old juvenile (parted from mother)
  # 2 = J3: 3-yr old juvenile
  # 3 = SA4: 4-yr old subadult
  # 4 = SA5: 5-yr old subadult 
  # 5 = A01: Female with 1 CoY
  # 6 = A02: Female with 2 CoYs
  # 7 = A11: Female with 1 yearling
  # 8 = A12: Female with 2 yearlings
  # 9 = AS1: Female who successfully raised 1 2-yr old
  # 10 = AS2: Female who successfully raised 2 2-yr old
  # 11 = LA: Lone Adult
  # 12 = D: Dead.
  
  # -------------------------------------------------+
  # Events (y):
  # 1 = J2: capture of a 2-yr old juvenile (parted from mother)
  # 2 = J3: capture of a 3-yr old juvenile
  # 3 = SA4: capture of a 4-yr old subadult
  # 4 = SA5: capture of a 5-yr old subadult (i.e. that does not have cubs at capture)
  # 5 = A01: capture of a female with 1 CoY
  # 6 = A02: capture of a female with 2 CoYs
  # 7 = A11: capture of a female with 1 yearling
  # 8 = A12: capture of a female with 2 yearlings
  # 9 = A21: capture of a female who successfully raised 1 2-yr old
  # 10 = A22: capture of a female who successfully raised 2 2-yr old
  # 11 = A: capture of a lone female
  # 12 = NA: Not captured 
  # -------------------------------------------------+
  
  # +++++++++++++++++++++++++ transition matrices ++++++++++++++++++++++++++++++
  
  
  
  # ~~~ a. E0 ----------------------------------
  E0[1, 1] <- 1
  E0[1, 2] <- 0
  E0[1, 3] <- 0
  E0[1, 4] <- 0
  E0[1, 5] <- 0
  E0[1, 6] <- 0
  E0[1, 7] <- 0
  E0[1, 8] <- 0
  E0[1, 9] <- 0
  E0[1, 10] <- 0
  E0[1, 11] <- 0
  E0[1, 12] <- 0
  
  E0[2, 1] <- 0
  E0[2, 2] <- 1
  E0[2, 3] <- 0
  E0[2, 4] <- 0
  E0[2, 5] <- 0
  E0[2, 6] <- 0
  E0[2, 7] <- 0
  E0[2, 8] <- 0
  E0[2, 9] <- 0
  E0[2, 10] <- 0
  E0[2, 11] <- 0
  E0[2, 12] <- 0
  
  E0[3, 1] <- 0
  E0[3, 2] <- 0
  E0[3, 3] <- 1
  E0[3, 4] <- 0
  E0[3, 5] <- 0
  E0[3, 6] <- 0
  E0[3, 7] <- 0
  E0[3, 8] <- 0
  E0[3, 9] <- 0
  E0[3, 10] <- 0
  E0[3, 11] <- 0
  E0[3, 12] <- 0
  
  E0[4, 1] <- 0
  E0[4, 2] <- 0
  E0[4, 3] <- 0
  E0[4, 4] <- 1
  E0[4, 5] <- 0
  E0[4, 6] <- 0
  E0[4, 7] <- 0
  E0[4, 8] <- 0
  E0[4, 9] <- 0
  E0[4, 10] <- 0
  E0[4, 11] <- 0
  E0[4, 12] <- 0
  
  E0[5, 1] <- 0
  E0[5, 2] <- 0
  E0[5, 3] <- 0
  E0[5, 4] <- 0
  E0[5, 5] <- 1
  E0[5, 6] <- 0
  E0[5, 7] <- 0
  E0[5, 8] <- 0
  E0[5, 9] <- 0
  E0[5, 10] <- 0
  E0[5, 11] <- 0
  E0[5, 12] <- 0
  
  E0[6, 1] <- 0
  E0[6, 2] <- 0
  E0[6, 3] <- 0
  E0[6, 4] <- 0
  E0[6, 5] <- 0
  E0[6, 6] <- 1
  E0[6, 7] <- 0
  E0[6, 8] <- 0
  E0[6, 9] <- 0
  E0[6, 10] <- 0
  E0[6, 11] <- 0
  E0[6, 12] <- 0
  
  E0[7, 1] <- 0
  E0[7, 2] <- 0
  E0[7, 3] <- 0
  E0[7, 4] <- 0
  E0[7, 5] <- 0
  E0[7, 6] <- 0
  E0[7, 7] <- 1
  E0[7, 8] <- 0
  E0[7, 9] <- 0
  E0[7, 10] <- 0
  E0[7, 11] <- 0
  E0[7, 12] <- 0
  
  E0[8, 1] <- 0
  E0[8, 2] <- 0
  E0[8, 3] <- 0
  E0[8, 4] <- 0
  E0[8, 5] <- 0
  E0[8, 6] <- 0
  E0[8, 7] <- 0
  E0[8, 8] <- 1
  E0[8, 9] <- 0
  E0[8, 10] <- 0
  E0[8, 11] <- 0
  E0[8, 12] <- 0
  
  E0[9, 1] <- 0
  E0[9, 2] <- 0
  E0[9, 3] <- 0
  E0[9, 4] <- 0
  E0[9, 5] <- 0
  E0[9, 6] <- 0
  E0[9, 7] <- 0
  E0[9, 8] <- 0
  E0[9, 9] <- 1-a
  E0[9, 10] <- 0
  E0[9, 11] <- a
  E0[9, 12] <- 0
  
  E0[10, 1] <- 0
  E0[10, 2] <- 0
  E0[10, 3] <- 0
  E0[10, 4] <- 0
  E0[10, 5] <- 0
  E0[10, 6] <- 0
  E0[10, 7] <- 0
  E0[10, 8] <- 0
  E0[10, 9] <- 2 * ( 1-a ) * a
  E0[10, 10] <- ( 1 - a )^2
  E0[10, 11] <- a^2
  E0[10, 12] <- 0
  
  E0[11, 1] <- 0
  E0[11, 2] <- 0
  E0[11, 3] <- 0
  E0[11, 4] <- 0
  E0[11, 5] <- 0
  E0[11, 6] <- 0
  E0[11, 7] <- 0
  E0[11, 8] <- 0
  E0[11, 9] <- 0
  E0[11, 10] <- 0
  E0[11, 11] <- 1
  E0[11, 12] <- 0
  
  E0[12, 1] <- 0
  E0[12, 2] <- 0
  E0[12, 3] <- 0
  E0[12, 4] <- 0
  E0[12, 5] <- 0
  E0[12, 6] <- 0
  E0[12, 7] <- 0
  E0[12, 8] <- 0
  E0[12, 9] <- 0
  E0[12, 10] <- 0
  E0[12, 11] <- 0
  E0[12, 12] <- 1
  
  # ~~~ b. E ----------------------------------
  
  E[1, 1] <- p
  E[1, 2] <- 0
  E[1, 3] <- 0
  E[1, 4] <- 0
  E[1, 5] <- 0
  E[1, 6] <- 0
  E[1, 7] <- 0
  E[1, 8] <- 0
  E[1, 9] <- 0
  E[1, 10] <- 0
  E[1, 11] <- 0
  E[1, 12] <- 1-p
  
  E[2, 1] <- 0
  E[2, 2] <- p
  E[2, 3] <- 0
  E[2, 4] <- 0
  E[2, 5] <- 0
  E[2, 6] <- 0
  E[2, 7] <- 0
  E[2, 8] <- 0
  E[2, 9] <- 0
  E[2, 10] <- 0
  E[2, 11] <- 0
  E[2, 12] <- 1 - p
  
  E[3, 1] <- 0
  E[3, 2] <- 0
  E[3, 3] <- p
  E[3, 4] <- 0
  E[3, 5] <- 0
  E[3, 6] <- 0
  E[3, 7] <- 0
  E[3, 8] <- 0
  E[3, 9] <- 0
  E[3, 10] <- 0
  E[3, 11] <- 0
  E[3, 12] <- 1 - p
  
  E[4, 1] <- 0
  E[4, 2] <- 0
  E[4, 3] <- 0
  E[4, 4] <- p
  E[4, 5] <- 0
  E[4, 6] <- 0
  E[4, 7] <- 0
  E[4, 8] <- 0
  E[4, 9] <- 0
  E[4, 10] <- 0
  E[4, 11] <- 0
  E[4, 12] <- 1 - p
  
  E[5, 1] <- 0
  E[5, 2] <- 0
  E[5, 3] <- 0
  E[5, 4] <- 0
  E[5, 5] <- p
  E[5, 6] <- 0
  E[5, 7] <- 0
  E[5, 8] <- 0
  E[5, 9] <- 0
  E[5, 10] <- 0
  E[5, 11] <- 0
  E[5, 12] <- 1 - p
  
  E[6, 1] <- 0
  E[6, 2] <- 0
  E[6, 3] <- 0
  E[6, 4] <- 0
  E[6, 5] <- 0
  E[6, 6] <- p
  E[6, 7] <- 0
  E[6, 8] <- 0
  E[6, 9] <- 0
  E[6, 10] <- 0
  E[6, 11] <- 0
  E[6, 12] <- 1 - p
  
  E[7, 1] <- 0
  E[7, 2] <- 0
  E[7, 3] <- 0
  E[7, 4] <- 0
  E[7, 5] <- 0
  E[7, 6] <- 0
  E[7, 7] <- p
  E[7, 8] <- 0
  E[7, 9] <- 0
  E[7, 10] <- 0
  E[7, 11] <- 0
  E[7, 12] <- 1 - p
  
  E[8, 1] <- 0
  E[8, 2] <- 0
  E[8, 3] <- 0
  E[8, 4] <- 0
  E[8, 5] <- 0
  E[8, 6] <- 0
  E[8, 7] <- 0
  E[8, 8] <- p
  E[8, 9] <- 0
  E[8, 10] <- 0
  E[8, 11] <- 0
  E[8, 12] <- 1 - p
  
  E[9, 1] <- 0
  E[9, 2] <- 0
  E[9, 3] <- 0
  E[9, 4] <- 0
  E[9, 5] <- 0
  E[9, 6] <- 0
  E[9, 7] <- 0
  E[9, 8] <- 0
  E[9, 9] <- (1 - a) * p
  E[9, 10] <- 0
  E[9, 11] <- a * p
  E[9, 12] <- 1 - p
  
  E[10, 1] <- 0
  E[10, 2] <- 0
  E[10, 3] <- 0
  E[10, 4] <- 0
  E[10, 5] <- 0
  E[10, 6] <- 0
  E[10, 7] <- 0
  E[10, 8] <- 0
  E[10, 9] <- 2 * (1 - a) * a * p
  E[10, 10] <- (1 - a)^2 * p
  E[10, 11] <- a^2 * p
  E[10, 12] <- 1 - p
  
  E[11, 1] <- 0
  E[11, 2] <- 0
  E[11, 3] <- 0
  E[11, 4] <- 0
  E[11, 5] <- 0
  E[11, 6] <- 0
  E[11, 7] <- 0
  E[11, 8] <- 0
  E[11, 9] <- 0
  E[11, 10] <- 0
  E[11, 11] <- p
  E[11, 12] <- 1 - p
  
  E[12, 1] <- 0
  E[12, 2] <- 0
  E[12, 3] <- 0
  E[12, 4] <- 0
  E[12, 5] <- 0
  E[12, 6] <- 0
  E[12, 7] <- 0
  E[12, 8] <- 0
  E[12, 9] <- 0
  E[12, 10] <- 0
  E[12, 11] <- 0
  E[12, 12] <- 1
  
  # ~~~ c. Regressions on parameters -------------------
  
  for (i in 1:n_ind) {
    # if (f[i] == n_year) next
    for (t in f[i]:(n_occasions - 1)) {
      
      # Probability of juvenile, subadult and adult survival
      phi[i, t] <- beta_phi[1] * AGE_2_4[i, t] +
        beta_phi[2] * (AGE_YOUNG[i, t] + AGE_MIDDLE[i, t]) +
        beta_phi[3] * AGE_OLD[i, t]
      
      # Probability of breeding
      beta[i, t] <- beta_beta[1] * AGE_4[i, t] +
        beta_beta[2] * AGE_YOUNG[i, t] +
        beta_beta[3] * AGE_MIDDLE[i, t] +
        beta_beta[4] * AGE_OLD[i, t]
      
      # litter size
      gamma[i, t] <- beta_gamma[1] * (AGE_4[i, t] + AGE_YOUNG[i, t]) +
        beta_gamma[2] * AGE_MIDDLE[i, t] +
        beta_gamma[3] * AGE_OLD[i, t]
      
      # ~~~ d. S ---------------------------------------------------------------
      
      S[1, 1, i, t] <- 0
      S[1, 2, i, t] <- phi[i, t]
      S[1, 3, i, t] <- 0
      S[1, 4, i, t] <- 0
      S[1, 5, i, t] <- 0
      S[1, 6, i, t] <- 0
      S[1, 7, i, t] <- 0
      S[1, 8, i, t] <- 0
      S[1, 9, i, t] <- 0
      S[1, 10, i, t] <- 0
      S[1, 11, i, t] <- 0
      S[1, 12, i, t] <- 1-phi[i, t]
      
      S[2, 1, i, t] <- 0
      S[2, 2, i, t] <- 0
      S[2, 3, i, t] <- phi[i, t]
      S[2, 4, i, t] <- 0
      S[2, 5, i, t] <- 0
      S[2, 6, i, t] <- 0
      S[2, 7, i, t] <- 0
      S[2, 8, i, t] <- 0
      S[2, 9, i, t] <- 0
      S[2, 10, i, t] <- 0
      S[2, 11, i, t] <- 0
      S[2, 12, i, t] <- 1-phi[i, t]
      
      S[3, 1, i, t] <- 0
      S[3, 2, i, t] <- 0
      S[3, 3, i, t] <- 0
      S[3, 4, i, t] <- phi[i, t] * (1- beta[i, t])
      S[3, 5, i, t] <- phi[i, t] * beta[i, t] * (1 - gamma[i, t])
      S[3, 6, i, t] <- phi[i, t] * beta[i, t] * gamma[i, t]
      S[3, 7, i, t] <- 0
      S[3, 8, i, t] <- 0
      S[3, 9, i, t] <- 0
      S[3, 10, i, t] <- 0
      S[3, 11, i, t] <- 0
      S[3, 12, i, t] <- 1-phi[i, t]
      
      S[4, 1, i, t] <- 0
      S[4, 2, i, t] <- 0
      S[4, 3, i, t] <- 0
      S[4, 4, i, t] <- 0
      S[4, 5, i, t] <- phi[i, t] * beta[i, t] * (1 - gamma[i, t])
      S[4, 6, i, t] <- phi[i, t] * beta[i, t] * gamma[i, t]
      S[4, 7, i, t] <- 0
      S[4, 8, i, t] <- 0
      S[4, 9, i, t] <- 0
      S[4, 10, i, t] <- 0
      S[4, 11, i, t] <- phi[i, t] * (1 - beta[i, t])
      S[4, 12, i, t] <- 1-phi[i, t]
      
      S[5, 1, i, t] <- 0
      S[5, 2, i, t] <- 0
      S[5, 3, i, t] <- 0
      S[5, 4, i, t] <- 0
      S[5, 5, i, t] <- phi[i, t] * (1-s0) * beta[i, t] * (1 - gamma[i, t])
      S[5, 6, i, t] <- phi[i, t] * (1-s0) * beta[i, t] * gamma[i, t]
      S[5, 7, i, t] <- phi[i, t] * s0
      S[5, 8, i, t] <- 0
      S[5, 9, i, t] <- 0
      S[5, 10, i, t] <- 0
      S[5, 11, i, t] <- phi[i, t] * (1-s0) * (1 - beta[i, t])
      S[5, 12, i, t] <- 1-phi[i, t]
      
      S[6, 1, i, t] <- 0
      S[6, 2, i, t] <- 0
      S[6, 3, i, t] <- 0
      S[6, 4, i, t] <- 0
      S[6, 5, i, t] <- phi[i, t] * (1-s0)^2 * beta[i, t] * (1 - gamma[i, t])
      S[6, 6, i, t] <- phi[i, t] * (1-s0)^2 * beta[i, t] * gamma[i, t]
      S[6, 7, i, t] <- phi[i, t] * 2 * s0 * (1-s0)
      S[6, 8, i, t] <- phi[i, t] * s0^2
      S[6, 9, i, t] <- 0
      S[6, 10, i, t] <- 0
      S[6, 11, i, t] <- phi[i, t] * (1-s0)^2 * (1 - beta[i, t])
      S[6, 12, i, t] <- 1-phi[i, t]
      
      S[7, 1, i, t] <- 0
      S[7, 2, i, t] <- 0
      S[7, 3, i, t] <- 0
      S[7, 4, i, t] <- 0
      S[7, 5, i, t] <- phi[i, t] * (1-s1) * beta[i, t] * (1 - gamma[i, t])
      S[7, 6, i, t] <- phi[i, t] * (1-s1) * beta[i, t] * gamma[i, t]
      S[7, 7, i, t] <- 0
      S[7, 8, i, t] <- 0
      S[7, 9, i, t] <- phi[i, t] * s1
      S[7, 10, i, t] <- 0
      S[7, 11, i, t] <- phi[i, t] * (1-s1) * (1 - beta[i, t])
      S[7, 12, i, t] <- 1-phi[i, t]
      
      S[8, 1, i, t] <- 0
      S[8, 2, i, t] <- 0
      S[8, 3, i, t] <- 0
      S[8, 4, i, t] <- 0
      S[8, 5, i, t] <- phi[i, t] * (1-s1)^2 * beta[i, t] * (1 - gamma[i, t])
      S[8, 6, i, t] <- phi[i, t] * (1-s1)^2 * beta[i, t] * gamma[i, t]
      S[8, 7, i, t] <- 0
      S[8, 8, i, t] <- 0
      S[8, 9, i, t] <- phi[i, t] * 2 * s1 * (1-s1)
      S[8, 10, i, t] <- phi[i, t] * s1^2
      S[8, 11, i, t] <- phi[i, t] * (1-s1)^2 * (1 - beta[i, t])
      S[8, 12, i, t] <- 1-phi[i, t]
      
      S[9, 1, i, t] <- 0
      S[9, 2, i, t] <- 0
      S[9, 3, i, t] <- 0
      S[9, 4, i, t] <- 0
      S[9, 5, i, t] <- phi[i, t] * beta[i, t] * (1 - gamma[i, t])
      S[9, 6, i, t] <- phi[i, t] * beta[i, t] * gamma[i, t]
      S[9, 7, i, t] <- 0
      S[9, 8, i, t] <- 0
      S[9, 9, i, t] <- 0
      S[9, 10, i, t] <- 0
      S[9, 11, i, t] <- phi[i, t] * (1 - beta[i, t])
      S[9, 12, i, t] <- 1-phi[i, t]
      
      S[10, 1, i, t] <- 0
      S[10, 2, i, t] <- 0
      S[10, 3, i, t] <- 0
      S[10, 4, i, t] <- 0
      S[10, 5, i, t] <- phi[i, t] * beta[i, t] * (1 - gamma[i, t])
      S[10, 6, i, t] <- phi[i, t] * beta[i, t] * gamma[i, t]
      S[10, 7, i, t] <- 0
      S[10, 8, i, t] <- 0
      S[10, 9, i, t] <- 0
      S[10, 10, i, t] <- 0
      S[10, 11, i, t] <- phi[i, t] * (1 - beta[i, t])
      S[10, 12, i, t] <- 1-phi[i, t]
      
      S[11, 1, i, t] <- 0
      S[11, 2, i, t] <- 0
      S[11, 3, i, t] <- 0
      S[11, 4, i, t] <- 0
      S[11, 5, i, t] <- phi[i, t] * beta[i, t] * (1 - gamma[i, t])
      S[11, 6, i, t] <- phi[i, t] * beta[i, t] * gamma[i, t]
      S[11, 7, i, t] <- 0
      S[11, 8, i, t] <- 0
      S[11, 9, i, t] <- 0
      S[11, 10, i, t] <- 0
      S[11, 11, i, t] <- phi[i, t] * (1 - beta[i, t])
      S[11, 12, i, t] <- 1-phi[i, t]
      
      S[12, 1, i, t] <- 0
      S[12, 2, i, t] <- 0
      S[12, 3, i, t] <- 0
      S[12, 4, i, t] <- 0
      S[12, 5, i, t] <- 0
      S[12, 6, i, t] <- 0
      S[12, 7, i, t] <- 0
      S[12, 8, i, t] <- 0
      S[12, 9, i, t] <- 0
      S[12, 10, i, t] <- 0
      S[12, 11, i, t] <- 0
      S[12, 12, i, t] <- 1
      
    }
  }
  
  # ++++++++++++++++++++++++++++++++ priors ++++++++++++++++++++++++++++++++++++
  
  #-------- State process --------#
  for (k in 1:3) {
    beta_phi[k] ~ dunif(0, 1)   # juvenile, subadult and adult survival
    beta_gamma[k] ~ dunif(0, 1) # litter size
  }
  s0 ~ dunif(0, 1)   # cub survival
  s1 ~ dunif(0, 1)   # yearling survival
  
  for (k in 1:4) {
    beta_beta[k] ~ dunif(0, 1)  # breeding probability
  }
  
  #----- observation process -----#
  p ~ dunif(0, 1)
  
  
  # +++++++++++++++++++++++++++ model & likelihood +++++++++++++++++++++++++++++
  
  # for (i in 1:n_ind){
  #   # latent state at first capture
  #   z[i, f[i]] ~ dcat(S0[1:12])
  #   y[i, f[i]] ~ dcat(E0[z[i, f[i]], 1:12, i])
  #   
  #   for (t in ( f[i]+1 ):n_occasions) {
  #     
  #     #-------- State process --------#
  #     # draw z(t) given z(t-1)
  #     z[i, t] ~ dcat( S[z[i, t-1], 1:12, i, t-1] )
  #     
  #     #----- Observation process -----#
  #     # draw y(t) given z(t)
  #     y[i, t] ~ dcat( E[z[i, t], 1:12, i, t-1] )
  #   }
  # }
  
  for(i in 1:n_ind){
    
    # First capture
    for (k in 1:n_states) {
      zeta[i, f[i], k] <- S0[k] * E0[k, CH[i, f[i]]]
    }
    
    # Subsequent captures
    for(t in f[i] : (n_occasions-1) ){
      for (k in 1:n_states) {
        zeta[i, (t+1), k] <- inprod(zeta[i, t, 1:n_states],
                                    S[1:n_states, k, i, t]) *
          E[k, CH[i, t+1]]
      }
    }
    lik[i] <- sum(zeta[i, n_occasions, 1:n_states])
    ones[i] ~ dbin(prob = lik[i], size = freq[i])
  }
})





# B. Generate datasets =========================================================

# ~ 1. Define parameters -------------------------------------------------------

{
  new_per_year <- 20
  n_occasions <- 8
  
  n_states <- 12
  n_states_simul <- 14
  n_events <- 12
  
  # Survival
  beta_phi <- NULL
  beta_phi[1] <- 0.90    # intercept 2-4yr 
  beta_phi[2] <- 0.96    # intercept 5-15yr 
  beta_phi[3] <- 0.85    # intercept >15yr 
  
  s0 <- 0.5  
  s1 <- 0.75  
  
  # Age dependent reproductive rates 
  beta_eta <- beta_beta <- beta_gamma <-  NULL
  
  # eta
  beta_eta[1] <- 0.25    # intercept SA4 
  beta_eta[2] <- 0.65  # intercept young 
  beta_eta[3] <- 0.75  # intercept middle aged 
  beta_eta[4] <- 0.70  # intercept old  
  
  # beta
  beta_beta[1] <- 0.65   # intercept 4-9yr 
  beta_beta[2] <- 0.76   # intercept middle aged 
  beta_beta[3] <- 0.67   # intercept old  
  # gamma
  beta_gamma[1] <- 0.65  # intercept 4-9yr 
  beta_gamma[2] <- 0.58  # intercept middle aged 
  beta_gamma[3] <- 0.55  # intercept old  
  
  # For females who successfully raised 2-yr-old cubs, the probability of having parted
  a <- 0.5
  
  f <- rep(c(1:n_occasions), each = new_per_year)
  
  # Recapture probability 
  values_p <- 0.05 # c(0.05, 0.10, 0.25, 0.5)
  
  # Probability of being equipped with a GPS at capture
  values_pGPS <- c(0)
  
  # Duration of collars,
  lambda_GPS <- 0.70
  
  # Observed proportions of each state (state is inferred by me) at first capture
  #                   1      2    3     4     5     6      7     8     9     10    11     12    13   14 
  #                  J2     J3   SA4  NDSA5 FDSA5  A01    A02   A11   A12   AS1   AS2    LNDA   FDA   D 
  f_state_prop <- c(0.07, 0.07, 0.07, 0.08, 0.04, 0.080, 0.150, 0.06, 0.03, 0.02, 0.01, 0.170, 0.150, 0)
  sum(f_state_prop)
  
  # Proportion of each age for each state at first capture (matrix with states in
  # lines and age in columns).
  n_trunc <- 10000  # Need a big n so that all values are represented
  AGE_f <- matrix(data = c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  # J2
                           0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  # J3
                           0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  # SA4
                           0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  # NDSA5
                           0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  # FDSA5
                           0, 0, 0, unname(table(round(truncnorm::rtruncnorm(n = n_trunc, a = 4, b = 24, mean = 10, sd = 5)))/n_trunc),     # A01
                           0, 0, 0, unname(table(round(truncnorm::rtruncnorm(n = n_trunc, a = 4, b = 24, mean = 10, sd = 5)))/n_trunc),     # A02
                           0, 0, 0, 0, unname(table(round(truncnorm::rtruncnorm(n = n_trunc, a = 5, b = 24, mean = 11, sd = 5)))/n_trunc),     # A11
                           0, 0, 0, 0, unname(table(round(truncnorm::rtruncnorm(n = n_trunc, a = 5, b = 24, mean = 11, sd = 5)))/n_trunc),     # A12                         
                           0, 0, 0, 0, 0, unname(table(round(truncnorm::rtruncnorm(n = n_trunc, a = 6, b = 24, mean = 12, sd = 5)))/n_trunc),  # AS1
                           0, 0, 0, 0, 0, unname(table(round(truncnorm::rtruncnorm(n = n_trunc, a = 6, b = 24, mean = 12, sd = 5)))/n_trunc),  # AS2
                           0, 0, 0, 0, unname(table(round(truncnorm::rtruncnorm(n = n_trunc, a = 5, b = 24, mean = 13, sd = 5)))/n_trunc),  # LNDA
                           0, 0, 0, 0, unname(table(round(truncnorm::rtruncnorm(n = n_trunc, a = 5, b = 24, mean = 13, sd = 5)))/n_trunc),  # FDA
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),  # D
                  nrow = n_states_simul, ncol = 24, byrow = TRUE)
  apply(AGE_f, 1, sum)
}



# ~ 2. Generate capture histories ----------------------------------------------

for (k in 1:length(values_p)) {
  for (l in 1:length(values_pGPS)) {
    
    #                                     no               yes
    COLLARING <- matrix(data = c(         1,                0,          # J2    
                                          1,                0,          # J3   
                                          1 - values_pGPS[l],  values_pGPS[l],   # SA4  
                                          1 - values_pGPS[l],  values_pGPS[l],   # NDSA5  
                                          1 - values_pGPS[l],  values_pGPS[l],   # FDSA5   
                                          1 - values_pGPS[l],  values_pGPS[l],   # A01  
                                          1 - values_pGPS[l],  values_pGPS[l],   # A02  
                                          1 - values_pGPS[l],  values_pGPS[l],   # A11
                                          1 - values_pGPS[l],  values_pGPS[l],   # A12
                                          1 - values_pGPS[l],  values_pGPS[l],   # AS1  
                                          1 - values_pGPS[l],  values_pGPS[l],   # AS2  
                                          1 - values_pGPS[l],  values_pGPS[l],   # LNDA    
                                          1 - values_pGPS[l],  values_pGPS[l],   # FDA
                                          1,             0),         # D
                        nrow = n_states_simul, ncol = 2, byrow = T)
    
    for (n in n_sim_start:n_sim) {
      print(paste0("simulation ", n, ":"))
      
      # Arrays to store the individual- and time-specific transition and observation matrices
      S <- array(data = NA, dim = c(n_states_simul, n_states_simul, (new_per_year*n_occasions - new_per_year), n_occasions))
      E <- E_p <- array(data = NA, dim = c(n_states_simul, n_events, (new_per_year*n_occasions - new_per_year), n_occasions))
      
      # Matrices for vital rates
      phi <- beta <- gamma <- eta <- p_realized <- matrix(NA, nrow = (new_per_year*n_occasions - new_per_year), 
                                                          ncol = n_occasions)
      # Matrices for data
      CH <- Z <- AGE <- matrix(NA, nrow = (new_per_year*n_occasions - new_per_year), 
                               ncol = n_occasions) 
      
      KHI <- matrix(0, nrow = (new_per_year*n_occasions - new_per_year), 
                    ncol = n_occasions) 
      f_event <- f_state <-  NULL
      
      # ~~~ a. First captures --------------------------------------------------------
      
      set.seed(n)
      for (i in 1:(new_per_year*n_occasions - new_per_year)) {
        Z[i, f[i]] <- f_state[i] <- sample(x = 1:n_states_simul,
                                           size = 1,
                                           prob = f_state_prop, replace = T)
        
        # Age at first capture
        AGE[i, f[i]] <- which(rmultinom(n = 1, size = 1,
                                        prob = AGE_f[Z[i, f[i]], ]) == 1) + 1 # +1 because J2 (state = 1) are 2yr old
        
        #               |------------------- physical only --------------------|| None |
        #               J2  J3  SA4 SA5 A01 A02 A11 A12    A21        A22    A     NA
        #                                5                             10                  
        E0 <- matrix(data =
                       c(1,  0,  0,  0,  0,  0,  0,  0,     0,         0,    0,     0,  # J2
                         0,  1,  0,  0,  0,  0,  0,  0,     0,         0,    0,     0,  # J3
                         0,  0,  1,  0,  0,  0,  0,  0,     0,         0,    0,     0,  # SA4
                         0,  0,  0,  1,  0,  0,  0,  0,     0,         0,    0,     0,  # NDSA5
                         0,  0,  0,  1,  0,  0,  0,  0,     0,         0,    0,     0,  # FDSA5  # 5
                         0,  0,  0,  0,  1,  0,  0,  0,     0,         0,    0,     0,  # A01
                         0,  0,  0,  0,  0,  1,  0,  0,     0,         0,    0,     0,  # A02
                         0,  0,  0,  0,  0,  0,  1,  0,     0,         0,    0,     0,  # A11    
                         0,  0,  0,  0,  0,  0,  0,  1,     0,         0,    0,     0,  # A12
                         0,  0,  0,  0,  0,  0,  0,  0,   (1-a),       0,    a,     0,  # AS1    # 10
                         0,  0,  0,  0,  0,  0,  0,  0, 2*(1-a)*a, (1-a)^2, a^2,    0,  # AS2
                         0,  0,  0,  0,  0,  0,  0,  0,     0,         0,    1,     0,  # LNDA
                         0,  0,  0,  0,  0,  0,  0,  0,     0,         0,    1,     0,  # FDA    
                         0,  0,  0,  0,  0,  0,  0,  0,     0,         0,    0,     1),  # D
                     nrow = 14, ncol = 12, byrow = TRUE)
        apply(E0, 1, sum)
        
        CH[i, f[i]] <- f_event[i] <- which(rmultinom(n = 1, size = 1, 
                                                     prob = E0[Z[i, f[i]], ]) == 1)
        # GPS collaring
        collared <- which(rmultinom(n = 1, size = 1, 
                                    prob = COLLARING[Z[i, f[i]], ]) == 1) - 1
        
        if (collared == 1) { # If individual has a collar, for how many years ?
          collared_years <- rpois(n = 1, lambda = lambda_GPS) + 1
          q <- 1
          while (q <= collared_years &  
                 (f[i] + q) <= n_occasions) {
            KHI[i, f[i] + q] <- 1
            q <- q + 1
          }
        }
        
        # Compute matrices with age
        AGE_4 <- AGE_5 <- AGE_YOUNG <- AGE_MIDDLE <- AGE_OLD <- AGE_2_4 <- 
          matrix(NA, nrow = (new_per_year*n_occasions - new_per_year), ncol = n_occasions) 
        for (i in 1:((new_per_year*n_occasions - new_per_year))) {
          AGE_4[i, f[i]] <- ifelse(AGE[i, f[i]] == 4, 1, 0)
          AGE_5[i, f[i]] <- ifelse(AGE[i, f[i]] == 5, 1, 0)
          AGE_YOUNG[i, f[i]] <- ifelse(AGE[i, f[i]] %in% 5:9, 1, 0)
          AGE_MIDDLE[i, f[i]] <- ifelse(AGE[i, f[i]] %in% 10:15, 1, 0)
          AGE_OLD[i, f[i]] <- ifelse(AGE[i, f[i]] %in% 16:50, 1, 0)
          AGE_2_4[i, f[i]] <- ifelse(AGE[i, f[i]] %in% 2:4, 1, 0)
          
          if (f[i] == n_occasions) next
          
          for (t in (f[i]+1):n_occasions) {
            AGE[i, t] <- AGE[i, t-1] + 1 
            AGE_4[i, t] <- ifelse(AGE[i, t] == 4, 1, 0)
            AGE_5[i, t] <- ifelse(AGE[i, t] == 5, 1, 0)
            AGE_YOUNG[i, t] <- ifelse(AGE[i, t] %in% 5:9, 1, 0)
            AGE_MIDDLE[i, t] <- ifelse(AGE[i, t] %in% 10:15, 1, 0)
            AGE_OLD[i, t] <- ifelse(AGE[i, t] %in% 16:50, 1, 0)
            AGE_2_4[i, t] <- ifelse(AGE[i, t] %in% 2:4, 1, 0)
          }
        }
      }  
      
      
      # ~~~ b. Subsequent captures ---------------------------------------------------
      
      
      set.seed(l)
      for (i in 1:(new_per_year*n_occasions - new_per_year)) {
        if (f[i] == n_occasions) next
        
        for (t in (f[i]:(n_occasions-1))) {
          # Here, note that at each time step t, I calculate state at t+1. Otherwise,
          # I would have to specify [i, t-1] instead of [i, t] when computing age-dependent
          # vital rates.
          
          phi[i, t] <- beta_phi[1] * AGE_2_4[i, t] +
            beta_phi[2] * (AGE_YOUNG[i, t] + AGE_MIDDLE[i, t]) +
            beta_phi[3] * AGE_OLD[i, t]
          
          eta[i, t] <- beta_eta[1] * AGE_4[i, t] +        # eta for lone females
            beta_eta[2] * AGE_YOUNG[i, t] +
            beta_eta[3] * AGE_MIDDLE[i, t] +
            beta_eta[4] * AGE_OLD[i, t]
          
          # Breeding success
          beta[i, t] <- beta_beta[1] * (AGE_4[i, t] + AGE_YOUNG[i, t]) +
            beta_beta[2] * AGE_MIDDLE[i, t] +
            beta_beta[3] * AGE_OLD[i, t]
          # Outcome 1 cub at emergence
          gamma[i, t] <- beta_gamma[1] * (AGE_4[i, t] + AGE_YOUNG[i, t]) +
            beta_gamma[2] * AGE_MIDDLE[i, t] +
            beta_gamma[3] * AGE_OLD[i, t]

          S[1, 1, i, t] <- 0
          S[1, 2, i, t] <- phi[i, t]
          S[1, 3, i, t] <- 0
          S[1, 4, i, t] <- 0
          S[1, 5, i, t] <- 0
          S[1, 6, i, t] <- 0
          S[1, 7, i, t] <- 0
          S[1, 8, i, t] <- 0
          S[1, 9, i, t] <- 0
          S[1, 10, i, t] <- 0
          S[1, 11, i, t] <- 0
          S[1, 12, i, t] <- 0
          S[1, 13, i, t] <- 0
          S[1, 14, i, t] <- 1-phi[i, t]
          
          S[2, 1, i, t] <- 0
          S[2, 2, i, t] <- 0
          S[2, 3, i, t] <- phi[i, t]
          S[2, 4, i, t] <- 0
          S[2, 5, i, t] <- 0
          S[2, 6, i, t] <- 0
          S[2, 7, i, t] <- 0
          S[2, 8, i, t] <- 0
          S[2, 9, i, t] <- 0
          S[2, 10, i, t] <- 0
          S[2, 11, i, t] <- 0
          S[2, 12, i, t] <- 0
          S[2, 13, i, t] <- 0
          S[2, 14, i, t] <- 1-phi[i, t]
          
          S[3, 1, i, t] <- 0
          S[3, 2, i, t] <- 0
          S[3, 3, i, t] <- 0
          S[3, 4, i, t] <- phi[i, t] * (1- eta[i, t])
          S[3, 5, i, t] <- phi[i, t] * eta[i, t] * (1-beta[i, t])
          S[3, 6, i, t] <- phi[i, t] * eta[i, t] * beta[i, t] * (1 - gamma[i, t])
          S[3, 7, i, t] <- phi[i, t] * eta[i, t] * beta[i, t] * gamma[i, t]
          S[3, 8, i, t] <- 0
          S[3, 9, i, t] <- 0
          S[3, 10, i, t] <- 0
          S[3, 11, i, t] <- 0
          S[3, 12, i, t] <- 0
          S[3, 13, i, t] <- 0
          S[3, 14, i, t] <- 1-phi[i, t]
          
          S[4, 1, i, t] <- 0
          S[4, 2, i, t] <- 0
          S[4, 3, i, t] <- 0
          S[4, 4, i, t] <- 0
          S[4, 5, i, t] <- 0
          S[4, 6, i, t] <- phi[i, t] * eta[i, t] * beta[i, t] * (1 - gamma[i, t])
          S[4, 7, i, t] <- phi[i, t] * eta[i, t] * beta[i, t] * gamma[i, t]
          S[4, 8, i, t] <- 0
          S[4, 9, i, t] <- 0
          S[4, 10, i, t] <- 0
          S[4, 11, i, t] <- 0
          S[4, 12, i, t] <- phi[i, t] * (1 - eta[i, t])
          S[4, 13, i, t] <- phi[i, t] * eta[i, t] * (1-beta[i, t])
          S[4, 14, i, t] <- 1-phi[i, t]
          
          S[5, 1, i, t] <- 0
          S[5, 2, i, t] <- 0
          S[5, 3, i, t] <- 0
          S[5, 4, i, t] <- 0
          S[5, 5, i, t] <- 0
          S[5, 6, i, t] <- phi[i, t] * eta[i, t] * beta[i, t] * (1 - gamma[i, t])
          S[5, 7, i, t] <- phi[i, t] * eta[i, t] * beta[i, t] * gamma[i, t]
          S[5, 8, i, t] <- 0
          S[5, 9, i, t] <- 0
          S[5, 10, i, t] <- 0
          S[5, 11, i, t] <- 0
          S[5, 12, i, t] <- phi[i, t] * (1 - eta[i, t])
          S[5, 13, i, t] <- phi[i, t] * eta[i, t] * (1-beta[i, t])
          S[5, 14, i, t] <- 1-phi[i, t]
          
          S[6, 1, i, t] <- 0
          S[6, 2, i, t] <- 0
          S[6, 3, i, t] <- 0
          S[6, 4, i, t] <- 0
          S[6, 5, i, t] <- 0
          S[6, 6, i, t] <- phi[i, t] * (1-s0) * eta[i, t] * beta[i, t] * (1 - gamma[i, t])
          S[6, 7, i, t] <- phi[i, t] * (1-s0) * eta[i, t] * beta[i, t] * gamma[i, t]
          S[6, 8, i, t] <- phi[i, t] * s0
          S[6, 9, i, t] <- 0
          S[6, 10, i, t] <- 0
          S[6, 11, i, t] <- 0
          S[6, 12, i, t] <- phi[i, t] * (1-s0) * (1-eta[i, t])
          S[6, 13, i, t] <- phi[i, t] * (1-s0) * eta[i, t] * (1-beta[i, t])
          S[6, 14, i, t] <- 1-phi[i, t]
          
          S[7, 1, i, t] <- 0
          S[7, 2, i, t] <- 0
          S[7, 3, i, t] <- 0
          S[7, 4, i, t] <- 0
          S[7, 5, i, t] <- 0
          S[7, 6, i, t] <- phi[i, t] * (1-s0)^2 * eta[i, t] * beta[i, t] * (1 - gamma[i, t])
          S[7, 7, i, t] <- phi[i, t] * (1-s0)^2 * eta[i, t] * beta[i, t] * gamma[i, t]
          S[7, 8, i, t] <- phi[i, t] * 2 * s0 * (1-s0)
          S[7, 9, i, t] <- phi[i, t] * s0^2
          S[7, 10, i, t] <- 0
          S[7, 11, i, t] <- 0
          S[7, 12, i, t] <- phi[i, t] * (1-s0)^2 * (1-eta[i, t])
          S[7, 13, i, t] <- phi[i, t] * (1-s0)^2 * eta[i, t] * (1-beta[i, t])
          S[7, 14, i, t] <- 1-phi[i, t]
          
          S[8, 1, i, t] <- 0
          S[8, 2, i, t] <- 0
          S[8, 3, i, t] <- 0
          S[8, 4, i, t] <- 0
          S[8, 5, i, t] <- 0
          S[8, 6, i, t] <- phi[i, t] * (1-s1) * eta[i, t] * beta[i, t] * (1 - gamma[i, t])
          S[8, 7, i, t] <- phi[i, t] * (1-s1) * eta[i, t] * beta[i, t] * gamma[i, t]
          S[8, 8, i, t] <- 0
          S[8, 9, i, t] <- 0
          S[8, 10, i, t] <- phi[i, t] * s1
          S[8, 11, i, t] <- 0
          S[8, 12, i, t] <- phi[i, t] * (1-s1) * (1-eta[i, t])
          S[8, 13, i, t] <- phi[i, t] * (1-s1) * eta[i, t] * (1-beta[i, t])
          S[8, 14, i, t] <- 1-phi[i, t]
          
          S[9, 1, i, t] <- 0
          S[9, 2, i, t] <- 0
          S[9, 3, i, t] <- 0
          S[9, 4, i, t] <- 0
          S[9, 5, i, t] <- 0
          S[9, 6, i, t] <- phi[i, t] * (1-s1)^2 * eta[i, t] * beta[i, t] * (1 - gamma[i, t])
          S[9, 7, i, t] <- phi[i, t] * (1-s1)^2 * eta[i, t] * beta[i, t] * gamma[i, t]
          S[9, 8, i, t] <- 0
          S[9, 9, i, t] <- 0
          S[9, 10, i, t] <- phi[i, t] * 2 * s1 * (1-s1)
          S[9, 11, i, t] <- phi[i, t] * s1^2
          S[9, 12, i, t] <- phi[i, t] * (1-s1)^2 * (1-eta[i, t])
          S[9, 13, i, t] <- phi[i, t] * (1-s1)^2 * eta[i, t] * (1-beta[i, t])
          S[9, 14, i, t] <- 1-phi[i, t]
          
          S[10, 1, i, t] <- 0
          S[10, 2, i, t] <- 0
          S[10, 3, i, t] <- 0
          S[10, 4, i, t] <- 0
          S[10, 5, i, t] <- 0
          S[10, 6, i, t] <- phi[i, t] * eta[i, t] * beta[i, t] * (1 - gamma[i, t])
          S[10, 7, i, t] <- phi[i, t] * eta[i, t] * beta[i, t] * gamma[i, t]
          S[10, 8, i, t] <- 0
          S[10, 9, i, t] <- 0
          S[10, 10, i, t] <- 0
          S[10, 11, i, t] <- 0
          S[10, 12, i, t] <- phi[i, t] * (1 - eta[i, t])
          S[10, 13, i, t] <- phi[i, t] * eta[i, t] * (1-beta[i, t])
          S[10, 14, i, t] <- 1-phi[i, t]
          
          S[11, 1, i, t] <- 0
          S[11, 2, i, t] <- 0
          S[11, 3, i, t] <- 0
          S[11, 4, i, t] <- 0
          S[11, 5, i, t] <- 0
          S[11, 6, i, t] <- phi[i, t] * eta[i, t] * beta[i, t] * (1 - gamma[i, t])
          S[11, 7, i, t] <- phi[i, t] * eta[i, t] * beta[i, t] * gamma[i, t]
          S[11, 8, i, t] <- 0
          S[11, 9, i, t] <- 0
          S[11, 10, i, t] <- 0
          S[11, 11, i, t] <- 0
          S[11, 12, i, t] <- phi[i, t] * (1 - eta[i, t])
          S[11, 13, i, t] <- phi[i, t] * eta[i, t] * (1-beta[i, t])
          S[11, 14, i, t] <- 1-phi[i, t]
          
          S[12, 1, i, t] <- 0
          S[12, 2, i, t] <- 0
          S[12, 3, i, t] <- 0
          S[12, 4, i, t] <- 0
          S[12, 5, i, t] <- 0
          S[12, 6, i, t] <- phi[i, t] * eta[i, t] * beta[i, t] * (1 - gamma[i, t])
          S[12, 7, i, t] <- phi[i, t] * eta[i, t] * beta[i, t] * gamma[i, t]
          S[12, 8, i, t] <- 0
          S[12, 9, i, t] <- 0
          S[12, 10, i, t] <- 0
          S[12, 11, i, t] <- 0
          S[12, 12, i, t] <- phi[i, t] * (1 - eta[i, t])
          S[12, 13, i, t] <- phi[i, t] * eta[i, t] * (1-beta[i, t])
          S[12, 14, i, t] <- 1-phi[i, t]
          
          S[13, 1, i, t] <- 0
          S[13, 2, i, t] <- 0
          S[13, 3, i, t] <- 0
          S[13, 4, i, t] <- 0
          S[13, 5, i, t] <- 0
          S[13, 6, i, t] <- phi[i, t] * eta[i, t] * beta[i, t] * (1 - gamma[i, t])
          S[13, 7, i, t] <- phi[i, t] * eta[i, t] * beta[i, t] * gamma[i, t]
          S[13, 8, i, t] <- 0
          S[13, 9, i, t] <- 0
          S[13, 10, i, t] <- 0
          S[13, 11, i, t] <- 0
          S[13, 12, i, t] <- phi[i, t] * (1 - eta[i, t])
          S[13, 13, i, t] <- phi[i, t] * eta[i, t] * (1-beta[i, t])
          S[13, 14, i, t] <- 1-phi[i, t]
          
          S[14, 1, i, t] <- 0
          S[14, 2, i, t] <- 0
          S[14, 3, i, t] <- 0
          S[14, 4, i, t] <- 0
          S[14, 5, i, t] <- 0
          S[14, 6, i, t] <- 0
          S[14, 7, i, t] <- 0
          S[14, 8, i, t] <- 0
          S[14, 9, i, t] <- 0
          S[14, 10, i, t] <- 0
          S[14, 11, i, t] <- 0
          S[14, 12, i, t] <- 0
          S[14, 13, i, t] <- 0
          S[14, 14, i, t] <- 1
          
          Z[i, t+1] <- which(rmultinom(n = 1, size = 1,
                                       prob = S[Z[i, t], , i, t]) == 1)
          
          p_realized[i, t+1] <- rbinom(n = 1, size = 1, prob = values_p[k])
          
          # Obs matrix (14 states but two pairs of states or indistinguishable because
          # absence of GPS-GLS data)
          E[1, 1, i, t] <- p_realized[i, t+1]
          E[1, 2, i, t] <- 0
          E[1, 3, i, t] <- 0
          E[1, 4, i, t] <- 0
          E[1, 5, i, t] <- 0
          E[1, 6, i, t] <- 0
          E[1, 7, i, t] <- 0
          E[1, 8, i, t] <- 0
          E[1, 9, i, t] <- 0
          E[1, 10, i, t] <- 0
          E[1, 11, i, t] <- 0
          E[1, 12, i, t] <- 1-p_realized[i, t+1]
          
          E[2, 1, i, t] <- 0
          E[2, 2, i, t] <- p_realized[i, t+1]
          E[2, 3, i, t] <- 0
          E[2, 4, i, t] <- 0
          E[2, 5, i, t] <- 0
          E[2, 6, i, t] <- 0
          E[2, 7, i, t] <- 0
          E[2, 8, i, t] <- 0
          E[2, 9, i, t] <- 0
          E[2, 10, i, t] <- 0
          E[2, 11, i, t] <- 0
          E[2, 12, i, t] <- 1 - p_realized[i, t+1]
          
          E[3, 1, i, t] <- 0
          E[3, 2, i, t] <- 0
          E[3, 3, i, t] <- p_realized[i, t+1]
          E[3, 4, i, t] <- 0
          E[3, 5, i, t] <- 0
          E[3, 6, i, t] <- 0
          E[3, 7, i, t] <- 0
          E[3, 8, i, t] <- 0
          E[3, 9, i, t] <- 0
          E[3, 10, i, t] <- 0
          E[3, 11, i, t] <- 0
          E[3, 12, i, t] <- 1 - p_realized[i, t+1]
          
          E[4, 1, i, t] <- 0
          E[4, 2, i, t] <- 0
          E[4, 3, i, t] <- 0
          E[4, 4, i, t] <- p_realized[i, t+1]
          E[4, 5, i, t] <- 0
          E[4, 6, i, t] <- 0
          E[4, 7, i, t] <- 0
          E[4, 8, i, t] <- 0
          E[4, 9, i, t] <- 0
          E[4, 10, i, t] <- 0
          E[4, 11, i, t] <- 0
          E[4, 12, i, t] <- 1 - p_realized[i, t+1]
          
          E[5, 1, i, t] <- 0
          E[5, 2, i, t] <- 0
          E[5, 3, i, t] <- 0
          E[5, 4, i, t] <- p_realized[i, t+1]
          E[5, 5, i, t] <- 0
          E[5, 6, i, t] <- 0
          E[5, 7, i, t] <- 0
          E[5, 8, i, t] <- 0
          E[5, 9, i, t] <- 0
          E[5, 10, i, t] <- 0
          E[5, 11, i, t] <- 0
          E[5, 12, i, t] <- 1 - p_realized[i, t+1]
          
          E[6, 1, i, t] <- 0
          E[6, 2, i, t] <- 0
          E[6, 3, i, t] <- 0
          E[6, 4, i, t] <- 0
          E[6, 5, i, t] <- p_realized[i, t+1]
          E[6, 6, i, t] <- 0
          E[6, 7, i, t] <- 0
          E[6, 8, i, t] <- 0
          E[6, 9, i, t] <- 0
          E[6, 10, i, t] <- 0
          E[6, 11, i, t] <- 0
          E[6, 12, i, t] <- 1 - p_realized[i, t+1]
          
          E[7, 1, i, t] <- 0
          E[7, 2, i, t] <- 0
          E[7, 3, i, t] <- 0
          E[7, 4, i, t] <- 0
          E[7, 5, i, t] <- 0
          E[7, 6, i, t] <- p_realized[i, t+1]
          E[7, 7, i, t] <- 0
          E[7, 8, i, t] <- 0
          E[7, 9, i, t] <- 0
          E[7, 10, i, t] <- 0
          E[7, 11, i, t] <- 0
          E[7, 12, i, t] <- 1 - p_realized[i, t+1]
          
          E[8, 1, i, t] <- 0
          E[8, 2, i, t] <- 0
          E[8, 3, i, t] <- 0
          E[8, 4, i, t] <- 0
          E[8, 5, i, t] <- 0
          E[8, 6, i, t] <- 0
          E[8, 7, i, t] <- p_realized[i, t+1]
          E[8, 8, i, t] <- 0
          E[8, 9, i, t] <- 0
          E[8, 10, i, t] <- 0
          E[8, 11, i, t] <- 0
          E[8, 12, i, t] <- 1 - p_realized[i, t+1]
          
          E[9, 1, i, t] <- 0
          E[9, 2, i, t] <- 0
          E[9, 3, i, t] <- 0
          E[9, 4, i, t] <- 0
          E[9, 5, i, t] <- 0
          E[9, 6, i, t] <- 0
          E[9, 7, i, t] <- 0
          E[9, 8, i, t] <- p_realized[i, t+1]
          E[9, 9, i, t] <- 0
          E[9, 10, i, t] <- 0
          E[9, 11, i, t] <- 0
          E[9, 12, i, t] <- 1 - p_realized[i, t+1]
          
          E[10, 1, i, t] <- 0
          E[10, 2, i, t] <- 0
          E[10, 3, i, t] <- 0
          E[10, 4, i, t] <- 0
          E[10, 5, i, t] <- 0
          E[10, 6, i, t] <- 0
          E[10, 7, i, t] <- 0
          E[10, 8, i, t] <- 0
          E[10, 9, i, t] <- (1 - a) * p_realized[i, t+1]
          E[10, 10, i, t] <- 0
          E[10, 11, i, t] <- a * p_realized[i, t+1]
          E[10, 12, i, t] <- 1 - p_realized[i, t+1]
          
          E[11, 1, i, t] <- 0
          E[11, 2, i, t] <- 0
          E[11, 3, i, t] <- 0
          E[11, 4, i, t] <- 0
          E[11, 5, i, t] <- 0
          E[11, 6, i, t] <- 0
          E[11, 7, i, t] <- 0
          E[11, 8, i, t] <- 0
          E[11, 9, i, t] <- 2 * (1 - a) * a * p_realized[i, t+1]
          E[11, 10, i, t] <- (1 - a)^2 * p_realized[i, t+1]
          E[11, 11, i, t] <- a^2 * p_realized[i, t+1]
          E[11, 12, i, t] <- 1 - p_realized[i, t+1]
          
          E[12, 1, i, t] <- 0
          E[12, 2, i, t] <- 0
          E[12, 3, i, t] <- 0
          E[12, 4, i, t] <- 0
          E[12, 5, i, t] <- 0
          E[12, 6, i, t] <- 0
          E[12, 7, i, t] <- 0
          E[12, 8, i, t] <- 0
          E[12, 9, i, t] <- 0
          E[12, 10, i, t] <- 0
          E[12, 11, i, t] <- p_realized[i, t+1]
          E[12, 12, i, t] <- 1 - p_realized[i, t+1]
          
          E[13, 1, i, t] <- 0
          E[13, 2, i, t] <- 0
          E[13, 3, i, t] <- 0
          E[13, 4, i, t] <- 0
          E[13, 5, i, t] <- 0
          E[13, 6, i, t] <- 0
          E[13, 7, i, t] <- 0
          E[13, 8, i, t] <- 0
          E[13, 9, i, t] <- 0
          E[13, 10, i, t] <- 0
          E[13, 11, i, t] <- p_realized[i, t+1]
          E[13, 12, i, t] <- 1 - p_realized[i, t+1]
          
          E[14, 1, i, t] <- 0
          E[14, 2, i, t] <- 0
          E[14, 3, i, t] <- 0
          E[14, 4, i, t] <- 0
          E[14, 5, i, t] <- 0
          E[14, 6, i, t] <- 0
          E[14, 7, i, t] <- 0
          E[14, 8, i, t] <- 0
          E[14, 9, i, t] <- 0
          E[14, 10, i, t] <- 0
          E[14, 11, i, t] <- 0
          E[14, 12, i, t] <- 1
          
          CH[i, t+1] <-  which(rmultinom(n = 1, size = 1,
                                         E[Z[i, t+1], , i, t]) == 1)
        }
        
        # GPS collaring ++++++++++++++++++++++++++++++++++++++
        if (p_realized[i, t+1] == 1) {  # If physically captured
          
          collared <- which(rmultinom(n = 1, size = 1, 
                                      prob = COLLARING[Z[i, f[i]], ]) == 1) - 1
          
          if (collared == 1) { # If individual has a collar, for how many years ?
            collared_years <- rpois(n = 1, lambda = lambda_GPS) + 1
            q <- 1
            while (q <= collared_years &  
                   (t + q + 1) <= n_occasions) {
              KHI[i, t + q + 1] <- 1
              q <- q + 1
            }
            
          }
        }
      }
      
      # ~ 3. Check data --------------------------------------------------------------
      
      # Let's check that all the transitions are possible. 
      S <- matrix(data = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,   # J2
                           0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,   # J3
                           0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1,   # SA4
                           0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1,   # NDSA5
                           0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1,   # FDSA5  5
                           0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1,   # A01
                           0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1,   # A02
                           0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1,   # A11    
                           0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1,   # A12
                           0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1,   # AS1    10
                           0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1,   # AS2
                           0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1,   # LNDA
                           0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1,   # FDA    
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),  # D
                  nrow = n_states_simul, ncol = n_states_simul, byrow = TRUE)
      check_Z <- CH ; check_Z[] <- NA
      for (i in 1:dim(CH)[1]) {
        if (f[i] == n_occasions) next
        for (t in f[i]:(n_occasions - 1)) {
          check_Z[i, t] <- ifelse(S[Z[i, t], Z[i, t+1]] == 0, 1, 0)
        }
      }
      print(sum(check_Z, na.rm = T) == 0)
      
      # Let's now check that the states are compatible with the observations
      E <- matrix(data = c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,    # J2
                           0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,    # J3
                           0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1,    # SA4
                           0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,    # NDSA5
                           0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,    # FDSA5  5
                           0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,    # A01
                           0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1,    # A02
                           0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1,    # A11    
                           0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1,    # A12
                           0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1,    # AS1    10
                           0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,    # AS2
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,    # LNDA
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,    # FDA    
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),   # D
                  nrow = n_states_simul, ncol = n_events, byrow = TRUE)
      
      check_E <- CH ; check_Z[] <- NA
      for (i in 1:dim(CH)[1]) {
        for (t in f[i]:n_occasions) {
          check_E[i, t] <- ifelse(E[Z[i, t], CH[i, t]] == 0, 1, 0)
          if (check_E[i, t] == 1) {
            print(paste0("i = ", i, " ; t = ", t))
            print(paste0("state: ", Z[i, t], ": event: ",  CH[i, t]))
          }
        }
      }
      print(sum(check_E, na.rm = T) == 0)
      
      # Change the values in Z to match the states that exist in the model (for the initial values of Z)
      for (i in 1:dim(CH)[1]) {
        f_state[i] <- ifelse(f_state[i] %in% 1:4, f_state[i],
                          ifelse(f_state[i] %in% 5:12, f_state[i] - 1,
                                 ifelse(f_state[i] %in% 13:14, f_state[i] - 2, NA)))
        for (t in f[i]:n_occasions) {
          Z[i, t] <- ifelse(Z[i, t] %in% 1:4, Z[i, t],
                            ifelse(Z[i, t] %in% 5:12, Z[i, t] - 1,
                                   ifelse(Z[i, t] %in% 13:14, Z[i, t] - 2, NA)))

        }
      }
      
      # Let's check again that all the transitions are possible with only twelve states
      S <- matrix(data = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,   # J2
                           0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1,   # J3
                           0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1,   # SA4
                           0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1,   # SA5
                           0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1,   # A01  5
                           0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1,   # A02
                           0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1,   # A11    
                           0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1,   # A12
                           0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1,   # AS1    10
                           0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1,   # AS2
                           0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1,   # A
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),  # D
                  nrow = n_states, ncol = n_states, byrow = TRUE)
      check_Z <- CH ; check_Z[] <- NA
      for (i in 1:dim(CH)[1]) {
        if (f[i] == n_occasions) next
        for (t in f[i]:(n_occasions - 1)) {
          check_Z[i, t] <- ifelse(S[Z[i, t], Z[i, t+1]] == 0, 1, 0)
        }
      }
      print(sum(check_Z, na.rm = T) == 0)
      
      # Let's now check again that the states are compatible with the observations
      E <- matrix(data = c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,    # J2
                           0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,    # J3
                           0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1,    # SA4
                           0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,    # SA5
                           0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,    # A01
                           0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1,    # A02
                           0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1,    # A11    
                           0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1,    # A12
                           0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1,    # AS1    10
                           0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,    # AS2
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,    # A
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),   # D
                  nrow = n_states, ncol = n_events, byrow = TRUE)
      
      check_E <- CH ; check_Z[] <- NA
      for (i in 1:dim(CH)[1]) {
        for (t in f[i]:n_occasions) {
          check_E[i, t] <- ifelse(E[Z[i, t], CH[i, t]] == 0, 1, 0)
          if (check_E[i, t] == 1) {
            print(paste0("i = ", i, " ; t = ", t))
            print(paste0("state: ", Z[i, t], ": event: ",  CH[i, t]))
          }
        }
      }
      print(sum(check_E, na.rm = T) == 0)
      
      # ~ 4. Bundle data ------------------------------------------------------------- 
      
      # Get simulated proportion of each state at first capture
      prop_raw <- table(f_state)/length(f_state)
      prop <- NULL
      for (o in 1:n_states) {
        x <- prop_raw[as.character(o)]
        prop[o] <- ifelse(!is.na(x), x, 0)
      }
      
      dat <- list(ones = rep(1, times = nrow(CH)))
      my.constants <- list(CH = CH,
                           freq = rep(1, times = nrow(CH)),
                           f = f,
                           n_ind = dim(CH)[1],
                           n_occasions = dim(CH)[2],
                           n_states = n_states,
                           S0 = prop,
                           # KHI = KHI,
                           a = a,
                           AGE_4 = AGE_4,
                           # AGE_5 = AGE_5,
                           AGE_YOUNG = AGE_YOUNG,
                           AGE_MIDDLE = AGE_MIDDLE,
                           AGE_OLD = AGE_OLD,
                           AGE_2_4 = AGE_2_4)
      
      # Save
      name <- paste0("03_outputs/simulations/simulated_datasets/albatross/no_biologging/",
                     "p_", values_p[k],
                     "_pGPS_", values_pGPS[l],
                     "_dataset_", n, ".RData")
      save(CR_model_biologging_sim, my.constants, dat, Z, # zeta,
           file = name)
      
    }
  }
}


end_global <- Sys.time() 
print(paste0("global end: ", end_global)) 
print(paste0("global duration: ", end_global - start_global))                 
 


