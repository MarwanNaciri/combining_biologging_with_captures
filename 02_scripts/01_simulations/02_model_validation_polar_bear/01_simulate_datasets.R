#==============================================================================#
#                                                                              #
#              Simulation of datasets for to validate the CR model             #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(nimble)

# The purpose of this script is to generate a dataset of capture histories that 
# is comparable to the observed dataset, in terms of:
# - values of demogrpahic parameters (both vital rates and recapture probability)
# - proportion of individuals that are equipped with a GPS or a GLS
# - lifespan of GPS/GLS
# - number of individuals




# How many datasets?
n_sim <- 100


N <- N_p <- N_3 <- N_GPS <- N_GLS <- N_remote <- NULL

# A. Build the model ===========================================================

CR_model <- nimbleCode({

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
    E0[1, 13] <- 0
    E0[1, 14] <- 0
    E0[1, 15] <- 0
    E0[1, 16] <- 0
    E0[1, 17] <- 0
    E0[1, 18] <- 0
    E0[1, 19] <- 0

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
    E0[2, 13] <- 0
    E0[2, 14] <- 0
    E0[2, 15] <- 0
    E0[2, 16] <- 0
    E0[2, 17] <- 0
    E0[2, 18] <- 0
    E0[2, 19] <- 0

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
    E0[3, 13] <- 0
    E0[3, 14] <- 0
    E0[3, 15] <- 0
    E0[3, 16] <- 0
    E0[3, 17] <- 0
    E0[3, 18] <- 0
    E0[3, 19] <- 0

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
    E0[4, 13] <- 0
    E0[4, 14] <- 0
    E0[4, 15] <- 0
    E0[4, 16] <- 0
    E0[4, 17] <- 0
    E0[4, 18] <- 0
    E0[4, 19] <- 0

    E0[5, 1] <- 0
    E0[5, 2] <- 0
    E0[5, 3] <- 0
    E0[5, 4] <- 1
    E0[5, 5] <- 0
    E0[5, 6] <- 0
    E0[5, 7] <- 0
    E0[5, 8] <- 0
    E0[5, 9] <- 0
    E0[5, 10] <- 0
    E0[5, 11] <- 0
    E0[5, 12] <- 0
    E0[5, 13] <- 0
    E0[5, 14] <- 0
    E0[5, 15] <- 0
    E0[5, 16] <- 0
    E0[5, 17] <- 0
    E0[5, 18] <- 0
    E0[5, 19] <- 0

    E0[6, 1] <- 0
    E0[6, 2] <- 0
    E0[6, 3] <- 0
    E0[6, 4] <- 0
    E0[6, 5] <- 1
    E0[6, 6] <- 0
    E0[6, 7] <- 0
    E0[6, 8] <- 0
    E0[6, 9] <- 0
    E0[6, 10] <- 0
    E0[6, 11] <- 0
    E0[6, 12] <- 0
    E0[6, 13] <- 0
    E0[6, 14] <- 0
    E0[6, 15] <- 0
    E0[6, 16] <- 0
    E0[6, 17] <- 0
    E0[6, 18] <- 0
    E0[6, 19] <- 0

    E0[7, 1] <- 0
    E0[7, 2] <- 0
    E0[7, 3] <- 0
    E0[7, 4] <- 0
    E0[7, 5] <- 0
    E0[7, 6] <- 1
    E0[7, 7] <- 0
    E0[7, 8] <- 0
    E0[7, 9] <- 0
    E0[7, 10] <- 0
    E0[7, 11] <- 0
    E0[7, 12] <- 0
    E0[7, 13] <- 0
    E0[7, 14] <- 0
    E0[7, 15] <- 0
    E0[7, 16] <- 0
    E0[7, 17] <- 0
    E0[7, 18] <- 0
    E0[7, 19] <- 0

    E0[8, 1] <- 0
    E0[8, 2] <- 0
    E0[8, 3] <- 0
    E0[8, 4] <- 0
    E0[8, 5] <- 0
    E0[8, 6] <- 0
    E0[8, 7] <- 1
    E0[8, 8] <- 0
    E0[8, 9] <- 0
    E0[8, 10] <- 0
    E0[8, 11] <- 0
    E0[8, 12] <- 0
    E0[8, 13] <- 0
    E0[8, 14] <- 0
    E0[8, 15] <- 0
    E0[8, 16] <- 0
    E0[8, 17] <- 0
    E0[8, 18] <- 0
    E0[8, 19] <- 0

    E0[9, 1] <- 0
    E0[9, 2] <- 0
    E0[9, 3] <- 0
    E0[9, 4] <- 0
    E0[9, 5] <- 0
    E0[9, 6] <- 0
    E0[9, 7] <- 0
    E0[9, 8] <- 1
    E0[9, 9] <- 0
    E0[9, 10] <- 0
    E0[9, 11] <- 0
    E0[9, 12] <- 0
    E0[9, 13] <- 0
    E0[9, 14] <- 0
    E0[9, 15] <- 0
    E0[9, 16] <- 0
    E0[9, 17] <- 0
    E0[9, 18] <- 0
    E0[9, 19] <- 0

    E0[10, 1] <- 0
    E0[10, 2] <- 0
    E0[10, 3] <- 0
    E0[10, 4] <- 0
    E0[10, 5] <- 0
    E0[10, 6] <- 0
    E0[10, 7] <- 0
    E0[10, 8] <- 0
    E0[10, 9] <- 1-a
    E0[10, 10] <- 0
    E0[10, 11] <- a
    E0[10, 12] <- 0
    E0[10, 13] <- 0
    E0[10, 14] <- 0
    E0[10, 15] <- 0
    E0[10, 16] <- 0
    E0[10, 17] <- 0
    E0[10, 18] <- 0
    E0[10, 19] <- 0

    E0[11, 1] <- 0
    E0[11, 2] <- 0
    E0[11, 3] <- 0
    E0[11, 4] <- 0
    E0[11, 5] <- 0
    E0[11, 6] <- 0
    E0[11, 7] <- 0
    E0[11, 8] <- 0
    E0[11, 9] <- 2 * ( 1-a ) * a
    E0[11, 10] <- ( 1 - a )^2
    E0[11, 11] <- a^2
    E0[11, 12] <- 0
    E0[11, 13] <- 0
    E0[11, 14] <- 0
    E0[11, 15] <- 0
    E0[11, 16] <- 0
    E0[11, 17] <- 0
    E0[11, 18] <- 0
    E0[11, 19] <- 0

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
    E0[12, 11] <- 1
    E0[12, 12] <- 0
    E0[12, 13] <- 0
    E0[12, 14] <- 0
    E0[12, 15] <- 0
    E0[12, 16] <- 0
    E0[12, 17] <- 0
    E0[12, 18] <- 0
    E0[12, 19] <- 0

    E0[13, 1] <- 0
    E0[13, 2] <- 0
    E0[13, 3] <- 0
    E0[13, 4] <- 0
    E0[13, 5] <- 0
    E0[13, 6] <- 0
    E0[13, 7] <- 0
    E0[13, 8] <- 0
    E0[13, 9] <- 0
    E0[13, 10] <- 0
    E0[13, 11] <- 1
    E0[13, 12] <- 0
    E0[13, 13] <- 0
    E0[13, 14] <- 0
    E0[13, 15] <- 0
    E0[13, 16] <- 0
    E0[13, 17] <- 0
    E0[13, 18] <- 0
    E0[13, 19] <- 0

    E0[14, 1] <- 0
    E0[14, 2] <- 0
    E0[14, 3] <- 0
    E0[14, 4] <- 0
    E0[14, 5] <- 0
    E0[14, 6] <- 0
    E0[14, 7] <- 0
    E0[14, 8] <- 0
    E0[14, 9] <- 0
    E0[14, 10] <- 0
    E0[14, 11] <- 0
    E0[14, 12] <- 0
    E0[14, 13] <- 0
    E0[14, 14] <- 0
    E0[14, 15] <- 0
    E0[14, 16] <- 0
    E0[14, 17] <- 0
    E0[14, 18] <- 0
    E0[14, 19] <- 1

    for (i in 1:n_ind) {
      
    for (t in f[i]:(n_occasions - 1)) {

      # ~~~ b. Regressions on parameters -------------------
      
      # Probability of juvenile, subadult and adult survival
      logit(phi[i, t]) <- beta_phi[1] * AGE_2_4[i, t] +
        beta_phi[2] * (AGE_YOUNG[i, t] + AGE_MIDDLE[i, t]) +
        beta_phi[3] * AGE_OLD[i, t] +
        beta_phi[4] * (AGE_2_4[i, t] + AGE_OLD[i, t]) * sea_ice_s[t]

      # Probability of coy survival
      logit(s01[i, t]) <- beta_s0[1] * AGE_YOUNG_1[i, t] +
        beta_s0[2] * AGE_MIDDLE_1[i, t] +
        beta_s0[3] * AGE_OLD_1[i, t] +
        beta_s0[4] +
        beta_s0[5] * sea_ice_s[t]
      
      logit(s02[i, t]) <- beta_s0[1] * AGE_YOUNG_1[i, t] +
        beta_s0[2] * AGE_MIDDLE_1[i, t] +
        beta_s0[3] * AGE_OLD_1[i, t] +
        beta_s0[5] * sea_ice_s[t]
      
      # Probability of yearling survival
      logit(s1[i, t]) <- beta_s1[1] * AGE_YOUNG_2[i, t] +
        beta_s1[2] * AGE_MIDDLE_2[i, t] +
        beta_s1[3] * AGE_OLD_2[i, t] +
        beta_s1[4] * sea_ice_s[t]

      # Probability of denning
      logit(eta[1, i, t]) <- beta_eta[1] * AGE_4[i, t] +   # Lone female
        beta_eta[2] * AGE_YOUNG[i, t] +
        beta_eta[3] * AGE_MIDDLE[i, t] +
        beta_eta[4] * AGE_OLD[i, t] +
        beta_eta[9] * sea_ice_s[t]
      
      logit(eta[2, i, t]) <- beta_eta[2] * AGE_YOUNG[i, t] + # FDA
        beta_eta[3] * AGE_MIDDLE[i, t] +
        beta_eta[4] * AGE_OLD[i, t] +
        beta_eta[5] +
        beta_eta[9] * sea_ice_s[t]
      
      logit(eta[3, i, t]) <- beta_eta[2] * AGE_YOUNG[i, t] + # A0-
        beta_eta[3] * AGE_MIDDLE[i, t] +
        beta_eta[4] * AGE_OLD[i, t] +
        beta_eta[6] +
        beta_eta[9] * sea_ice_s[t]
      
      logit(eta[4, i, t]) <- beta_eta[2] * AGE_YOUNG[i, t] + # A1-
        beta_eta[3] * AGE_MIDDLE[i, t] +
        beta_eta[4] * AGE_OLD[i, t] +
        beta_eta[7] +
        beta_eta[9] * sea_ice_s[t]
      
      logit(eta[5, i, t]) <- beta_eta[2] * AGE_YOUNG[i, t] + # AS
        beta_eta[3] * AGE_MIDDLE[i, t] +
        beta_eta[4] * AGE_OLD[i, t] +
        beta_eta[8] +
        beta_eta[9] * sea_ice_s[t]
      
      # Early litter survival
      logit(beta[i, t]) <- beta_beta[1] * (AGE_4[i, t] + AGE_YOUNG[i, t]) +
        beta_beta[2] * AGE_MIDDLE[i, t] +
        beta_beta[3] * AGE_OLD[i, t] +
        beta_beta[4] * sea_ice_s[t]
      
      # twinning prob
      logit(gamma[i, t]) <- beta_gamma[1] * (AGE_4[i, t] + AGE_YOUNG[i, t]) +
        beta_gamma[2] * AGE_MIDDLE[i, t] +
        beta_gamma[3] * AGE_OLD[i, t] +
        beta_gamma[4] * sea_ice_s[t]

      
      # ~~~ c. S ---------------------------------------------------------------

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
      S[3, 4, i, t] <- phi[i, t] * (1- eta[1, i, t])
      S[3, 5, i, t] <- phi[i, t] * eta[1, i, t] * (1-beta[i, t])
      S[3, 6, i, t] <- phi[i, t] * eta[1, i, t] * beta[i, t] * (1 - gamma[i, t])
      S[3, 7, i, t] <- phi[i, t] * eta[1, i, t] * beta[i, t] * gamma[i, t]
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
      S[4, 6, i, t] <- phi[i, t] * eta[1, i, t] * beta[i, t] * (1 - gamma[i, t])
      S[4, 7, i, t] <- phi[i, t] * eta[1, i, t] * beta[i, t] * gamma[i, t]
      S[4, 8, i, t] <- 0
      S[4, 9, i, t] <- 0
      S[4, 10, i, t] <- 0
      S[4, 11, i, t] <- 0
      S[4, 12, i, t] <- phi[i, t] * (1 - eta[1, i, t])
      S[4, 13, i, t] <- phi[i, t] * eta[1, i, t] * (1-beta[i, t])
      S[4, 14, i, t] <- 1-phi[i, t]

      S[5, 1, i, t] <- 0
      S[5, 2, i, t] <- 0
      S[5, 3, i, t] <- 0
      S[5, 4, i, t] <- 0
      S[5, 5, i, t] <- 0
      S[5, 6, i, t] <- phi[i, t] * eta[2, i, t] * beta[i, t] * (1 - gamma[i, t])
      S[5, 7, i, t] <- phi[i, t] * eta[2, i, t] * beta[i, t] * gamma[i, t]
      S[5, 8, i, t] <- 0
      S[5, 9, i, t] <- 0
      S[5, 10, i, t] <- 0
      S[5, 11, i, t] <- 0
      S[5, 12, i, t] <- phi[i, t] * (1 - eta[2, i, t])
      S[5, 13, i, t] <- phi[i, t] * eta[2, i, t] * (1-beta[i, t])
      S[5, 14, i, t] <- 1-phi[i, t]

      S[6, 1, i, t] <- 0
      S[6, 2, i, t] <- 0
      S[6, 3, i, t] <- 0
      S[6, 4, i, t] <- 0
      S[6, 5, i, t] <- 0
      S[6, 6, i, t] <- phi[i, t] * eta[3, i, t] * (1-s01[i, t]) * beta[i, t] * (1 - gamma[i, t])
      S[6, 7, i, t] <- phi[i, t] * eta[3, i, t] * (1-s01[i, t]) * beta[i, t] * gamma[i, t]
      S[6, 8, i, t] <- phi[i, t] * s01[i, t]
      S[6, 9, i, t] <- 0
      S[6, 10, i, t] <- 0
      S[6, 11, i, t] <- 0
      S[6, 12, i, t] <- phi[i, t] * (1-s01[i, t]) * (1-eta[3, i, t])
      S[6, 13, i, t] <- phi[i, t] * eta[3, i, t] * (1-s01[i, t]) * (1-beta[i, t])
      S[6, 14, i, t] <- 1-phi[i, t]

      S[7, 1, i, t] <- 0
      S[7, 2, i, t] <- 0
      S[7, 3, i, t] <- 0
      S[7, 4, i, t] <- 0
      S[7, 5, i, t] <- 0
      S[7, 6, i, t] <- phi[i, t] * eta[3, i, t] * (1-s02[i, t])^2 * beta[i, t] * (1 - gamma[i, t])
      S[7, 7, i, t] <- phi[i, t] * eta[3, i, t] * (1-s02[i, t])^2 * beta[i, t] * gamma[i, t]
      S[7, 8, i, t] <- phi[i, t] * 2 * s02[i, t] * (1-s02[i, t])
      S[7, 9, i, t] <- phi[i, t] * s02[i, t]^2
      S[7, 10, i, t] <- 0
      S[7, 11, i, t] <- 0
      S[7, 12, i, t] <- phi[i, t] * (1-s02[i, t])^2 * (1-eta[3, i, t])
      S[7, 13, i, t] <- phi[i, t] * eta[3, i, t] * (1-s02[i, t])^2 * (1-beta[i, t])
      S[7, 14, i, t] <- 1-phi[i, t]

      S[8, 1, i, t] <- 0
      S[8, 2, i, t] <- 0
      S[8, 3, i, t] <- 0
      S[8, 4, i, t] <- 0
      S[8, 5, i, t] <- 0
      S[8, 6, i, t] <- phi[i, t] * eta[4, i, t] * (1-s1[i, t]) * beta[i, t] * (1 - gamma[i, t])
      S[8, 7, i, t] <- phi[i, t] * eta[4, i, t] * (1-s1[i, t]) * beta[i, t] * gamma[i, t]
      S[8, 8, i, t] <- 0
      S[8, 9, i, t] <- 0
      S[8, 10, i, t] <- phi[i, t] * s1[i, t]
      S[8, 11, i, t] <- 0
      S[8, 12, i, t] <- phi[i, t] * (1-s1[i, t]) * (1-eta[4, i, t])
      S[8, 13, i, t] <- phi[i, t] * eta[4, i, t] * (1-s1[i, t]) * (1-beta[i, t])
      S[8, 14, i, t] <- 1-phi[i, t]

      S[9, 1, i, t] <- 0
      S[9, 2, i, t] <- 0
      S[9, 3, i, t] <- 0
      S[9, 4, i, t] <- 0
      S[9, 5, i, t] <- 0
      S[9, 6, i, t] <- phi[i, t] * eta[4, i, t] * (1-s1[i, t])^2 * beta[i, t] * (1 - gamma[i, t])   
      S[9, 7, i, t] <- phi[i, t] * eta[4, i, t] * (1-s1[i, t])^2 * beta[i, t] * gamma[i, t]
      S[9, 8, i, t] <- 0
      S[9, 9, i, t] <- 0
      S[9, 10, i, t] <- phi[i, t] * 2 * s1[i, t] * (1-s1[i, t])
      S[9, 11, i, t] <- phi[i, t] * s1[i, t]^2
      S[9, 12, i, t] <- phi[i, t] * (1-s1[i, t])^2 * (1-eta[4, i, t])
      S[9, 13, i, t] <- phi[i, t] * eta[4, i, t] * (1-s1[i, t])^2 * (1-beta[i, t])
      S[9, 14, i, t] <- 1-phi[i, t]

      S[10, 1, i, t] <- 0
      S[10, 2, i, t] <- 0
      S[10, 3, i, t] <- 0
      S[10, 4, i, t] <- 0
      S[10, 5, i, t] <- 0
      S[10, 6, i, t] <- phi[i, t] * eta[5, i, t] * beta[i, t] * (1 - gamma[i, t])
      S[10, 7, i, t] <- phi[i, t] * eta[5, i, t] * beta[i, t] * gamma[i, t]
      S[10, 8, i, t] <- 0
      S[10, 9, i, t] <- 0
      S[10, 10, i, t] <- 0
      S[10, 11, i, t] <- 0
      S[10, 12, i, t] <- phi[i, t] * (1 - eta[5, i, t])
      S[10, 13, i, t] <- phi[i, t] * eta[5, i, t] * (1-beta[i, t])
      S[10, 14, i, t] <- 1-phi[i, t]

      S[11, 1, i, t] <- 0
      S[11, 2, i, t] <- 0
      S[11, 3, i, t] <- 0
      S[11, 4, i, t] <- 0
      S[11, 5, i, t] <- 0
      S[11, 6, i, t] <- phi[i, t] * eta[5, i, t] * beta[i, t] * (1 - gamma[i, t])
      S[11, 7, i, t] <- phi[i, t] * eta[5, i, t] * beta[i, t] * gamma[i, t]
      S[11, 8, i, t] <- 0
      S[11, 9, i, t] <- 0
      S[11, 10, i, t] <- 0
      S[11, 11, i, t] <- 0
      S[11, 12, i, t] <- phi[i, t] * (1 - eta[5, i, t])
      S[11, 13, i, t] <- phi[i, t] * eta[5, i, t] * (1-beta[i, t])
      S[11, 14, i, t] <- 1-phi[i, t]

      S[12, 1, i, t] <- 0
      S[12, 2, i, t] <- 0
      S[12, 3, i, t] <- 0
      S[12, 4, i, t] <- 0
      S[12, 5, i, t] <- 0
      S[12, 6, i, t] <- phi[i, t] * eta[1, i, t] * beta[i, t] * (1 - gamma[i, t])
      S[12, 7, i, t] <- phi[i, t] * eta[1, i, t] * beta[i, t] * gamma[i, t]
      S[12, 8, i, t] <- 0
      S[12, 9, i, t] <- 0
      S[12, 10, i, t] <- 0
      S[12, 11, i, t] <- 0
      S[12, 12, i, t] <- phi[i, t] * (1 - eta[1, i, t])
      S[12, 13, i, t] <- phi[i, t] * eta[1, i, t] * (1-beta[i, t])
      S[12, 14, i, t] <- 1-phi[i, t]

      S[13, 1, i, t] <- 0
      S[13, 2, i, t] <- 0
      S[13, 3, i, t] <- 0
      S[13, 4, i, t] <- 0
      S[13, 5, i, t] <- 0
      S[13, 6, i, t] <- phi[i, t] * eta[2, i, t] * beta[i, t] * (1 - gamma[i, t])
      S[13, 7, i, t] <- phi[i, t] * eta[2, i, t] * beta[i, t] * gamma[i, t]
      S[13, 8, i, t] <- 0
      S[13, 9, i, t] <- 0
      S[13, 10, i, t] <- 0
      S[13, 11, i, t] <- 0
      S[13, 12, i, t] <- phi[i, t] * (1 - eta[2, i, t])
      S[13, 13, i, t] <- phi[i, t] * eta[2, i, t] * (1-beta[i, t])
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

      # ~~~ d. E ----------------------------------

      E[1, 1, i, t] <- p
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
      E[1, 12, i, t] <- 0
      E[1, 13, i, t] <- 0
      E[1, 14, i, t] <- 0
      E[1, 15, i, t] <- 0
      E[1, 16, i, t] <- 0
      E[1, 17, i, t] <- 0
      E[1, 18, i, t] <- 0
      E[1, 19, i, t] <- 1-p

      E[2, 1, i, t] <- 0
      E[2, 2, i, t] <- KHI[i, t+1] + (1 - KHI[i, t+1]) * p
      E[2, 3, i, t] <- 0
      E[2, 4, i, t] <- 0
      E[2, 5, i, t] <- 0
      E[2, 6, i, t] <- 0
      E[2, 7, i, t] <- 0
      E[2, 8, i, t] <- 0
      E[2, 9, i, t] <- 0
      E[2, 10, i, t] <- 0
      E[2, 11, i, t] <- 0
      E[2, 12, i, t] <- 0
      E[2, 13, i, t] <- 0
      E[2, 14, i, t] <- 0
      E[2, 15, i, t] <- 0
      E[2, 16, i, t] <- 0
      E[2, 17, i, t] <- 0
      E[2, 18, i, t] <- 0
      E[2, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p)

      E[3, 1, i, t] <- 0
      E[3, 2, i, t] <- 0
      E[3, 3, i, t] <- KHI[i, t+1] + (1 - KHI[i, t+1]) * p
      E[3, 4, i, t] <- 0
      E[3, 5, i, t] <- 0
      E[3, 6, i, t] <- 0
      E[3, 7, i, t] <- 0
      E[3, 8, i, t] <- 0
      E[3, 9, i, t] <- 0
      E[3, 10, i, t] <- 0
      E[3, 11, i, t] <- 0
      E[3, 12, i, t] <- 0
      E[3, 13, i, t] <- 0
      E[3, 14, i, t] <- 0
      E[3, 15, i, t] <- 0
      E[3, 16, i, t] <- 0
      E[3, 17, i, t] <- 0
      E[3, 18, i, t] <- 0
      E[3, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p)

      E[4, 1, i, t] <- 0
      E[4, 2, i, t] <- 0
      E[4, 3, i, t] <- 0
      E[4, 4, i, t] <- (1 - KHI[i, t+1]) * p
      E[4, 5, i, t] <- 0
      E[4, 6, i, t] <- 0
      E[4, 7, i, t] <- 0
      E[4, 8, i, t] <- 0
      E[4, 9, i, t] <- 0
      E[4, 10, i, t] <- 0
      E[4, 11, i, t] <- 0
      E[4, 12, i, t] <- KHI[i, t+1] * (1-p) + KHI[i, t+1] * p
      E[4, 13, i, t] <- 0
      E[4, 14, i, t] <- 0
      E[4, 15, i, t] <- 0
      E[4, 16, i, t] <- 0
      E[4, 17, i, t] <- 0
      E[4, 18, i, t] <- 0
      E[4, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p)

      E[5, 1, i, t] <- 0
      E[5, 2, i, t] <- 0
      E[5, 3, i, t] <- 0
      E[5, 4, i, t] <- (1 - KHI[i, t+1]) * p
      E[5, 5, i, t] <- 0
      E[5, 6, i, t] <- 0
      E[5, 7, i, t] <- 0
      E[5, 8, i, t] <- 0
      E[5, 9, i, t] <- 0
      E[5, 10, i, t] <- 0
      E[5, 11, i, t] <- 0
      E[5, 12, i, t] <- 0
      E[5, 13, i, t] <- KHI[i, t+1] * (1 - p)
      E[5, 14, i, t] <- 0
      E[5, 15, i, t] <- 0
      E[5, 16, i, t] <- KHI[i, t+1] * p
      E[5, 17, i, t] <- 0
      E[5, 18, i, t] <- 0
      E[5, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p)

      E[6, 1, i, t] <- 0
      E[6, 2, i, t] <- 0
      E[6, 3, i, t] <- 0
      E[6, 4, i, t] <- 0
      E[6, 5, i, t] <- p
      E[6, 6, i, t] <- 0
      E[6, 7, i, t] <- 0
      E[6, 8, i, t] <- 0
      E[6, 9, i, t] <- 0
      E[6, 10, i, t] <- 0
      E[6, 11, i, t] <- 0
      E[6, 12, i, t] <- 0
      E[6, 13, i, t] <- KHI[i, t+1] * (1 - p) * AGE_5[i, t+1]
      E[6, 14, i, t] <- 0
      E[6, 15, i, t] <- KHI[i, t+1] * (1 - p) * (1 - AGE_5[i, t+1])
      E[6, 16, i, t] <- 0
      E[6, 17, i, t] <- 0
      E[6, 18, i, t] <- 0
      E[6, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p)

      E[7, 1, i, t] <- 0
      E[7, 2, i, t] <- 0
      E[7, 3, i, t] <- 0
      E[7, 4, i, t] <- 0
      E[7, 5, i, t] <- 0
      E[7, 6, i, t] <- p
      E[7, 7, i, t] <- 0
      E[7, 8, i, t] <- 0
      E[7, 9, i, t] <- 0
      E[7, 10, i, t] <- 0
      E[7, 11, i, t] <- 0
      E[7, 12, i, t] <- 0
      E[7, 13, i, t] <- KHI[i, t+1] * (1 - p) * AGE_5[i, t+1]
      E[7, 14, i, t] <- 0
      E[7, 15, i, t] <- KHI[i, t+1] * (1 - p) * (1 - AGE_5[i, t+1])
      E[7, 16, i, t] <- 0
      E[7, 17, i, t] <- 0
      E[7, 18, i, t] <- 0
      E[7, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p)

      E[8, 1, i, t] <- 0
      E[8, 2, i, t] <- 0
      E[8, 3, i, t] <- 0
      E[8, 4, i, t] <- 0
      E[8, 5, i, t] <- 0
      E[8, 6, i, t] <- 0
      E[8, 7, i, t] <- p
      E[8, 8, i, t] <- 0
      E[8, 9, i, t] <- 0
      E[8, 10, i, t] <- 0
      E[8, 11, i, t] <- 0
      E[8, 12, i, t] <- 0
      E[8, 13, i, t] <- 0
      E[8, 14, i, t] <- KHI[i, t+1] * (1 - p)
      E[8, 15, i, t] <- 0
      E[8, 16, i, t] <- 0
      E[8, 17, i, t] <- 0
      E[8, 18, i, t] <- 0
      E[8, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p)

      E[9, 1, i, t] <- 0
      E[9, 2, i, t] <- 0
      E[9, 3, i, t] <- 0
      E[9, 4, i, t] <- 0
      E[9, 5, i, t] <- 0
      E[9, 6, i, t] <- 0
      E[9, 7, i, t] <- 0
      E[9, 8, i, t] <- p
      E[9, 9, i, t] <- 0
      E[9, 10, i, t] <- 0
      E[9, 11, i, t] <- 0
      E[9, 12, i, t] <- 0
      E[9, 13, i, t] <- 0
      E[9, 14, i, t] <- KHI[i, t+1] * (1 - p)
      E[9, 15, i, t] <- 0
      E[9, 16, i, t] <- 0
      E[9, 17, i, t] <- 0
      E[9, 18, i, t] <- 0
      E[9, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p)

      E[10, 1, i, t] <- 0
      E[10, 2, i, t] <- 0
      E[10, 3, i, t] <- 0
      E[10, 4, i, t] <- 0
      E[10, 5, i, t] <- 0
      E[10, 6, i, t] <- 0
      E[10, 7, i, t] <- 0
      E[10, 8, i, t] <- 0
      E[10, 9, i, t] <- (1 - a) * p
      E[10, 10, i, t] <- 0
      E[10, 11, i, t] <- a * (1 - KHI[i, t+1]) * p
      E[10, 12, i, t] <- 0
      E[10, 13, i, t] <- 0
      E[10, 14, i, t] <- KHI[i, t+1] * (1 - p)
      E[10, 15, i, t] <- 0
      E[10, 16, i, t] <- 0
      E[10, 17, i, t] <- a * KHI[i, t+1] * p
      E[10, 18, i, t] <- 0
      E[10, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p)

      E[11, 1, i, t] <- 0
      E[11, 2, i, t] <- 0
      E[11, 3, i, t] <- 0
      E[11, 4, i, t] <- 0
      E[11, 5, i, t] <- 0
      E[11, 6, i, t] <- 0
      E[11, 7, i, t] <- 0
      E[11, 8, i, t] <- 0
      E[11, 9, i, t] <- 2 * (1 - a) * a * p
      E[11, 10, i, t] <- (1 - a)^2 * p
      E[11, 11, i, t] <- a^2 * (1 - KHI[i, t+1]) * p
      E[11, 12, i, t] <- 0
      E[11, 13, i, t] <- 0
      E[11, 14, i, t] <- KHI[i, t+1] * (1 - p)
      E[11, 15, i, t] <- 0
      E[11, 16, i, t] <- 0
      E[11, 17, i, t] <- a^2 * KHI[i, t+1] * p
      E[11, 18, i, t] <- 0
      E[11, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p)

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
      E[12, 11, i, t] <- (1 - KHI[i, t+1]) * p
      E[12, 12, i, t] <- 0
      E[12, 13, i, t] <- 0
      E[12, 14, i, t] <- KHI[i, t+1] * (1 - p)
      E[12, 15, i, t] <- 0
      E[12, 16, i, t] <- 0
      E[12, 17, i, t] <- KHI[i, t+1] * p
      E[12, 18, i, t] <- 0
      E[12, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p)

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
      E[13, 11, i, t] <- (1 - KHI[i, t+1]) * p
      E[13, 12, i, t] <- 0
      E[13, 13, i, t] <- 0
      E[13, 14, i, t] <- 0
      E[13, 15, i, t] <- KHI[i, t+1] * (1 - p)
      E[13, 16, i, t] <- 0
      E[13, 17, i, t] <- 0
      E[13, 18, i, t] <- KHI[i, t+1] * p
      E[13, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p)

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
      E[14, 12, i, t] <- 0
      E[14, 13, i, t] <- 0
      E[14, 14, i, t] <- 0
      E[14, 15, i, t] <- 0
      E[14, 16, i, t] <- 0
      E[14, 17, i, t] <- 0
      E[14, 18, i, t] <- 0
      E[14, 19, i, t] <- 1
    }
  }

  # ++++++++++++++++++++++++++++++++ priors ++++++++++++++++++++++++++++++++++++

  #-------- State process --------#
  for (k in 1:3) {
    beta_phi[k] ~ dnorm(0, sd = 1.5)  # juvenile, subadult, and adult survival
    beta_s0[k] ~ dnorm(0, sd = 1.5)   # cub survival
    beta_s1[k] ~ dnorm(0, sd = 1.5)   # yearling survival
    beta_beta[k] ~ dnorm(0, sd = 1.5) # early litter survival
    beta_gamma[k] ~ dnorm(0, sd = 1.5) # twinning
  }
  beta_phi[4] ~ dnorm(0, sd = 3)   # cub survival (additive parameter)
  beta_s0[4] ~ dnorm(0, sd = 3)   # cub survival (additive parameter)
  beta_s0[5] ~ dnorm(0, sd = 3)   # cub survival (additive parameter)
  beta_s1[4] ~ dnorm(0, sd = 3)   # cub survival (additive parameter)
  beta_beta[4] ~ dnorm(0, sd = 3)   # early litter survival (additive parameter)
  beta_gamma[4] ~ dnorm(0, sd = 3)   # twinning (additive parameter)
  
  for (k in 1:4) {
    beta_eta[k] ~ dnorm(0, sd = 1.5)  # denning probability
  }
  for (k in 5:9) {
    beta_eta[k] ~ dnorm(0, sd = 3)  # denning probability (additive parameters)
  }
  


  #----- observation process -----#
  p ~ dunif(0, 1)


  # +++++++++++++++++++++++++++ model & likelihood +++++++++++++++++++++++++++++

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
          E[k, CH[i, t+1], i, t]
      }
    }
    lik[i] <- sum(zeta[i, n_occasions, 1:n_states])
    ones[i] ~ dbin(prob = lik[i], size = freq[i])
  }
})


# B. Generate datasets =========================================================

# ~ 1. Define parameters -------------------------------------------------------

{
new_per_year <- 25
n_occasions <- 12

n_states <- 14
not_captured <- 19

# Parameters of the linear regression of the probability of having parted from mother
# for a 2-yr old against day of capture
intercept_alpha <- -3.046745
slope_alpha <- 0.03542667
a <- plogis(intercept_alpha + slope_alpha*105)  # All captures are assumed to occur on the same date (day of year 105)

# Observed proportions of each state (state is inferred by me) at first capture
#                   1      2      3      4      5      6      7      8      9      10    11      12    13   14    
f_state_prop <- c(0.098, 0.052, 0.045, 0.035, 0.029, 0.082, 0.176, 0.068, 0.027, 0.021, 0.006, 0.188, 0.173, 0)


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
                nrow = n_states, ncol = 24, byrow = TRUE)
apply(AGE_f, 1, sum)


# Probability of being collared (matrix with state in line and collaring outcome
# - yes or no - in column)
# Here I use a probability of being tagged that seems reasonable given Jon told me
# he deploys 10-20 collars each year. If we say that ~ 50 bears are captured each 
# year, of which at most 25 are independent females, then 25*0.4 = 10 collars. In 
# the simulated datasets, more females than that are captured each year (25 new marks
# + recaptures). If we consider that the same proportion of female captured is 
# equipped with collars in the simulated datasets (more captures per year), then I 
# can keep 0.4. 
prob_GPS_A <- 0.4 
prob_GPS_SA4 <- 0.4 

#                                 no            yes
COLLARING <- matrix(data = c(      1,               0,       # J2    
                                   1,               0,       # J3   
                              1-prob_GPS_SA4, prob_GPS_SA4,  # SA4  
                              1- prob_GPS_A,   prob_GPS_A,   # NDSA5  
                              1- prob_GPS_A,   prob_GPS_A,   # FDSA5   
                              1- prob_GPS_A,   prob_GPS_A,   # A01  
                              1- prob_GPS_A,   prob_GPS_A,   # A02  
                              1- prob_GPS_A,   prob_GPS_A,   # A11
                              1- prob_GPS_A,   prob_GPS_A,   # A12
                              1- prob_GPS_A,   prob_GPS_A,   # AS1  
                              1- prob_GPS_A,   prob_GPS_A,   # AS2  
                              1- prob_GPS_A,   prob_GPS_A,   # LNDA    
                              1- prob_GPS_A,   prob_GPS_A,   # FDA
                                   1,              0),         # D
                    nrow = 14, ncol = 2, byrow = T)

# Duration of collars,
#                     1 year                          # 4 years
collar_duration <- c(0.5384615, 0.2955466, 0.1174089, 0.0485830)


# GLS tagging and duration
# Here I use a probability of being equipped with a GLS tags that reflects the 
# fact that almost all females are equipped since 2013 or so. In some years, they
# run out of GLS tags, so not all females are equipped. Also, some tags fail entirely
# or the data can't be accessed. Thus, I use:
prob_GLS_tagging <- 0.8 
# A tag as 1 chance out of three of lasting 3, 4 or 5 years (validated by Jon)
#                 3yrs 4yrs 5yrs
GLS_duration <- c(1/3, 1/3, 1/3)

# Survival
beta_phi <- NULL
beta_phi[1] <- logit(0.92)    # intercept 2-4yr 
beta_phi[2] <- logit(0.96)   # intercept 5-15yr 
beta_phi[3] <- logit(0.85)  # intercept >15yr 
beta_phi[4] <- -0.2         # SI

# Age dependent reproductive rates 
beta_s0 <- beta_s1 <- beta_eta <- beta_beta <- beta_gamma <-  NULL

# s0
beta_s0[1] = logit(0.35)  # intercept young 
beta_s0[2] = logit(0.55)  # intercept middle aged 
beta_s0[3] = logit(0.45)  # intercept old 
beta_s0[4] = -0.3    # effect of being singleton (i.e. ref is twin)
beta_s0[5] = -0.2   # SI
# s1
beta_s1[1] = logit(0.75)  # intercept young 
beta_s1[2] = logit(0.8)   # intercept middle aged 
beta_s1[3] = logit(0.85)   # intercept old 
beta_s1[4] = -0.2   # SI

# eta
beta_eta[1] = logit(0.25) # intercept SA4 
beta_eta[2] = logit(0.70)  # intercept young 
beta_eta[3] = logit(0.80)  # intercept middle aged 
beta_eta[4] = logit(0.75)   # intercept old 
beta_eta[5] = -0.25 # cost of reproduction FDA
beta_eta[6] = -1.3  # cost of reproduction A0-
beta_eta[7] = -0.5  # cost of reproduction A1-
beta_eta[8] = 0.35  # cost of reproduction (individual quality) AS
beta_eta[9] = -0.2  # SI

# beta
beta_beta[1] = logit(0.65)   # intercept 4-9yr 
beta_beta[2] = logit(0.75)   # intercept middle aged 
beta_beta[3] = logit(0.65)   # intercept old  
beta_beta[4] = -0.4          # SI
# gamma
beta_gamma[1] = logit(0.65)  # intercept 4-9yr 
beta_gamma[2] = logit(0.75)  # intercept middle aged 
beta_gamma[3] = logit(0.70)  # intercept old  
beta_gamma[4] = -0.3  # SI


# Recapture probability 
p <- 0.25 

f <- rep(c(1:n_occasions), each = new_per_year)

set.seed(1)
sea_ice_s <- as.vector(scale(rnorm(n = n_occasions-1, sd = 1)))
}


# ~ 2. Generate capture histories ----------------------------------------------

for (l in 1:n_sim) {
  print(paste0("simulation ", l, ":"))
  
  # Arrays to store the individual- and time-specific transition and observation matrices
  S <- array(data = NA, dim = c(14, 14, (new_per_year*n_occasions - new_per_year), n_occasions))
  E <- E_p <- array(data = NA, dim = c(14, 19, (new_per_year*n_occasions - new_per_year), n_occasions))
  
  # Matrices for vital rates
  phi <-  s01 <- s02 <- s1 <- beta <- gamma <- p_realized <- matrix(NA, nrow = (new_per_year*n_occasions - new_per_year), 
                                                                                        ncol = n_occasions)
  eta <- array(data = NA, dim = c(5, (new_per_year*n_occasions - new_per_year), n_occasions))
  
  # Matrices for data
  CH <- Z <- AGE <- matrix(NA, nrow = (new_per_year*n_occasions - new_per_year), 
                                           ncol = n_occasions) 
  
  PI <- RHO <- TAU <- KHI <- matrix(0, nrow = (new_per_year*n_occasions - new_per_year), 
                                    ncol = n_occasions) 
  f_event <- f_state <-  NULL
  
  collaring_save <- matrix(0, nrow = (new_per_year*n_occasions - new_per_year), 
                           ncol = n_occasions) 
  # ~~~ a. First captures --------------------------------------------------------
  
  set.seed(l)
  for (i in 1:(new_per_year*n_occasions - new_per_year)) {
    Z[i, f[i]] <- f_state[i] <- sample(x = 1:14,
                                       size = 1,
                                       prob = f_state_prop, replace = T)
    
    # Age at first capture
    AGE[i, f[i]] <- which(rmultinom(n = 1, size = 1,
                                    prob = AGE_f[Z[i, f[i]], ]) == 1) + 1 # +1 because J2 (state = 1) are 2yr old
    
    #               |-------------------- physical only ------------------------||------- GPS only ------||----- Both -----|| None |
    #               J2  J3  SA4 SA5 A01  A02  A11  A12     A21       A22     A     NDSA5 DSA5   NDA   DA   FDSA5 LNDA  FDA     NA
    E0 <- matrix(data =
                   c(1,  0,  0,  0,  0,   0,   0,   0,      0,        0,     0,      0,    0,    0,   0,     0,    0,   0,     0,    # J2
                     0,  1,  0,  0,  0,   0,   0,   0,      0,        0,     0,      0,    0,    0,   0,     0,    0,   0,     0,    # J3
                     0,  0,  1,  0,  0,   0,   0,   0,      0,        0,     0,      0,    0,    0,   0,     0,    0,   0,     0,    # SA4
                     0,  0,  0,  1,  0,   0,   0,   0,      0,        0,     0,      0,    0,    0,   0,     0,    0,   0,     0,    # NDSA5
                     0,  0,  0,  1,  0,   0,   0,   0,      0,        0,     0,      0,    0,    0,   0,     0,    0,   0,     0,    # FDSA5  # 5
                     0,  0,  0,  0,  1,   0,   0,   0,      0,        0,     0,      0,    0,    0,   0,     0,    0,   0,     0,    # A01
                     0,  0,  0,  0,  0,   1,   0,   0,      0,        0,     0,      0,    0,    0,   0,     0,    0,   0,     0,    # A02
                     0,  0,  0,  0,  0,   0,   1,   0,      0,        0,     0,      0,    0,    0,   0,     0,    0,   0,     0,    # A11    # 10
                     0,  0,  0,  0,  0,   0,   0,   1,      0,        0,     0,      0,    0,    0,   0,     0,    0,   0,     0,    # A12
                     0,  0,  0,  0,  0,   0,   0,   0,    (1-a),      0,     a,      0,    0,    0,   0,     0,    0,   0,     0,    # AS1
                     0,  0,  0,  0,  0,   0,   0,   0,  2*(1-a)*a, (1-a)^2, a^2,     0,    0,    0,   0,     0,    0,   0,     0,    # AS2
                     0,  0,  0,  0,  0,   0,   0,   0,      0,        0,     1,      0,    0,    0,   0,     0,    0,   0,     0,    # LNDA
                     0,  0,  0,  0,  0,   0,   0,   0,      0,        0,     1,      0,    0,    0,   0,     0,    0,   0,     0,    # FDA    # 15
                     0,  0,  0,  0,  0,   0,   0,   0,      0,        0,     0,      0,    0,    0,   0,     0,    0,   0,     1),   # D
                 nrow = 14, ncol = 19, byrow = TRUE)
    
    CH[i, f[i]] <- f_event[i] <- which(rmultinom(n = 1, size = 1, 
                                                                  prob = E0[Z[i, f[i]], ]) == 1)
    # GPS collaring
    collared <- collaring_save[i, f[i]] <- which(rmultinom(n = 1, size = 1, 
                                prob = COLLARING[Z[i, f[i]], ]) == 1) - 1
    if (collared == 1) { # If individual has a collar, for how many years ?
      collared_years <- sample(x = 1:4, size = 1, prob = collar_duration)
      if (collared_years >= 1) {
        if (f[i]+1 <= n_occasions) {
          PI[i, f[i]+1] <- 1
        } 
      }
      if (collared_years >= 2) {
        if (f[i]+2 <= n_occasions) {
          PI[i, f[i]+2] <- 1
        }
      }
      if (collared_years >= 3) {
        if (f[i]+3 <= n_occasions) {
          PI[i, f[i]+3] <- 1
        }
      }
    }
    
    # GLS tagging
    # Will the individual be tagged, and how many years will the tag record ?
    tag <- rbinom(n = 1, size = 1, prob = prob_GLS_tagging) * 
      sample(x = c(3, 4, 5), size = 1, prob = GLS_duration)

    if (tag >= 3) {
      if (f[i]+1 <= n_occasions) {
        RHO[i, f[i]+1] <- 1
      } 
      if (f[i]+2 <= n_occasions) {
        RHO[i, f[i]+2] <- 1
      } 
      if (f[i]+3 <= n_occasions) {
        RHO[i, f[i]+3] <- 1
      } 
    }
    if (tag >= 4) {
      if (f[i]+4 <= n_occasions) {
        RHO[i, f[i]+4] <- 1
      }
    }
    if (tag == 5) {
      if (f[i]+5 <= n_occasions) {
        RHO[i, f[i]+5] <- 1
      }
    }
  }
  
  # Compute matrices with age
  AGE_2_4 <- AGE_4 <- AGE_5 <- AGE_YOUNG <- AGE_MIDDLE <- AGE_OLD <- 
    AGE_YOUNG_1 <- AGE_MIDDLE_1 <- AGE_OLD_1 <- AGE_YOUNG_2 <- AGE_MIDDLE_2 <- AGE_OLD_2 <-    
    matrix(NA, nrow = (new_per_year*n_occasions - new_per_year), ncol = n_occasions) 
  for (i in 1:((new_per_year*n_occasions - new_per_year))) {
    AGE_2_4[i, f[i]] <- ifelse(AGE[i, f[i]] %in% 2:4, 1, 0)
    AGE_4[i, f[i]] <- ifelse(AGE[i, f[i]] == 4, 1, 0)
    AGE_5[i, f[i]] <- ifelse(AGE[i, f[i]] == 5, 1, 0)
    AGE_YOUNG[i, f[i]] <- ifelse(AGE[i, f[i]] %in% 5:8, 1, 0)
    AGE_MIDDLE[i, f[i]] <- ifelse(AGE[i, f[i]] %in% 9:14, 1, 0)
    AGE_OLD[i, f[i]] <- ifelse(AGE[i, f[i]] %in% 15:75, 1, 0)
    AGE_YOUNG_1[i, f[i]] <- ifelse(AGE[i, f[i]] %in% 5:9, 1, 0)
    AGE_MIDDLE_1[i, f[i]] <- ifelse(AGE[i, f[i]] %in% 10:15, 1, 0)
    AGE_OLD_1[i, f[i]] <- ifelse(AGE[i, f[i]] %in% 16:75, 1, 0)
    AGE_YOUNG_2[i, f[i]] <- ifelse(AGE[i, f[i]] %in% 5:10, 1, 0)
    AGE_MIDDLE_2[i, f[i]] <- ifelse(AGE[i, f[i]] %in% 11:16, 1, 0)
    AGE_OLD_2[i, f[i]] <- ifelse(AGE[i, f[i]] %in% 17:75, 1, 0)
    if (f[i] == n_occasions) next
    for (t in (f[i]+1):n_occasions) {
      AGE[i, t] <- AGE[i, t-1] +1
      AGE_2_4[i, t] <- ifelse(AGE[i, t] %in% 2:4, 1, 0)
      AGE_4[i, t] <- ifelse(AGE[i, t] == 4, 1, 0)
      AGE_5[i, t] <- ifelse(AGE[i, t] == 5, 1, 0)
      AGE_YOUNG[i, t] <- ifelse(AGE[i, t] %in% 5:8, 1, 0)
      AGE_MIDDLE[i, t] <- ifelse(AGE[i, t] %in% 9:14, 1, 0)
      AGE_OLD[i, t] <- ifelse(AGE[i, t] %in% 15:75, 1, 0)
      AGE_YOUNG_1[i, t] <- ifelse(AGE[i, t] %in% 5:9, 1, 0)
      AGE_MIDDLE_1[i, t] <- ifelse(AGE[i, t] %in% 10:15, 1, 0)
      AGE_OLD_1[i, t] <- ifelse(AGE[i, t] %in% 16:75, 1, 0)
      AGE_YOUNG_2[i, t] <- ifelse(AGE[i, t] %in% 5:10, 1, 0)
      AGE_MIDDLE_2[i, t] <- ifelse(AGE[i, t] %in% 11:16, 1, 0)
      AGE_OLD_2[i, t] <- ifelse(AGE[i, t] %in% 17:75, 1, 0)
    }
  }
  
  # ~~~ b. Subsequent captures ---------------------------------------------------
  
  set.seed(l)
  for (i in 1:(new_per_year*n_occasions - new_per_year)) {
    if (f[i] == n_occasions) next
    
    for (t in f[i]:(n_occasions-1)) {
      
      # ~~~~~ Physical ---------------------------------------------------------
      phi[i, t] <- plogis(beta_phi[1] * AGE_2_4[i, t] +
                            beta_phi[2] * (AGE_YOUNG[i, t] + AGE_MIDDLE[i, t]) +
                            beta_phi[3] * AGE_OLD[i, t] +
                            beta_phi[4] * (AGE_2_4[i, t] + AGE_OLD[i, t]) * sea_ice_s[t])
      
      # Probability of coy survival
      s01[i, t] <- plogis(beta_s0[1] * AGE_YOUNG_1[i, t] +
                            beta_s0[2] * AGE_MIDDLE_1[i, t] +
                            beta_s0[3] * AGE_OLD_1[i, t] +
                            beta_s0[4] +
                            beta_s0[5] * sea_ice_s[t])
      
      s02[i, t] <- plogis(beta_s0[1] * AGE_YOUNG_1[i, t] +
                            beta_s0[2] * AGE_MIDDLE_1[i, t] +
                            beta_s0[3] * AGE_OLD_1[i, t] +
                            beta_s0[5] * sea_ice_s[t])

      # Probability of yearling survival
      s1[i, t] <- plogis(beta_s1[1] * AGE_YOUNG_2[i, t] +
                           beta_s1[2] * AGE_MIDDLE_2[i, t] +
                           beta_s1[3] * AGE_OLD_2[i, t] +
                           beta_s1[4] * sea_ice_s[t])
      
      # Probability of denning
      eta[1, i, t] <- plogis(beta_eta[1] * AGE_4[i, t] +   # Lone female
                               beta_eta[2] * AGE_YOUNG[i, t] +
                               beta_eta[3] * AGE_MIDDLE[i, t] +
                               beta_eta[4] * AGE_OLD[i, t] +
                               beta_eta[9] * sea_ice_s[t])
      
      eta[2, i, t] <- plogis(beta_eta[2] * AGE_YOUNG[i, t] + # FDA
                               beta_eta[3] * AGE_MIDDLE[i, t] +
                               beta_eta[4] * AGE_OLD[i, t] +
                               beta_eta[5] +
                               beta_eta[9] * sea_ice_s[t])
      
      eta[3, i, t] <- plogis(beta_eta[2] * AGE_YOUNG[i, t] + # A0-
                               beta_eta[3] * AGE_MIDDLE[i, t] +
                               beta_eta[4] * AGE_OLD[i, t] +
                               beta_eta[6] +
                               beta_eta[9] * sea_ice_s[t])
      
      eta[4, i, t] <- plogis(beta_eta[2] * AGE_YOUNG[i, t] + # A1-
                               beta_eta[3] * AGE_MIDDLE[i, t] +
                               beta_eta[4] * AGE_OLD[i, t] +
                               beta_eta[7] +
                               beta_eta[9] * sea_ice_s[t])
      
      eta[5, i, t] <- plogis(beta_eta[2] * AGE_YOUNG[i, t] + # AS
                               beta_eta[3] * AGE_MIDDLE[i, t] +
                               beta_eta[4] * AGE_OLD[i, t] +
                               beta_eta[8] +
                               beta_eta[9] * sea_ice_s[t])
      
      # Early litter survival
      beta[i, t] <- plogis(beta_beta[1] * (AGE_4[i, t] + AGE_YOUNG[i, t]) +
        beta_beta[2] * AGE_MIDDLE[i, t] +
        beta_beta[3] * AGE_OLD[i, t] +
        beta_beta[4] * sea_ice_s[t])
      
      # twinning prob
      gamma[i, t] <- plogis(beta_gamma[1] * (AGE_4[i, t] + AGE_YOUNG[i, t]) +
        beta_gamma[2] * AGE_MIDDLE[i, t] +
        beta_gamma[3] * AGE_OLD[i, t] +
        beta_gamma[4] * sea_ice_s[t])

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
      S[3, 4, i, t] <- phi[i, t] * (1- eta[1, i, t])
      S[3, 5, i, t] <- phi[i, t] * eta[1, i, t] * (1-beta[i, t])
      S[3, 6, i, t] <- phi[i, t] * eta[1, i, t] * beta[i, t] * (1 - gamma[i, t])
      S[3, 7, i, t] <- phi[i, t] * eta[1, i, t] * beta[i, t] * gamma[i, t]
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
      S[4, 6, i, t] <- phi[i, t] * eta[1, i, t] * beta[i, t] * (1 - gamma[i, t])
      S[4, 7, i, t] <- phi[i, t] * eta[1, i, t] * beta[i, t] * gamma[i, t]
      S[4, 8, i, t] <- 0
      S[4, 9, i, t] <- 0
      S[4, 10, i, t] <- 0
      S[4, 11, i, t] <- 0
      S[4, 12, i, t] <- phi[i, t] * (1 - eta[1, i, t])
      S[4, 13, i, t] <- phi[i, t] * eta[1, i, t] * (1-beta[i, t])
      S[4, 14, i, t] <- 1-phi[i, t]
      
      S[5, 1, i, t] <- 0
      S[5, 2, i, t] <- 0
      S[5, 3, i, t] <- 0
      S[5, 4, i, t] <- 0
      S[5, 5, i, t] <- 0
      S[5, 6, i, t] <- phi[i, t] * eta[2, i, t] * beta[i, t] * (1 - gamma[i, t])
      S[5, 7, i, t] <- phi[i, t] * eta[2, i, t] * beta[i, t] * gamma[i, t]
      S[5, 8, i, t] <- 0
      S[5, 9, i, t] <- 0
      S[5, 10, i, t] <- 0
      S[5, 11, i, t] <- 0
      S[5, 12, i, t] <- phi[i, t] * (1 - eta[2, i, t])
      S[5, 13, i, t] <- phi[i, t] * eta[2, i, t] * (1-beta[i, t])
      S[5, 14, i, t] <- 1-phi[i, t]
      
      S[6, 1, i, t] <- 0
      S[6, 2, i, t] <- 0
      S[6, 3, i, t] <- 0
      S[6, 4, i, t] <- 0
      S[6, 5, i, t] <- 0
      S[6, 6, i, t] <- phi[i, t] * eta[3, i, t] * (1-s01[i, t]) * beta[i, t] * (1 - gamma[i, t])
      S[6, 7, i, t] <- phi[i, t] * eta[3, i, t] * (1-s01[i, t]) * beta[i, t] * gamma[i, t]
      S[6, 8, i, t] <- phi[i, t] * s01[i, t]
      S[6, 9, i, t] <- 0
      S[6, 10, i, t] <- 0
      S[6, 11, i, t] <- 0
      S[6, 12, i, t] <- phi[i, t] * (1-s01[i, t]) * (1-eta[3, i, t])
      S[6, 13, i, t] <- phi[i, t] * eta[3, i, t] * (1-s01[i, t]) * (1-beta[i, t])
      S[6, 14, i, t] <- 1-phi[i, t]
      
      S[7, 1, i, t] <- 0
      S[7, 2, i, t] <- 0
      S[7, 3, i, t] <- 0
      S[7, 4, i, t] <- 0
      S[7, 5, i, t] <- 0
      S[7, 6, i, t] <- phi[i, t] * eta[3, i, t] * (1-s02[i, t])^2 * beta[i, t] * (1 - gamma[i, t])
      S[7, 7, i, t] <- phi[i, t] * eta[3, i, t] * (1-s02[i, t])^2 * beta[i, t] * gamma[i, t]
      S[7, 8, i, t] <- phi[i, t] * 2 * s02[i, t] * (1-s02[i, t])
      S[7, 9, i, t] <- phi[i, t] * s02[i, t]^2
      S[7, 10, i, t] <- 0
      S[7, 11, i, t] <- 0
      S[7, 12, i, t] <- phi[i, t] * (1-s02[i, t])^2 * (1-eta[3, i, t])
      S[7, 13, i, t] <- phi[i, t] * eta[3, i, t] * (1-s02[i, t])^2 * (1-beta[i, t])
      S[7, 14, i, t] <- 1-phi[i, t]
      
      S[8, 1, i, t] <- 0
      S[8, 2, i, t] <- 0
      S[8, 3, i, t] <- 0
      S[8, 4, i, t] <- 0
      S[8, 5, i, t] <- 0
      S[8, 6, i, t] <- phi[i, t] * eta[4, i, t] * (1-s1[i, t]) * beta[i, t] * (1 - gamma[i, t])
      S[8, 7, i, t] <- phi[i, t] * eta[4, i, t] * (1-s1[i, t]) * beta[i, t] * gamma[i, t]
      S[8, 8, i, t] <- 0
      S[8, 9, i, t] <- 0
      S[8, 10, i, t] <- phi[i, t] * s1[i, t]
      S[8, 11, i, t] <- 0
      S[8, 12, i, t] <- phi[i, t] * (1-s1[i, t]) * (1-eta[4, i, t])
      S[8, 13, i, t] <- phi[i, t] * eta[4, i, t] * (1-s1[i, t]) * (1-beta[i, t])
      S[8, 14, i, t] <- 1-phi[i, t]
      
      S[9, 1, i, t] <- 0
      S[9, 2, i, t] <- 0
      S[9, 3, i, t] <- 0
      S[9, 4, i, t] <- 0
      S[9, 5, i, t] <- 0
      S[9, 6, i, t] <- phi[i, t] * eta[4, i, t] * (1-s1[i, t])^2 * beta[i, t] * (1 - gamma[i, t])   
      S[9, 7, i, t] <- phi[i, t] * eta[4, i, t] * (1-s1[i, t])^2 * beta[i, t] * gamma[i, t]
      S[9, 8, i, t] <- 0
      S[9, 9, i, t] <- 0
      S[9, 10, i, t] <- phi[i, t] * 2 * s1[i, t] * (1-s1[i, t])
      S[9, 11, i, t] <- phi[i, t] * s1[i, t]^2
      S[9, 12, i, t] <- phi[i, t] * (1-s1[i, t])^2 * (1-eta[4, i, t])
      S[9, 13, i, t] <- phi[i, t] * eta[4, i, t] * (1-s1[i, t])^2 * (1-beta[i, t])
      S[9, 14, i, t] <- 1-phi[i, t]
      
      S[10, 1, i, t] <- 0
      S[10, 2, i, t] <- 0
      S[10, 3, i, t] <- 0
      S[10, 4, i, t] <- 0
      S[10, 5, i, t] <- 0
      S[10, 6, i, t] <- phi[i, t] * eta[5, i, t] * beta[i, t] * (1 - gamma[i, t])
      S[10, 7, i, t] <- phi[i, t] * eta[5, i, t] * beta[i, t] * gamma[i, t]
      S[10, 8, i, t] <- 0
      S[10, 9, i, t] <- 0
      S[10, 10, i, t] <- 0
      S[10, 11, i, t] <- 0
      S[10, 12, i, t] <- phi[i, t] * (1 - eta[5, i, t])
      S[10, 13, i, t] <- phi[i, t] * eta[5, i, t] * (1-beta[i, t])
      S[10, 14, i, t] <- 1-phi[i, t]
      
      S[11, 1, i, t] <- 0
      S[11, 2, i, t] <- 0
      S[11, 3, i, t] <- 0
      S[11, 4, i, t] <- 0
      S[11, 5, i, t] <- 0
      S[11, 6, i, t] <- phi[i, t] * eta[5, i, t] * beta[i, t] * (1 - gamma[i, t])
      S[11, 7, i, t] <- phi[i, t] * eta[5, i, t] * beta[i, t] * gamma[i, t]
      S[11, 8, i, t] <- 0
      S[11, 9, i, t] <- 0
      S[11, 10, i, t] <- 0
      S[11, 11, i, t] <- 0
      S[11, 12, i, t] <- phi[i, t] * (1 - eta[5, i, t])
      S[11, 13, i, t] <- phi[i, t] * eta[5, i, t] * (1-beta[i, t])
      S[11, 14, i, t] <- 1-phi[i, t]
      
      S[12, 1, i, t] <- 0
      S[12, 2, i, t] <- 0
      S[12, 3, i, t] <- 0
      S[12, 4, i, t] <- 0
      S[12, 5, i, t] <- 0
      S[12, 6, i, t] <- phi[i, t] * eta[1, i, t] * beta[i, t] * (1 - gamma[i, t])
      S[12, 7, i, t] <- phi[i, t] * eta[1, i, t] * beta[i, t] * gamma[i, t]
      S[12, 8, i, t] <- 0
      S[12, 9, i, t] <- 0
      S[12, 10, i, t] <- 0
      S[12, 11, i, t] <- 0
      S[12, 12, i, t] <- phi[i, t] * (1 - eta[1, i, t])
      S[12, 13, i, t] <- phi[i, t] * eta[1, i, t] * (1-beta[i, t])
      S[12, 14, i, t] <- 1-phi[i, t]
      
      S[13, 1, i, t] <- 0
      S[13, 2, i, t] <- 0
      S[13, 3, i, t] <- 0
      S[13, 4, i, t] <- 0
      S[13, 5, i, t] <- 0
      S[13, 6, i, t] <- phi[i, t] * eta[2, i, t] * beta[i, t] * (1 - gamma[i, t])
      S[13, 7, i, t] <- phi[i, t] * eta[2, i, t] * beta[i, t] * gamma[i, t]
      S[13, 8, i, t] <- 0
      S[13, 9, i, t] <- 0
      S[13, 10, i, t] <- 0
      S[13, 11, i, t] <- 0
      S[13, 12, i, t] <- phi[i, t] * (1 - eta[2, i, t])
      S[13, 13, i, t] <- phi[i, t] * eta[2, i, t] * (1-beta[i, t])
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
      
      p_realized[i, t+1] <- rbinom(n = 1, size = 1, prob = p)
      
      
      E_p[1, 1, i, t] <- p_realized[i, t+1]
      E_p[1, 2, i, t] <- 0
      E_p[1, 3, i, t] <- 0
      E_p[1, 4, i, t] <- 0
      E_p[1, 5, i, t] <- 0
      E_p[1, 6, i, t] <- 0
      E_p[1, 7, i, t] <- 0
      E_p[1, 8, i, t] <- 0
      E_p[1, 9, i, t] <- 0
      E_p[1, 10, i, t] <- 0
      E_p[1, 11, i, t] <- 0
      E_p[1, 12, i, t] <- 0
      E_p[1, 13, i, t] <- 0
      E_p[1, 14, i, t] <- 0
      E_p[1, 15, i, t] <- 0
      E_p[1, 16, i, t] <- 0
      E_p[1, 17, i, t] <- 0
      E_p[1, 18, i, t] <- 0
      E_p[1, 19, i, t] <- 1-p_realized[i, t+1]
      
      E_p[2, 1, i, t] <- 0
      E_p[2, 2, i, t] <- p_realized[i, t+1]
      E_p[2, 3, i, t] <- 0
      E_p[2, 4, i, t] <- 0
      E_p[2, 5, i, t] <- 0
      E_p[2, 6, i, t] <- 0
      E_p[2, 7, i, t] <- 0
      E_p[2, 8, i, t] <- 0
      E_p[2, 9, i, t] <- 0
      E_p[2, 10, i, t] <- 0
      E_p[2, 11, i, t] <- 0
      E_p[2, 12, i, t] <- 0
      E_p[2, 13, i, t] <- 0
      E_p[2, 14, i, t] <- 0
      E_p[2, 15, i, t] <- 0
      E_p[2, 16, i, t] <- 0
      E_p[2, 17, i, t] <- 0
      E_p[2, 18, i, t] <- 0
      E_p[2, 19, i, t] <- (1 - p_realized[i, t+1])
      
      E_p[3, 1, i, t] <- 0
      E_p[3, 2, i, t] <- 0
      E_p[3, 3, i, t] <- p_realized[i, t+1]
      E_p[3, 4, i, t] <- 0
      E_p[3, 5, i, t] <- 0
      E_p[3, 6, i, t] <- 0
      E_p[3, 7, i, t] <- 0
      E_p[3, 8, i, t] <- 0
      E_p[3, 9, i, t] <- 0
      E_p[3, 10, i, t] <- 0
      E_p[3, 11, i, t] <- 0
      E_p[3, 12, i, t] <- 0
      E_p[3, 13, i, t] <- 0
      E_p[3, 14, i, t] <- 0
      E_p[3, 15, i, t] <- 0
      E_p[3, 16, i, t] <- 0
      E_p[3, 17, i, t] <- 0
      E_p[3, 18, i, t] <- 0
      E_p[3, 19, i, t] <- (1 - p_realized[i, t+1])
      
      E_p[4, 1, i, t] <- 0
      E_p[4, 2, i, t] <- 0
      E_p[4, 3, i, t] <- 0
      E_p[4, 4, i, t] <- p_realized[i, t+1]
      E_p[4, 5, i, t] <- 0
      E_p[4, 6, i, t] <- 0
      E_p[4, 7, i, t] <- 0
      E_p[4, 8, i, t] <- 0
      E_p[4, 9, i, t] <- 0
      E_p[4, 10, i, t] <- 0
      E_p[4, 11, i, t] <- 0
      E_p[4, 12, i, t] <- 0
      E_p[4, 13, i, t] <- 0
      E_p[4, 14, i, t] <- 0
      E_p[4, 15, i, t] <- 0
      E_p[4, 16, i, t] <- 0
      E_p[4, 17, i, t] <- 0
      E_p[4, 18, i, t] <- 0
      E_p[4, 19, i, t] <- (1 - p_realized[i, t+1])
      
      E_p[5, 1, i, t] <- 0
      E_p[5, 2, i, t] <- 0
      E_p[5, 3, i, t] <- 0
      E_p[5, 4, i, t] <- p_realized[i, t+1]
      E_p[5, 5, i, t] <- 0
      E_p[5, 6, i, t] <- 0
      E_p[5, 7, i, t] <- 0
      E_p[5, 8, i, t] <- 0
      E_p[5, 9, i, t] <- 0
      E_p[5, 10, i, t] <- 0
      E_p[5, 11, i, t] <- 0
      E_p[5, 12, i, t] <- 0
      E_p[5, 13, i, t] <- 0
      E_p[5, 14, i, t] <- 0
      E_p[5, 15, i, t] <- 0
      E_p[5, 16, i, t] <- 0
      E_p[5, 17, i, t] <- 0
      E_p[5, 18, i, t] <- 0
      E_p[5, 19, i, t] <- (1 - p_realized[i, t+1])
      
      E_p[6, 1, i, t] <- 0
      E_p[6, 2, i, t] <- 0
      E_p[6, 3, i, t] <- 0
      E_p[6, 4, i, t] <- 0
      E_p[6, 5, i, t] <- p_realized[i, t+1]
      E_p[6, 6, i, t] <- 0
      E_p[6, 7, i, t] <- 0
      E_p[6, 8, i, t] <- 0
      E_p[6, 9, i, t] <- 0
      E_p[6, 10, i, t] <- 0
      E_p[6, 11, i, t] <- 0
      E_p[6, 12, i, t] <- 0
      E_p[6, 13, i, t] <- 0
      E_p[6, 14, i, t] <- 0
      E_p[6, 15, i, t] <- 0
      E_p[6, 16, i, t] <- 0
      E_p[6, 17, i, t] <- 0
      E_p[6, 18, i, t] <- 0
      E_p[6, 19, i, t] <- (1 - p_realized[i, t+1])
      
      E_p[7, 1, i, t] <- 0
      E_p[7, 2, i, t] <- 0
      E_p[7, 3, i, t] <- 0
      E_p[7, 4, i, t] <- 0
      E_p[7, 5, i, t] <- 0
      E_p[7, 6, i, t] <- p_realized[i, t+1]
      E_p[7, 7, i, t] <- 0
      E_p[7, 8, i, t] <- 0
      E_p[7, 9, i, t] <- 0
      E_p[7, 10, i, t] <- 0
      E_p[7, 11, i, t] <- 0
      E_p[7, 12, i, t] <- 0
      E_p[7, 13, i, t] <- 0
      E_p[7, 14, i, t] <- 0
      E_p[7, 15, i, t] <- 0
      E_p[7, 16, i, t] <- 0
      E_p[7, 17, i, t] <- 0
      E_p[7, 18, i, t] <- 0
      E_p[7, 19, i, t] <- (1 - p_realized[i, t+1])
      
      E_p[8, 1, i, t] <- 0
      E_p[8, 2, i, t] <- 0
      E_p[8, 3, i, t] <- 0
      E_p[8, 4, i, t] <- 0
      E_p[8, 5, i, t] <- 0
      E_p[8, 6, i, t] <- 0
      E_p[8, 7, i, t] <- p_realized[i, t+1]
      E_p[8, 8, i, t] <- 0
      E_p[8, 9, i, t] <- 0
      E_p[8, 10, i, t] <- 0
      E_p[8, 11, i, t] <- 0
      E_p[8, 12, i, t] <- 0
      E_p[8, 13, i, t] <- 0
      E_p[8, 14, i, t] <- 0
      E_p[8, 15, i, t] <- 0
      E_p[8, 16, i, t] <- 0
      E_p[8, 17, i, t] <- 0
      E_p[8, 18, i, t] <- 0
      E_p[8, 19, i, t] <- (1 - p_realized[i, t+1])
      
      E_p[9, 1, i, t] <- 0
      E_p[9, 2, i, t] <- 0
      E_p[9, 3, i, t] <- 0
      E_p[9, 4, i, t] <- 0
      E_p[9, 5, i, t] <- 0
      E_p[9, 6, i, t] <- 0
      E_p[9, 7, i, t] <- 0
      E_p[9, 8, i, t] <- p_realized[i, t+1]
      E_p[9, 9, i, t] <- 0
      E_p[9, 10, i, t] <- 0
      E_p[9, 11, i, t] <- 0
      E_p[9, 12, i, t] <- 0
      E_p[9, 13, i, t] <- 0
      E_p[9, 14, i, t] <- 0
      E_p[9, 15, i, t] <- 0
      E_p[9, 16, i, t] <- 0
      E_p[9, 17, i, t] <- 0
      E_p[9, 18, i, t] <- 0
      E_p[9, 19, i, t] <- (1 - p_realized[i, t+1])
      
      E_p[10, 1, i, t] <- 0
      E_p[10, 2, i, t] <- 0
      E_p[10, 3, i, t] <- 0
      E_p[10, 4, i, t] <- 0
      E_p[10, 5, i, t] <- 0
      E_p[10, 6, i, t] <- 0
      E_p[10, 7, i, t] <- 0
      E_p[10, 8, i, t] <- 0
      E_p[10, 9, i, t] <- (1 - a) * p_realized[i, t+1]
      E_p[10, 10, i, t] <- 0
      E_p[10, 11, i, t] <- a * p_realized[i, t+1]
      E_p[10, 12, i, t] <- 0
      E_p[10, 13, i, t] <- 0
      E_p[10, 14, i, t] <- 0
      E_p[10, 15, i, t] <- 0
      E_p[10, 16, i, t] <- 0
      E_p[10, 17, i, t] <- 0
      E_p[10, 18, i, t] <- 0
      E_p[10, 19, i, t] <- (1 - p_realized[i, t+1])
      
      E_p[11, 1, i, t] <- 0
      E_p[11, 2, i, t] <- 0
      E_p[11, 3, i, t] <- 0
      E_p[11, 4, i, t] <- 0
      E_p[11, 5, i, t] <- 0
      E_p[11, 6, i, t] <- 0
      E_p[11, 7, i, t] <- 0
      E_p[11, 8, i, t] <- 0
      E_p[11, 9, i, t] <- 2 * (1 - a) * a * p_realized[i, t+1]
      E_p[11, 10, i, t] <- (1 - a)^2 * p_realized[i, t+1]
      E_p[11, 11, i, t] <- a^2 * p_realized[i, t+1]
      E_p[11, 12, i, t] <- 0
      E_p[11, 13, i, t] <- 0
      E_p[11, 14, i, t] <- 0
      E_p[11, 15, i, t] <- 0
      E_p[11, 16, i, t] <- 0
      E_p[11, 17, i, t] <- 0
      E_p[11, 18, i, t] <- 0
      E_p[11, 19, i, t] <- (1 - p_realized[i, t+1])
      
      E_p[12, 1, i, t] <- 0
      E_p[12, 2, i, t] <- 0
      E_p[12, 3, i, t] <- 0
      E_p[12, 4, i, t] <- 0
      E_p[12, 5, i, t] <- 0
      E_p[12, 6, i, t] <- 0
      E_p[12, 7, i, t] <- 0
      E_p[12, 8, i, t] <- 0
      E_p[12, 9, i, t] <- 0
      E_p[12, 10, i, t] <- 0
      E_p[12, 11, i, t] <- p_realized[i, t+1]
      E_p[12, 12, i, t] <- 0
      E_p[12, 13, i, t] <- 0
      E_p[12, 14, i, t] <- 0
      E_p[12, 15, i, t] <- 0
      E_p[12, 16, i, t] <- 0
      E_p[12, 17, i, t] <- 0
      E_p[12, 18, i, t] <- 0
      E_p[12, 19, i, t] <- (1 - p_realized[i, t+1])
      
      E_p[13, 1, i, t] <- 0
      E_p[13, 2, i, t] <- 0
      E_p[13, 3, i, t] <- 0
      E_p[13, 4, i, t] <- 0
      E_p[13, 5, i, t] <- 0
      E_p[13, 6, i, t] <- 0
      E_p[13, 7, i, t] <- 0
      E_p[13, 8, i, t] <- 0
      E_p[13, 9, i, t] <- 0
      E_p[13, 10, i, t] <- 0
      E_p[13, 11, i, t] <- p_realized[i, t+1]
      E_p[13, 12, i, t] <- 0
      E_p[13, 13, i, t] <- 0
      E_p[13, 14, i, t] <- 0
      E_p[13, 15, i, t] <- 0
      E_p[13, 16, i, t] <- 0
      E_p[13, 17, i, t] <- 0
      E_p[13, 18, i, t] <- 0
      E_p[13, 19, i, t] <- (1 - p_realized[i, t+1])
      
      E_p[14, 1, i, t] <- 0
      E_p[14, 2, i, t] <- 0
      E_p[14, 3, i, t] <- 0
      E_p[14, 4, i, t] <- 0
      E_p[14, 5, i, t] <- 0
      E_p[14, 6, i, t] <- 0
      E_p[14, 7, i, t] <- 0
      E_p[14, 8, i, t] <- 0
      E_p[14, 9, i, t] <- 0
      E_p[14, 10, i, t] <- 0
      E_p[14, 11, i, t] <- 0
      E_p[14, 12, i, t] <- 0
      E_p[14, 13, i, t] <- 0
      E_p[14, 14, i, t] <- 0
      E_p[14, 15, i, t] <- 0
      E_p[14, 16, i, t] <- 0
      E_p[14, 17, i, t] <- 0
      E_p[14, 18, i, t] <- 0
      E_p[14, 19, i, t] <- 1
      
      Z[i, t+1] <- which(rmultinom(n = 1, size = 1,
                                   prob = S[Z[i, t], , i, t]) == 1)
  
      CH[i, t+1] <-  which(rmultinom(n = 1, size = 1,
                                     E_p[Z[i, t+1], , i, t]) == 1)
      
      # ~~~~~ GPS/GLS ----------------------------------------------------------
      
      # Tag deployment
      if (CH[i, t+1] != not_captured) {  # If (physically) captured
        
        # GPS collaring
        collared <- collaring_save[i, t+1] <- which(rmultinom(n = 1, size = 1, 
                                                              prob = COLLARING[Z[i, t+1], ]) == 1) - 1
        if (collared == 1) { # If individual receives a collar, for how many years ?
          collared_years <- sample(x = 1:4, size = 1, prob = collar_duration)
          if (collared_years >= 1) {
            if (t+2 <= n_occasions) {
              PI[i, t+2] <- 1
            } 
          }
          if (collared_years >= 2) {
            if (t+3 <= n_occasions) {
              PI[i, t+3] <- 1
            }
          }
          if (collared_years >= 3) {
            if (t+4 <= n_occasions) {
              PI[i, t+4] <- 1
            }
          }
        }
        
        # GLS tagging
        # Will the individual be tagged, and how many years will the tag record ?
        tag <- rbinom(n = 1, size = 1, prob = prob_GLS_tagging) * 
          sample(x = c(3, 4, 5), size = 1, prob = GLS_duration)
        
        if (tag >= 3) {
          if (t+2 <= n_occasions) {
            RHO[i, t+2] <- 1
          } 
          if (t+3 <= n_occasions) {
            RHO[i, t+3] <- 1
          } 
          if (t+4 <= n_occasions) {
            RHO[i, t+4] <- 1
          } 
        }
        if (tag >= 4) {
          if (t+5 <= n_occasions) {
            RHO[i, t+5] <- 1
          }
        }
        if (tag == 5) {
          if (t+6 <= n_occasions) {
            RHO[i, t+6] <- 1
          }
        }
      }
    }
  }
  
  CH_p <- CH

  # GLS tag data is available only if the bear is recaptured. Let's take that into account
  for (i in 1:(new_per_year*n_occasions - new_per_year)) {
    if (f[i] == n_occasions) next

    for (t in (f[i]+1):n_occasions) {

      # If the bear wore a GLS tag and was recaptured (necessarily physically since
      # remote captures not taken into account yet) at some point
      if (RHO[i, t] == 1 && sum(CH[i, t:n_occasions] != not_captured) > 0) {
        TAU[i, t] <- 1
      }
    }
  }
  KHI <- PI + TAU
  KHI[KHI == 2] <- 1


  # Update the observation matrix according to the availability of GPS/GLS data
  for (i in 1:(new_per_year*n_occasions - new_per_year)) {
    if (f[i] == n_occasions) next

    for (t in f[i]:(n_occasions-1)) {

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
      E[1, 12, i, t] <- 0
      E[1, 13, i, t] <- 0
      E[1, 14, i, t] <- 0
      E[1, 15, i, t] <- 0
      E[1, 16, i, t] <- 0
      E[1, 17, i, t] <- 0
      E[1, 18, i, t] <- 0
      E[1, 19, i, t] <- 1-p_realized[i, t+1]

      E[2, 1, i, t] <- 0
      E[2, 2, i, t] <- KHI[i, t+1] + (1 - KHI[i, t+1]) * p_realized[i, t+1]
      E[2, 3, i, t] <- 0
      E[2, 4, i, t] <- 0
      E[2, 5, i, t] <- 0
      E[2, 6, i, t] <- 0
      E[2, 7, i, t] <- 0
      E[2, 8, i, t] <- 0
      E[2, 9, i, t] <- 0
      E[2, 10, i, t] <- 0
      E[2, 11, i, t] <- 0
      E[2, 12, i, t] <- 0
      E[2, 13, i, t] <- 0
      E[2, 14, i, t] <- 0
      E[2, 15, i, t] <- 0
      E[2, 16, i, t] <- 0
      E[2, 17, i, t] <- 0
      E[2, 18, i, t] <- 0
      E[2, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p_realized[i, t+1])

      E[3, 1, i, t] <- 0
      E[3, 2, i, t] <- 0
      E[3, 3, i, t] <- KHI[i, t+1] + (1 - KHI[i, t+1]) * p_realized[i, t+1]
      E[3, 4, i, t] <- 0
      E[3, 5, i, t] <- 0
      E[3, 6, i, t] <- 0
      E[3, 7, i, t] <- 0
      E[3, 8, i, t] <- 0
      E[3, 9, i, t] <- 0
      E[3, 10, i, t] <- 0
      E[3, 11, i, t] <- 0
      E[3, 12, i, t] <- 0
      E[3, 13, i, t] <- 0
      E[3, 14, i, t] <- 0
      E[3, 15, i, t] <- 0
      E[3, 16, i, t] <- 0
      E[3, 17, i, t] <- 0
      E[3, 18, i, t] <- 0
      E[3, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p_realized[i, t+1])

      E[4, 1, i, t] <- 0
      E[4, 2, i, t] <- 0
      E[4, 3, i, t] <- 0
      E[4, 4, i, t] <- (1 - KHI[i, t+1]) * p_realized[i, t+1]
      E[4, 5, i, t] <- 0
      E[4, 6, i, t] <- 0
      E[4, 7, i, t] <- 0
      E[4, 8, i, t] <- 0
      E[4, 9, i, t] <- 0
      E[4, 10, i, t] <- 0
      E[4, 11, i, t] <- 0
      E[4, 12, i, t] <- KHI[i, t+1] * (1-p_realized[i, t+1]) + KHI[i, t+1] * p_realized[i, t+1]
      E[4, 13, i, t] <- 0
      E[4, 14, i, t] <- 0
      E[4, 15, i, t] <- 0
      E[4, 16, i, t] <- 0
      E[4, 17, i, t] <- 0
      E[4, 18, i, t] <- 0
      E[4, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p_realized[i, t+1])

      E[5, 1, i, t] <- 0
      E[5, 2, i, t] <- 0
      E[5, 3, i, t] <- 0
      E[5, 4, i, t] <- (1 - KHI[i, t+1]) * p_realized[i, t+1]
      E[5, 5, i, t] <- 0
      E[5, 6, i, t] <- 0
      E[5, 7, i, t] <- 0
      E[5, 8, i, t] <- 0
      E[5, 9, i, t] <- 0
      E[5, 10, i, t] <- 0
      E[5, 11, i, t] <- 0
      E[5, 12, i, t] <- 0
      E[5, 13, i, t] <- KHI[i, t+1] * (1 - p_realized[i, t+1])
      E[5, 14, i, t] <- 0
      E[5, 15, i, t] <- 0
      E[5, 16, i, t] <- KHI[i, t+1] * p_realized[i, t+1]
      E[5, 17, i, t] <- 0
      E[5, 18, i, t] <- 0
      E[5, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p_realized[i, t+1])

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
      E[6, 12, i, t] <- 0
      E[6, 13, i, t] <- KHI[i, t+1] * (1 - p_realized[i, t+1]) * AGE_5[i, t+1]
      E[6, 14, i, t] <- 0
      E[6, 15, i, t] <- KHI[i, t+1] * (1 - p_realized[i, t+1]) * (1 - AGE_5[i, t+1])
      E[6, 16, i, t] <- 0
      E[6, 17, i, t] <- 0
      E[6, 18, i, t] <- 0
      E[6, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p_realized[i, t+1])

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
      E[7, 12, i, t] <- 0
      E[7, 13, i, t] <- KHI[i, t+1] * (1 - p_realized[i, t+1]) * AGE_5[i, t+1]
      E[7, 14, i, t] <- 0
      E[7, 15, i, t] <- KHI[i, t+1] * (1 - p_realized[i, t+1]) * (1 - AGE_5[i, t+1])
      E[7, 16, i, t] <- 0
      E[7, 17, i, t] <- 0
      E[7, 18, i, t] <- 0
      E[7, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p_realized[i, t+1])

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
      E[8, 12, i, t] <- 0
      E[8, 13, i, t] <- 0
      E[8, 14, i, t] <- KHI[i, t+1] * (1 - p_realized[i, t+1])
      E[8, 15, i, t] <- 0
      E[8, 16, i, t] <- 0
      E[8, 17, i, t] <- 0
      E[8, 18, i, t] <- 0
      E[8, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p_realized[i, t+1])

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
      E[9, 12, i, t] <- 0
      E[9, 13, i, t] <- 0
      E[9, 14, i, t] <- KHI[i, t+1] * (1 - p_realized[i, t+1])
      E[9, 15, i, t] <- 0
      E[9, 16, i, t] <- 0
      E[9, 17, i, t] <- 0
      E[9, 18, i, t] <- 0
      E[9, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p_realized[i, t+1])

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
      E[10, 11, i, t] <- a * (1 - KHI[i, t+1]) * p_realized[i, t+1]
      E[10, 12, i, t] <- 0
      E[10, 13, i, t] <- 0
      E[10, 14, i, t] <- KHI[i, t+1] * (1 - p_realized[i, t+1])
      E[10, 15, i, t] <- 0
      E[10, 16, i, t] <- 0
      E[10, 17, i, t] <- a * KHI[i, t+1] * p_realized[i, t+1]
      E[10, 18, i, t] <- 0
      E[10, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p_realized[i, t+1])

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
      E[11, 11, i, t] <- a^2 * (1 - KHI[i, t+1]) * p_realized[i, t+1]
      E[11, 12, i, t] <- 0
      E[11, 13, i, t] <- 0
      E[11, 14, i, t] <- KHI[i, t+1] * (1 - p_realized[i, t+1])
      E[11, 15, i, t] <- 0
      E[11, 16, i, t] <- 0
      E[11, 17, i, t] <- a^2 * KHI[i, t+1] * p_realized[i, t+1]
      E[11, 18, i, t] <- 0
      E[11, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p_realized[i, t+1])

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
      E[12, 11, i, t] <- (1 - KHI[i, t+1]) * p_realized[i, t+1]
      E[12, 12, i, t] <- 0
      E[12, 13, i, t] <- 0
      E[12, 14, i, t] <- KHI[i, t+1] * (1 - p_realized[i, t+1])
      E[12, 15, i, t] <- 0
      E[12, 16, i, t] <- 0
      E[12, 17, i, t] <- KHI[i, t+1] * p_realized[i, t+1]
      E[12, 18, i, t] <- 0
      E[12, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p_realized[i, t+1])

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
      E[13, 11, i, t] <- (1 - KHI[i, t+1]) * p_realized[i, t+1]
      E[13, 12, i, t] <- 0
      E[13, 13, i, t] <- 0
      E[13, 14, i, t] <- 0
      E[13, 15, i, t] <- KHI[i, t+1] * (1 - p_realized[i, t+1])
      E[13, 16, i, t] <- 0
      E[13, 17, i, t] <- 0
      E[13, 18, i, t] <- KHI[i, t+1] * p_realized[i, t+1]
      E[13, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p_realized[i, t+1])

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
      E[14, 12, i, t] <- 0
      E[14, 13, i, t] <- 0
      E[14, 14, i, t] <- 0
      E[14, 15, i, t] <- 0
      E[14, 16, i, t] <- 0
      E[14, 17, i, t] <- 0
      E[14, 18, i, t] <- 0
      E[14, 19, i, t] <- 1


    }
  }

  for (i in 1:(new_per_year*n_occasions - new_per_year)) {

    for (t in ( f[i]+1 ):n_occasions) {

      #----- Observation process -----#
      # draw y(t) given z(t) using the observation matrix including remote observation
      CH[i, t] <- which(rmultinom(n = 1, size = 1,
                                  E[Z[i, t], , i, t-1]) == 1)    }
  }
  
  
  # Compute zeta to give the model as initial values
  zeta <- array(data = NA, dim = c((new_per_year*n_occasions - new_per_year), 
                                   n_occasions, n_states))
  lik <- NULL
  for(i in 1:(new_per_year*n_occasions - new_per_year)){
    
    # First capture
    for (j in 1:n_states) {
      zeta[i, f[i], j] <- f_state_prop[j] * E0[j, CH[i, f[i]]]
    }
    
    # Subsequent captures
    for(t in f[i] : (n_occasions-1) ){
      for (j in 1:n_states) {
        zeta[i, (t+1), j] <- inprod(
          zeta[i, t, 1:n_states],  S[1:n_states, j, i, t])  *
          E[j, CH[i, t+1], i, t]
      }
    }
    lik[i] <- sum(zeta[i, n_occasions, 1:n_states])
  }


  
  # ~ 3. Bundle data -------------------------------------------------------------

  # Get simulated proportion of each state at first capture
  prop_raw <- table(f_state)/length(f_state)
  prop <- NULL
  for (k in 1:n_states) {
    x <- prop_raw[as.character(k)]
    prop[k] <- ifelse(!is.na(x), x, 0)
  }

  dat <- list(ones = rep(1, times = nrow(CH)))

  my.constants <- list(CH = CH,
                       freq = rep(1, times = nrow(CH)),
                       f = f,
                       n_ind = nrow(CH),
                       n_occasions = ncol(CH),
                       S0 = f_state_prop,
                       KHI = KHI,
                       a = a,
                       AGE_4 = AGE_4,
                       AGE_5 = AGE_5,
                       AGE_YOUNG = AGE_YOUNG,
                       AGE_MIDDLE = AGE_MIDDLE,
                       AGE_OLD = AGE_OLD,
                       AGE_YOUNG_1 = AGE_YOUNG_1,
                       AGE_MIDDLE_1 = AGE_MIDDLE_1,
                       AGE_OLD_1 = AGE_OLD_1,
                       AGE_YOUNG_2 = AGE_YOUNG_2,
                       AGE_MIDDLE_2 = AGE_MIDDLE_2,
                       AGE_OLD_2 = AGE_OLD_2,
                       AGE_2_4 = AGE_2_4,
                       sea_ice_s = sea_ice_s)
  
  # ~ 4. Check data --------------------------------------------------------------
  
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
              nrow = n_states, ncol = n_states, byrow = TRUE)
  check_Z <- CH ; check_Z[] <- NA
  for (i in 1:nrow(CH)) {
    if (f[i] == n_occasions) next
    for (t in f[i]:(n_occasions - 1)) {
      check_Z[i, t] <- ifelse(S[Z[i, t], Z[i, t+1]] == 0, 1, 0)
    }
  }
  print(sum(check_Z, na.rm = T) == 0)
  
  # Let's now check that the states are compatible with the observations
  E <- matrix(data = c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,    # J2
                       0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,    # J3
                       0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,    # SA4
                       0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,    # NDSA5
                       0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1,    # FDSA5  5
                       0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1,    # A01
                       0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1,    # A02
                       0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1,    # A11    
                       0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1,    # A12
                       0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1,    # AS1    10
                       0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1,    # AS2
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1,    # LNDA
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1,    # FDA    
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),   # D
              nrow = n_states, ncol = not_captured, byrow = TRUE)
  
  check_Z <- CH ; check_Z[] <- NA
  for (i in 1:nrow(CH)) {
    for (t in f[i]:n_occasions) {
      check_Z[i, t] <- ifelse(E[Z[i, t], CH[i, t]] == 0, 1, 0)
      if (check_Z[i, t] == 1) {
        print(paste0("i = ", i, " ; t = ", t))
        print(paste0("state: ", Z[i, t], ": event: ",  CH[i, t]))
      }
    }
  }
  print(sum(check_Z, na.rm = T) == 0)
  
  N[l] <- sum(CH %in% 1:18)                # Number of captures
  N_p[l] <- sum(CH %in% c(1:11, 16:18))    # Number of physical captures
  N_3[l] <- sum(CH %in% 3:18)              # Number of captures of inividuals old 
  # enough to be equipped
  N_remote[l] <- sum(CH %in% 12:18)        # Number of remote captures that bring
  # additional information
  PI[which(Z == n_states)] <- 0
  TAU[which(Z == n_states)] <- 0
  N_GPS[l] <- sum(PI == 1)                 # Number of GPS remote captures
  N_GLS[l] <- sum(TAU == 1)                # Number of GLS remote captures
  
  name <- paste0("03_outputs/simulations/simulated_datasets/polar_bear_model_validation_all_covariates/validation_dataset_", 
                 l, ".RData")
  save(CR_model, my.constants, dat, Z,
       file = name)
  
}

mean(N)         # ~970 captures (remote or physical)
mean(N_3)       # ~915 captures (remote or physical) of females old enough to be equipped
mean(N_p)       # ~575 physical captures
mean(N_GPS)     # ~250 female-years with GPS data available
mean(N_GLS)     # ~510 female-years with GLS data available
mean(N_remote)  # ~585 captures where remote data added information

mean(N_GPS/N_3)  # ~27% of captures of females old enough to be equipped for which GPS data was available
mean(N_GLS/N_3)  # ~55% of captures of females old enough to be equipped for which GLS data was available

mean(N_remote/N_3)  # ~55% of captures of individuals old enough to be equipped for which remote 
                    # data added information


