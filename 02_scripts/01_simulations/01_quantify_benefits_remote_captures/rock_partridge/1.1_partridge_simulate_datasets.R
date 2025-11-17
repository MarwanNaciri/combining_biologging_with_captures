#==============================================================================#
#                                                                              #
#                   Simulate datasets partridge example                        #
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
  # phi: survival probability
  # beta: breeding probability
  # eta: early offspring survival
  # gamma: probability of a large litter/brood

  # -------------------------------------------------+
  # States (z):
  # 1 = NB: Non-breeder
  # 2 = FB: Failed-breeder 
  # 3 = SB1: Successful breeder: litter/brood size = 1
  # 4 = SB2: Successful breeder: litter/brood size = 2
  # 5 = D: Dead
  
  # -------------------------------------------------+
  # Events (y):
  # 1 = L: Capture of a lone individual 
  # 2 = SB1: Capture of an individual with litter/brood size = 1
  # 3 = SB2: Capture of an individual with litter/brood size = 2
  # 4 = NB: Remote capture of a non-breeder
  # 5 = B: Remote capture of a breeder
  # 6 = FB: Physical + remote capture of lone individual who bred 
  # 7 = NA: Not captured/observed
  
  # +++++++++++++++++++++++++++++++ matrices +++++++++++++++++++++++++++++++++++

  for (t in 1:(n_occasions - 1)) {
    
    # S ++++++++++++++++
    logit(phi[1, t]) <- a1_phi + a3_phi * time[t]
    logit(phi[2, t]) <- a2_phi + a3_phi * time[t]
    
    S[1, 1, t] <- 0
    S[1, 2, t] <- phi[1, t] * (1 - beta[1])
    S[1, 3, t] <- phi[1, t] * beta[1] * (1 - eta[1])
    S[1, 4, t] <- phi[1, t] * beta[1] * eta[1] * (1 - gamma)
    S[1, 5, t] <- phi[1, t] * beta[1] * eta[1] * gamma
    S[1, 6, t] <-  1-phi[1, t]
    
    S[2, 1, t] <- 0
    S[2, 2, t] <- phi[2, t] * (1 - beta[2])
    S[2, 3, t] <- phi[2, t] * beta[2] * (1 - eta[2])
    S[2, 4, t] <- phi[2, t] * beta[2] * eta[2] * (1 - gamma)
    S[2, 5, t] <- phi[2, t] * beta[2] * eta[2] * gamma
    S[2, 6, t] <-  1-phi[2, t]
    
    S[3, 1, t] <- 0
    S[3, 2, t] <- phi[2, t] * (1 - beta[2])
    S[3, 3, t] <- phi[2, t] * beta[2] * (1 - eta[2])
    S[3, 4, t] <- phi[2, t] * beta[2] * eta[2] * (1 - gamma)
    S[3, 5, t] <- phi[2, t] * beta[2] * eta[2] * gamma
    S[3, 6, t] <-  1-phi[2, t]
    
    S[4, 1, t] <- 0
    S[4, 2, t] <- phi[2, t] * (1 - beta[3])
    S[4, 3, t] <- phi[2, t] * beta[3] * (1 - eta[3])
    S[4, 4, t] <- phi[2, t] * beta[3] * eta[3] * (1 - gamma)
    S[4, 5, t] <- phi[2, t] * beta[3] * eta[3] * gamma
    S[4, 6, t] <-  1-phi[2, t]
    
    S[5, 1, t] <- 0
    S[5, 2, t] <- phi[2, t] * (1 - beta[3])
    S[5, 3, t] <- phi[2, t] * beta[3] * (1 - eta[3])
    S[5, 4, t] <- phi[2, t] * beta[3] * eta[3] * (1 - gamma)
    S[5, 5, t] <- phi[2, t] * beta[3] * eta[3] * gamma
    S[5, 6, t] <-  1-phi[2, t]
    
    S[6, 1, t] <- 0
    S[6, 2, t] <- 0
    S[6, 3, t] <- 0
    S[6, 4, t] <- 0
    S[6, 5, t] <- 0
    S[6, 6, t] <- 1
    
    for (i in 1:n_ind) {
      
      # E ++++++++++++++++
      
      E[1, 1, i, t] <- p 
      E[1, 2, i, t] <- 0
      E[1, 3, i, t] <- 0
      E[1, 4, i, t] <- 0
      E[1, 5, i, t] <- 0
      E[1, 6, i, t] <- 0
      E[1, 7, i, t] <- 0
      E[1, 8, i, t] <- 1 - p
      
      E[2, 1, i, t] <- 0
      E[2, 2, i, t] <- p * (1- KHI[i, t+1])  
      E[2, 3, i, t] <- 0
      E[2, 4, i, t] <- 0
      E[2, 5, i, t] <- KHI[i, t+1]
      E[2, 6, i, t] <- 0
      E[2, 7, i, t] <- 0
      E[2, 8, i, t] <- (1 - p) * (1 - KHI[i, t+1])
      
      E[3, 1, i, t] <- 0
      E[3, 2, i, t] <- p * (1- KHI[i, t+1]) 
      E[3, 3, i, t] <- 0
      E[3, 4, i, t] <- 0
      E[3, 5, i, t] <- 0
      E[3, 6, i, t] <- (1 - p) * KHI[i, t+1]
      E[3, 7, i, t] <- p * KHI[i, t+1]
      E[3, 8, i, t] <- (1 - p) * (1 - KHI[i, t+1])
      
      E[4, 1, i, t] <- 0
      E[4, 2, i, t] <- 0
      E[4, 3, i, t] <- p
      E[4, 4, i, t] <- 0
      E[4, 5, i, t] <- 0
      E[4, 6, i, t] <- (1 - p) * KHI[i, t+1]
      E[4, 7, i, t] <- 0
      E[4, 8, i, t] <- (1 - p) * (1 - KHI[i, t+1])
      
      E[5, 1, i, t] <- 0
      E[5, 2, i, t] <- 0
      E[5, 3, i, t] <- 0
      E[5, 4, i, t] <- p
      E[5, 5, i, t] <- 0
      E[5, 6, i, t] <- (1 - p) * KHI[i, t+1]
      E[5, 7, i, t] <- 0
      E[5, 8, i, t] <- (1 - p) * (1 - KHI[i, t+1])
      
      E[6, 1, i, t] <- 0 
      E[6, 2, i, t] <- 0
      E[6, 3, i, t] <- 0
      E[6, 4, i, t] <- 0
      E[6, 5, i, t] <- 0
      E[6, 6, i, t] <- 0
      E[6, 7, i, t] <- 0
      E[6, 8, i, t] <- 1
    }
  }

  # E0 ++++++++++++++++
  
  E0[1, 1] <- 1 
  E0[1, 2] <- 0
  E0[1, 3] <- 0
  E0[1, 4] <- 0
  E0[1, 5] <- 0
  E0[1, 6] <- 0
  E0[1, 7] <- 0
  E0[1, 8] <- 0
  
  E0[2, 1] <- 0
  E0[2, 2] <- 1
  E0[2, 3] <- 0
  E0[2, 4] <- 0
  E0[2, 5] <- 0
  E0[2, 6] <- 0
  E0[2, 7] <- 0
  E0[2, 8] <- 0
  
  E0[3, 1] <- 0
  E0[3, 2] <- 1
  E0[3, 3] <- 0
  E0[3, 4] <- 0
  E0[3, 5] <- 0
  E0[3, 6] <- 0
  E0[3, 7] <- 0
  E0[3, 8] <- 0
  
  E0[4, 1] <- 0
  E0[4, 2] <- 0
  E0[4, 3] <- 1
  E0[4, 4] <- 0
  E0[4, 5] <- 0
  E0[4, 6] <- 0
  E0[4, 7] <- 0
  E0[4, 8] <- 0
  
  E0[5, 1] <- 0 
  E0[5, 2] <- 0
  E0[5, 3] <- 0
  E0[5, 4] <- 4
  E0[5, 5] <- 0
  E0[5, 6] <- 0
  E0[5, 7] <- 0
  E0[5, 8] <- 0
  
  E0[6, 1] <- 0 
  E0[6, 2] <- 0
  E0[6, 3] <- 0
  E0[6, 4] <- 0
  E0[6, 5] <- 0
  E0[6, 6] <- 0
  E0[6, 7] <- 0
  E0[6, 8] <- 1
  
  # ++++++++++++++++++++++++++++++++ priors ++++++++++++++++++++++++++++++++++++

  #-------- State process --------#
  a1_phi ~ dnorm(0, 1.5)   # survival (J)
  a2_phi ~ dnorm(0, 1.5)   # survival (Adults)
  a3_phi ~ dnorm(0, 3)   # survival: time effect
  beta[1] ~ dunif(0, 1)  # breeding probability (J)
  beta[2] ~ dunif(0, 1)  # breeding probability (NB/FB)
  beta[3] ~ dunif(0, 1)  # breeding probability (SB)
  eta[1] ~ dunif(0, 1)  # early offspring survival probability (J)
  eta[2] ~ dunif(0, 1)  # early offspring survival probability (NB/FB)
  eta[3] ~ dunif(0, 1)  # early offspring survival probability (SB)
  gamma ~ dunif(0, 1)  # probability of a large litter/brood
  
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
                                    S[1:n_states, k, t]) *
          E[k, CH[i, t+1], i, t]
      }
    }
    lik[i] <- sum(zeta[i, n_occasions, 1:n_states])
    ones[i] ~ dbin(prob = lik[i], size = freq[i])
  }
})


# # Matrices
# #                     J          NB               FB                  SB1                    SB2          D
# S <- matrix(data = c( 0,    phi*(1-beta), phi*beta*(1-eta), phi*beta*eta*(1-gamma), phi*beta*eta*gamma, 1-phi,   # J
#                       0,    phi*(1-beta), phi*beta*(1-eta), phi*beta*eta*(1-gamma), phi*beta*eta*gamma, 1-phi,   # NB
#                       0,    phi*(1-beta), phi*beta*(1-eta), phi*beta*eta*(1-gamma), phi*beta*eta*gamma, 1-phi,   # FB
#                       0,    phi*(1-beta), phi*beta*(1-eta), phi*beta*eta*(1-gamma), phi*beta*eta*gamma, 1-phi,   # SB1
#                       0,    phi*(1-beta), phi*beta*(1-eta), phi*beta*eta*(1-gamma), phi*beta*eta*gamma, 1-phi,   # SB2
#                       0,         0,              0,                    0,                     0,          1),    # D
#             nrow = 6, ncol = 6, byrow = TRUE)
# 
# p <- pi <- 0.5
# #                    J       L        SB1  SB2  NB      B       FB           D
# E <- matrix(data = c(p,     0,        0,   0,   0,      0,      0,          1-p,       # J
#                      0,  p*(1-pi),    0,   0,   pi,     0,      0,      (1-p)*(1-pi),  # NB
#                      0,  p*(1-pi),    0,   0,   0,  (1-p)*pi,  p*pi,    (1-p)*(1-pi),  # FB
#                      0,     0,        p,   0,   0,  (1-p)*pi,   0,      (1-p)*(1-pi),  # SB1
#                      0,     0,        0,   p,   0,  (1-p)*pi,   0,      (1-p)*(1-pi),  # SB2
#                      0,     0,        0,   0,   0,      0,      0,           1),       # D
#             nrow = 6, ncol = 8, byrow = TRUE)
# apply(E, 1, sum)
# 
# #                     J    L   SB1  SB2  NB    B    FB    D
# E0 <- matrix(data = c(1,   0,   0,   0,   0,   0,    0,   0,  # J
#                       0,   1,   0,   0,   0,   0,    0,   0,  # NB
#                       0,   1,   0,   0,   0,   0,    0,   0,  # FB
#                       0,   0,   1,   0,   0,   0,    0,   0,  # SB1
#                       0,   0,   0,   1,   0,   0,    0,   0,  # SB2
#                       0,   0,   0,   0,   0,   0,    0,   1),  # D
#             nrow = 6, ncol = 8, byrow = TRUE)
# apply(E, 1, sum)



# B. Generate datasets =========================================================

# ~ 1. Define parameters -------------------------------------------------------

{
  new_per_year <- 15
  n_occasions <- 8
  
  n_states <- 6
  n_events <- 8
  physical <- c(1, 2, 3, 4, 7)
  
  f <- rep(c(1:n_occasions), each = new_per_year)
  
  beta <- eta <- NULL
  a1_phi <- logit(0.35)  # Survival  
  a2_phi <- logit(0.45)
  a3_phi <- -0.4  # time effect   
  time <- as.vector(scale(1:8))
  beta[1] <- 0.4  # Breeding probability (J) 
  beta[2] <- 0.7  # Breeding probability (NB/FB) 
  beta[3] <- 0.9  # Breeding probability (SB) 
  eta[1] <- 0.4  # Early offspring survival probability (J) 
  eta[2] <- 0.4  # Early offspring survival probability (NB/FB) 
  eta[3] <- 0.6  # Early offspring survival probability (SB) 
  gamma <- 0.4 # Probability of large litter/brood
  
  # Recapture probability 
  values_p <- c(0.05, 0.10, 0.25, 0.5)
  
  # Probability of being equipped with a GPS at capture
  values_pGPS <- c(0.10, 0.25, 0.5, 0.75, 1)
  
  # Duration of collars,
  lambda_GPS <- 0.70
  
  # Error rate in assignment of breeding status based on remote capture (only false
  # negative i.e. missed breeding events)
  values_pError <- 0 # c(0, 0.05, 0.10, 0.20)
  
  # Proportions of each state at first capture
  f_state_prop <- c(0.4,  # j
                    0.15,  # NB
                    0.1,  # FB
                    0.15,  # SB1
                    0.1,  # SB2
                    0)    # D
  
}



# ~ 2. Generate capture histories ----------------------------------------------

for (k in 1:length(values_p)) {
  for (l in 1:length(values_pGPS)) {
    for (m in 1:length(values_pError)) {
      for (n in n_sim_start:n_sim) {
        print(paste0("simulation ", n, ":"))
        
        E <- array(data = NA, 
                   dim = c(n_states, n_events, (new_per_year*n_occasions - new_per_year), n_occasions))
        
        #                     J    L   SB1  SB2  NB    B    FB    D
        E0 <- matrix(data = c(1,   0,   0,   0,   0,   0,    0,   0,  # J
                              0,   1,   0,   0,   0,   0,    0,   0,  # NB
                              0,   1,   0,   0,   0,   0,    0,   0,  # FB
                              0,   0,   1,   0,   0,   0,    0,   0,  # SB1
                              0,   0,   0,   1,   0,   0,    0,   0,  # SB2
                              0,   0,   0,   0,   0,   0,    0,   1),  # D 
                     nrow = n_states, ncol = n_events, byrow = TRUE)
        
        S <- array(data = NA, dim = c(n_states, n_states, n_occasions))

        # Matrices for data
        CH <- Z <- matrix(NA, nrow = (new_per_year*n_occasions - new_per_year), 
                          ncol = n_occasions) 
        phi <- matrix(NA, nrow = 2, 
                      ncol = n_occasions) 
        KHI <- matrix(0, nrow = (new_per_year*n_occasions - new_per_year), 
                      ncol = n_occasions) 
        
        f_event <- f_state <-  NULL
        
        # ~~~ a. First captures --------------------------------------------------------
        
        set.seed(n)
        for (i in 1:(new_per_year*n_occasions - new_per_year)) {
          Z[i, f[i]] <- f_state[i] <- sample(x = 1:n_states,
                                             size = 1,
                                             prob = f_state_prop, replace = T)
          
          CH[i, f[i]] <- f_event[i] <- which(rmultinom(n = 1, size = 1, 
                                                       prob = E0[Z[i, f[i]], ]) == 1)
          
          # GPS collaring
          collared <- rbinom(n = 1, size = 1, prob = values_pGPS[l])
          
          if (collared == 1) { # If individual has a collar, for how many years ?
            collared_years <- rpois(n = 1, lambda = lambda_GPS) + 1
            q <- 1
            while (q <= collared_years &  
                   (f[i] + q) <= n_occasions) {
              KHI[i, f[i] + q] <- 1
              q <- q + 1
            }
          }
        }  
        
        # ~~~ b. Subsequent captures ---------------------------------------------------
        
        set.seed(n)
        for (t in 1:(n_occasions-1)) {
          
          phi[1, t] <- plogis(a1_phi + a3_phi * time[t])
          phi[2, t] <- plogis(a2_phi + a3_phi * time[t])
          
          S[1, 1, t] <- 0
          S[1, 2, t] <- phi[1, t] * (1 - beta[1])
          S[1, 3, t] <- phi[1, t] * beta[1] * (1 - eta[1])
          S[1, 4, t] <- phi[1, t] * beta[1] * eta[1] * (1 - gamma)
          S[1, 5, t] <- phi[1, t] * beta[1] * eta[1] * gamma
          S[1, 6, t] <-  1-phi[1, t]
          
          S[2, 1, t] <- 0
          S[2, 2, t] <- phi[2, t] * (1 - beta[2])
          S[2, 3, t] <- phi[2, t] * beta[2] * (1 - eta[2])
          S[2, 4, t] <- phi[2, t] * beta[2] * eta[2] * (1 - gamma)
          S[2, 5, t] <- phi[2, t] * beta[2] * eta[2] * gamma
          S[2, 6, t] <-  1-phi[2, t]
          
          S[3, 1, t] <- 0
          S[3, 2, t] <- phi[2, t] * (1 - beta[2])
          S[3, 3, t] <- phi[2, t] * beta[2] * (1 - eta[2])
          S[3, 4, t] <- phi[2, t] * beta[2] * eta[2] * (1 - gamma)
          S[3, 5, t] <- phi[2, t] * beta[2] * eta[2] * gamma
          S[3, 6, t] <-  1-phi[2, t]
          
          S[4, 1, t] <- 0
          S[4, 2, t] <- phi[2, t] * (1 - beta[3])
          S[4, 3, t] <- phi[2, t] * beta[3] * (1 - eta[3])
          S[4, 4, t] <- phi[2, t] * beta[3] * eta[3] * (1 - gamma)
          S[4, 5, t] <- phi[2, t] * beta[3] * eta[3] * gamma
          S[4, 6, t] <-  1-phi[2, t]
          
          S[5, 1, t] <- 0
          S[5, 2, t] <- phi[2, t] * (1 - beta[3])
          S[5, 3, t] <- phi[2, t] * beta[3] * (1 - eta[3])
          S[5, 4, t] <- phi[2, t] * beta[3] * eta[3] * (1 - gamma)
          S[5, 5, t] <- phi[2, t] * beta[3] * eta[3] * gamma
          S[5, 6, t] <-  1-phi[2, t]
          
          S[6, 1, t] <- 0
          S[6, 2, t] <- 0
          S[6, 3, t] <- 0
          S[6, 4, t] <- 0
          S[6, 5, t] <- 0
          S[6, 6, t] <- 1
        }
        
        for (i in 1:(new_per_year*n_occasions - new_per_year)) {
          if (f[i] == n_occasions) next
          
          for (t in f[i]:(n_occasions-1)) {
            
            Z[i, t+1] <- which(rmultinom(n = 1, size = 1,
                                         prob = S[Z[i, t], , t]) == 1)
            
            E[1, 1, i, t] <- values_p[k][k] 
            E[1, 2, i, t] <- 0
            E[1, 3, i, t] <- 0
            E[1, 4, i, t] <- 0
            E[1, 5, i, t] <- 0
            E[1, 6, i, t] <- 0
            E[1, 7, i, t] <- 0
            E[1, 8, i, t] <- 1 - values_p[k]
            
            E[2, 1, i, t] <- 0
            E[2, 2, i, t] <- values_p[k] * (1- KHI[i, t+1])  
            E[2, 3, i, t] <- 0
            E[2, 4, i, t] <- 0
            E[2, 5, i, t] <- KHI[i, t+1]
            E[2, 6, i, t] <- 0
            E[2, 7, i, t] <- 0
            E[2, 8, i, t] <- (1 - values_p[k]) * (1 - KHI[i, t+1])
            
            E[3, 1, i, t] <- 0
            E[3, 2, i, t] <- values_p[k] * (1- KHI[i, t+1]) 
            E[3, 3, i, t] <- 0
            E[3, 4, i, t] <- 0
            E[3, 5, i, t] <- 0
            E[3, 6, i, t] <- (1 - values_p[k]) * KHI[i, t+1]
            E[3, 7, i, t] <- values_p[k] * KHI[i, t+1]
            E[3, 8, i, t] <- (1 - values_p[k]) * (1 - KHI[i, t+1])
            
            E[4, 1, i, t] <- 0
            E[4, 2, i, t] <- 0
            E[4, 3, i, t] <- values_p[k]
            E[4, 4, i, t] <- 0
            E[4, 5, i, t] <- 0
            E[4, 6, i, t] <- (1 - values_p[k]) * KHI[i, t+1]
            E[4, 7, i, t] <- 0
            E[4, 8, i, t] <- (1 - values_p[k]) * (1 - KHI[i, t+1])
            
            E[5, 1, i, t] <- 0
            E[5, 2, i, t] <- 0
            E[5, 3, i, t] <- 0
            E[5, 4, i, t] <- values_p[k]
            E[5, 5, i, t] <- 0
            E[5, 6, i, t] <- (1 - values_p[k]) * KHI[i, t+1]
            E[5, 7, i, t] <- 0
            E[5, 8, i, t] <- (1 - values_p[k]) * (1 - KHI[i, t+1])
            
            E[6, 1, i, t] <- 0 
            E[6, 2, i, t] <- 0
            E[6, 3, i, t] <- 0
            E[6, 4, i, t] <- 0
            E[6, 5, i, t] <- 0
            E[6, 6, i, t] <- 0
            E[6, 7, i, t] <- 0
            E[6, 8, i, t] <- 1
            
          
            CH[i, t+1] <-  which(rmultinom(n = 1, size = 1,
                                           E[Z[i, t+1], , i, t]) == 1)
            
            # GPS collaring ++++++++++++++++++++++++++++++++++++++
            if (CH[i, t+1] %in% physical) {  # If physically captured
              
              collared <- rbinom(n = 1, size = 1, prob = values_pGPS[l])
              
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
                zeta[i, t, 1:n_states],  S[1:n_states, j, t])  *
                E[j, CH[i, t+1], i, t]
            }
          }
          lik[i] <- sum(zeta[i, n_occasions, 1:n_states])
          # ones[i] ~ dbin(prob = lik[i], size = freq[i])
        }
        
        # ~ 3. Bundle data -------------------------------------------------------------
        
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
                             time = time,
                             S0 = prop,
                             KHI = KHI)
  
        
        # ~ 4. Check data --------------------------------------------------------------
        
        # Let's check that all the transitions are possible. 
        S <- matrix(data = c(0, 1, 1, 1, 1, 1,  # K 
                             0, 1, 1, 1, 1, 1,  # NB 
                             0, 1, 1, 1, 1, 1,  # FB    
                             0, 1, 1, 1, 1, 1,  # SB    
                             0, 1, 1, 1, 1, 1,  # SB    
                             0, 0, 0, 0, 0, 1), # D    
                    nrow = n_states, ncol = n_states, byrow = TRUE)
        
        check_Z <- CH ; check_Z[] <- NA
        for (i in 1:dim(CH)[1]) {
          if (f[i] == n_occasions) next
          for (t in f[i]:(n_occasions - 1)) {
            check_Z[i, t] <- ifelse(S[Z[i, t], Z[i, t+1]] == 0, 1, 0)
          }
        }
        print(sum(check_Z, na.rm = T) == 0)
        
        # Let's now check that the states are compatible with the observations
        #                   J    L   SB1  SB2  NB    B    FB    D
        E <- matrix(data = c(1,   0,   0,   0,   0,   0,    0,   1,  # J
                             0,   1,   0,   0,   1,   0,    0,   1,  # NB
                             0,   1,   0,   0,   0,   1,    1,   1,  # FB
                             0,   0,   1,   0,   0,   1,    0,   1,  # SB1
                             0,   0,   0,   1,   0,   1,    0,   1,  # SB2
                             0,   0,   0,   0,   0,   0,    0,   1),  # D
                    nrow = n_states, ncol = n_events, byrow = TRUE)
        
        check_Z <- CH ; check_Z[] <- NA
        for (i in 1:dim(CH)[1]) {
          for (t in f[i]:n_occasions) {
            check_Z[i, t] <- ifelse(E[Z[i, t], CH[i, t]] == 0, 1, 0)
            if (check_Z[i, t] == 1) {
              print(paste0("i = ", i, " ; t = ", t))
              print(paste0("state: ", Z[i, t], ": event: ",  CH[i, t]))
            }
          }
        }
        print(sum(check_Z, na.rm = T) == 0)
        
        
        # Save
        name <- paste0("03_outputs/simulations/simulated_datasets/rock_partridge/with_biologging/",
                       "p_", values_p[k],
                       "_pGPS_", values_pGPS[l],
                       "_pError_", values_pError[m],
                       "_dataset_", n, ".RData")
        save(CR_model_biologging_sim, my.constants, dat, Z, zeta,
             file = name)
        
      }
    }  # pError
  }  # pGPS
}  # p

end_global <- Sys.time() 
print(paste0("global end: ", end_global)) 
print(paste0("global duration: ", end_global - start_global))                 
