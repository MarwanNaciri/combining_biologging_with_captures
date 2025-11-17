#==============================================================================#
#                                                                              #
#                    Simulate datasets albatross example                       #
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
  # ps1: probability to leave the colony 
  # psi2: probability to come back to the colony
  # beta: breeding probability
  # eta: probability of fledging the chick
  # beta: probability of breeding outcome assignement.
  
  # -------------------------------------------------+
  # States (z):
  # 1 = C: Chick (at ringing)    
  # 2 = PrBaway: Pre-Breeder that is away from the colony (not available for capture)
  # 3 = PrB: Pre-Breeder that is in the colony (available for capture)
  # 4 = NB: Non Breeder   
  # 5 = FB: Failed Breeder (chick died before fledging)
  # 6 = SB: Succsessful Breeder (chick fledged)
  # 7 = SAB: Adult in sabbatical (already bred once in its life)
  # 8 = D: Dead     
  
  # -------------------------------------------------+
  # Events (y):
  # 1 = C: Capture of chick
  # 2 = PrBaway: Remote capture of a Pre-Breeder that is away from the colony 
  # 3 = PrB: Physical or remote capture of a Pre-Breeder that is in the colony 
  # 4 = NB: Physical or remote capture of a non-breeder
  # 5 = SB: Physical or remote capture of a successful breeder
  # 6 = B: Physcial capture of a breeder (outcome of reproduction unknown)
  # 7 = SAB: Remote capture of an adult in a sabbatical year (away from the colony) 
  # 8 = NA: Not captured/observed
  
  # +++++++++++++++++++++++++++++++ matrices +++++++++++++++++++++++++++++++++++
  
  # E ++++++++++++++++
  
  # C
  E[1, 1] <- 1
  E[1, 2] <- 0
  E[1, 3] <- 0
  E[1, 4] <- 0
  E[1, 5] <- 0
  E[1, 6] <- 0
  E[1, 7] <- 0
  
  # PrBaway
  E[2, 1] <- 0
  E[2, 2] <- 0
  E[2, 3] <- 0
  E[2, 4] <- 0
  E[2, 5] <- 0
  E[2, 6] <- 0
  E[2, 7] <- 1
  
  # PrB
  E[3, 1] <- 0
  E[3, 2] <- p
  E[3, 3] <- 0
  E[3, 4] <- 0
  E[3, 5] <- 0
  E[3, 6] <- 0
  E[3, 7] <- 1-p
  
  # NB
  E[4, 1] <- 0
  E[4, 2] <- 0
  E[4, 3] <- p
  E[4, 4] <- 0
  E[4, 5] <- 0
  E[4, 6] <- 0
  E[4, 7] <- 1-p
  
  # FB
  E[5, 1] <- 0
  E[5, 2] <- 0
  E[5, 3] <- 0
  E[5, 4] <- p*theta 
  E[5, 5] <- 0
  E[5, 6] <- p*(1-theta)
  E[5, 7] <- 1-p
  
  # SB
  E[6, 1] <- 0
  E[6, 2] <- 0
  E[6, 3] <- 0
  E[6, 4] <- 0 
  E[6, 5] <- p*theta
  E[6, 6] <- p*(1-theta)
  E[6, 7] <- 1-p
  
  # SAB
  E[7, 1] <- 0
  E[7, 2] <- 0
  E[7, 3] <- 0
  E[7, 4] <- 0
  E[7, 5] <- 0
  E[7, 6] <- 0
  E[7, 7] <- 1
  
  # D
  E[8, 1] <- 0
  E[8, 2] <- 0
  E[8, 3] <- 0
  E[8, 4] <- 0
  E[8, 5] <- 0
  E[8, 6] <- 0
  E[8, 7] <- 1
  
  for (i in 1:n_ind) {
    
    for (t in f[i]:(n_occasions - 1)) {
      
      phi[i, t] <- (AGE_0[i, t] * b_phi[1] + # 0 yr old (actually 0.9 -> 1.9 yr old)
                      AGE_1[i, t] * b_phi[1] + # 1 yr old
                      AGE_2_5[i, t] * b_phi[2] + # 3-5 yr old
                      AGE_6_plus[i, t] * b_phi[3]) # 6+ yr old 
      
      # Psi1: probability of departure
      # Pre-Breeder
      psi1[1, i, t] <- (AGE_2_5[i, t] * b_psi1[1] + # 2-5 yr old
                          AGE_6_plus[i, t] * b_psi1[2])  # 6+ yr old 
      # Non-breeder
      psi1[2, i, t] <- AGE_6_plus[i, t] * b_psi1[2]
      # Failed-breeder
      psi1[3, i, t] <- AGE_6_plus[i, t] * b_psi1[3]
      # Successful-breeder
      psi1[4, i, t] <- AGE_6_plus[i, t] * b_psi1[4]
      
      
      # Psi2: probability of returning
      psi2[i, t] <- (AGE_2_5[i, t] * b_psi2[1] + # 2-5 yr old
                     AGE_6_plus[i, t] * b_psi2[2]) # 6+ yr old 
      
      
      # Eta: probability of breeding (egg-laying)
      # Pre-breeder (away of not)
      eta[1, i, t] <- (AGE_6[i, t] * b_eta[1] +
                         AGE_7_8[i, t] * b_eta[2] +
                         AGE_9_plus[i, t] * b_eta[3])
      # Adults (not away)
      eta[2, i, t] <- (AGE_7_8[i, t] + AGE_9_plus[i, t]) * b_eta[4]
      # Adults (away)
      eta[3, i, t] <- (AGE_7_8[i, t] + AGE_9_plus[i, t]) * b_eta[5]
      
      # Beta: probability of nest success
      beta[i, t] <- AGE_6_plus[i, t] * a1_beta
      
      # C
      S[1, 1, i, t] <- 0
      S[1, 2, i, t] <- phi[i, t] 
      S[1, 3, i, t] <- 0
      S[1, 4, i, t] <- 0
      S[1, 5, i, t] <- 0
      S[1, 6, i, t] <- 0
      S[1, 7, i, t] <- 0
      S[1, 8, i, t] <- 1-phi[i, t]
      
      # PrBaway
      S[2, 1, i, t] <- 0
      S[2, 2, i, t] <- phi[i, t]*(1-psi2[i, t])
      S[2, 3, i, t] <- phi[i, t]*psi2[i, t]*(1-eta[1, i, t])
      S[2, 4, i, t] <- 0
      S[2, 5, i, t] <- phi[i, t]*psi2[i, t]*eta[1, i, t]*(1-beta[i, t])
      S[2, 6, i, t] <- phi[i, t]*psi2[i, t]*eta[1, i, t]*beta[i, t]
      S[2, 7, i, t] <- 0
      S[2, 8, i, t] <- 1-phi[i, t]
      
      # PrB
      S[3, 1, i, t] <- 0
      S[3, 2, i, t] <- phi[i, t]*psi1[1, i, t]
      S[3, 3, i, t] <- phi[i, t]*(1-psi1[1, i, t])*(1-eta[1, i, t])
      S[3, 4, i, t] <- 0
      S[3, 5, i, t] <- phi[i, t]*(1-psi1[1, i, t])*eta[1, i, t]*(1-beta[i, t])
      S[3, 6, i, t] <- phi[i, t]*(1-psi1[1, i, t])*eta[1, i, t]*beta[i, t]
      S[3, 7, i, t] <- 0
      S[3, 8, i, t] <- 1-phi[i, t]
      
      # NB
      S[4, 1, i, t] <- 0
      S[4, 2, i, t] <- 0
      S[4, 3, i, t] <- 0
      S[4, 4, i, t] <- phi[i, t]*(1-psi1[2, i, t])*(1-eta[2, i, t])
      S[4, 5, i, t] <- phi[i, t]*(1-psi1[2, i, t])*eta[2, i, t]*(1-beta[i, t])
      S[4, 6, i, t] <- phi[i, t]*(1-psi1[2, i, t])*eta[2, i, t]*beta[i, t]
      S[4, 7, i, t] <- phi[i, t]*psi1[2, i, t]
      S[4, 8, i, t] <- 1-phi[i, t]
      
      # FB
      S[5, 1, i, t] <- 0
      S[5, 2, i, t] <- 0
      S[5, 3, i, t] <- 0
      S[5, 4, i, t] <- phi[i, t]*(1-psi1[3, i, t])*(1-eta[2, i, t])
      S[5, 5, i, t] <- phi[i, t]*(1-psi1[3, i, t])*eta[2, i, t]*(1-beta[i, t])
      S[5, 6, i, t] <- phi[i, t]*(1-psi1[3, i, t])*eta[2, i, t]*beta[i, t]
      S[5, 7, i, t] <- phi[i, t]*psi1[3, i, t]
      S[5, 8, i, t] <- 1-phi[i, t]
      
      # SB
      S[6, 1, i, t] <- 0
      S[6, 2, i, t] <- 0
      S[6, 3, i, t] <- 0
      S[6, 4, i, t] <- phi[i, t]*(1-psi1[4, i, t])
      S[6, 5, i, t] <- 0
      S[6, 6, i, t] <- 0
      S[6, 7, i, t] <- phi[i, t]*psi1[4, i, t]
      S[6, 8, i, t] <- 1-phi[i, t]
      
      # SAB
      S[7, 1, i, t] <- 0
      S[7, 2, i, t] <- 0
      S[7, 3, i, t] <- 0
      S[7, 4, i, t] <- phi[i, t]*psi2[i, t]*(1-eta[3, i, t])
      S[7, 5, i, t] <- phi[i, t]*psi2[i, t]*eta[3, i, t]*(1-beta[i, t])
      S[7, 6, i, t] <- phi[i, t]*psi2[i, t]*eta[3, i, t]*beta[i, t]
      S[7, 7, i, t] <- phi[i, t]*(1-psi2[i, t])
      S[7, 8, i, t] <- 1-phi[i, t]
      
      S[8, 1, i, t] <- 0
      S[8, 2, i, t] <- 0 
      S[8, 3, i, t] <- 0
      S[8, 4, i, t] <- 0
      S[8, 5, i, t] <- 0
      S[8, 6, i, t] <- 0
      S[8, 7, i, t] <- 0
      S[8, 8, i, t] <- 1
      
    }
  }
  
  # ++++++++++++++++++++++++++++++++ priors ++++++++++++++++++++++++++++++++++++
  
  #-------- State process --------#
  for (k in 1:3) {
    b_phi[k] ~ dunif(0, 1)   # survival
  }
  for (k in 1:4) {
    b_psi1[k] ~ dunif(0, 1)   # emigration probability (sabbatical) 
  }
  for (k in 1:2) {
    b_psi2[k] ~ dunif(0, 1)   # immigration probability (return)
  }
  for (k in 1:5) {
  b_eta[k] ~ dunif(0, 1)    # breeding probability 
  }
  
  a1_beta ~ dunif(0, 1)  # nest success
  
  #----- observation process -----#
  p ~ dunif(0, 1)
  theta ~ dunif(0, 1)  # prob of state assignment
  
  
  # +++++++++++++++++++++++++++ model & likelihood +++++++++++++++++++++++++++++
  
  # for (i in 1:n_ind){
  #   # latent state at first capture (chick)
  #   z[i, f[i]] <- 1
  #   
  #   for (t in ( f[i]+1 ):n_occasions) {
  #     
  #     #-------- State process --------#
  #     # draw z(t) given z(t-1)
  #     z[i, t] ~ dcat( S[z[i, t-1], 1:n_states, i, t-1] )
  #     
  #     #----- Observation process -----#
  #     # draw y(t) given z(t)
  #     y[i, t] ~ dcat( E[z[i, t], 1:n_events, i, t-1] )
  #   }
  # }
  # 
  
  for(i in 1:n_ind){

    # First capture
    zeta[i, f[i], 1] <- 1
    for (k in 2:n_states) {
      zeta[i, f[i], k] <- 0
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
  new_per_year <- 15
  n_occasions <- 15
  
  n_states <- 8
  n_events <- 7
  
  f <- rep(c(1:n_occasions), each = new_per_year)
  
  b_phi <- b_psi1 <- b_psi2 <- b_eta <- NULL
  b_phi[1] <- 0.7 # 0 yr old (actually 0.9 yr old)
  b_phi[2] <- 0.8 # 1 yr old
  b_phi[3] <- 0.88 # 2-5 yr old
  b_phi[4] <- 0.95 # 6+ yr old

  b_psi1[1] <- 0.2 # 2-5 yr old
  b_psi1[2] <- 0.2 # 6 yr+ non-breeder/pre-breeder
  b_psi1[3] <- 0.3 # 6 yr+ failed breeder 
  b_psi1[4] <- 0.8 # 6 yr+ successful breeder 

  b_psi2[1] <- 0.4 # 2-5 yr old
  b_psi2[2] <- 0.75 # 6+ yr old

  b_eta[1] <- 0.4 # 6 yr old
  b_eta[2] <- 0.7 # 7-8 yr old
  b_eta[3] <- 0.8 # 9+ yr old
  b_eta[4] <- 0.3 # Adults (not-away)
  b_eta[5] <- 0.8 # Adults (away)
  
  a1_beta <- 0.75 # nest success
  
  theta <- 0.6
  
  # Recapture probability 
  values_p <- c(0.05, 0.1, 0.25, 0.5)
  
  # Probability of being equipped with a GPS at capture
  values_pGPS <- 0 
  
  # Duration of collars,
  lambda_GPS <- 0.70
  
  # Error rate in assignment of breeding status based on remote capture (only false
  # negative i.e. missed breeding events)
  values_pError <- 0 # c(0, 0.05, 0.10, 0.20)

}

# Age matrices
AGE <- AGE_0 <- AGE_1 <- AGE_2_5 <- AGE_6 <- AGE_6_plus <- AGE_7_8 <- 
  AGE_9_plus <- matrix(NA, nrow = (new_per_year*n_occasions - new_per_year), 
                       ncol = n_occasions)
for (i in 1:(new_per_year*n_occasions - new_per_year)) {
  AGE[i, f[i]] <- age <- AGE_1[i, f[i]] <- AGE_2_5[i, f[i]] <- 
    AGE_6[i, f[i]] <- AGE_6_plus[i, f[i]] <- AGE_7_8[i, f[i]] <- AGE_9_plus[i, f[i]] <- 0
  AGE_0[i, f[i]] <- 1
  for (t in (f[i]+1):n_occasions){
    age <- age + 1
    AGE[i, t] <- age
    AGE_0[i, t] <- 0
    AGE_1[i, t] <- ifelse(AGE[i, t] == 1, 1, 0)
    AGE_2_5[i, t] <- ifelse(AGE[i, t] %in% c(2:5), 1, 0)
    AGE_6[i, t] <- ifelse(AGE[i, t] == 6, 1, 0)
    AGE_6_plus[i, t] <- ifelse(AGE[i, t] >= 6, 1, 0)
    AGE_7_8[i, t] <- ifelse(AGE[i, t] %in% c(7:8), 1, 0)
    AGE_9_plus[i, t] <- ifelse(AGE[i, t] >= 9, 1, 0)
  }
}


# ~ 2. Generate capture histories ----------------------------------------------

for (k in 1:length(values_p)) {
  for (l in 1:length(values_pGPS)) {
    for (m in 1:length(values_pError)) {
      for (n in n_sim_start:n_sim) {
        print(paste0("simulation ", n, ":"))
        
        S <- array(data = NA, 
                   dim = c(n_states, n_states, (new_per_year*n_occasions - new_per_year), n_occasions))
        
        E <- matrix(data = NA,
                    nrow = n_states, ncol = n_events, byrow = TRUE)
        
        phi <- psi2 <- beta <- matrix(NA, nrow = (new_per_year*n_occasions - new_per_year), 
                      ncol = n_occasions)
        psi1 <- array(data = NA, 
                   dim = c(4, (new_per_year*n_occasions - new_per_year), 
                           n_occasions))
        eta <- array(data = NA, 
                     dim = c(3, (new_per_year*n_occasions - new_per_year), 
                             n_occasions))
        
        # Matrices for data
        CH <- Z <-  matrix(NA, nrow = (new_per_year*n_occasions - new_per_year), 
                              ncol = n_occasions) 
        KHI <- p_realized <- matrix(0, nrow = (new_per_year*n_occasions - new_per_year), 
                                    ncol = n_occasions) 
        
        # ~~~ a. First captures --------------------------------------------------------
        
        set.seed(n)
        for (i in 1:(new_per_year*n_occasions - new_per_year)) {
          Z[i, f[i]] <- 1
          
          CH[i, f[i]] <- 1 # ringed as chicks
          
          
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
        
        for (i in 1:(new_per_year*n_occasions - new_per_year)) {
          if (f[i] == n_occasions) next
          
          for (t in f[i]:(n_occasions-1)) {
            
            phi[i, t] <- (AGE_0[i, t] * b_phi[1] + # 0 yr old (actually 0.9 -> 1.9 yr old)
                            AGE_1[i, t] * b_phi[1] + # 1 yr old
                            AGE_2_5[i, t] * b_phi[2] + # 3-5 yr old
                            AGE_6_plus[i, t] * b_phi[3]) # 6+ yr old 
            
            # Psi1: probability of departure
            # Pre-Breeder
            psi1[1, i, t] <- (AGE_2_5[i, t] * b_psi1[1] + # 3-5 yr old
                                AGE_6_plus[i, t] * b_psi1[2])  # 6+ yr old 
            # Non-breeder
            psi1[2, i, t] <- AGE_6_plus[i, t] * b_psi1[2]
            # Failed-breeder
            psi1[3, i, t] <- AGE_6_plus[i, t] * b_psi1[3]
            # Successful-breeder
            psi1[4, i, t] <- AGE_6_plus[i, t] * b_psi1[4]
            
            
            # Psi2: probability of returning
            psi2[i, t] <- (AGE_2_5[i, t] * b_psi2[1] + # 2-5 yr old
                             AGE_6_plus[i, t] * b_psi2[2]) # 6+ yr old 
            
            
            # Eta: probability of breeding (egg-laying)
            # Pre-breeder (away of not)
            eta[1, i, t] <- (AGE_6[i, t] * b_eta[1] +
                               AGE_7_8[i, t] * b_eta[2] +
                               AGE_9_plus[i, t] * b_eta[3])
            # Adults (not away)
            eta[2, i, t] <- (AGE_7_8[i, t] + AGE_9_plus[i, t]) * b_eta[4]
            # Adults (away)
            eta[3, i, t] <- (AGE_7_8[i, t] + AGE_9_plus[i, t]) * b_eta[5]
            
            # Beta: probability of nest success
            beta[i, t] <- AGE_6_plus[i, t] * a1_beta

            # C
            S[1, 1, i, t] <- 0
            S[1, 2, i, t] <- phi[i, t] 
            S[1, 3, i, t] <- 0
            S[1, 4, i, t] <- 0
            S[1, 5, i, t] <- 0
            S[1, 6, i, t] <- 0
            S[1, 7, i, t] <- 0
            S[1, 8, i, t] <- 1-phi[i, t]
            
            # PrBaway
            S[2, 1, i, t] <- 0
            S[2, 2, i, t] <- phi[i, t]*(1-psi2[i, t])
            S[2, 3, i, t] <- phi[i, t]*psi2[i, t]*(1-eta[1, i, t])
            S[2, 4, i, t] <- 0
            S[2, 5, i, t] <- phi[i, t]*psi2[i, t]*eta[1, i, t]*(1-beta[i, t])
            S[2, 6, i, t] <- phi[i, t]*psi2[i, t]*eta[1, i, t]*beta[i, t]
            S[2, 7, i, t] <- 0
            S[2, 8, i, t] <- 1-phi[i, t]
            
            # PrB
            S[3, 1, i, t] <- 0
            S[3, 2, i, t] <- phi[i, t]*psi1[1, i, t]
            S[3, 3, i, t] <- phi[i, t]*(1-psi1[1, i, t])*(1-eta[1, i, t])
            S[3, 4, i, t] <- 0
            S[3, 5, i, t] <- phi[i, t]*(1-psi1[1, i, t])*eta[1, i, t]*(1-beta[i, t])
            S[3, 6, i, t] <- phi[i, t]*(1-psi1[1, i, t])*eta[1, i, t]*beta[i, t]
            S[3, 7, i, t] <- 0
            S[3, 8, i, t] <- 1-phi[i, t]
            
            # NB
            S[4, 1, i, t] <- 0
            S[4, 2, i, t] <- 0
            S[4, 3, i, t] <- 0
            S[4, 4, i, t] <- phi[i, t]*(1-psi1[2, i, t])*(1-eta[2, i, t])
            S[4, 5, i, t] <- phi[i, t]*(1-psi1[2, i, t])*eta[2, i, t]*(1-beta[i, t])
            S[4, 6, i, t] <- phi[i, t]*(1-psi1[2, i, t])*eta[2, i, t]*beta[i, t]
            S[4, 7, i, t] <- phi[i, t]*psi1[2, i, t]
            S[4, 8, i, t] <- 1-phi[i, t]
            
            # FB
            S[5, 1, i, t] <- 0
            S[5, 2, i, t] <- 0
            S[5, 3, i, t] <- 0
            S[5, 4, i, t] <- phi[i, t]*(1-psi1[3, i, t])*(1-eta[2, i, t])
            S[5, 5, i, t] <- phi[i, t]*(1-psi1[3, i, t])*eta[2, i, t]*(1-beta[i, t])
            S[5, 6, i, t] <- phi[i, t]*(1-psi1[3, i, t])*eta[2, i, t]*beta[i, t]
            S[5, 7, i, t] <- phi[i, t]*psi1[3, i, t]
            S[5, 8, i, t] <- 1-phi[i, t]
            
            # SB
            S[6, 1, i, t] <- 0
            S[6, 2, i, t] <- 0
            S[6, 3, i, t] <- 0
            S[6, 4, i, t] <- phi[i, t]*(1-psi1[4, i, t])
            S[6, 5, i, t] <- 0
            S[6, 6, i, t] <- 0
            S[6, 7, i, t] <- phi[i, t]*psi1[4, i, t]
            S[6, 8, i, t] <- 1-phi[i, t]
            
            # SAB
            S[7, 1, i, t] <- 0
            S[7, 2, i, t] <- 0
            S[7, 3, i, t] <- 0
            S[7, 4, i, t] <- phi[i, t]*psi2[i, t]*(1-eta[3, i, t])
            S[7, 5, i, t] <- phi[i, t]*psi2[i, t]*eta[3, i, t]*(1-beta[i, t])
            S[7, 6, i, t] <- phi[i, t]*psi2[i, t]*eta[3, i, t]*beta[i, t]
            S[7, 7, i, t] <- phi[i, t]*(1-psi2[i, t])
            S[7, 8, i, t] <- 1-phi[i, t]
            
            S[8, 1, i, t] <- 0
            S[8, 2, i, t] <- 0 
            S[8, 3, i, t] <- 0
            S[8, 4, i, t] <- 0
            S[8, 5, i, t] <- 0
            S[8, 6, i, t] <- 0
            S[8, 7, i, t] <- 0
            S[8, 8, i, t] <- 1
            
            
            Z[i, t+1] <- which(rmultinom(n = 1, size = 1,
                                         prob = S[Z[i, t], , i, t]) == 1)
            
            p_realized[i, t+1] <- rbinom(n = 1, size = 1, prob = values_p[k])
            
            # C
            E[1, 1] <- 1
            E[1, 2] <- 0
            E[1, 3] <- 0
            E[1, 4] <- 0
            E[1, 5] <- 0
            E[1, 6] <- 0
            E[1, 7] <- 0
            
            # PrBaway
            E[2, 1] <- 0
            E[2, 2] <- 0
            E[2, 3] <- 0
            E[2, 4] <- 0
            E[2, 5] <- 0
            E[2, 6] <- 0
            E[2, 7] <- 1
            
            # PrB
            E[3, 1] <- 0
            E[3, 2] <- p_realized[i, t+1]
            E[3, 3] <- 0
            E[3, 4] <- 0
            E[3, 5] <- 0
            E[3, 6] <- 0
            E[3, 7] <- 1-p_realized[i, t+1]
            
            # NB
            E[4, 1] <- 0
            E[4, 2] <- 0
            E[4, 3] <- p_realized[i, t+1]
            E[4, 4] <- 0
            E[4, 5] <- 0
            E[4, 6] <- 0
            E[4, 7] <- 1-p_realized[i, t+1]
            
            # FB
            E[5, 1] <- 0
            E[5, 2] <- 0
            E[5, 3] <- 0
            E[5, 4] <- p_realized[i, t+1]*theta 
            E[5, 5] <- 0
            E[5, 6] <- p_realized[i, t+1]*(1-theta)
            E[5, 7] <- 1-p_realized[i, t+1]
            
            # SB
            E[6, 1] <- 0
            E[6, 2] <- 0
            E[6, 3] <- 0
            E[6, 4] <- 0 
            E[6, 5] <- p_realized[i, t+1]*theta
            E[6, 6] <- p_realized[i, t+1]*(1-theta)
            E[6, 7] <- 1-p_realized[i, t+1]
            
            # SAB
            E[7, 1] <- 0
            E[7, 2] <- 0
            E[7, 3] <- 0
            E[7, 4] <- 0
            E[7, 5] <- 0
            E[7, 6] <- 0
            E[7, 7] <- 1
            
            # D
            E[8, 1] <- 0
            E[8, 2] <- 0
            E[8, 3] <- 0
            E[8, 4] <- 0
            E[8, 5] <- 0
            E[8, 6] <- 0
            E[8, 7] <- 1
            
            CH[i, t+1] <-  which(rmultinom(n = 1, size = 1,
                                           E[Z[i, t+1], ]) == 1)
            
            
            # GPS collaring ++++++++++++++++++++++++++++++++++++++
            if (p_realized[i, t+1] == 1) {  # If physically captured
              
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
          zeta[i, f[i], ] <- c(1, rep(0, times = n_states-1))
    
          # Subsequent captures
          for(t in f[i] : (n_occasions-1) ){
            for (j in 1:n_states) {
              zeta[i, (t+1), j] <- inprod(
                zeta[i, t, 1:n_states],  S[1:n_states, j, i, t])  *
                E[j, CH[i, t+1]]
            }
          }
          lik[i] <- sum(zeta[i, n_occasions, 1:n_states])
          # ones[i] ~ dbin(prob = lik[i], size = freq[i])
        }
        
        # ~ 3. Bundle data -------------------------------------------------------------
        
        dat <- list(ones = rep(1, times = nrow(CH)))
        my.constants <- list(CH = CH,
                             freq = rep(1, times = nrow(CH)),
                             f = f,
                             n_ind = dim(CH)[1],
                             n_occasions = dim(CH)[2],
                             n_states = n_states,
                             KHI = KHI,
                             AGE_0 = AGE_0,
                             AGE_1 = AGE_1,
                             AGE_2_5 = AGE_2_5,
                             AGE_6 = AGE_6,
                             AGE_6_plus = AGE_6_plus,
                             AGE_7_8 = AGE_7_8,
                             AGE_9_plus = AGE_9_plus)
        
        # dat <- list(y = CH)
        # my.constants <- list(f = f,
        #                      n_ind = dim(CH)[1],
        #                      n_occasions = dim(CH)[2],
        #                      n_states = n_states,
        #                      n_events = n_events,
        #                      KHI = KHI,
        #                      AGE_1 = AGE_1,
        #                      AGE_2 = AGE_2,
        #                      AGE_2_5 = AGE_2_5,
        #                      AGE_6 = AGE_6,
        #                      AGE_6_plus = AGE_6_plus,
        #                      AGE_7_8 = AGE_7_8,
        #                      AGE_9_plus = AGE_9_plus)
                             
        
        # ~ 4. Check data --------------------------------------------------------------
        
        # Let's check that all the transitions are possible. 
        S <- matrix(data = c(0, 1, 0, 0, 0, 0, 0, 1,
                             0, 1, 1, 0, 1, 1, 0, 1,
                             0, 1, 1, 0, 1, 1, 0, 1,
                             0, 0, 0, 1, 1, 1, 1, 1,
                             0, 0, 0, 1, 1, 1, 1, 1,
                             0, 0, 0, 1, 0, 0, 1, 1,
                             0, 0, 0, 1, 1, 1, 1, 1,
                             0, 0, 0, 0, 0, 0, 0, 1),  
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
        #                    
        E <- matrix(data = c(1,  0,  0,  0,  0,  0,  0,   # C    
                             0,  0,  0,  0,  0,  0,  1,   # PrBaway 
                             0,  1,  0,  0,  0,  0,  1,  # PrB
                             0,  0,  1,  0,  0,  0,  1,  # NB    
                             0,  0,  0,  1,  0,  1,  1,  # FB     5
                             0,  0,  0,  0,  1,  1,  1,  # SB    
                             0,  0,  0,  0,  0,  0,  1,   # SAB
                             0,  0,  0,  0,  0,  0,  1),   # D  
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
        name <- paste0("03_outputs/simulations/simulated_datasets/albatross/no_biologging/",
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

