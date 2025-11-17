library(tidyverse)
library(nimble)


## 1. Load and process data ----------------------------------------------------

data_CR_events <- read_csv("01_inputs/CR_events.csv", show_col_types = F) %>%
  filter(year < 1994)

# Remove individuals marked on the last occasion
nbr_captures <- table(data_CR_events$ID_NR)
to_remove <- NULL 
for (k in 1:length(nbr_captures)) {
  if (nbr_captures[k] == 1) {
    if (data_CR_events$year[data_CR_events$ID_NR == names(nbr_captures[k])] == max(data_CR_events$year)) {
      to_remove <- c(to_remove, names(nbr_captures[k]))
    }
  }
}
data_CR_events <- data_CR_events %>% 
  filter(!ID_NR %in% to_remove)

n_states <- 12
not_captured <- 12

CH <- matrix(data = NA, 
             nrow = length(unique(data_CR_events$ID_NR)), 
             ncol = length(min(data_CR_events$year):max(data_CR_events$year)),
             dimnames = list(unique(data_CR_events$ID_NR),
                             min(data_CR_events$year):max(data_CR_events$year)))

for (k in 1:nrow(data_CR_events)) {
  CH[which(rownames(CH) == data_CR_events$ID_NR[k]), 
     which(colnames(CH) == as.character(data_CR_events$year[k]))] <- data_CR_events$event[k]
}
CH[is.na(CH)] <- not_captured 



# First & last capture +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
f <- f_event <- l <- NULL
for (i in 1:dim(CH)[1]) {
  f[i] <- min(data_CR_events$year[as.character(data_CR_events$ID_NR) == rownames(CH)[i]]) - 
    min(data_CR_events$year) + 1
  f_event[i] <- CH[i, f[i]]
  
  l[i] <- ifelse(is.na(unique(data_CR_events$year_death[as.character(data_CR_events$ID_NR) == rownames(CH)[i]])),
                 dim(CH)[2], 
                 data_CR_events$year_death[as.character(data_CR_events$ID_NR) == rownames(CH)[i]] - min(data_CR_events$year) + 1)
  if (f[i] > 1) {
    CH[i, 1:(f[i]-1)] <- NA
  }
}
# l[l > dim(CH)[2]] <- dim(CH)[2]



# DOY ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DOY <- matrix(data = NA, 
              nrow = length(unique(data_CR_events$ID_NR)), 
              ncol = length(min(data_CR_events$year):max(data_CR_events$year)),
              dimnames = list(unique(data_CR_events$ID_NR),
                              min(data_CR_events$year):max(data_CR_events$year)))
for (k in 1:nrow(data_CR_events)) {
  i <- which(rownames(CH) == data_CR_events$ID_NR[k])
  t <- which(colnames(CH) == as.character(data_CR_events$year[k]))
  
  if (t == f[i]) next
  DOY[i, t-1] <- data_CR_events$day_number[k]
}


# Here I don't directly compute the ALPAH matrix from the DOY matrix because the 
# indexation by column is not the same 
ALPHA <- matrix(data = NA, 
                nrow = length(unique(data_CR_events$ID_NR)), 
                ncol = length(min(data_CR_events$year):max(data_CR_events$year)),
                dimnames = list(unique(data_CR_events$ID_NR),
                                min(data_CR_events$year):max(data_CR_events$year)))
for (k in 1:nrow(data_CR_events)) {
  i <- which(rownames(CH) == data_CR_events$ID_NR[k])
  t <- which(colnames(CH) == as.character(data_CR_events$year[k]))
  
  ALPHA[i, t] <- plogis(-3.046745 + 0.03542667*data_CR_events$day_number[k]) 
}
ALPHA[is.na(ALPHA)] <- 0

DOY <- (DOY - mean(DOY, na.rm = T))/sd(DOY, na.rm = T)
DOY[is.na(DOY)] <- 0


# Age ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AGE <- matrix(data = NA, nrow = length(unique(data_CR_events$ID_NR)), 
              ncol = length(min(data_CR_events$year):max(data_CR_events$year)),
              dimnames = list(unique(data_CR_events$ID_NR),
                              min(data_CR_events$year):max(data_CR_events$year)))

AGE_4 <- AGE_5 <- AGE_YOUNG <- AGE_MIDDLE <- AGE_OLD <- AGE_2_4 <-
  AGE_YOUNG_1 <- AGE_MIDDLE_1 <- AGE_OLD_1 <- AGE_YOUNG_2 <- AGE_MIDDLE_2 <- AGE_OLD_2 <-
  matrix(data = 0, nrow = length(unique(data_CR_events$ID_NR)), 
         ncol = length(min(data_CR_events$year):max(data_CR_events$year)),
         dimnames = list(unique(data_CR_events$ID_NR),
                         min(data_CR_events$year):max(data_CR_events$year)))
for (i in 1:dim(CH)[1]) {
  AGE[i, f[i]] <- data_CR_events$age_for_analyses[data_CR_events$ID_NR == rownames(CH)[i]][1]
  
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
  if (f[i] == dim(CH)[2]) next
  for (t in (f[i]+1): dim(CH)[2]) {
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

# Year 2020 without captures ++++++++++++++++++++++++++++++++++++++++++++++++++
# I'll use a vector that takes the value one if the year is not 2020 and 0 if 2020
# Then I'll use this vector to multiply the probability of capture. This way, it
# is null in 2020 but not in other years.
year_not_2020 <- ifelse(as.numeric(colnames(CH)[1:( dim(CH)[2]-1)])+1 != 2020, 1, 0)


# Years 1989, 90, and 91 with no REcaptures of A0 females ++++++++++++++++++++++
year_post_1992 <- ifelse((as.numeric(colnames(CH)[1:( dim(CH)[2]-1)])+1) %in% 1989:1991, 0, 1) 



# Space-use strategy (=ecotype) ++++++++++++++++++++++++++++++++++++++++++++++++
SpPel <- SpUnk <- NULL
for (i in 1:nrow(CH)) {
  SpPel[i] <- ifelse(data_CR_events$ecotype[data_CR_events$ID_NR == rownames(CH)[i]][1] %in% "O", 
                     1, 0)
  SpUnk[i] <- ifelse(is.na(data_CR_events$ecotype[data_CR_events$ID_NR == rownames(CH)[i]][1]), 
                     1, 0)
}





sea_ice_data <- read_csv("01_inputs/sea_ice_metrics.csv", show_col_types = F)
sea_ice <- NULL
for (t in 1:(ncol(CH)-1)) {
  sea_ice[t] <- sea_ice_data$ice_free_days[which(sea_ice_data$year == as.numeric(colnames(CH)[t]))]   
}
sea_ice_s <- as.vector(scale(sea_ice))



# 2. Compute proportion of each state at first capture -------------------------

# With the marginalized likelihood, we don't need initial values for the latent 
# state. However, to compute the proportion of each state at first capture, we can
# use the proposed initial value for the first capture of each individual
z_inits <- NULL
for (i in 1:nrow(CH)) {
  obs_i <- CH[i, ]
  z_inits_i <- rep(11, times = ncol(CH))
  
  # First, the "constrained" observations
  J2 <- which(obs_i %in% 1)
  J3 <- which(obs_i %in% 2)
  SA4 <- which(obs_i %in% 3)
  SA5 <- which(obs_i %in% 4)
  A01 <- which(obs_i %in% 5)
  A02 <- which(obs_i %in% 6)
  A11 <- which(obs_i %in% 7)
  A12 <- which(obs_i %in% 8)
  A21 <- which(obs_i %in% 9)
  A22 <- which(obs_i %in% 10)
  
  # If captured as J2
  if (length(J2) > 0) {
    z_inits_i[J2] <- 1        # State J2 at t
    if (J2 < l[i]) {
      z_inits_i[J2 + 1] <- 2     # State J3 at t+1
      if ((J2 + 1) < l[i]) {
        z_inits_i[J2 + 2] <- 3     # State SA4 at t+2
        if ((J2 + 2) < l[i]) {
          z_inits_i[J2 + 3] <- 4     # State SA5 at t+3
        }
      }
    }
  }
  
  # If captured as J3 
  if (length(J3) > 0) {
    z_inits_i[J3] <- 2       # state J3 at t
    if (J3 < l[i]) {
      z_inits_i[J3 + 1] <- 3     # state SA4 at t+1
      if ((J3 + 1) < l[i]) {
        z_inits_i[J3 + 2] <- 4     # State SA5 at t+3
      }
    }
  }
  
  # If captured as SA4
  if (length(SA4) > 0) {
    z_inits_i[SA4] <- 3     # state SA4 at t
    if (SA4 < l[i]) {
      z_inits_i[SA4 + 1] <- 4     # state SA5 at t+1
    }
  }
  
  # If captured as SA5
  if (length(SA5) > 0) {
    z_inits_i[SA5] <- 4     # state SA5 at t
  }
  
  # If captured as A21
  if (length(A21) == 1) {    # If captured only once as A21 
    z_inits_i[A21] <- 9          # state AS1 at t 
    if (A21 > f[i]) {
      z_inits_i[A21 - 1] <- 7       # state A11 at t-1
      if (A21 > (f[i] + 1)) {
        z_inits_i[A21 - 2] <- 6       # state A02 at t-2
      }
    }  
  } else {
    if (length(A21) >= 2) {  # If captured more than once as A21
      z_inits_i[A21[1]] <- 9          # state AS1 at t 
      if (A21[1] > f[i]) {
        z_inits_i[A21[1] - 1] <- 7       # state A11 at t-1
        if (A21[1] > (f[i] + 1)) {
          z_inits_i[A21[1] - 2] <- 6        # state A02 at t-2
        }
      } 
      for (t in 2:length(A21)) {
        z_inits_i[A21[t]] <- 9        # state AS1 at t
        z_inits_i[A21[t] - 1] <- 7     # state A11 at t-1
        z_inits_i[A21[t] - 2] <- 6     # state A02 at t-2
      }
    }
  }
  
  # If captured as A22
  if (length(A22) == 1) {    # If captured only once as A22 
    z_inits_i[A22] <- 10          # state AS2 at t 
    if (A22 > f[i]) {
      z_inits_i[A22 - 1] <- 8       # state A12 at t-1
      if (A22 > (f[i] + 1)) {
        z_inits_i[A22 - 2] <- 6       # state A02 at t-2
      }
    }  
  } else {
    if (length(A22) >= 2) {  # If captured more than once as A22
      z_inits_i[A22[1]] <- 10          # state AS2 at t 
      if (A22[1] > f[i]) {
        z_inits_i[A22[1] - 1] <- 8       # state A12 at t-1
        if (A22[1] > (f[i] + 1)) {
          z_inits_i[A22[1] - 2] <- 6        # state A02 at t-2
        }
      } 
      for (t in 2:length(A22)) {
        z_inits_i[A22[t]] <- 10        # state AS2 at t
        z_inits_i[A22[t] - 1] <- 8    # state A12 at t-1
        z_inits_i[A22[t] - 2] <- 6     # state A02 at t-2
      }
    }
  }
  
  # If captured as A11
  if (length(A11) == 1) {      # If captured only once as A11 
    z_inits_i[A11] <- 7            # state A11 at t
    if (A11 > f[i]) {
      z_inits_i[A11 - 1] <- 6       # state A02 at t-1 
    }
  } else {
    if (length(A11) >= 2) {    # If captured more than once as A11 
      z_inits_i[A11[1]] <- 7        # state A11 at t
      if (A11[1] > f[i]) {
        z_inits_i[A11[1] - 1] <- 6     # state A02 at t-1 
      }
      for (t in 2:length(A11)) {
        z_inits_i[A11[t]] <- 7       # state A11 at t
        z_inits_i[A11[t] - 1] <- 6       # state A02 at t-1 
      }
    }
  }
  
  # If captured as A12
  if (length(A12) == 1) {      # If captured only once as A12
    z_inits_i[A12] <- 8          # state A12 at t 
    if (A12 > f[i]) {
      z_inits_i[A12 - 1] <- 6       # state A02 at t-1 
    }
  } else {
    if (length(A12) >= 2) {    # If captured more than once as A12
      z_inits_i[A12[1]] <- 8        # state A12 at t
      if (A12[1] > f[i]) {
        z_inits_i[A12[1] - 1] <- 6     # state A02 at t-1 
      }
      for (t in 2:length(A12)) {
        z_inits_i[A12[t]] <- 8       # state A12 at t
        z_inits_i[A12[t] - 1] <- 6       # state A02 at t-1 
      }
    }
  }
  
  # If captured as A01
  if (length(A01) == 1) {    # If captured only once as A01 
    z_inits_i[A01] <- 5   # state A01
  } else {    
    if (length(A01) >= 2) {  # If captured more than once as A01
      z_inits_i[A01[1]] <- 5   # state A01
      for (t in 2:length(A01)) {
        z_inits_i[A01[t]] <- 5      # state A01
      }
    }
  }
  
  # If captured as A02
  if (length(A02) == 1) {    # If captured only once as A02 
    z_inits_i[A02] <- 6   # state A02
  } else {    
    if (length(A02) >= 2) {  # If captured more than once as A02
      z_inits_i[A02[1]] <- 6  # state A02
      for (t in 2:length(A02)) {
        z_inits_i[A02[t]] <- 6      # state A02
      }
    }
  }
  
  # fill remaining cells
  if (f[i] > 1) {
    z_inits_i[1:(f[i]-1)] <- NA
  }
  if (l[i] < ncol(CH)) {
    z_inits_i[(l[i]+1):ncol(CH)] <- NA
  }  
  
  z_inits <- rbind(z_inits, z_inits_i)
}
rownames(z_inits) <- rownames(CH) # 1:dim(z_inits)[1]
colnames(z_inits) <- colnames(CH)
table(z_inits)


# Let's check that all the transitions are possible. 
S <- matrix(data = c(0, 1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  # J2    
                     0, 0,  1, 0, 0, 0, 0, 0, 0, 0, 0, 1,  # J3   
                     0, 0,  0, 1, 1, 1, 0, 0, 0, 0, 0, 1,  # SA4 
                     0, 0,  0, 0, 1, 1, 0, 0, 0, 0, 1, 1,  # SA5 
                     0, 0,  0, 0, 1, 1, 1, 0, 0, 0, 1, 1,  # A01   5  
                     0, 0,  0, 0, 1, 1, 1, 1, 0, 0, 1, 1,  # A02  
                     0, 0,  0, 0, 1, 1, 0, 0, 1, 0, 1, 1,  # A11    
                     0, 0,  0, 0, 1, 1, 0, 0, 1, 1, 1, 1,  # A12   
                     0, 0,  0, 0, 1, 1, 0, 0, 0, 0, 1, 1,  # AS1
                     0, 0,  0, 0, 1, 1, 0, 0, 0, 0, 1, 1,  # AS2   10
                     0, 0,  0, 0, 1, 1, 0, 0, 0, 0, 1, 1,  # A
                     0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1), # D    
            nrow = 12, ncol = 12, byrow = TRUE)

check_z_inits <- CH ; check_z_inits[] <- NA
for (i in 1:nrow(CH)) {
  if (f[i] == ncol(CH)) next
  for (t in f[i]:min(ncol(CH)-1, l[i]-1)) {
    check_z_inits[i, t] <- ifelse(S[z_inits[i, t], z_inits[i, t+1]] == 0, 1, 0)
  }
}
sum(check_z_inits, na.rm = T) == 0
# All good

# Let's now check that the states are compatbile with the observations
E <- matrix(data = c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,    # J2
                     0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,    # J3
                     0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1,    # SA4
                     0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,    # SA5  
                     0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,    # A01
                     0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1,    # A02
                     0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1,    # A11    
                     0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1,    # A12
                     0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1,    # AS1
                     0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,    # AS2
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,    # FA    
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),   # D
            nrow = 12, ncol = 12, byrow = TRUE)

check_z_inits <- CH ; check_z_inits[] <- NA
for (i in 1:nrow(CH)) {
  for (t in f[i]:min(ncol(CH), l[i])) {
    check_z_inits[i, t] <- ifelse(E[z_inits[i, t], CH[i, t]] == 0, 1, 0)
    if (check_z_inits[i, t] == 1) {
      print(paste0("state: ", z_inits[i, t], ": event: ",  CH[i, t]))
    }
  }
}
sum(check_z_inits, na.rm = T) == 0
# All good


# Proportion of each state at first capture ++++++++++++++++++++++++++
# (must be derived from the data)
f_state <- NULL
for (i in 1:nrow(CH)) {
  f_state[i] <- z_inits[i, f[i]]
}
prop_raw <- table(f_state)/length(f_state)
prop <- NULL
for (k in 1:n_states) {
  x <- prop_raw[as.character(k)]
  prop[k] <- ifelse(!is.na(x), x, 0)
}
names(prop) <- 1:n_states



# 3. Bundle data ---------------------------------------------------------------

dat <- list(ones = rep(1, times = nrow(CH)),
            constraint = rep(1, times = 3))

my.constants <- list(CH = CH,
                     f = f,
                     l = l,
                     S0 = prop,
                     n_ind = nrow(CH),
                     n_intervals = ncol(CH)-1,
                     n_states = n_states,
                     year_not_2020 = year_not_2020,
                     year_post_1992 = year_post_1992,
                     
                     AGE_2_4 = AGE_2_4,
                     AGE_4 = AGE_4,
                     AGE_YOUNG = AGE_YOUNG,
                     AGE_MIDDLE = AGE_MIDDLE,
                     AGE_OLD = AGE_OLD,
                     AGE_YOUNG_1 = AGE_YOUNG_1,
                     AGE_MIDDLE_1 = AGE_MIDDLE_1,
                     AGE_OLD_1 = AGE_OLD_1,
                     AGE_YOUNG_2 = AGE_YOUNG_2,
                     AGE_MIDDLE_2 = AGE_MIDDLE_2,
                     AGE_OLD_2 = AGE_OLD_2,
                     
                     DOY = DOY,
                     ALPHA = ALPHA,
                     SpUnk = SpUnk,
                     SpPel = SpPel,
                     sea_ice_s = sea_ice_s,
                     time_s = as.vector(scale(as.numeric(colnames(CH)))),
                     
                     freq = rep(1, times = nrow(CH)))


inits <- function() list(beta_phi = c(2, 3.18, 1.6, 0, 0) + rnorm(n = 5, mean = 0, sd = 0.1),
                         beta_s0  = c(-0.24, 0.66, 0.08, 0, 0, 0) + rnorm(n = 6, mean = 0, sd = 0.05), 
                         beta_s1  = c(0.30, 0.85, 1, 0) + rnorm(n = 4, mean = 0, sd = 0.05), 
                         
                         beta_beta  = c(-1.45, 0.6, 0.7, 0.60, 0, -2, -0.8, 0.4, 0) + rnorm(n = 9, mean = 0, sd = 0.5),
                         
                         beta_gamma = c(0.95, 0.8, 0.7, 0, 0, 0) + rnorm(n = 6, mean = 0, sd = 0.5),
                         
                         beta_p = c(-1.35, -2.75, -1.6, -0.1) + rnorm(n = 4, mean = 0, sd = 0.05),
                         omega = rnorm(ncol(CH)-1, mean = 0, sd = 1),
                         sigma_omega = runif(1, min = 0, max = 5))

inits_values <- list(inits(), inits())

# Parameters monitored 
params <- c("beta_phi",
            "beta_s0", "beta_s1",
            "beta_beta", "beta_gamma", 
            "beta_p",  "omega", "sigma_omega")


# 4. Build the model -----------------------------------------------------------

# Marginalized formulation of the CR model to speed up computation and more
# effectively sample the posterior distribution.

CR_model <- nimbleCode({
  # +++++++++++++++++++++++++++++++ matrices +++++++++++++++++++++++++++++++++++
  
  for (i in 1:n_ind) {
    
    # ~~~ a. E0 ----------------------------------
    E0[1, 1, i] <- 1
    E0[1, 2, i] <- 0
    E0[1, 3, i] <- 0
    E0[1, 4, i] <- 0
    E0[1, 5, i] <- 0
    E0[1, 6, i] <- 0
    E0[1, 7, i] <- 0
    E0[1, 8, i] <- 0
    E0[1, 9, i] <- 0
    E0[1, 10, i] <- 0
    E0[1, 11, i] <- 0
    E0[1, 12, i] <- 0
    
    E0[2, 1, i] <- 0
    E0[2, 2, i] <- 1
    E0[2, 3, i] <- 0
    E0[2, 4, i] <- 0
    E0[2, 5, i] <- 0
    E0[2, 6, i] <- 0
    E0[2, 7, i] <- 0
    E0[2, 8, i] <- 0
    E0[2, 9, i] <- 0
    E0[2, 10, i] <- 0
    E0[2, 11, i] <- 0
    E0[2, 12, i] <- 0
    
    E0[3, 1, i] <- 0
    E0[3, 2, i] <- 0
    E0[3, 3, i] <- 1
    E0[3, 4, i] <- 0
    E0[3, 5, i] <- 0
    E0[3, 6, i] <- 0
    E0[3, 7, i] <- 0
    E0[3, 8, i] <- 0
    E0[3, 9, i] <- 0
    E0[3, 10, i] <- 0
    E0[3, 11, i] <- 0
    E0[3, 12, i] <- 0
    
    E0[4, 1, i] <- 0
    E0[4, 2, i] <- 0
    E0[4, 3, i] <- 0
    E0[4, 4, i] <- 1
    E0[4, 5, i] <- 0
    E0[4, 6, i] <- 0
    E0[4, 7, i] <- 0
    E0[4, 8, i] <- 0
    E0[4, 9, i] <- 0
    E0[4, 10, i] <- 0
    E0[4, 11, i] <- 0
    E0[4, 12, i] <- 0
    
    E0[5, 1, i] <- 0
    E0[5, 2, i] <- 0
    E0[5, 3, i] <- 0
    E0[5, 4, i] <- 0
    E0[5, 5, i] <- 1
    E0[5, 6, i] <- 0
    E0[5, 7, i] <- 0
    E0[5, 8, i] <- 0
    E0[5, 9, i] <- 0
    E0[5, 10, i] <- 0
    E0[5, 11, i] <- 0
    E0[5, 12, i] <- 0
    
    E0[6, 1, i] <- 0
    E0[6, 2, i] <- 0
    E0[6, 3, i] <- 0
    E0[6, 4, i] <- 0
    E0[6, 5, i] <- 0
    E0[6, 6, i] <- 1
    E0[6, 7, i] <- 0
    E0[6, 8, i] <- 0
    E0[6, 9, i] <- 0
    E0[6, 10, i] <- 0
    E0[6, 11, i] <- 0
    E0[6, 12, i] <- 0
    
    E0[7, 1, i] <- 0
    E0[7, 2, i] <- 0
    E0[7, 3, i] <- 0
    E0[7, 4, i] <- 0
    E0[7, 5, i] <- 0
    E0[7, 6, i] <- 0
    E0[7, 7, i] <- 1
    E0[7, 8, i] <- 0
    E0[7, 9, i] <- 0
    E0[7, 10, i] <- 0
    E0[7, 11, i] <- 0
    E0[7, 12, i] <- 0
    
    E0[8, 1, i] <- 0
    E0[8, 2, i] <- 0
    E0[8, 3, i] <- 0
    E0[8, 4, i] <- 0
    E0[8, 5, i] <- 0
    E0[8, 6, i] <- 0
    E0[8, 7, i] <- 0
    E0[8, 8, i] <- 1
    E0[8, 9, i] <- 0
    E0[8, 10, i] <- 0
    E0[8, 11, i] <- 0
    E0[8, 12, i] <- 0
    
    E0[9, 1, i] <- 0
    E0[9, 2, i] <- 0
    E0[9, 3, i] <- 0
    E0[9, 4, i] <- 0
    E0[9, 5, i] <- 0
    E0[9, 6, i] <- 0
    E0[9, 7, i] <- 0
    E0[9, 8, i] <- 0
    E0[9, 9, i] <- 1-ALPHA[i, f[i]]
    E0[9, 10, i] <- 0
    E0[9, 11, i] <- ALPHA[i, f[i]]
    E0[9, 12, i] <- 0
    
    E0[10, 1, i] <- 0
    E0[10, 2, i] <- 0
    E0[10, 3, i] <- 0
    E0[10, 4, i] <- 0
    E0[10, 5, i] <- 0
    E0[10, 6, i] <- 0
    E0[10, 7, i] <- 0
    E0[10, 8, i] <- 0
    E0[10, 9, i] <- 2 * ( 1-ALPHA[i, f[i]] ) * ALPHA[i, f[i]]
    E0[10, 10, i] <- ( 1 - ALPHA[i, f[i]] )^2
    E0[10, 11, i] <- ALPHA[i, f[i]]^2
    E0[10, 12, i] <- 0
    
    E0[11, 1, i] <- 0
    E0[11, 2, i] <- 0
    E0[11, 3, i] <- 0
    E0[11, 4, i] <- 0
    E0[11, 5, i] <- 0
    E0[11, 6, i] <- 0
    E0[11, 7, i] <- 0
    E0[11, 8, i] <- 0
    E0[11, 9, i] <- 0
    E0[11, 10, i] <- 0
    E0[11, 11, i] <- 1
    E0[11, 12, i] <- 0
    
    E0[12, 1, i] <- 0
    E0[12, 2, i] <- 0
    E0[12, 3, i] <- 0
    E0[12, 4, i] <- 0
    E0[12, 5, i] <- 0
    E0[12, 6, i] <- 0
    E0[12, 7, i] <- 0
    E0[12, 8, i] <- 0
    E0[12, 9, i] <- 0
    E0[12, 10, i] <- 0
    E0[12, 11, i] <- 0
    E0[12, 12, i] <- 1
    
    # ~~~ b. Regressions on parameters -------------------
    
    for (t in f[i]:(l[i] - 1)) { 
      
      # Recapture probability
      # All females except A0 females
      p[i, t+1, 1] <- ilogit(beta_p[1] +   
                               beta_p[2] * SpUnk[i] +
                               beta_p[3] * SpPel[i] +
                               beta_p[4] * SpPel[i] * time_s[t] + 
                               omega[t]) * year_not_2020[t]
      # A0 females
      p[i, t+1, 2] <- ilogit(beta_p[1] +   
                               beta_p[2] * SpUnk[i] +
                               beta_p[3] * SpPel[i] +
                               beta_p[4] * SpPel[i] * time_s[t] + 
                               omega[t]) * year_not_2020[t] * year_post_1992[t]
      
      # Probability of survival
      logit(phi[i, t]) <- beta_phi[1] * AGE_2_4[i, t] +
        beta_phi[2] * (AGE_YOUNG[i, t] + AGE_MIDDLE[i, t]) +
        beta_phi[3] * AGE_OLD[i, t] +
        beta_phi[4] * (SpUnk[i] + SpPel[i]) +
        beta_phi[5] * (AGE_2_4[i, t] + AGE_OLD[i, t]) * sea_ice_s[t]
      
      
      # Probability of coy survival
      logit(s01[i, t]) <- beta_s0[1] * AGE_YOUNG_1[i, t] +
        beta_s0[2] * AGE_MIDDLE_1[i, t] +
        beta_s0[3] * AGE_OLD_1[i, t] +
        beta_s0[4] +
        beta_s0[5] * (SpUnk[i] + SpPel[i]) +
        beta_s0[6] * sea_ice_s[t]
      
      logit(s02[i, t]) <- beta_s0[1] * AGE_YOUNG_1[i, t] +
        beta_s0[2] * AGE_MIDDLE_1[i, t] +
        beta_s0[3] * AGE_OLD_1[i, t] +
        beta_s0[5] * (SpUnk[i] + SpPel[i]) +
        beta_s0[6] * sea_ice_s[t]
      
      # Probability of yearling survival
      logit(s1[i, t]) <- beta_s1[1] * AGE_YOUNG_2[i, t] +
        beta_s1[2] * AGE_MIDDLE_2[i, t] +
        beta_s1[3] * AGE_OLD_2[i, t] +
        beta_s1[4] * sea_ice_s[t]
      
      # Probability of breeding
      logit(beta[i, t, 1]) <- beta_beta[1] * AGE_4[i, t] +   # Lone female
        beta_beta[2] * AGE_YOUNG[i, t] +
        beta_beta[3] * AGE_MIDDLE[i, t] +
        beta_beta[4] * AGE_OLD[i, t] +
        beta_beta[5] * DOY[i, t] +
        beta_beta[9] * sea_ice_s[t]
      
      logit(beta[i, t, 2]) <- beta_beta[2] * AGE_YOUNG[i, t] + # A0-
        beta_beta[3] * AGE_MIDDLE[i, t] +
        beta_beta[4] * AGE_OLD[i, t] +
        beta_beta[5] * DOY[i, t] +
        beta_beta[6] +
        beta_beta[9] * sea_ice_s[t]
      
      logit(beta[i, t, 3]) <- beta_beta[2] * AGE_YOUNG[i, t] + # A1-
        beta_beta[3] * AGE_MIDDLE[i, t] +
        beta_beta[4] * AGE_OLD[i, t] +
        beta_beta[5] * DOY[i, t] +
        beta_beta[7] +
        beta_beta[9] * sea_ice_s[t]
      
      logit(beta[i, t, 4]) <- beta_beta[2] * AGE_YOUNG[i, t] + # AS
        beta_beta[3] * AGE_MIDDLE[i, t] +
        beta_beta[4] * AGE_OLD[i, t] +
        beta_beta[5] * DOY[i, t] +
        beta_beta[8] +
        beta_beta[9] * sea_ice_s[t]
      
      # Twinning prob
      logit(gamma[i, t]) <- beta_gamma[1] * (AGE_4[i, t] + AGE_YOUNG[i, t]) +
        beta_gamma[2] * AGE_MIDDLE[i, t] +
        beta_gamma[3] * AGE_OLD[i, t] +
        beta_gamma[4] * (AGE_4[i, t] + AGE_YOUNG[i, t]) * DOY[i, t] +
        beta_gamma[5] * (SpUnk[i] + SpPel[i]) +
        beta_gamma[6] * sea_ice_s[t]
      
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
      S[3, 4, i, t] <- phi[i, t] * (1- beta[i, t, 1])
      S[3, 5, i, t] <- phi[i, t] * beta[i, t, 1] * (1 - gamma[i, t])
      S[3, 6, i, t] <- phi[i, t] * beta[i, t, 1] * gamma[i, t]
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
      S[4, 5, i, t] <- phi[i, t] * beta[i, t, 1] * (1 - gamma[i, t])
      S[4, 6, i, t] <- phi[i, t] * beta[i, t, 1] * gamma[i, t]
      S[4, 7, i, t] <- 0
      S[4, 8, i, t] <- 0
      S[4, 9, i, t] <- 0
      S[4, 10, i, t] <- 0
      S[4, 11, i, t] <- phi[i, t] * (1 - beta[i, t, 1])
      S[4, 12, i, t] <- 1-phi[i, t]
      
      S[5, 1, i, t] <- 0
      S[5, 2, i, t] <- 0
      S[5, 3, i, t] <- 0
      S[5, 4, i, t] <- 0
      S[5, 5, i, t] <- phi[i, t] * (1-s01[i, t]) * beta[i, t, 2] * (1 - gamma[i, t])
      S[5, 6, i, t] <- phi[i, t] * (1-s01[i, t]) * beta[i, t, 2] * gamma[i, t]
      S[5, 7, i, t] <- phi[i, t] * s01[i, t]
      S[5, 8, i, t] <- 0
      S[5, 9, i, t] <- 0
      S[5, 10, i, t] <- 0
      S[5, 11, i, t] <- phi[i, t] * (1-s01[i, t]) * (1 - beta[i, t, 2])
      S[5, 12, i, t] <- 1-phi[i, t]
      
      S[6, 1, i, t] <- 0
      S[6, 2, i, t] <- 0
      S[6, 3, i, t] <- 0
      S[6, 4, i, t] <- 0
      S[6, 5, i, t] <- phi[i, t] * (1-s02[i, t])^2 * beta[i, t, 2] * (1 - gamma[i, t])
      S[6, 6, i, t] <- phi[i, t] * (1-s02[i, t])^2 * beta[i, t, 2] * gamma[i, t]
      S[6, 7, i, t] <- phi[i, t] * 2 * s02[i, t] * (1-s02[i, t])
      S[6, 8, i, t] <- phi[i, t] * s02[i, t]^2
      S[6, 9, i, t] <- 0
      S[6, 10, i, t] <- 0
      S[6, 11, i, t] <- phi[i, t] * (1-s02[i, t])^2 * (1 - beta[i, t, 2])
      S[6, 12, i, t] <- 1-phi[i, t]
      
      S[7, 1, i, t] <- 0
      S[7, 2, i, t] <- 0
      S[7, 3, i, t] <- 0
      S[7, 4, i, t] <- 0
      S[7, 5, i, t] <- phi[i, t] * (1-s1[i, t]) * beta[i, t, 3] * (1 - gamma[i, t])
      S[7, 6, i, t] <- phi[i, t] * (1-s1[i, t]) * beta[i, t, 3] * gamma[i, t]
      S[7, 7, i, t] <- 0
      S[7, 8, i, t] <- 0
      S[7, 9, i, t] <- phi[i, t] * s1[i, t]
      S[7, 10, i, t] <- 0
      S[7, 11, i, t] <- phi[i, t] * (1-s1[i, t]) * (1 - beta[i, t, 3])
      S[7, 12, i, t] <- 1-phi[i, t]
      
      S[8, 1, i, t] <- 0
      S[8, 2, i, t] <- 0
      S[8, 3, i, t] <- 0
      S[8, 4, i, t] <- 0
      S[8, 5, i, t] <- phi[i, t] * (1-s1[i, t])^2 * beta[i, t, 3] * (1 - gamma[i, t])
      S[8, 6, i, t] <- phi[i, t] * (1-s1[i, t])^2 * beta[i, t, 3] * gamma[i, t]
      S[8, 7, i, t] <- 0
      S[8, 8, i, t] <- 0
      S[8, 9, i, t] <- phi[i, t] * 2 * s1[i, t] * (1-s1[i, t])
      S[8, 10, i, t] <- phi[i, t] * s1[i, t]^2
      S[8, 11, i, t] <- phi[i, t] * (1-s1[i, t])^2 * (1 - beta[i, t, 3])
      S[8, 12, i, t] <- 1-phi[i, t]
      
      S[9, 1, i, t] <- 0
      S[9, 2, i, t] <- 0
      S[9, 3, i, t] <- 0
      S[9, 4, i, t] <- 0
      S[9, 5, i, t] <- phi[i, t] * beta[i, t, 4] * (1 - gamma[i, t])
      S[9, 6, i, t] <- phi[i, t] * beta[i, t, 4] * gamma[i, t]
      S[9, 7, i, t] <- 0
      S[9, 8, i, t] <- 0
      S[9, 9, i, t] <- 0
      S[9, 10, i, t] <- 0
      S[9, 11, i, t] <- phi[i, t] * (1 - beta[i, t, 4])
      S[9, 12, i, t] <- 1-phi[i, t]
      
      S[10, 1, i, t] <- 0
      S[10, 2, i, t] <- 0
      S[10, 3, i, t] <- 0
      S[10, 4, i, t] <- 0
      S[10, 5, i, t] <- phi[i, t] * beta[i, t, 4] * (1 - gamma[i, t])
      S[10, 6, i, t] <- phi[i, t] * beta[i, t, 4] * gamma[i, t]
      S[10, 7, i, t] <- 0
      S[10, 8, i, t] <- 0
      S[10, 9, i, t] <- 0
      S[10, 10, i, t] <- 0
      S[10, 11, i, t] <- phi[i, t] * (1 - beta[i, t, 4])
      S[10, 12, i, t] <- 1-phi[i, t]
      
      S[11, 1, i, t] <- 0
      S[11, 2, i, t] <- 0
      S[11, 3, i, t] <- 0
      S[11, 4, i, t] <- 0
      S[11, 5, i, t] <- phi[i, t] * beta[i, t, 1] * (1 - gamma[i, t])
      S[11, 6, i, t] <- phi[i, t] * beta[i, t, 1] * gamma[i, t]
      S[11, 7, i, t] <- 0
      S[11, 8, i, t] <- 0
      S[11, 9, i, t] <- 0
      S[11, 10, i, t] <- 0
      S[11, 11, i, t] <- phi[i, t] * (1 - beta[i, t, 1])
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
      
      # ~~~ d. E ----------------------------------
      
      E[1, 1, i, t] <- p[i, t+1, 1]
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
      E[1, 12, i, t] <- 1-p[i, t+1, 1]
      
      E[2, 1, i, t] <- 0
      E[2, 2, i, t] <- p[i, t+1, 1]
      E[2, 3, i, t] <- 0
      E[2, 4, i, t] <- 0
      E[2, 5, i, t] <- 0
      E[2, 6, i, t] <- 0
      E[2, 7, i, t] <- 0
      E[2, 8, i, t] <- 0
      E[2, 9, i, t] <- 0
      E[2, 10, i, t] <- 0
      E[2, 11, i, t] <- 0
      E[2, 12, i, t] <- 1 - p[i, t+1, 1]
      
      E[3, 1, i, t] <- 0
      E[3, 2, i, t] <- 0
      E[3, 3, i, t] <- p[i, t+1, 1]
      E[3, 4, i, t] <- 0
      E[3, 5, i, t] <- 0
      E[3, 6, i, t] <- 0
      E[3, 7, i, t] <- 0
      E[3, 8, i, t] <- 0
      E[3, 9, i, t] <- 0
      E[3, 10, i, t] <- 0
      E[3, 11, i, t] <- 0
      E[3, 12, i, t] <- 1 - p[i, t+1, 1]
      
      E[4, 1, i, t] <- 0
      E[4, 2, i, t] <- 0
      E[4, 3, i, t] <- 0
      E[4, 4, i, t] <- p[i, t+1, 1]
      E[4, 5, i, t] <- 0
      E[4, 6, i, t] <- 0
      E[4, 7, i, t] <- 0
      E[4, 8, i, t] <- 0
      E[4, 9, i, t] <- 0
      E[4, 10, i, t] <- 0
      E[4, 11, i, t] <- 0
      E[4, 12, i, t] <- 1 - p[i, t+1, 1]
      
      E[5, 1, i, t] <- 0
      E[5, 2, i, t] <- 0
      E[5, 3, i, t] <- 0
      E[5, 4, i, t] <- 0
      E[5, 5, i, t] <- p[i, t+1, 2]
      E[5, 6, i, t] <- 0
      E[5, 7, i, t] <- 0
      E[5, 8, i, t] <- 0
      E[5, 9, i, t] <- 0
      E[5, 10, i, t] <- 0
      E[5, 11, i, t] <- 0
      E[5, 12, i, t] <- 1 - p[i, t+1, 2]
      
      E[6, 1, i, t] <- 0
      E[6, 2, i, t] <- 0
      E[6, 3, i, t] <- 0
      E[6, 4, i, t] <- 0
      E[6, 5, i, t] <- 0
      E[6, 6, i, t] <- p[i, t+1, 2]
      E[6, 7, i, t] <- 0
      E[6, 8, i, t] <- 0
      E[6, 9, i, t] <- 0
      E[6, 10, i, t] <- 0
      E[6, 11, i, t] <- 0
      E[6, 12, i, t] <- 1 - p[i, t+1, 2]
      
      E[7, 1, i, t] <- 0
      E[7, 2, i, t] <- 0
      E[7, 3, i, t] <- 0
      E[7, 4, i, t] <- 0
      E[7, 5, i, t] <- 0
      E[7, 6, i, t] <- 0
      E[7, 7, i, t] <- p[i, t+1, 1]
      E[7, 8, i, t] <- 0
      E[7, 9, i, t] <- 0
      E[7, 10, i, t] <- 0
      E[7, 11, i, t] <- 0
      E[7, 12, i, t] <- 1 - p[i, t+1, 1]
      
      E[8, 1, i, t] <- 0
      E[8, 2, i, t] <- 0
      E[8, 3, i, t] <- 0
      E[8, 4, i, t] <- 0
      E[8, 5, i, t] <- 0
      E[8, 6, i, t] <- 0
      E[8, 7, i, t] <- 0
      E[8, 8, i, t] <- p[i, t+1, 1]
      E[8, 9, i, t] <- 0
      E[8, 10, i, t] <- 0
      E[8, 11, i, t] <- 0
      E[8, 12, i, t] <- 1 - p[i, t+1, 1]
      
      E[9, 1, i, t] <- 0
      E[9, 2, i, t] <- 0
      E[9, 3, i, t] <- 0
      E[9, 4, i, t] <- 0
      E[9, 5, i, t] <- 0
      E[9, 6, i, t] <- 0
      E[9, 7, i, t] <- 0
      E[9, 8, i, t] <- 0
      E[9, 9, i, t] <- (1 - ALPHA[i, t+1]) * p[i, t+1, 1]
      E[9, 10, i, t] <- 0
      E[9, 11, i, t] <- ALPHA[i, t+1] * p[i, t+1, 1]
      E[9, 12, i, t] <- 1 - p[i, t+1, 1]
      
      E[10, 1, i, t] <- 0
      E[10, 2, i, t] <- 0
      E[10, 3, i, t] <- 0
      E[10, 4, i, t] <- 0
      E[10, 5, i, t] <- 0
      E[10, 6, i, t] <- 0
      E[10, 7, i, t] <- 0
      E[10, 8, i, t] <- 0
      E[10, 9, i, t] <- 2 * (1 - ALPHA[i, t+1]) * ALPHA[i, t+1] * p[i, t+1, 1]
      E[10, 10, i, t] <- (1 - ALPHA[i, t+1])^2 * p[i, t+1, 1]
      E[10, 11, i, t] <- ALPHA[i, t+1]^2 * p[i, t+1, 1]
      E[10, 12, i, t] <- 1 - p[i, t+1, 1]
      
      E[11, 1, i, t] <- 0
      E[11, 2, i, t] <- 0
      E[11, 3, i, t] <- 0
      E[11, 4, i, t] <- 0
      E[11, 5, i, t] <- 0
      E[11, 6, i, t] <- 0
      E[11, 7, i, t] <- 0
      E[11, 8, i, t] <- 0
      E[11, 9, i, t] <- 0
      E[11, 10, i, t] <- 0
      E[11, 11, i, t] <- p[i, t+1, 1]
      E[11, 12, i, t] <- 1 - p[i, t+1, 1]
      
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
      E[12, 11, i, t] <- 0
      E[12, 12, i, t] <- 1
    }
  }
  
  # ++++++++++++++++++++++++++++++++ priors ++++++++++++++++++++++++++++++++++++
  
  #-------- State process --------#
  for (k in 1:3) {
    beta_phi[k] ~ dnorm(0, sd = 1.5)  # survival
  }
  beta_phi[4] ~ dnorm(0, sd = 3)  # survival
  beta_phi[5] ~ dnorm(0, sd = 3)  # survival
  
  for (k in 1:3) {
    beta_s0[k] ~ dnorm(0, sd = 1.5)   # cub survival
  }
  for (k in 4:6) {
    beta_s0[k] ~ dnorm(0, sd = 3)   # cub survival (additive)
  }
  
  for (k in 1:3) {
    beta_s1[k] ~ dnorm(0, sd = 1.5)   # yearling survival
  }
  beta_s1[4] ~ dnorm(0, sd = 3)
  
  for (k in 1:4) {
    beta_beta[k] ~ dnorm(0, sd = 1.5)  # breeding probability (intercepts)
  }
  for (k in 5:9) {
    beta_beta[k] ~ dnorm(0, sd = 3)  # breeding probability (additive effects)
  }
  
  for (k in 1:3) {
    beta_gamma[k] ~ dnorm(0, sd = 1.5)   # litter size
  }
  beta_gamma[4] ~ dnorm(0, sd = 3) # litter size (effect of date capture on young females)
  beta_gamma[5] ~ dnorm(0, sd = 3)
  beta_gamma[6] ~ dnorm(0, sd = 3)
  
  # Constraint
  constraint[1] ~ dconstraint(beta_s1[1] > beta_s0[1])
  constraint[2] ~ dconstraint(beta_s1[2] > beta_s0[2])
  constraint[3] ~ dconstraint(beta_s1[3] > beta_s0[3])
  
  #----- observation process -----#
  for (k in 1:4) {
    beta_p[k] ~ dnorm(0, sd = 1.5)
  }
  # Yearly random effects on recapture
  for (t in 1:n_intervals) {
    omega[t] ~ dnorm(0, sd = sigma_omega)
  }
  sigma_omega ~ dunif(0, 10)
  
  
  # +++++++++++++++++++++++++++ model & likelihood +++++++++++++++++++++++++++++
  
  for(i in 1:n_ind){
    
    # First capture
    for (k in 1:n_states) {
      zeta[i, f[i], k] <- S0[k] * E0[k, CH[i, f[i]], i]
    }
    
    # Subsequent captures
    for(t in f[i] : (l[i]-1) ){
      for (k in 1:n_states) {
        zeta[i, (t+1), k] <- inprod(zeta[i, t, 1:n_states], 
                                    S[1:n_states, k, i, t]) * 
          E[k, CH[i, t+1], i, t]
      }
    }
    
    lik[i] <- sum(zeta[i, l[i], 1:n_states])
    ones[i] ~ dbin(prob = lik[i], size = freq[i])
  }
})



# 5. Run the model -------------------------------------------------------------

# Need to run that twice (in parallel) to get two chains. ~18h per chain on an 
# HPC cluster


start <- Sys.time()  # 
# Create R model
inits_values <- inits()
CR_model_R <- nimbleModel(code = CR_model,
                      constants = my.constants,
                      data = dat,
                      inits = inits_values,
                      calculate = FALSE)

# Compile model (in C++)
CR_model_C <- compileNimble(CR_model_R,
                        showCompilerOutput = FALSE)

# Configure of MCMC
MCMC_conf <- configureMCMC(CR_model_R,
                           monitors = params)


# Compile MCMC
MCMC_R <- buildMCMC(MCMC_conf)
MCMC_C <- compileNimble(MCMC_R, project = CR_model_R)

fit_without_remote_1 <- runMCMC(mcmc = MCMC_C, 
                                nchains = 1,
                                niter = 18000,
                                nburnin = 5000,
                                thin = 5,
                                inits = inits())

save(fit_without_remote_1,
     file = "03_outputs/fit_without_remote_1.RData")
end <- Sys.time() ; end - start

