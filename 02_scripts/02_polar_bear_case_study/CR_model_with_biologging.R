library(tidyverse)
library(nimble)


## 1. Load and process data ----------------------------------------------------

data_CR_events <-  read_csv("01_inputs/CR_GPS_GLS_events.csv", show_col_types = F) %>%
  mutate(date = ymd(paste0(year, "-01-01")) + day_number - 1,
         date_ref = ymd(paste0(year-1, "-01-01")),
         doy_capture = as.numeric(date - date_ref + 1)) 

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


n_states <- 14
not_captured <- 19
physical_capture <- c(1:11, 16:18)
GPS_GLS_capture <- 12:18

CH  <- matrix(data = NA, 
              nrow = length(unique(data_CR_events$ID_NR)), 
              ncol = length(min(data_CR_events$year):max(data_CR_events$year)),
              dimnames = list(unique(data_CR_events$ID_NR),
                              min(data_CR_events$year):max(data_CR_events$year)))

for (k in 1:nrow(data_CR_events)) {
  i <- which(rownames(CH) == data_CR_events$ID_NR[k])
  t <- which(colnames(CH) == as.character(data_CR_events$year[k]))
  
  CH[i, t] <- data_CR_events$event[k]
}
CH[is.na(CH)] <- not_captured 


# Score whether GPS/GLS capture ++++++++++++++++++++++++++++++++++++++++++++++++
KHI <- matrix(data = NA, 
              nrow = length(unique(data_CR_events$ID_NR)), 
              ncol = length(min(data_CR_events$year):max(data_CR_events$year)),
              dimnames = list(unique(data_CR_events$ID_NR),
                              min(data_CR_events$year):max(data_CR_events$year))) 
for (i in 1:dim(CH)[1]) {
  for (j in 1:dim(CH)[2]) {
    KHI[i, j] <- ifelse(CH[i, j] %in% GPS_GLS_capture, 1, 0)
  }
}


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
l[l > dim(CH)[2]] <- dim(CH)[2]

# Physical capture at first capture?
P0 <- ifelse(f_event %in% physical_capture, 1, 0)


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
  DOY[i, t-1] <- data_CR_events$doy_capture[k]
}



# ALPHA ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# The probability that a female who has successfully raised cubs to age 2.5 is 
# alone at capture (i.e., 2-year-old cubs have parted from the mother) depending
# on the date of capture is computed using a GLM on all spring captures of 
# two-year-old cubs. To do so, we create a dummy variable:
# 0 = two-year-old cub is still with mother alone
# 1 = two-year-old cub is alone

two_yr_old_cubs <- read_csv("01_inputs/two_year_old_cubs.csv", show_col_types = F) 

logistic <- glm(parted ~ day_number, data = two_yr_old_cubs, family = "binomial")

ALPHA <- matrix(data = NA, 
                nrow = length(unique(data_CR_events$ID_NR)), 
                ncol = length(min(data_CR_events$year):max(data_CR_events$year)),
                dimnames = list(unique(data_CR_events$ID_NR),
                                min(data_CR_events$year):max(data_CR_events$year)))
for (k in 1:nrow(data_CR_events)) {
  i <- which(rownames(CH) == data_CR_events$ID_NR[k])
  t <- which(colnames(CH) == as.character(data_CR_events$year[k]))
  
  ALPHA[i, t] <- plogis(logistic$coefficients[1] + logistic$coefficients[2]*data_CR_events$day_number[k]) 
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
# The indexing is a bit confusing but by trial and error, I've reach this indexing.



# Years 1989, 90, and 91 with no REcaptures of A0 females +++++++++++++++++++++
year_post_1992 <- ifelse((as.numeric(colnames(CH)[1:( dim(CH)[2]-1)])+1) %in% 1989:1991, 0, 1) 
# The indexing is a bit confusing but by trial and error, I've reach this indexing.


# Ecotype ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SpPel <- SpUnk <- NULL
for (i in 1:dim(CH)[1]) {
  SpPel[i] <- ifelse(data_CR_events$ecotype[data_CR_events$ID_NR == rownames(CH)[i]][1] %in% "O", 
                     1, 0)
  SpUnk[i] <- ifelse(is.na(data_CR_events$ecotype[data_CR_events$ID_NR == rownames(CH)[i]][1]), 
                     1, 0)
}



# Sea ice availability
sea_ice_data <- read_csv("01_inputs/sea_ice_metrics.csv", show_col_types = F)
sea_ice <- NULL
for (t in 1:(ncol(CH)-1)) {
  sea_ice[t] <- sea_ice_data$ice_free_days[which(sea_ice_data$year == as.numeric(colnames(CH)[t]))]   # here I remove the '-1'
}
sea_ice_s <- as.vector(scale(sea_ice))


# 2. Compute initial state -----------------------------------------------------

# With the marginalized likelihood, we don't need initial values for the latent 
# state. However, to compute the proportion of each state at first capture, we can
# use the proposed initial value for the first capture of each individual
z_inits <- NULL
for (i in 1:dim(CH)[1]) {
  obs_i <- CH[i, ]
  z_inits_i <- rep(13, times = dim(CH)[2])
  
  # First, the "constrained" observations
  J2 <- which(obs_i == 1)
  J3 <- which(obs_i == 2)
  SA4 <- which(obs_i == 3)
  
  SA5 <- which(obs_i == 4)
  NDSA5 <- which(obs_i == 12)
  DSA5 <- which(obs_i == 13)
  FDSA5 <- which(obs_i == 16)
  
  A01 <- which(obs_i == 5)
  A02 <- which(obs_i == 6)
  
  A11 <- which(obs_i == 7)
  A12 <- which(obs_i == 8)
  A21 <- which(obs_i == 9)
  A22 <- which(obs_i == 10)
  
  A <- which(obs_i == 11)
  NDA <- which(obs_i == 14)
  DA <- which(obs_i == 15)
  LNDA <- which(obs_i == 17)
  FDA <- which(obs_i == 18)
  
  # If captured as J2
  if (length(J2) > 0) {
    z_inits_i[J2] <- 1        # State J2 at t
    if (J2 < l[i]) {
      z_inits_i[J2 + 1] <- 2     # State J3 at t+1
      if ((J2 + 1) < l[i]) {
        z_inits_i[J2 + 2] <- 3     # State SA4 at t+2
        if ((J2 + 2) < l[i]) {
          z_inits_i[J2 + 3] <- 4     # State NDSA5 at t+3
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
        z_inits_i[J3 + 2] <- 4     # State NDSA5 at t+3
      }
    }
  }
  # If captured as SA4
  if (length(SA4) > 0) {
    z_inits_i[SA4] <- 3     # state SA4 at t
    if (SA4 < l[i]) {
      z_inits_i[SA4 + 1] <- 4     # state NDSA5 at t+1
    }
  }
  # If captured as SA5
  if (length(SA5) > 0) {
    z_inits_i[SA5] <- sample(4:5, size = 1)   # state NDSA5 or FDSA5
  }
  # If captured as NDSA5
  if (length(NDSA5) > 0) {
    z_inits_i[NDSA5] <- 4   # state NDSA5
  }
  # If captured as DSA5
  if (length(DSA5) > 0) {
    z_inits_i[DSA5] <- sample(5:7, size = 1)  # state FDSA5, A01 or A02
  }
  # If captured as FDSA5
  if (length(FDSA5) > 0) {
    z_inits_i[FDSA5] <- 5  # state FDSA5
  }
  # If captured as NDA                           # I must fill this cells of the matrix
  if (length(NDA) > 0) {                         # now so that the values may be overwritten
    z_inits_i[NDA] <- 12  # state LNDA           # if the females were actually accompanied by
  }                                              # yearlings  
  # If captured as DA
  if (length(DA) > 0) {
    z_inits_i[DA] <- sample(c(6, 7, 13), size = 1)  # state A01, A02 or FDA
  }
  # If captured as A                             # Same idea for the cells corresponding
  if (length(A) > 0) {                           # to the event A
    z_inits_i[A] <- sample(12:13, size = 1)  # state LNDA or FDA
  }
  # If captured as A21
  if (length(A21) == 1) {    # If captured only once as A21 
    z_inits_i[A21] <- 10          # state AS1 at t 
    if (A21 > f[i]) {
      z_inits_i[A21 - 1] <- 8       # state A11 at t-1
      if (A21 > (f[i] + 1)) {
        z_inits_i[A21 - 2] <- 7       # state A02 at t-2
      }
    }  
  } else {
    if (length(A21) >= 2) {  # If captured more than once as A21
      z_inits_i[A21[1]] <- 10          # state AS1 at t 
      if (A21[1] > f[i]) {
        z_inits_i[A21[1] - 1] <- 8       # state A11 at t-1
        if (A21[1] > (f[i] + 1)) {
          z_inits_i[A21[1] - 2] <- 7        # state A02 at t-2
        }
      } 
      for (t in 2:length(A21)) {
        z_inits_i[A21[t]] <- 10        # state AS1 at t
        z_inits_i[A21[t] - 1] <- 8    # state A11 at t-1
        z_inits_i[A21[t] - 2] <- 7     # state A02 at t-2
      }
    }
  }
  # If captured as A22
  if (length(A22) == 1) {    # If captured only once as A22 
    z_inits_i[A22] <- 11          # state AS2 at t 
    if (A22 > f[i]) {
      z_inits_i[A22 - 1] <- 9       # state A12 at t-1
      if (A22 > (f[i] + 1)) {
        z_inits_i[A22 - 2] <- 7       # state A02 at t-2
      }
    }  
  } else {
    if (length(A22) >= 2) {  # If captured more than once as A22
      z_inits_i[A22[1]] <- 11          # state AS2 at t 
      if (A22[1] > f[i]) {
        z_inits_i[A22[1] - 1] <- 9       # state A12 at t-1
        if (A22[1] > (f[i] + 1)) {
          z_inits_i[A22[1] - 2] <- 7        # state A02 at t-2
        }
      } 
      for (t in 2:length(A22)) {
        z_inits_i[A22[t]] <- 11        # state AS2 at t
        z_inits_i[A22[t] - 1] <- 9    # state A12 at t-1
        z_inits_i[A22[t] - 2] <- 7     # state A02 at t-2
      }
    }
  }
  # If captured as A11
  if (length(A11) == 1) {      # If captured only once as A11 
    z_inits_i[A11] <- 8            # state A11 at t
    if (A11 > f[i]) {
      z_inits_i[A11 - 1] <- 7       # state A02 at t-1 
    }
  } else {
    if (length(A11) >= 2) {    # If captured more than once as A11 
      z_inits_i[A11[1]] <- 8        # state A11 at t
      if (A11[1] > f[i]) {
        z_inits_i[A11[1] - 1] <- 7     # state A02 at t-1 
      }
      for (t in 2:length(A11)) {
        z_inits_i[A11[t]] <- 8       # state A11 at t
        z_inits_i[A11[t] - 1] <- 7       # state A02 at t-1 
      }
    }
  }
  # If captured as A12
  if (length(A12) == 1) {      # If captured only once as A12
    z_inits_i[A12] <- 9          # state A12 at t 
    if (A12 > f[i]) {
      z_inits_i[A12 - 1] <- 7       # state A02 at t-1 
    }
  } else {
    if (length(A12) >= 2) {    # If captured more than once as A12
      z_inits_i[A12[1]] <- 9        # state A12 at t
      if (A12[1] > f[i]) {
        z_inits_i[A12[1] - 1] <- 7     # state A02 at t-1 
      }
      for (t in 2:length(A12)) {
        z_inits_i[A12[t]] <- 9       # state A12 at t
        z_inits_i[A12[t] - 1] <- 7       # state A02 at t-1 
      }
    }
  }
  # If captured as A01
  for (t in 1:length(A01)) {
    z_inits_i[A01[t]] <- 6      # state A01
  }
  # If captured as A02
  for (t in 1:length(A02)) {
    z_inits_i[A02[t]] <- 7      # state A02
  }
  # If captured as LNDA
  if (length(LNDA) > 0) {
    z_inits_i[LNDA] <- 12  # state LNDA
  }
  # If captured as FDA
  if (length(FDA) > 0) {
    z_inits_i[FDA] <- 13  # state FDA
  }
  
  # fill remaining cells
  if (f[i] > 1) {
    z_inits_i[1:(f[i]-1)] <- NA
  }
  # if (l[i] < dim(CH)[2]) {
  #   z_inits_i[(l[i]+1):dim(CH)[2]] <- NA
  # }  
  
  z_inits <- rbind(z_inits, z_inits_i)
}
rownames(z_inits) <- rownames(CH) # 1:dim(z_inits)[1]
colnames(z_inits) <- colnames(CH)
table(z_inits)


# Proportion of each state at first capture ++++++++++++++++++++++++++
# (must be derived from the data)
f_state <- NULL
for (i in 1:dim(CH)[1]) {
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
                     freq = rep(1, times = nrow(CH)),
                     f = f,
                     l = l,
                     S0 = prop,
                     n_CH = nrow(CH),
                     n_occasions = ncol(CH),
                     n_states = n_states,
                     year_not_2020 = year_not_2020,
                     year_post_1992 = year_post_1992,
                     P0 = P0,
                     KHI = KHI,
                     
                     AGE_2_4 = AGE_2_4,
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
                     
                     DOY = DOY,
                     ALPHA = ALPHA,
                     SpUnk = SpUnk,
                     SpPel = SpPel,
                     sea_ice_s = sea_ice_s,
                     time_s = as.vector(scale(as.numeric(colnames(CH)))))

inits <- function() list(beta_phi = c(2, 3.18, 1.6, -0.35, 0) + rnorm(n = 5, mean = 0, sd = 0.1),
                         
                         
                         beta_s0  = rnorm(n = 6, mean = 0, sd = 0.1), 
                         beta_s1  = c(1, 1, 1, 0) + rnorm(n = 4, mean = 0, sd = 0.1), 
                         
                         beta_eta  = c(-1.45, 0.6, 0.7, 0.60, -0.6, -2, -0.8, 0.4, 0) + rnorm(n = 9, mean = 0, sd = 0.1),
                         
                         beta_beta = c(1, 1, 1, 0, 0) + rnorm(n = 5, mean = 0, sd = 0.5),
                         beta_gamma = c(0.75, 0.75, 0.75, 0, 0, 0) + rnorm(n = 6, mean = 0, sd = 0.5),
                         
                         beta_p = c(-1.35, -2.75, -1.6, -0.1) + rnorm(n = 4, mean = 0, sd = 0.05),
                         omega = rnorm(dim(CH)[2]-1, mean = 0, sd = 1),
                         sigma_omega = runif(1, min = 0, max = 5))

# Parameters monitored 
params <- c("beta_phi",
            "beta_s0", "beta_s1",
            "beta_eta", "beta_beta", "beta_gamma", 
            "beta_p",  "omega", "sigma_omega")



# 4. Build the model -----------------------------------------------------------

CR_model_biologging <- nimbleCode({
  
  # +++++++++++++++++++++++++++++++ matrices +++++++++++++++++++++++++++++++++++
  
  for (i in 1:n_CH) {
    
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
    E0[1, 13, i] <- 0
    E0[1, 14, i] <- 0
    E0[1, 15, i] <- 0
    E0[1, 16, i] <- 0
    E0[1, 17, i] <- 0
    E0[1, 18, i] <- 0
    E0[1, 19, i] <- 0
    
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
    E0[2, 13, i] <- 0
    E0[2, 14, i] <- 0
    E0[2, 15, i] <- 0
    E0[2, 16, i] <- 0
    E0[2, 17, i] <- 0
    E0[2, 18, i] <- 0
    E0[2, 19, i] <- 0
    
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
    E0[3, 13, i] <- 0
    E0[3, 14, i] <- 0
    E0[3, 15, i] <- 0
    E0[3, 16, i] <- 0
    E0[3, 17, i] <- 0
    E0[3, 18, i] <- 0
    E0[3, 19, i] <- 0
    
    E0[4, 1, i] <- 0
    E0[4, 2, i] <- 0
    E0[4, 3, i] <- 0
    E0[4, 4, i] <- ( 1-KHI[i, f[i]] ) * P0[i]
    E0[4, 5, i] <- 0
    E0[4, 6, i] <- 0
    E0[4, 7, i] <- 0
    E0[4, 8, i] <- 0
    E0[4, 9, i] <- 0
    E0[4, 10, i] <- 0
    E0[4, 11, i] <- 0
    E0[4, 12, i] <- KHI[i, f[i]]
    E0[4, 13, i] <- 0
    E0[4, 14, i] <- 0
    E0[4, 15, i] <- 0
    E0[4, 16, i] <- 0
    E0[4, 17, i] <- 0
    E0[4, 18, i] <- 0
    E0[4, 19, i] <- 0
    
    E0[5, 1, i] <- 0
    E0[5, 2, i] <- 0
    E0[5, 3, i] <- 0
    E0[5, 4, i] <- ( 1-KHI[i, f[i]] ) * P0[i]
    E0[5, 5, i] <- 0
    E0[5, 6, i] <- 0
    E0[5, 7, i] <- 0
    E0[5, 8, i] <- 0
    E0[5, 9, i] <- 0
    E0[5, 10, i] <- 0
    E0[5, 11, i] <- 0
    E0[5, 12, i] <- 0
    E0[5, 13, i] <- KHI[i, f[i]] * ( 1-P0[i] )
    E0[5, 14, i] <- 0
    E0[5, 15, i] <- 0
    E0[5, 16, i] <- KHI[i, f[i]] * P0[i]
    E0[5, 17, i] <- 0
    E0[5, 18, i] <- 0
    E0[5, 19, i] <- 0
    
    E0[6, 1, i] <- 0
    E0[6, 2, i] <- 0
    E0[6, 3, i] <- 0
    E0[6, 4, i] <- 0
    E0[6, 5, i] <- P0[i]
    E0[6, 6, i] <- 0
    E0[6, 7, i] <- 0
    E0[6, 8, i] <- 0
    E0[6, 9, i] <- 0
    E0[6, 10, i] <- 0
    E0[6, 11, i] <- 0
    E0[6, 12, i] <- 0
    E0[6, 13, i] <- KHI[i, f[i]] * ( 1-P0[i] ) * AGE_5[i, f[i]]
    E0[6, 14, i] <- 0
    E0[6, 15, i] <- KHI[i, f[i]] * ( 1-P0[i] ) * ( 1-AGE_5[i, f[i]] )
    E0[6, 16, i] <- 0
    E0[6, 17, i] <- 0
    E0[6, 18, i] <- 0
    E0[6, 19, i] <- 0
    
    E0[7, 1, i] <- 0
    E0[7, 2, i] <- 0
    E0[7, 3, i] <- 0
    E0[7, 4, i] <- 0
    E0[7, 5, i] <- 0
    E0[7, 6, i] <- P0[i]
    E0[7, 7, i] <- 0
    E0[7, 8, i] <- 0
    E0[7, 9, i] <- 0
    E0[7, 10, i] <- 0
    E0[7, 11, i] <- 0
    E0[7, 12, i] <- 0
    E0[7, 13, i] <- KHI[i, f[i]] * ( 1-P0[i] ) * AGE_5[i, f[i]]
    E0[7, 14, i] <- 0
    E0[7, 15, i] <- KHI[i, f[i]] * ( 1-P0[i] ) * ( 1-AGE_5[i, f[i]] )
    E0[7, 16, i] <- 0
    E0[7, 17, i] <- 0
    E0[7, 18, i] <- 0
    E0[7, 19, i] <- 0
    
    E0[8, 1, i] <- 0
    E0[8, 2, i] <- 0
    E0[8, 3, i] <- 0
    E0[8, 4, i] <- 0
    E0[8, 5, i] <- 0
    E0[8, 6, i] <- 0
    E0[8, 7, i] <- P0[i]
    E0[8, 8, i] <- 0
    E0[8, 9, i] <- 0
    E0[8, 10, i] <- 0
    E0[8, 11, i] <- 0
    E0[8, 12, i] <- 0
    E0[8, 13, i] <- 0
    E0[8, 14, i] <- KHI[i, f[i]] * ( 1-P0[i] )
    E0[8, 15, i] <- 0
    E0[8, 16, i] <- 0
    E0[8, 17, i] <- 0
    E0[8, 18, i] <- 0
    E0[8, 19, i] <- 0
    
    E0[9, 1, i] <- 0
    E0[9, 2, i] <- 0
    E0[9, 3, i] <- 0
    E0[9, 4, i] <- 0
    E0[9, 5, i] <- 0
    E0[9, 6, i] <- 0
    E0[9, 7, i] <- 0
    E0[9, 8, i] <- P0[i]
    E0[9, 9, i] <- 0
    E0[9, 10, i] <- 0
    E0[9, 11, i] <- 0
    E0[9, 12, i] <- 0
    E0[9, 13, i] <- 0
    E0[9, 14, i] <- KHI[i, f[i]] * ( 1-P0[i] )
    E0[9, 15, i] <- 0
    E0[9, 16, i] <- 0
    E0[9, 17, i] <- 0
    E0[9, 18, i] <- 0
    E0[9, 19, i] <- 0
    
    E0[10, 1, i] <- 0
    E0[10, 2, i] <- 0
    E0[10, 3, i] <- 0
    E0[10, 4, i] <- 0
    E0[10, 5, i] <- 0
    E0[10, 6, i] <- 0
    E0[10, 7, i] <- 0
    E0[10, 8, i] <- 0
    E0[10, 9, i] <- ( 1-ALPHA[i, f[i]] ) * P0[i]
    E0[10, 10, i] <- 0
    E0[10, 11, i] <- ALPHA[i, f[i]] * ( 1-KHI[i, f[i]] ) * P0[i]
    E0[10, 12, i] <- 0
    E0[10, 13, i] <- 0
    E0[10, 14, i] <- KHI[i, f[i]] * ( 1-P0[i] )
    E0[10, 15, i] <- 0
    E0[10, 16, i] <- 0
    E0[10, 17, i] <- ALPHA[i, f[i]] * KHI[i, f[i]] * P0[i]
    E0[10, 18, i] <- 0
    E0[10, 19, i] <- 0
    
    E0[11, 1, i] <- 0
    E0[11, 2, i] <- 0
    E0[11, 3, i] <- 0
    E0[11, 4, i] <- 0
    E0[11, 5, i] <- 0
    E0[11, 6, i] <- 0
    E0[11, 7, i] <- 0
    E0[11, 8, i] <- 0
    E0[11, 9, i] <- 2 * ( 1-ALPHA[i, f[i]] ) * ALPHA[i, f[i]] * P0[i]
    E0[11, 10, i] <- ( 1-ALPHA[i, f[i]] )^2 * P0[i]
    E0[11, 11, i] <- ALPHA[i, f[i]]^2 * (1-KHI[i, f[i]]) * P0[i] 
    E0[11, 12, i] <- 0
    E0[11, 13, i] <- 0
    E0[11, 14, i] <- KHI[i, f[i]] * ( 1-P0[i] )
    E0[11, 15, i] <- 0
    E0[11, 16, i] <- 0
    E0[11, 17, i] <- ALPHA[i, f[i]]^2 * KHI[i, f[i]] * P0[i] 
    E0[11, 18, i] <- 0
    E0[11, 19, i] <- 0
    
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
    E0[12, 11, i] <- 1-KHI[i, f[i]]
    E0[12, 12, i] <- 0
    E0[12, 13, i] <- 0
    E0[12, 14, i] <- KHI[i, f[i]] * ( 1-P0[i] )
    E0[12, 15, i] <- 0
    E0[12, 16, i] <- 0
    E0[12, 17, i] <- KHI[i, f[i]] * P0[i] 
    E0[12, 18, i] <- 0
    E0[12, 19, i] <- 0
    
    E0[13, 1, i] <- 0
    E0[13, 2, i] <- 0
    E0[13, 3, i] <- 0
    E0[13, 4, i] <- 0
    E0[13, 5, i] <- 0
    E0[13, 6, i] <- 0
    E0[13, 7, i] <- 0
    E0[13, 8, i] <- 0
    E0[13, 9, i] <- 0
    E0[13, 10, i] <- 0
    E0[13, 11, i] <- 1-KHI[i, f[i]]
    E0[13, 12, i] <- 0
    E0[13, 13, i] <- 0
    E0[13, 14, i] <- 0
    E0[13, 15, i] <- KHI[i, f[i]] * ( 1-P0[i] )
    E0[13, 16, i] <- 0
    E0[13, 17, i] <- 0
    E0[13, 18, i] <- KHI[i, f[i]] * P0[i] 
    E0[13, 19, i] <- 0
    
    E0[14, 1, i] <- 0
    E0[14, 2, i] <- 0
    E0[14, 3, i] <- 0
    E0[14, 4, i] <- 0
    E0[14, 5, i] <- 0
    E0[14, 6, i] <- 0
    E0[14, 7, i] <- 0
    E0[14, 8, i] <- 0
    E0[14, 9, i] <- 0
    E0[14, 10, i] <- 0
    E0[14, 11, i] <- 0
    E0[14, 12, i] <- 0
    E0[14, 13, i] <- 0
    E0[14, 14, i] <- 0
    E0[14, 15, i] <- 0
    E0[14, 16, i] <- 0
    E0[14, 17, i] <- 0
    E0[14, 18, i] <- 0
    E0[14, 19, i] <- 1
    
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
        beta_s0[4] * (SpUnk[i] + SpPel[i]) +
        beta_s0[5] +
        beta_s0[6] * sea_ice_s[t]
      
      logit(s02[i, t]) <- beta_s0[1] * AGE_YOUNG_1[i, t] +
        beta_s0[2] * AGE_MIDDLE_1[i, t] +
        beta_s0[3] * AGE_OLD_1[i, t] +
        beta_s0[4] * (SpUnk[i] + SpPel[i]) +
        beta_s0[6] * sea_ice_s[t]
      
      # Probability of yearling survival
      logit(s1[i, t]) <- beta_s1[1] * AGE_YOUNG_2[i, t] +
        beta_s1[2] * AGE_MIDDLE_2[i, t] +
        beta_s1[3] * AGE_OLD_2[i, t] +
        beta_s1[4] * sea_ice_s[t]
      
      # Denning probability
      logit(eta[i, t, 1]) <- beta_eta[1] * AGE_4[i, t] +   # Lone female
        beta_eta[2] * AGE_YOUNG[i, t] +
        beta_eta[3] * AGE_MIDDLE[i, t] +
        beta_eta[4] * AGE_OLD[i, t] +
        beta_eta[9] * sea_ice_s[t]
      
      logit(eta[i, t, 2]) <- beta_eta[2] * AGE_YOUNG[i, t] + # FDA
        beta_eta[3] * AGE_MIDDLE[i, t] +
        beta_eta[4] * AGE_OLD[i, t] +
        beta_eta[5] +
        beta_eta[9] * sea_ice_s[t]
      
      logit(eta[i, t, 3]) <- beta_eta[2] * AGE_YOUNG[i, t] + # A0-
        beta_eta[3] * AGE_MIDDLE[i, t] +
        beta_eta[4] * AGE_OLD[i, t] +
        beta_eta[6] +
        beta_eta[9] * sea_ice_s[t]
      
      logit(eta[i, t, 4]) <- beta_eta[2] * AGE_YOUNG[i, t] + # A1-
        beta_eta[3] * AGE_MIDDLE[i, t] +
        beta_eta[4] * AGE_OLD[i, t] +
        beta_eta[7] +
        beta_eta[9] * sea_ice_s[t]
      
      logit(eta[i, t, 5]) <- beta_eta[2] * AGE_YOUNG[i, t] + # AS
        beta_eta[3] * AGE_MIDDLE[i, t] +
        beta_eta[4] * AGE_OLD[i, t] +
        beta_eta[8] +
        beta_eta[9] * sea_ice_s[t]
      
      # Early litter survival
      logit(beta[i, t]) <- beta_beta[1] * (AGE_4[i, t] + AGE_YOUNG[i, t]) +
        beta_beta[2] * AGE_MIDDLE[i, t] +
        beta_beta[3] * AGE_OLD[i, t] +
        beta_beta[4] * DOY[i, t] +
        beta_beta[5] * sea_ice_s[t]
      
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
      S[3, 4, i, t] <- phi[i, t] * (1- eta[i, t, 1])
      S[3, 5, i, t] <- phi[i, t] * eta[i, t, 1] * (1 - beta[i, t])
      S[3, 6, i, t] <- phi[i, t] * eta[i, t, 1] * beta[i, t] * (1 - gamma[i, t])
      S[3, 7, i, t] <- phi[i, t] * eta[i, t, 1] * beta[i, t] * gamma[i, t]
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
      S[4, 6, i, t] <- phi[i, t] * eta[i, t, 1] * beta[i, t] * (1 - gamma[i, t])
      S[4, 7, i, t] <- phi[i, t] * eta[i, t, 1] * beta[i, t] * gamma[i, t]
      S[4, 8, i, t] <- 0
      S[4, 9, i, t] <- 0
      S[4, 10, i, t] <- 0
      S[4, 11, i, t] <- 0
      S[4, 12, i, t] <- phi[i, t] * (1 - eta[i, t, 1])
      S[4, 13, i, t] <- phi[i, t] * eta[i, t, 1] * (1 - beta[i, t])
      S[4, 14, i, t] <- 1-phi[i, t] 
      
      S[5, 1, i, t] <- 0
      S[5, 2, i, t] <- 0
      S[5, 3, i, t] <- 0
      S[5, 4, i, t] <- 0
      S[5, 5, i, t] <- 0
      S[5, 6, i, t] <- phi[i, t] * eta[i, t, 2] * beta[i, t] * (1 - gamma[i, t])
      S[5, 7, i, t] <- phi[i, t] * eta[i, t, 2] * beta[i, t] * gamma[i, t]
      S[5, 8, i, t] <- 0
      S[5, 9, i, t] <- 0
      S[5, 10, i, t] <- 0
      S[5, 11, i, t] <- 0
      S[5, 12, i, t] <- phi[i, t] * (1 - eta[i, t, 2])
      S[5, 13, i, t] <- phi[i, t] * eta[i, t, 2] * (1 - beta[i, t])
      S[5, 14, i, t] <- 1-phi[i, t] 
      
      S[6, 1, i, t] <- 0
      S[6, 2, i, t] <- 0
      S[6, 3, i, t] <- 0
      S[6, 4, i, t] <- 0
      S[6, 5, i, t] <- 0
      S[6, 6, i, t] <- phi[i, t] * eta[i, t, 3] * (1-s01[i, t]) * beta[i, t] * (1 - gamma[i, t])
      S[6, 7, i, t] <- phi[i, t] * eta[i, t, 3] * (1-s01[i, t]) * beta[i, t] * gamma[i, t]
      S[6, 8, i, t] <- phi[i, t] * s01[i, t]
      S[6, 9, i, t] <- 0
      S[6, 10, i, t] <- 0
      S[6, 11, i, t] <- 0
      S[6, 12, i, t] <- phi[i, t] * (1-s01[i, t]) * (1-eta[i, t, 3])
      S[6, 13, i, t] <- phi[i, t] * eta[i, t, 3] * (1-s01[i, t]) * (1 - beta[i, t])
      S[6, 14, i, t] <- 1-phi[i, t]
      
      S[7, 1, i, t] <- 0
      S[7, 2, i, t] <- 0
      S[7, 3, i, t] <- 0
      S[7, 4, i, t] <- 0
      S[7, 5, i, t] <- 0
      S[7, 6, i, t] <- phi[i, t] * eta[i, t, 3] * (1-s02[i, t])^2 * beta[i, t] * (1 - gamma[i, t])
      S[7, 7, i, t] <- phi[i, t] * eta[i, t, 3] * (1-s02[i, t])^2 * beta[i, t] * gamma[i, t]
      S[7, 8, i, t] <- phi[i, t] * 2 * s02[i, t] * (1-s02[i, t])
      S[7, 9, i, t] <- phi[i, t] * s02[i, t]^2
      S[7, 10, i, t] <- 0
      S[7, 11, i, t] <- 0
      S[7, 12, i, t] <- phi[i, t] * (1-s02[i, t])^2 * (1-eta[i, t, 3])
      S[7, 13, i, t] <- phi[i, t] * (1-s02[i, t])^2 * eta[i, t, 3] * (1 - beta[i, t])
      S[7, 14, i, t] <- 1-phi[i, t]
      
      S[8, 1, i, t] <- 0
      S[8, 2, i, t] <- 0
      S[8, 3, i, t] <- 0
      S[8, 4, i, t] <- 0
      S[8, 5, i, t] <- 0
      S[8, 6, i, t] <- phi[i, t] * eta[i, t, 4] * (1-s1[i, t]) * beta[i, t] * (1 - gamma[i, t])
      S[8, 7, i, t] <- phi[i, t] * eta[i, t, 4] * (1-s1[i, t]) * beta[i, t] * gamma[i, t]
      S[8, 8, i, t] <- 0
      S[8, 9, i, t] <- 0
      S[8, 10, i, t] <- phi[i, t] * s1[i, t]
      S[8, 11, i, t] <- 0
      S[8, 12, i, t] <- phi[i, t] * (1-s1[i, t]) * (1-eta[i, t, 4])
      S[8, 13, i, t] <- phi[i, t] * (1-s1[i, t]) * eta[i, t, 4] * (1 - beta[i, t])
      S[8, 14, i, t] <- 1-phi[i, t]
      
      S[9, 1, i, t] <- 0
      S[9, 2, i, t] <- 0
      S[9, 3, i, t] <- 0
      S[9, 4, i, t] <- 0
      S[9, 5, i, t] <- 0
      S[9, 6, i, t] <- phi[i, t] * eta[i, t, 4] * (1-s1[i, t])^2 * beta[i, t] * (1 - gamma[i, t])
      S[9, 7, i, t] <- phi[i, t] * eta[i, t, 4] * (1-s1[i, t])^2 * beta[i, t] * gamma[i, t]
      S[9, 8, i, t] <- 0
      S[9, 9, i, t] <- 0
      S[9, 10, i, t] <- phi[i, t] * 2 * s1[i, t] * (1-s1[i, t])
      S[9, 11, i, t] <- phi[i, t] * s1[i, t]^2
      S[9, 12, i, t] <- phi[i, t] * (1-s1[i, t])^2 * (1-eta[i, t, 4])
      S[9, 13, i, t] <- phi[i, t] * (1-s1[i, t])^2 * eta[i, t, 4] * (1 - beta[i, t])
      S[9, 14, i, t] <- 1-phi[i, t]
      
      S[10, 1, i, t] <- 0
      S[10, 2, i, t] <- 0
      S[10, 3, i, t] <- 0
      S[10, 4, i, t] <- 0
      S[10, 5, i, t] <- 0
      S[10, 6, i, t] <- phi[i, t] * eta[i, t, 5] * beta[i, t] * (1 - gamma[i, t])
      S[10, 7, i, t] <- phi[i, t] * eta[i, t, 5] * beta[i, t] * gamma[i, t]
      S[10, 8, i, t] <- 0
      S[10, 9, i, t] <- 0
      S[10, 10, i, t] <- 0
      S[10, 11, i, t] <- 0
      S[10, 12, i, t] <- phi[i, t] * (1 - eta[i, t, 5])
      S[10, 13, i, t] <- phi[i, t] * eta[i, t, 5] * (1 - beta[i, t])
      S[10, 14, i, t] <- 1-phi[i, t] 
      
      S[11, 1, i, t] <- 0
      S[11, 2, i, t] <- 0
      S[11, 3, i, t] <- 0
      S[11, 4, i, t] <- 0
      S[11, 5, i, t] <- 0
      S[11, 6, i, t] <- phi[i, t] * eta[i, t, 5] * beta[i, t] * (1 - gamma[i, t])
      S[11, 7, i, t] <- phi[i, t] * eta[i, t, 5] * beta[i, t] * gamma[i, t]
      S[11, 8, i, t] <- 0
      S[11, 9, i, t] <- 0
      S[11, 10, i, t] <- 0
      S[11, 11, i, t] <- 0
      S[11, 12, i, t] <- phi[i, t] * (1 - eta[i, t, 5])
      S[11, 13, i, t] <- phi[i, t] * eta[i, t, 5] * (1 - beta[i, t])
      S[11, 14, i, t] <- 1-phi[i, t] 
      
      S[12, 1, i, t] <- 0
      S[12, 2, i, t] <- 0
      S[12, 3, i, t] <- 0
      S[12, 4, i, t] <- 0
      S[12, 5, i, t] <- 0
      S[12, 6, i, t] <- phi[i, t] * eta[i, t, 1] * beta[i, t] * (1 - gamma[i, t])
      S[12, 7, i, t] <- phi[i, t] * eta[i, t, 1] * beta[i, t] * gamma[i, t]
      S[12, 8, i, t] <- 0
      S[12, 9, i, t] <- 0
      S[12, 10, i, t] <- 0
      S[12, 11, i, t] <- 0
      S[12, 12, i, t] <- phi[i, t] * (1 - eta[i, t, 1])
      S[12, 13, i, t] <- phi[i, t] * eta[i, t, 1] * (1 - beta[i, t])
      S[12, 14, i, t] <- 1-phi[i, t] 
      
      S[13, 1, i, t] <- 0
      S[13, 2, i, t] <- 0
      S[13, 3, i, t] <- 0
      S[13, 4, i, t] <- 0
      S[13, 5, i, t] <- 0
      S[13, 6, i, t] <- phi[i, t] * eta[i, t, 2] * beta[i, t] * (1 - gamma[i, t])
      S[13, 7, i, t] <- phi[i, t] * eta[i, t, 2] * beta[i, t] * gamma[i, t]
      S[13, 8, i, t] <- 0
      S[13, 9, i, t] <- 0
      S[13, 10, i, t] <- 0
      S[13, 11, i, t] <- 0
      S[13, 12, i, t] <- phi[i, t] * (1 - eta[i, t, 2])
      S[13, 13, i, t] <- phi[i, t] * eta[i, t, 2] * (1 - beta[i, t])
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
      E[1, 12, i, t] <- 0
      E[1, 13, i, t] <- 0
      E[1, 14, i, t] <- 0
      E[1, 15, i, t] <- 0
      E[1, 16, i, t] <- 0
      E[1, 17, i, t] <- 0
      E[1, 18, i, t] <- 0
      E[1, 19, i, t] <- 1-p[i, t+1, 1]      
      
      E[2, 1, i, t] <- 0
      E[2, 2, i, t] <- KHI[i, t+1] + (1 - KHI[i, t+1]) * p[i, t+1, 1]
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
      E[2, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p[i, t+1, 1]) 
      
      E[3, 1, i, t] <- 0
      E[3, 2, i, t] <- 0
      E[3, 3, i, t] <- KHI[i, t+1] + (1 - KHI[i, t+1]) * p[i, t+1, 1]
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
      E[3, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p[i, t+1, 1])   
      
      E[4, 1, i, t] <- 0
      E[4, 2, i, t] <- 0
      E[4, 3, i, t] <- 0
      E[4, 4, i, t] <- (1 - KHI[i, t+1]) * p[i, t+1, 1]
      E[4, 5, i, t] <- 0
      E[4, 6, i, t] <- 0
      E[4, 7, i, t] <- 0
      E[4, 8, i, t] <- 0
      E[4, 9, i, t] <- 0
      E[4, 10, i, t] <- 0
      E[4, 11, i, t] <- 0
      E[4, 12, i, t] <- KHI[i, t+1] * (1-p[i, t+1, 1]) + KHI[i, t+1] * p[i, t+1, 1]
      E[4, 13, i, t] <- 0
      E[4, 14, i, t] <- 0
      E[4, 15, i, t] <- 0
      E[4, 16, i, t] <- 0
      E[4, 17, i, t] <- 0
      E[4, 18, i, t] <- 0
      E[4, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p[i, t+1, 1])   
      
      E[5, 1, i, t] <- 0
      E[5, 2, i, t] <- 0
      E[5, 3, i, t] <- 0
      E[5, 4, i, t] <- (1 - KHI[i, t+1]) * p[i, t+1, 2]
      E[5, 5, i, t] <- 0
      E[5, 6, i, t] <- 0
      E[5, 7, i, t] <- 0
      E[5, 8, i, t] <- 0
      E[5, 9, i, t] <- 0
      E[5, 10, i, t] <- 0
      E[5, 11, i, t] <- 0
      E[5, 12, i, t] <- 0
      E[5, 13, i, t] <- KHI[i, t+1] * (1 - p[i, t+1, 2])
      E[5, 14, i, t] <- 0
      E[5, 15, i, t] <- 0
      E[5, 16, i, t] <- KHI[i, t+1] * p[i, t+1, 2]
      E[5, 17, i, t] <- 0
      E[5, 18, i, t] <- 0
      E[5, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p[i, t+1, 2])   
      
      E[6, 1, i, t] <- 0
      E[6, 2, i, t] <- 0
      E[6, 3, i, t] <- 0
      E[6, 4, i, t] <- 0
      E[6, 5, i, t] <- p[i, t+1, 2]
      E[6, 6, i, t] <- 0
      E[6, 7, i, t] <- 0
      E[6, 8, i, t] <- 0
      E[6, 9, i, t] <- 0
      E[6, 10, i, t] <- 0
      E[6, 11, i, t] <- 0
      E[6, 12, i, t] <- 0
      E[6, 13, i, t] <- KHI[i, t+1] * (1 - p[i, t+1, 2]) * AGE_5[i, t+1]
      E[6, 14, i, t] <- 0
      E[6, 15, i, t] <- KHI[i, t+1] * (1 - p[i, t+1, 2]) * ( 1-AGE_5[i, t+1] )
      E[6, 16, i, t] <- 0
      E[6, 17, i, t] <- 0
      E[6, 18, i, t] <- 0
      E[6, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p[i, t+1, 2])   
      
      E[7, 1, i, t] <- 0
      E[7, 2, i, t] <- 0
      E[7, 3, i, t] <- 0
      E[7, 4, i, t] <- 0
      E[7, 5, i, t] <- 0
      E[7, 6, i, t] <- p[i, t+1, 1]
      E[7, 7, i, t] <- 0
      E[7, 8, i, t] <- 0
      E[7, 9, i, t] <- 0
      E[7, 10, i, t] <- 0
      E[7, 11, i, t] <- 0
      E[7, 12, i, t] <- 0
      E[7, 13, i, t] <- KHI[i, t+1] * (1 - p[i, t+1, 1]) * AGE_5[i, t+1]
      E[7, 14, i, t] <- 0
      E[7, 15, i, t] <- KHI[i, t+1] * (1 - p[i, t+1, 1]) * ( 1-AGE_5[i, t+1] )
      E[7, 16, i, t] <- 0
      E[7, 17, i, t] <- 0
      E[7, 18, i, t] <- 0
      E[7, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p[i, t+1, 1])   
      
      E[8, 1, i, t] <- 0
      E[8, 2, i, t] <- 0
      E[8, 3, i, t] <- 0
      E[8, 4, i, t] <- 0
      E[8, 5, i, t] <- 0
      E[8, 6, i, t] <- 0
      E[8, 7, i, t] <- p[i, t+1, 1]
      E[8, 8, i, t] <- 0
      E[8, 9, i, t] <- 0
      E[8, 10, i, t] <- 0
      E[8, 11, i, t] <- 0
      E[8, 12, i, t] <- 0
      E[8, 13, i, t] <- 0
      E[8, 14, i, t] <- KHI[i, t+1] * (1 - p[i, t+1, 1])
      E[8, 15, i, t] <- 0
      E[8, 16, i, t] <- 0
      E[8, 17, i, t] <- 0
      E[8, 18, i, t] <- 0
      E[8, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p[i, t+1, 1])   
      
      E[9, 1, i, t] <- 0
      E[9, 2, i, t] <- 0
      E[9, 3, i, t] <- 0
      E[9, 4, i, t] <- 0
      E[9, 5, i, t] <- 0
      E[9, 6, i, t] <- 0
      E[9, 7, i, t] <- 0
      E[9, 8, i, t] <- p[i, t+1, 1]
      E[9, 9, i, t] <- 0
      E[9, 10, i, t] <- 0
      E[9, 11, i, t] <- 0
      E[9, 12, i, t] <- 0
      E[9, 13, i, t] <- 0
      E[9, 14, i, t] <- KHI[i, t+1] * (1 - p[i, t+1, 1])
      E[9, 15, i, t] <- 0
      E[9, 16, i, t] <- 0
      E[9, 17, i, t] <- 0
      E[9, 18, i, t] <- 0
      E[9, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p[i, t+1, 1])   
      
      E[10, 1, i, t] <- 0
      E[10, 2, i, t] <- 0
      E[10, 3, i, t] <- 0
      E[10, 4, i, t] <- 0
      E[10, 5, i, t] <- 0
      E[10, 6, i, t] <- 0
      E[10, 7, i, t] <- 0
      E[10, 8, i, t] <- 0
      E[10, 9, i, t] <- (1 - ALPHA[i, t+1]) * p[i, t+1, 1]
      E[10, 10, i, t] <- 0
      E[10, 11, i, t] <- ALPHA[i, t+1] * (1 - KHI[i, t+1]) * p[i, t+1, 1]
      E[10, 12, i, t] <- 0
      E[10, 13, i, t] <- 0
      E[10, 14, i, t] <- KHI[i, t+1] * (1 - p[i, t+1, 1])
      E[10, 15, i, t] <- 0
      E[10, 16, i, t] <- 0
      E[10, 17, i, t] <- ALPHA[i, t+1] * KHI[i, t+1] * p[i, t+1, 1]
      E[10, 18, i, t] <- 0
      E[10, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p[i, t+1, 1])   
      
      E[11, 1, i, t] <- 0
      E[11, 2, i, t] <- 0
      E[11, 3, i, t] <- 0
      E[11, 4, i, t] <- 0
      E[11, 5, i, t] <- 0
      E[11, 6, i, t] <- 0
      E[11, 7, i, t] <- 0
      E[11, 8, i, t] <- 0
      E[11, 9, i, t] <- 2 * (1 - ALPHA[i, t+1]) * ALPHA[i, t+1] * p[i, t+1, 1]
      E[11, 10, i, t] <- (1 - ALPHA[i, t+1])^2 * p[i, t+1, 1]
      E[11, 11, i, t] <- ALPHA[i, t+1]^2 * (1 - KHI[i, t+1]) * p[i, t+1, 1]
      E[11, 12, i, t] <- 0
      E[11, 13, i, t] <- 0
      E[11, 14, i, t] <- KHI[i, t+1] * (1 - p[i, t+1, 1])
      E[11, 15, i, t] <- 0
      E[11, 16, i, t] <- 0
      E[11, 17, i, t] <- ALPHA[i, t+1]^2 * KHI[i, t+1] * p[i, t+1, 1]
      E[11, 18, i, t] <- 0
      E[11, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p[i, t+1, 1])   
      
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
      E[12, 11, i, t] <- (1 - KHI[i, t+1]) * p[i, t+1, 1]
      E[12, 12, i, t] <- 0
      E[12, 13, i, t] <- 0
      E[12, 14, i, t] <- KHI[i, t+1] * (1 - p[i, t+1, 1])
      E[12, 15, i, t] <- 0
      E[12, 16, i, t] <- 0
      E[12, 17, i, t] <- KHI[i, t+1] * p[i, t+1, 1]
      E[12, 18, i, t] <- 0
      E[12, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p[i, t+1, 1])   
      
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
      E[13, 11, i, t] <- (1 - KHI[i, t+1]) * p[i, t+1, 1]
      E[13, 12, i, t] <- 0
      E[13, 13, i, t] <- 0
      E[13, 14, i, t] <- 0
      E[13, 15, i, t] <- KHI[i, t+1] * (1 - p[i, t+1, 1])
      E[13, 16, i, t] <- 0
      E[13, 17, i, t] <- 0
      E[13, 18, i, t] <- KHI[i, t+1] * p[i, t+1, 1]
      E[13, 19, i, t] <- (1 - KHI[i, t+1]) * (1 - p[i, t+1, 1])   
      
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
    beta_phi[k] ~ dnorm(0, sd = 1.5)  # survival
  }
  beta_phi[4] ~ dnorm(0, sd = 3)  # survival
  beta_phi[5] ~ dnorm(0, sd = 3)  # survival
  
  for (k in 1:3) {
    beta_s0[k] ~ dnorm(0, sd = 1.5)   # cub survival
  }
  for (k in 4:6) {
    beta_s0[k] ~ dnorm(0, sd = 3)   # cub survival
  }
  
  for (k in 1:3) {
    beta_s1[k] ~ dnorm(0, sd = 1.5)   # yearling survival
  }
  beta_s1[4] ~ dnorm(0, sd = 3)
  
  for (k in 1:4) {
    beta_eta[k] ~ dnorm(0, sd = 1.5)  # denning probability (intercepts)
  }
  for (k in 5:9) {
    beta_eta[k] ~ dnorm(0, sd = 3)  # denning probability (additive effects)
  }
  
  for (k in 1:3) {
    beta_beta[k] ~ dnorm(0, sd = 1.5)   # early litter survival
  }
  beta_beta[4] ~ dnorm(0, sd = 3)
  beta_beta[5] ~ dnorm(0, sd = 3)
  
  for (k in 1:3) {
    beta_gamma[k] ~ dnorm(0, sd = 1.5)   # litter size
  }
  for (k in 4:6) {
    beta_gamma[k] ~ dnorm(0, sd = 3)   # litter size
  }
  
  # Constraint
  constraint[1] ~ dconstraint(beta_s1[1] > beta_s0[1])
  constraint[2] ~ dconstraint(beta_s1[2] > beta_s0[2])
  constraint[3] ~ dconstraint(beta_s1[3] > beta_s0[3])
  
  #----- observation process -----#
  beta_p[1] ~ dnorm(0, sd = 1.5)
  beta_p[2] ~ dnorm(0, sd = 1)
  for (k in 3:4) {
    beta_p[k] ~ dnorm(0, sd = 3)
  }
  # Yearly random effects on recapture
  for (t in 1:(n_occasions-1)) {
    omega[t] ~ dnorm(0, sd = sigma_omega)
  }
  sigma_omega ~ dunif(0, 15)
  
  # +++++++++++++++++++++++++++ model & likelihood +++++++++++++++++++++++++++++
  
  for(i in 1:n_CH){
    
    # First capture
    for (k in 1:n_states) {
      zeta[i, f[i], k] <- S0[k] * E0[k, CH[i, f[i]], i]
    }
    
    # Subsequent captures
    for(t in f[i] : (l[i]-1) ){
      for (k in 1:n_states) {
        zeta[i, t+1, k] <- inprod(zeta[i, t, 1:n_states],
                                  S[1:n_states, k, i, t]) *
          E[k, CH[i, t+1], i, t]
      }
    }
    lik[i] <- sum(zeta[i, l[i], 1:n_states]) 
    ones[i] ~ dbin(prob = lik[i], size = freq[i])
  }
  
  # for (i in 1:n_ind) {
  #   # latent state at first capture
  #   z[i, f[i]] ~ dcat(S0[1:14])
  #   y[i, f[i]] ~ dcat(E0[z[i, f[i]], 1:19, i])
  #   
  #   for (t in ( f[i]+1 ):l[i]) {
  #     
  #     #-------- State process --------#
  #     # draw z(t) given z(t-1)
  #     z[i, t] ~ dcat( S[z[i, t-1], 1:14, i, t-1] )
  #     
  #     #----- Observation process -----#
  #     # draw y(t) given z(t)
  #     y[i, t] ~ dcat( E[z[i, t], 1:19, i, t-1] )
  #   }
  # }
})


# 5. Run the model -------------------------------------------------------------

# Need to run that twice (in parallel) to get two chains. ~22h per chain on an 
# HPC cluster.

start <- Sys.time()   
# Create R model
inits_values <- inits()
CR_model_biologging_R <- nimbleModel(code = CR_model_biologging,
                               constants = my.constants,
                               data = dat,
                               inits = inits_values,
                               calculate = FALSE)

# Compile model (in C++)
CR_model_biologging_C <- compileNimble(CR_model_biologging_R,
                                 showCompilerOutput = FALSE)

# Configure of MCMC
MCMC_conf <- configureMCMC(CR_model_biologging_R,
                           monitors = params)

# Compile MCMC
MCMC_R <- buildMCMC(MCMC_conf)
MCMC_C <- compileNimble(MCMC_R, project = CR_model_biologging_R)

fit_remote_1 <- runMCMC(mcmc = MCMC_C,
                               nchains = 1,
                               niter = 25000,
                               nburnin = 10000,
                               thin = 50,
                               inits = inits())

save(fit_remote_1,
     file = "03_outputs/fit_remote_1.RData")
end <- Sys.time() ; end - start


