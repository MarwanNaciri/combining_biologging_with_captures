#==============================================================================#
#                                                                              #
#                          Plot supplementary figures                          #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(lubridate)
library(ggdist)
library(patchwork)
library(nimble) # for logit function
Sys.setenv(LANG = "en")
source("05_script/functions_for_models.R")

get_median_and_CI <- function(x, lower, upper) {
  output <- data.frame(y = median(x), 
                       ymin = quantile(x, probs = lower),
                       ymax = quantile(x, probs = upper))
  return(output)
}



# ~ Fig. S4. Improvements estimates rock partridge ------------------------------

# Real value
{beta <- eta <- delta <- NULL
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
delta[1] <- eta[1]*beta[1]
delta[2] <- eta[2]*beta[2]
delta[3] <- eta[3]*beta[3]
gamma <- 0.4 # Probability of large litter/brood
}

# ~~~ a. Load and process data -------------------------------------------------

path <- c("03_outputs/simulations/model_outputs/rock_partridge/")

list_files_all <- list.files(path)

p <- c(0.05, 0.1, 0.25, 0.5) 
p_GPS <- c(0, 0.1, 0.25, 0.5,  1) 
pError <- 0
cols <- c("phiJ", "phiA", "a2_phi", 
          "delta_1", "delta_2", "delta_3", 
          "gamma") 

real_values <- data.frame(parameter = factor(cols, levels = cols),
                          real_value = c(plogis(a1_phi), plogis(a2_phi), a3_phi,
                                         delta[1], delta[2], delta[3], 
                                         gamma))

df_all <- data.frame()
for (n in 1:length(p)) {
  for (m in 1:length(p_GPS)) {
    
    print(paste0("p = ", p[n]))
    print(paste0("pGPS = ", p_GPS[m]))
    
    scenario <- paste0("p_", p[n], "_pGPS_", p_GPS[m], "_")
    file_list <- list_files_all[grep(list_files_all, pattern = scenario)]
    
    model_fits <- list() 
    for (l in 1:length(file_list)) {
      load(paste0(path, file_list[l]))
      
      
      model_fits[[l]] <- rbind(get(substr(file_list[l], 1, nchar(file_list[l])-6))[[1]],
                               get(substr(file_list[l], 1, nchar(file_list[l])-6))[[2]])
      
      rm(list = substr(file_list[l], 1, nchar(file_list[l])-6))
    }
    
    #++++++++++++++++++++++++++++++ Without GPS ++++++++++++++++++++++++++++++++#
    if (p_GPS[m] == 0) {
      
      df <- as.data.frame(do.call(rbind, model_fits)) %>%
        janitor::clean_names() %>%
        mutate(simulation = rep(1:length(file_list), 
                                each = nrow(model_fits[[1]])),
               iteration = rep(1:nrow(model_fits[[1]]), 
                               times = length(file_list))) %>%
        mutate(phiJ = plogis(a1_phi),
               phiA = plogis(a2_phi)) %>%
        dplyr::select(all_of(cols), simulation) %>%
        pivot_longer(cols = all_of(cols),
                     names_to = "parameter") %>%
        group_by(parameter,
                 simulation) %>%
        summarize(median = median(value), sd = sd(value)) %>%
        group_by(parameter) %>%
        summarize(upper_95 = quantile(median, probs = 0.975),
                  lower_95 = quantile(median, probs = 0.025),
                  upper_50 = quantile(median, probs = 0.75),
                  lower_50 = quantile(median, probs = 0.25),
                  mean = mean(median)) %>%
        mutate(p = p[n],
               p_GPS = p_GPS[m]) 
      
      df_all <- rbind(df_all, df)
      
    } else { 
      
      #+++++++++++++++++++++++++++++ With GPS ++++++++++++++++++++++++++++++++#
      
      df <- as.data.frame(do.call(rbind, model_fits)) %>%
        janitor::clean_names() %>%
        mutate(simulation = rep(1:length(file_list), 
                                each = nrow(model_fits[[1]])),
               iteration = rep(1:nrow(model_fits[[1]]), 
                               times = length(file_list))) %>%
        mutate(phiJ = plogis(a1_phi),
               phiA = plogis(a2_phi),
               delta_1 = beta_1 * eta_1,
               delta_2 = beta_2 * eta_2,
               delta_3 = beta_3 * eta_3) %>%
        dplyr::select(all_of(cols), simulation) %>%
        pivot_longer(cols = all_of(cols),
                     names_to = "parameter") %>%
        group_by(parameter,
                 simulation) %>%
        summarize(median = median(value), sd = sd(value)) %>%
        group_by(parameter) %>%
        summarize(upper_95 = quantile(median, probs = 0.975),
                  lower_95 = quantile(median, probs = 0.025),
                  upper_50 = quantile(median, probs = 0.75),
                  lower_50 = quantile(median, probs = 0.25),
                  mean = mean(median)) %>%
        mutate(p = p[n],
               p_GPS = p_GPS[m]) 
      
      df_all <- rbind(df_all, df)
    }
  }
}


df_plot <- df_all %>%
  mutate(p_GPS = as.factor(p_GPS),
         p = as.factor(p))

df_plot_title <- c(expression(paste(italic("\u03C6")['J'])),
                   expression(paste(italic("\u03C6")['Ad'])),
                   expression(paste(italic("\u03C6")['time'])),

                   # delta
                   expression(paste(italic("\u03B4")['J'])),
                   expression(paste(italic("\u03B4")['NB-FB'])),
                   expression(paste(italic("\u03B4")['SB'])),

                   # gamma
                   expression(paste(italic("\u03B3"))))

ymin <- c(0, 0, -1.5, 0, 0, 0, 0)
ymax <- c(1, 1, 1.5, 1, 1, 1, 1)


# ~~~ b. Plot ------------------------------------------------------------------

for (k in 1:length(cols)) {
  df_plot_k <- df_plot %>%
    filter(parameter == cols[k])
  
  plot_title_k <- parse(text = df_plot_title[k])[[1]]
  
  assign(x = paste0("plot_", k),
         value =  ggplot(data = df_plot_k,
                         aes(x = p , y = mean, color = p_GPS)) + 
           geom_pointrange(aes(ymin = lower_95, ymax = upper_95),
                           position = position_dodge(width = 0.4)) +
           geom_pointrange(aes(ymin = lower_50, ymax = upper_50),
                           position = position_dodge(width = 0.4),
                           linewidth = 1) +
           geom_hline(yintercept = real_values$real_value[k]) +
           scale_color_grey(start = 0.8, end = 0.2) +
           scale_y_continuous(limits = c(ymin[k], ymax[k])) +
           theme_bw() +
           theme(plot.title = element_text(hjust = 0.5, size = 10),
                 legend.title = element_text(hjust = 0.5),
                 axis.title.y = element_blank()) +
           labs(x = expression(italic("p")), y = bquote(.(plot_title_k)), 
                title = bquote(.(plot_title_k)),
                color = expression(kappa))
  )
}


full_plot <- plot_1 + plot_2 + plot_3 + plot_4 + plot_5 + plot_6 + plot_7 +
  plot_layout(guides = 'collect')

ggsave(plot = full_plot, 
       filename = "03_outputs/Figure S4.png",
       units = 'cm', width = 18, height = 16, dpi = 600)



# ~ Fig. S5. Improvements estimates albatross ---------------------------------------

# Real value
{b_phi <- b_psi1 <- b_psi2 <- b_eta <- NULL
b_phi[1] <- 0.7 # 1 yr old
b_phi[2] <- 0.8 # 2 yr old
b_phi[3] <- 0.88 # 3-5 yr old
b_phi[4] <- 0.95 # 6+ yr old

b_psi1[1] <- 0.2 # 3-5 yr old
b_psi1[2] <- 0.2 # 6 yr+ non-breeder/pre-breeder
b_psi1[3] <- 0.3 # 6 yr+ failed breeder 
b_psi1[4] <- 0.8 # 6 yr+ successful breeder 

b_psi2[1] <- 0.4 # 3-5 yr old
b_psi2[2] <- 0.75 # 6+ yr old

b_eta[1] <- 0.4 # 6 yr old
b_eta[2] <- 0.7 # 7-8 yr old
b_eta[3] <- 0.8 # 9+ yr old
b_eta[4] <- 0.3 # Adults (not-away)
b_eta[5] <- 0.8 # Adults (away)

a1_beta <- 0.75 # nest success

theta <- 0.6
}


# ~~~ a. Load and process data -------------------------------------------------

path <- c("03_outputs/simulations/model_outputs/albatross/")

list_files_all <- list.files(path)


p <- c(0.1, 0.25, 0.5) 
pGPS <- c(0, 0.1, 0.25, 0.5, 1) 
pError <- 0
cols <- c("b_phi_1", "b_phi_2", "b_phi_3",
          "b_psi1_1", "b_psi1_2", "b_psi1_3", "b_psi1_4",
          "b_psi2_1", "b_psi2_2",
          "b_eta_1", "b_eta_2", "b_eta_3", "b_eta_4", "b_eta_5",
          "a1_beta") # , "theta") 
real_values <- data.frame(parameter = factor(cols, levels = cols),
                          real_value = c(b_phi[1]*b_phi[2], b_phi[3], b_phi[4],
                                         b_psi1[1], b_psi1[2], b_psi1[3], b_psi1[4],
                                         b_psi2[1], b_psi2[2],
                                         b_eta[1], b_eta[2], b_eta[3], b_eta[4], b_eta[5],
                                         a1_beta)) # , theta))

df_all <- data.frame()
for (n in 1:length(p)) {
  for (m in 1:length(pGPS)) {
    
    print(paste0("p = ", p[n]))
    print(paste0("pGPS = ", pGPS[m]))
    
    scenario <- paste0("p_", p[n], "_pGPS_", pGPS[m], "_")
    file_list <- list_files_all[grep(list_files_all, pattern = scenario)]
    
    model_fits <- list() 
    for (l in 1:length(file_list)) {
      load(paste0(path, file_list[l]))
      
      
      model_fits[[l]] <- rbind(get(substr(file_list[l], 1, nchar(file_list[l])-6))[[1]],
                               get(substr(file_list[l], 1, nchar(file_list[l])-6))[[2]])
      
      rm(list = substr(file_list[l], 1, nchar(file_list[l])-6))
    }
    
    #++++++++++++++++++++++++++++++ Without GPS ++++++++++++++++++++++++++++++++#
    if (pGPS[m] == 0) {
      
      df_1 <- as.data.frame(do.call(rbind, model_fits))
      
      df_2 <- df_1 %>%
        mutate(simulation = rep(1:length(file_list), 
                                each = nrow(model_fits[[1]]))) %>%
        janitor::clean_names() %>%
        dplyr::select(all_of(cols), simulation) %>%
        pivot_longer(cols = all_of(cols),
                     names_to = "parameter") %>%
        group_by(parameter,
                 simulation) %>%
        summarize(median = median(value), sd = sd(value)) %>%
        group_by(parameter) %>%
        summarize(upper_95 = quantile(median, probs = 0.975),
                  lower_95 = quantile(median, probs = 0.025),
                  upper_50 = quantile(median, probs = 0.75),
                  lower_50 = quantile(median, probs = 0.25),
                  mean = mean(median)) %>%
        mutate(p = p[n],
               p_GPS = pGPS[m]) 
      
      df_all <- rbind(df_all, df_2)
    } else { 
      
      #+++++++++++++++++++++++++++++ With GPS ++++++++++++++++++++++++++++++++#
      
      df_1 <- as.data.frame(do.call(rbind, model_fits))
      
      df_2 <- df_1 %>%
        # mutate(simulation = rep(1:length(file_list), 
        #                         each = nrow(model_fits[[1]])),
        #        iteration = rep(1:nrow(model_fits[[1]]), 
        #                        times = length(file_list))) %>%
        mutate(simulation = rep(1:length(file_list), 
                                each = nrow(model_fits[[1]]))) %>%
        janitor::clean_names() %>%
        mutate(b_phi_1 = b_phi_1 * b_phi_2,
               b_phi_2 = b_phi_3,
               b_phi_3 = b_phi_4) %>%
        dplyr::select(all_of(cols), simulation) %>%
        pivot_longer(cols = all_of(cols),
                     names_to = "parameter") %>%
        group_by(parameter,
                 simulation) %>%
        summarize(median = median(value), sd = sd(value)) %>%
        group_by(parameter) %>%
        summarize(upper_95 = quantile(median, probs = 0.975),
                  lower_95 = quantile(median, probs = 0.025),
                  upper_50 = quantile(median, probs = 0.75),
                  lower_50 = quantile(median, probs = 0.25),
                  mean = mean(median)) %>%
        mutate(p = p[n],
               p_GPS = pGPS[m]) 
      
      df_all <- rbind(df_all, df_2)
    }
    
  }
}

df_plot <- df_all %>%
  mutate(p_GPS = as.factor(p_GPS),
         p = as.factor(p))


df_plot_title <- c(expression(paste(italic("\u03C6")['1-2yr (biennal)'])),
                   expression(paste(italic("\u03C6")['3-5yr'])),
                   expression(paste(italic("\u03C6")['\u22656yr'])),
                   # Psi 1
                   expression(italic("\u03A8")['3-5yr']^1),
                   expression(paste(italic("\u03A8")['6yr, \u22656yr NB']^1)),
                   expression(paste(italic("\u03A8")['\u22656yr FB']^1)),
                   expression(paste(italic("\u03A8")['\u22656yr SB']^1)),
                   # Psi 2
                   expression(paste(italic("\u03A8")['3-5yr']^2)),
                   expression(paste(italic("\u03A8")['\u22656yr']^2)),
                   # Eta
                   expression(paste(italic("\u03B7")['6yr PrB'])),
                   expression(paste(italic("\u03B7")['7-8yr PrB'])),
                   expression(paste(italic("\u03B7")['\u22659yr PrB'])),
                   expression(paste(italic("\u03B7")['Adult'])),
                   expression(paste(italic("\u03B7")['Adult sabbatical'])),
                   # Beta
                   expression(paste(italic("\u03B2"))),
                   # Theta
                   expression(paste(italic("\u03B8"))))

# ~~~ b. Plot --------------------------------------------------------------------

ymin <- 0
ymax <- 1

for (k in 1:length(cols)) {
  df_plot_k <- df_plot %>%
    filter(parameter == cols[k])
  
  plot_title_k <- parse(text = df_plot_title[k])[[1]]
  assign(x = paste0("plot_", k),
         value =  ggplot(data = df_plot_k,
                         aes(x = p , y = mean, color = p_GPS)) + 
           geom_pointrange(aes(ymin = lower_95, ymax = upper_95),
                           position = position_dodge(width = 0.4)) +
           geom_pointrange(aes(ymin = lower_50, ymax = upper_50),
                           position = position_dodge(width = 0.4),
                           linewidth = 1) +
           geom_hline(yintercept = real_values$real_value[k]) +
           scale_color_grey(start = 0.8, end = 0.2) +
           scale_y_continuous(limits = c(ymin, ymax)) +
           theme_bw() +
           theme(plot.title = element_text(hjust = 0.5, size = 10),
                 legend.title = element_text(hjust = 0.5),
                 axis.title.y = element_blank()) +
           labs(x = expression(italic("p")), y = bquote(.(plot_title_k)), 
                title = bquote(.(plot_title_k)),
                color = expression(kappa))
  )
}



full_plot <- plot_1 + plot_2 + plot_3 + plot_4 + plot_5 + plot_6 + 
  plot_7 + plot_8 + plot_9 + plot_10 + plot_11 + plot_12 + 
  plot_13 + plot_14 + plot_15 + 
  plot_layout(guides = 'collect')

ggsave(plot = full_plot, 
       filename = "03_outputs/Figure S5.png",
       units = 'cm', width = 18, height = 20, dpi = 600)



# ~ Fig. S6. Improvements estimates polar bear --------------------------------------

# ~~~ a. Load and process ------------------------------------------------------

{# Survival
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
  
  # Compute the breeding probability (eta * beta)
  delta <- NULL
  delta[1] <- beta_eta[1] * beta_beta[1]
  delta[2] <- beta_eta[2] * beta_beta[1]
  delta[3] <- beta_eta[3] * beta_beta[2]
  delta[4] <- beta_eta[4] * beta_beta[3]
  # Recapture probability 
}


cols <- c("phi_1", "phi_2", "phi_3", 
          "delta_1", "delta_2", "delta_3", "delta_4",
          "gamma_1", "gamma_2", "gamma_3",
          "s0", "s1") 
real_values <- data.frame(parameter = factor(cols, levels = cols),
                          real_value = c(beta_phi[1],  beta_phi[2], beta_phi[3],
                                         delta[1], delta[2], delta[3],  delta[4],
                                         beta_gamma[1], beta_gamma[2], beta_gamma[3],
                                         s0, s1))


path_model_outputs <- "03_outputs/simulations/model_outputs/polar_bear/"
file_list_all <- list.files(paste0(path_model_outputs))


p <- c(0.05, 0.1, 0.25, 0.5) 
pGPS <- c(0, 0.1, 0.25, 0.5, 1)
df_all <- data.frame()
for (n in 1:length(p)) {
  
  # index_start <- ifelse(n %in% 1, 2, 1)
  for (m in 1:length(pGPS)) {
    
    print(paste0("p = ", p[n]))
    print(paste0("pGPS = ", pGPS[m]))
    
    scenario <- paste0("p_", p[n], "_pGPS_",pGPS[m], "_")
    file_list <- file_list_all[grep(file_list_all, pattern = scenario)]
    file_list <- file_list#[1:50]
    
    model_fits <- list() 
    for (l in 1:length(file_list)) {
      load(paste0(path_model_outputs, file_list[l]))
      
      
      model_fits[[l]] <- rbind(get(substr(file_list[l], 1, nchar(file_list[l])-6))[[1]],
                               get(substr(file_list[l], 1, nchar(file_list[l])-6))[[2]])
      
      rm(list = substr(file_list[l], 1, nchar(file_list[l])-6))
    }
    
    if (pGPS[m] == 0) {
      df_1 <- as.data.frame(do.call(rbind, model_fits)) %>%
        mutate(simulation = rep(1:length(file_list), 
                                each = nrow(model_fits[[1]])),
               iteration = rep(1:nrow(model_fits[[1]]), 
                               times = length(file_list))) %>%
        janitor::clean_names() %>%
        rename(phi_1 = beta_phi_1,
               phi_2 = beta_phi_2,
               phi_3 = beta_phi_3,
               delta_1 = beta_beta_1,
               delta_2 = beta_beta_2,
               delta_3 = beta_beta_3,
               delta_4 = beta_beta_4,
               gamma_1 = beta_gamma_1,
               gamma_2 = beta_gamma_2,
               gamma_3 = beta_gamma_3) %>%
        dplyr::select(all_of(cols), simulation) %>%
        pivot_longer(cols = all_of(cols),
                     names_to = "parameter") %>%
        group_by(parameter,
                 simulation) %>%
        summarize(median = median(value), sd = sd(value)) %>%
        group_by(parameter) %>%
        summarize(upper_95 = quantile(median, probs = 0.975),
                  lower_95 = quantile(median, probs = 0.025),
                  upper_50 = quantile(median, probs = 0.75),
                  lower_50 = quantile(median, probs = 0.25),
                  mean = mean(median)) %>%
        mutate(p = p[n],
               p_GPS = pGPS[m]) 
      
      df_all <- rbind(df_all, df_1)
      
    } else { # If pGPS = 0.5 or 1
      df_1 <- as.data.frame(do.call(rbind, model_fits)) %>%
        mutate(simulation = rep(1:length(file_list), 
                                each = nrow(model_fits[[1]])),
               iteration = rep(1:nrow(model_fits[[1]]), 
                               times = length(file_list))) %>%
        janitor::clean_names() %>%
        mutate(delta_1 = beta_eta_1 * beta_beta_1,
               delta_2 = beta_eta_2 * beta_beta_1,
               delta_3 = beta_eta_3 * beta_beta_2,
               delta_4 = beta_eta_4 * beta_beta_3) %>%
        rename(phi_1 = beta_phi_1,
               phi_2 = beta_phi_2,
               phi_3 = beta_phi_3,
               gamma_1 = beta_gamma_1,
               gamma_2 = beta_gamma_2,
               gamma_3 = beta_gamma_3) %>%
        dplyr::select(all_of(cols), simulation) %>%
        pivot_longer(cols = all_of(cols),
                     names_to = "parameter") %>%
        group_by(parameter,
                 simulation) %>%
        summarize(median = median(value), sd = sd(value)) %>%
        group_by(parameter) %>%
        summarize(upper_95 = quantile(median, probs = 0.975),
                  lower_95 = quantile(median, probs = 0.025),
                  upper_50 = quantile(median, probs = 0.75),
                  lower_50 = quantile(median, probs = 0.25),
                  mean = mean(median)) %>%
        mutate(p = p[n],
               p_GPS = pGPS[m]) 
      
      df_all <- rbind(df_all, df_1)
    }
  }
}

df_plot <- df_all %>%
  mutate(p_GPS = as.factor(p_GPS),
         p = as.factor(p))


df_plot_title <- c(expression(paste(italic("\u03C6")['2-4yr '])),
                   expression(paste(italic("\u03C6")['5-14yr'])),
                   expression(paste(italic("\u03C6")['\u226515yr'])),
                   # delta
                   expression(paste(italic("\u03B4")['4yr'])),
                   expression(paste(italic("\u03B4")['5-8yr'])),
                   expression(paste(italic("\u03B4")['9-14yr'])),
                   expression(paste(italic("\u03B4")['\u226515yr'])),
                   # gamma
                   expression(paste(italic("\u03B3")['4-8yr'])),
                   expression(paste(italic("\u03B3")['9-14yr'])),
                   expression(paste(italic("\u03B3")['\u226515yr'])),
                   # 
                   expression(paste(italic("s")[0])),
                   expression(paste(italic("s")[1])))

ymin <- rep(0, times = length(df_plot_title))
ymax <- rep(1, times = length(df_plot_title))


# ~~~ b. Plot --------------------------------------------------------------------

for (k in 1:length(cols)) {
  df_plot_k <- df_plot %>%
    filter(parameter == cols[k])
  
  plot_title_k <- parse(text = df_plot_title[k])[[1]]
  
  if (k %in% c(1, 4)) {   # If first column
    if (k %in% c(4)) {  # if last line
      assign(x = paste0("plot_", k),
             value =  ggplot(data = df_plot_k,
                             aes(x = p , y = mean, color = p_GPS)) + 
               geom_pointrange(aes(ymin = lower_95, ymax = upper_95),
                               position = position_dodge(width = 0.4)) +
               geom_pointrange(aes(ymin = lower_50, ymax = upper_50),
                               position = position_dodge(width = 0.4), 
                               linewidth = 1) +
               geom_hline(yintercept = real_values$real_value[k]) +
               scale_color_grey(start = 0.8, end = 0.2) +
               scale_y_continuous(limits = c(ymin[k], ymax[k])) +
               theme_bw() +
               theme(plot.title = element_text(hjust = 0.5, size = 10),
                     legend.title = element_text(hjust = 0.5),
                     axis.title.y = element_blank()) +
               labs(x = "p", y = bquote(.(plot_title_k)), 
                    title = bquote(.(plot_title_k)),
                    color = expression(kappa))
      )
    } else { # If first two line
      assign(x = paste0("plot_", k),
             value =  ggplot(data = df_plot_k,
                             aes(x = p , y = mean, color = p_GPS)) + 
               geom_pointrange(aes(ymin = lower_95, ymax = upper_95),
                               position = position_dodge(width = 0.4)) +
               geom_pointrange(aes(ymin = lower_50, ymax = upper_50),
                               position = position_dodge(width = 0.4), 
                               linewidth = 1) +
               geom_hline(yintercept = real_values$real_value[k]) +
               scale_color_grey(start = 0.8, end = 0.2) +
               scale_y_continuous(limits = c(ymin[k], ymax[k])) +
               theme_bw() +
               theme(plot.title = element_text(hjust = 0.5, size = 10),
                     legend.title = element_text(hjust = 0.5),
                     axis.title.y = element_blank()) +
               labs(x = "p", y = bquote(.(plot_title_k)), 
                    title = bquote(.(plot_title_k)),
                    color = expression(kappa))
      )
    }
  } else {   # If not first column
    if (k %in% c(5)) {  # if last line
      assign(x = paste0("plot_", k),
             value =  ggplot(data = df_plot_k,
                             aes(x = p , y = mean, color = p_GPS)) + 
               geom_pointrange(aes(ymin = lower_95, ymax = upper_95),
                               position = position_dodge(width = 0.4)) +
               geom_pointrange(aes(ymin = lower_50, ymax = upper_50),
                               position = position_dodge(width = 0.4), 
                               linewidth = 1) +
               geom_hline(yintercept = real_values$real_value[k]) +
               scale_color_grey(start = 0.8, end = 0.2) +
               scale_y_continuous(limits = c(ymin[k], ymax[k])) +
               theme_bw() +
               theme(plot.title = element_text(hjust = 0.5, size = 10),
                     legend.title = element_text(hjust = 0.5),
                     axis.title.y = element_blank()) +
               labs(x = "p", y = bquote(.(plot_title_k)), 
                    title = bquote(.(plot_title_k)),
                    color = expression(kappa))
      )
    } else {
      assign(x = paste0("plot_", k),
             value =  ggplot(data = df_plot_k,
                             aes(x = p , y = mean, color = p_GPS)) + 
               geom_pointrange(aes(ymin = lower_95, ymax = upper_95),
                               position = position_dodge(width = 0.4)) +
               geom_pointrange(aes(ymin = lower_50, ymax = upper_50),
                               position = position_dodge(width = 0.4), 
                               linewidth = 1) +
               geom_hline(yintercept = real_values$real_value[k]) +
               scale_color_grey(start = 0.8, end = 0.2) +
               scale_y_continuous(limits = c(ymin[k], ymax[k])) +
               theme_bw() +
               theme(plot.title = element_text(hjust = 0.5, size = 10),
                     legend.title = element_text(hjust = 0.5),
                     axis.title.y = element_blank()) +
               labs(x = expression(italic("p")), y = bquote(.(plot_title_k)), 
                    title = bquote(.(plot_title_k)),
                    color = expression(kappa))
      )
    }
  }
}

full_plot <- plot_1 + plot_2 + plot_3 + plot_4 + plot_5 + plot_6 + plot_7 + plot_8 +
  plot_9 + plot_10 + plot_11 + plot_12 +
  plot_layout(guides = 'collect')

ggsave(full_plot,
       filename = "03_outputs/Figure S6.png",
       width = 18, height = 20, units = "cm", dpi = 600)






# ~ Fig. S7: Improvement estimates roe deer increasing information content -------

# Real value
{beta <- eta <- delta <- NULL
a1_phi <- logit(0.4)  # Survival (Fawn)   
a2_phi <- logit(0.85)  # Survival 
a3_phi <- -0.4  # time effect   
time <- as.vector(scale(1:8))
beta[1] <- 0.3  # Breeding probability (Y) 
beta[2] <- 0.4  # Breeding probability (NB/FB) 
beta[3] <- 0.75  # Breeding probability (SB) 
eta[1] <- 0.3  # Early offspring survival probability (Y) 
eta[2] <- 0.4  # Early offspring survival probability (NB/FB) 
eta[3] <- 0.75  # Early offspring survival probability (SB) 
gamma <- 0.4 # Probability of large litter/brood
delta[1] <- beta[1] * eta[1]
delta[2] <- beta[2] * eta[2]
delta[3] <- beta[3] * eta[3]
gamma <- 0.4 # Probability of large litter/brood
}


# ~~~ a. Load and process data -------------------------------------------------

p <- 0.05
p_GPS <- c(0, 0.1, 0.25, 0.5, 1) 
p_Error <- 0

cols <- c("phiJ", "phiA", "a3_phi", 
          "delta_1", "delta_2", "delta_3", 
          "gamma") 

real_values <- data.frame(parameter = factor(cols, levels = cols),
                          real_value = c(plogis(a1_phi), plogis(a2_phi),  a3_phi, 
                                         delta[1], delta[2], delta[3], 
                                         gamma))
scenario <- c("baseline", "known_offspring_survival", "known_litter_size")
df_all <- data.frame()
for (n in 1:length(scenario)) {
  
  print(paste0("scenario = ", scenario[n]))
  
  path <- paste0("03_outputs/simulations/model_outputs/main_example/",
                 scenario[n], "/")
  list_files_all <- list.files(path)
  
  for (k in 1:length(p)) {
    
    # to accommodate the fact that I only have pGPS = 0 in the first scenario
    index_start <- ifelse(n == 1, 1, 2)
    for (m in index_start:length(p_GPS)) {
      
      parameter_values <- paste0("p_", p[k], "_pGPS_", p_GPS[m], "_")
      file_list <- list_files_all[grep(list_files_all, pattern = parameter_values)]
      
      model_fits <- list() 
      for (l in 1:length(file_list)) {
        load(paste0(path, file_list[l]))
        
        
        model_fits[[l]] <- rbind(get(substr(file_list[l], 1, nchar(file_list[l])-6))[[1]],
                                 get(substr(file_list[l], 1, nchar(file_list[l])-6))[[2]])
        
        rm(list = substr(file_list[l], 1, nchar(file_list[l])-6))
      }
      
      #++++++++++++++++++++++++++++++ Without GPS ++++++++++++++++++++++++++++++++#
      if (p_GPS[m] == 0) {
        
        df <- as.data.frame(do.call(rbind, model_fits)) %>%
          janitor::clean_names() %>%
          mutate(simulation = rep(1:length(file_list), 
                                  each = nrow(model_fits[[1]])),
                 iteration = rep(1:nrow(model_fits[[1]]), 
                                 times = length(file_list))) %>%
          mutate(phiJ = plogis(a1_phi),
                 phiA = plogis(a2_phi)) %>%
          dplyr::select(all_of(cols), simulation) %>%
          pivot_longer(cols = all_of(cols),
                       names_to = "parameter") %>%
          group_by(parameter,
                   simulation) %>%
          summarize(median = median(value), sd = sd(value)) %>%
          group_by(parameter) %>%
          summarize(upper_95 = quantile(median, probs = 0.975),
                    lower_95 = quantile(median, probs = 0.025),
                    upper_50 = quantile(median, probs = 0.75),
                    lower_50 = quantile(median, probs = 0.25),
                    mean = mean(median)) %>%
          mutate(p = p[k],
                 scenario = n,
                 p_GPS = p_GPS[m]) 
        
        df_all <- rbind(df_all, df)
        
      } else { 
        
        #+++++++++++++++++++++++++++++ With GPS ++++++++++++++++++++++++++++++++#
        
        df <- as.data.frame(do.call(rbind, model_fits)) %>%
          janitor::clean_names() %>%
          mutate(simulation = rep(1:length(file_list), 
                                  each = nrow(model_fits[[1]])),
                 iteration = rep(1:nrow(model_fits[[1]]), 
                                 times = length(file_list))) %>%
          mutate(phiJ = plogis(a1_phi),
                 phiA = plogis(a2_phi),
                 delta_1 = beta_1 * eta_1,
                 delta_2 = beta_2 * eta_2,
                 delta_3 = beta_3 * eta_3) %>%
          dplyr::select(all_of(cols), simulation) %>%
          pivot_longer(cols = all_of(cols),
                       names_to = "parameter") %>%
          group_by(parameter,
                   simulation) %>%
          summarize(median = median(value), sd = sd(value)) %>%
          group_by(parameter) %>%
          summarize(upper_95 = quantile(median, probs = 0.975),
                    lower_95 = quantile(median, probs = 0.025),
                    upper_50 = quantile(median, probs = 0.75),
                    lower_50 = quantile(median, probs = 0.25),
                    mean = mean(median)) %>%
          mutate(p = p[k],
                 scenario = n,
                 p_GPS = p_GPS[m]) 
        
        df_all <- rbind(df_all, df)
      }
    }
  }
}
df_plot <- df_all %>%
  mutate(p_GPS = as.factor(p_GPS),
         scenario = as.factor(scenario))

df_plot_title <- c(expression(paste(italic("\u03C6")['1'])),
                   expression(paste(italic("\u03C6")['2'])),
                   expression(paste(italic("\u03C6")['time'])),
                   
                   # delta
                   expression(paste(italic("\u03B4")['1'])),
                   expression(paste(italic("\u03B4")['2'])),
                   expression(paste(italic("\u03B4")['3'])),
                   
                   # gamma
                   expression(paste(italic("\u03B3"))))

ymin <- c(0, 0, -1.5, 0, 0, 0, 0)
ymax <- c(1, 1, 1.5, 1, 1, 1, 1)


# ~~~ b. Plot --------------------------------------------------------------------

for (k in 1:length(cols)) {
  df_plot_k <- df_plot %>%
    filter(parameter == cols[k])
  
  plot_title_k <- parse(text = df_plot_title[k])[[1]]
  
  assign(x = paste0("plot_", k),
         value =  ggplot(data = df_plot_k,
                         aes(x = scenario , y = mean, color = p_GPS)) + 
           geom_pointrange(aes(ymin = lower_95, ymax = upper_95),
                           position = position_dodge(width = 0.4)) +
           geom_pointrange(aes(ymin = lower_50, ymax = upper_50),
                           position = position_dodge(width = 0.4),
                           linewidth = 1) +
           geom_hline(yintercept = real_values$real_value[k]) +
           scale_color_grey(start = 0.8, end = 0.2) +
           scale_y_continuous(limits = c(ymin[k], ymax[k])) +
           theme_bw() +
           theme(plot.title = element_text(hjust = 0.5, size = 10),
                 legend.title = element_text(hjust = 0.5),
                 axis.title.y = element_blank()) +
           labs(x = "scenario", y = bquote(.(plot_title_k)), 
                title = bquote(.(plot_title_k)),
                color = expression(kappa))
  )
}

full_plot <- plot_1 + plot_2 + plot_3 + plot_4 + plot_5 + plot_6 + plot_7 +
  plot_layout(guides = 'collect')

ggsave(full_plot,
       filename = "03_outputs/Figure S7.png",
       width = 18, height = 20, units = "cm", dpi = 600)






# ~ Fig. S8. Performance model polar bear intercepts ---------------------------

{beta_phi <- NULL
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

p <- 0.25
}

path <- "07_results/model_outputs/CR_GPS_GLS/simulations_first_sub/0.0_base/"
file_list <- sort(list.files(path,
                             pattern = ".RData"))
n_sim <- length(file_list)

outputs_processed <- list() 
for (l in 1:n_sim) {
  load(paste0(path, file_list[l]))
  mcmc_output <- list(get(substr(file_list[l], 1, nchar(file_list[l])-6))[[1]],
                      get(substr(file_list[l], 1, nchar(file_list[l])-6))[[2]]) 
  outputs_processed[[l]] <- do.call(rbind, mcmc_output)
}


# First process data for plots themselves ++++++++++++++++++++++++++++++++++++++
cols_bt <- c("phi2_4", "phi5_19", "phi20", 
             "eta_1", "eta_2", "eta_3", "eta_4", 
             "beta_1", "beta_2", "beta_3", 
             "gamma_1", "gamma_2", "gamma_3",
             "s0_1", "s0_2", "s0_3", 
             "s1_1", "s1_2", "s1_3",
             "p") 

real_values_bt <- data.frame(parameter = factor(cols_bt, levels = cols_bt),
                             real_value = c(plogis(beta_phi[1]),  plogis(beta_phi[2]), plogis(beta_phi[3]),
                                            plogis(beta_eta[1]), plogis(beta_eta[2]), plogis(beta_eta[3]), plogis(beta_eta[4]), 
                                            plogis(beta_beta[1]), plogis(beta_beta[2]), plogis(beta_beta[3]),
                                            plogis(beta_gamma[1]), plogis(beta_gamma[2]), plogis(beta_gamma[3]),
                                            plogis(beta_s0[1]), plogis(beta_s0[2]), plogis(beta_s0[3]), 
                                            plogis(beta_s1[1]), plogis(beta_s1[2]), plogis(beta_s1[3]),
                                            p))

all_outputs_processed <- as.data.frame(cbind(do.call(rbind, outputs_processed), 
                                             rep(1:n_sim, each = 2*nrow(get(substr(file_list[l], 1, nchar(file_list[l])-6))$chain1)),
                                             rep(1:nrow(get(substr(file_list[l], 1, nchar(file_list[l])-6))$chain1), times = n_sim))) %>%
  rename(simulation = V32, iteration = V33) %>%
  janitor::clean_names() %>%
  mutate(phi2_4 = plogis(beta_phi_1),
         phi5_19 = plogis(beta_phi_2),
         phi20 = plogis(beta_phi_3),
         s0_1 = plogis(beta_s0_1),
         s0_2 = plogis(beta_s0_2),
         s0_3 = plogis(beta_s0_3),
         s1_1 = plogis(beta_s1_1),
         s1_2 = plogis(beta_s1_2),
         s1_3 = plogis(beta_s1_3),
         eta_1 = plogis(beta_eta_1),
         eta_2 = plogis(beta_eta_2),
         eta_3 = plogis(beta_eta_3),
         eta_4 = plogis(beta_eta_4),
         beta_1 = plogis(beta_beta_1),
         beta_2 = plogis(beta_beta_2),
         beta_3 = plogis(beta_beta_3),
         gamma_1 = plogis(beta_gamma_1),
         gamma_2 = plogis(beta_gamma_2),
         gamma_3 = plogis(beta_gamma_3)) %>%
  dplyr::select(-beta_phi_1, -beta_phi_2, -beta_phi_2, -beta_s0_1, -beta_s0_2, -beta_s0_3, 
                -beta_s1_1, -beta_s1_2, -beta_s1_3,
                -beta_eta_1, -beta_eta_2, -beta_eta_3, -beta_eta_4,
                -beta_beta_1, -beta_beta_2, -beta_beta_3,
                -beta_gamma_1, -beta_gamma_2, -beta_gamma_3) %>%
  pivot_longer(cols = all_of(cols_bt),
               names_to = "parameter")

df_caterpillar <- all_outputs_processed %>%
  group_by(simulation, parameter) %>%
  left_join(x = .,
            y = real_values_bt,
            by = "parameter") %>%
  mutate(parameter = factor(parameter, levels = cols_bt),
         to_remove = ifelse(iteration %% 2 == 0, 1, 0)) %>%
  filter(to_remove == 0) %>%
  dplyr::select(-to_remove)


# Process the data for the RB, RMSE and coverage values ++++++++++++++++++++++++
df_performance <- all_outputs_processed %>%
  # Calculate RB and RMSE
  group_by(simulation, parameter) %>%
  summarize(mean = mean(value),
            median = median(value)) %>%
  left_join(x = .,
            y = real_values_bt,
            by = "parameter") %>%
  mutate(diff = median - real_value,  
         RB = 100 * diff/real_value, # Calculate RB of each posterior for each scenario
         RMSE = sqrt(diff^2)) %>% # Calculate RMSE of each posterior for each scenario (= abs(diff))
  group_by(parameter) %>%
  summarize(RB = round( mean(RB), 1 ),
            RMSE = round( mean(RMSE), 3 )) %>% 
  # Add coverage
  left_join(x = .,
            y = all_outputs_processed %>%
              mutate(to_remove = ifelse(iteration %% 2 == 0, 1, 0)) %>%
              filter(to_remove == 0) %>%
              group_by(simulation, parameter) %>%
              summarize(lower = quantile(value, 0.025),
                        upper = quantile(value, 0.975)) %>%
              left_join(x = .,
                        y = real_values_bt,
                        by = "parameter") %>%
              mutate(contains = ifelse(real_value > lower & real_value < upper, 1, 0)) %>%
              ungroup() %>%
              count(parameter, contains) %>%
              group_by(parameter) %>%
              mutate(coverage = round(100*n / sum(n))) %>%    # Calculate the coverage
              filter(contains == 1),
            by = "parameter")


df_plot_title <- c(expression(paste(italic("\u03C6")['2-4yr'])),
                   expression(paste(italic("\u03C6")['5-14yr'])),
                   expression(paste(italic("\u03C6")['\u226515yr'])),
                   # Eta
                   expression(paste(italic("\u03B7")['4yr'])),
                   expression(paste(italic("\u03B7")['5-8yr'])),
                   expression(paste(italic("\u03B7")['9-14yr'])),
                   expression(paste(italic("\u03B7")['\u226515yr'])),
                   # Beta
                   expression(paste(italic("\u03B2")['4-8yr'])),
                   expression(paste(italic("\u03B2")['9-14yr'])),
                   expression(paste(italic("\u03B2")['\u226515yr'])),
                   # Gamma
                   expression(paste(italic("\u03B3")['4-8yr'])),
                   expression(paste(italic("\u03B3")['9-14yr'])),
                   expression(paste(italic("\u03B3")['\u226515yr'])),
                   # s0
                   expression(paste(italic("s")['0 4-8yr'])),
                   expression(paste(italic("s")['0 9-14yr'])),
                   expression(paste(italic("s")['0 \u226515yr'])),
                   # s1
                   expression(paste(italic("s")['1 4-8yr'])),
                   expression(paste(italic("s")['1 9-14yr'])),
                   expression(paste(italic("s")['1 \u226515yr'])),
                   # p
                   expression(italic("p")))

for (k in 1:length(cols_bt)) {
  df_plot_k <- df_caterpillar %>%
    filter(parameter == cols_bt[k])
  df_performance_k <- df_performance %>%
    filter(parameter == cols_bt[k])
  
  plot_title_k <- parse(text = df_plot_title[k])[[1]]
  
  assign(x = paste0("plot_", k),
         value =  ggplot(data = df_plot_k,
                         aes(x = simulation , y = value)) +
           geom_hline(aes(yintercept = real_value), data = df_plot_k) +
           stat_summary(fun.data = get_mean_and_CI,
                        fun.args = list(lower = 0.025, upper = 0.975),
                        geom = "pointrange", size = 0.05, linewidth = 0.2, alpha = 0.75) +
           theme_bw() +
           theme(axis.title = element_blank(), axis.text.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 plot.title = element_text(hjust = 0.5, size = 10),
                 plot.subtitle = element_text(hjust = 0.5, size = 7)) +
           scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
           coord_flip() +
           labs(title = bquote(.(plot_title_k)),
                subtitle = paste0("RB = ", df_performance_k$RB, "%\nRMSE = ",  df_performance_k$RMSE,
                                  "\ncoverage = ",  df_performance_k$coverage, "%"))
  )
}


full_plot <- plot_1 + plot_2 + plot_3 + plot_4 + plot_5 + plot_6 + plot_7 + plot_8 +
  plot_9 + plot_10 + plot_11 + plot_12 + plot_13 + plot_14 + plot_15 + plot_16 +
  plot_17 + plot_18 + plot_19 + plot_20 

ggsave(full_plot,
       filename = "03_outputs/Figure S8.png",
       width = 18, height = 20, units = "cm", dpi = 600)



# ~ Fig. S9. Performance model polar bear additive parameters -----------------------

{beta_phi <- NULL
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
}


path <- "03_outputs/simulations/model_outputs/polar_bear_all_covariates/"
file_list <- sort(list.files(path,
                             pattern = ".RData"))
n_sim <- length(file_list)

outputs_processed <- list() 
for (l in 1:n_sim) {
  load(paste0(path, file_list[l]))
  mcmc_output <- list(get(substr(file_list[l], 1, nchar(file_list[l])-6))[[1]],
                      get(substr(file_list[l], 1, nchar(file_list[l])-6))[[2]]) 
  outputs_processed[[l]] <- do.call(rbind, mcmc_output)
}


# ~~~ a. Additive intercepts -----------------------------------------------------

# First process data for plots themselves ++++++++++++++++++++++++++++++++++++++
cols_bt <- c("beta_eta_5", "beta_eta_6", "beta_eta_7", "beta_eta_8",
             "beta_s0_4") 

real_values_bt <- data.frame(parameter = factor(cols_bt, levels = cols_bt),
                             real_value = c(beta_eta[5],  beta_eta[6], beta_eta[7], beta_eta[8], 
                                            beta_s0[4]))


df_caterpillar <- as.data.frame(cbind(do.call(rbind, outputs_processed), 
                                      rep(1:n_sim, each = 2*nrow(get(substr(file_list[l], 1, nchar(file_list[l])-6))$chain1)),
                                      rep(1:nrow(get(substr(file_list[l], 1, nchar(file_list[l])-6))$chain1), times = n_sim))) %>%
  rename(simulation = V32, iteration = V33) %>%
  janitor::clean_names() %>%
  dplyr::select(all_of(c(c("simulation", "iteration"),
                         cols_bt))) %>%
  pivot_longer(cols = all_of(cols_bt),
               names_to = "parameter") %>%
  group_by(simulation, parameter) %>%
  left_join(x = .,
            y = real_values_bt,
            by = "parameter") %>%
  mutate(parameter = factor(parameter, levels = cols_bt)) 



# First the raw RB, RMSE and coverage values +++++++++++++++++++++++++++++++++++
cols_bt <- c("beta_eta_5", "beta_eta_6", "beta_eta_7", "beta_eta_8",
             "beta_s0_4") 
df_performance_raw <- df_caterpillar %>%
  group_by(simulation, parameter) %>%
  summarize(mean = mean(value),
            median = median(value)) %>%
  left_join(x = .,
            y = real_values_bt,
            by = "parameter") %>%
  mutate(diff = median - real_value,  
         RB = 100 * diff/real_value, # Calculate RB of each posterior for each scenario
         RMSE = sqrt(diff^2)) %>% # Calculate RMSE of each posterior for each scenario (= abs(diff))
  group_by(parameter) %>%
  summarize(RB = round( mean(RB), 1 ),
            RMSE = round( mean(RMSE), 3 ),
            abs_RB = abs(RB)) %>% 
  # Add coverage
  left_join(x = .,
            y = df_caterpillar %>%
              group_by(simulation, parameter) %>%
              summarize(lower = quantile(value, 0.025),
                        upper = quantile(value, 0.975)) %>%
              left_join(x = .,
                        y = real_values_bt,
                        by = "parameter") %>%
              mutate(contains = ifelse(real_value > lower & real_value < upper, 1, 0)) %>%
              ungroup() %>%
              count(parameter, contains) %>%
              group_by(parameter) %>%
              mutate(coverage = round(100*n / sum(n))) %>%    # Calculate the coverage
              filter(contains == 1),
            by = "parameter") %>% 
  dplyr::select(parameter, RB, RMSE, coverage)



# Now the RB, RMSE and coverage values when added to intercepts ++++++++++++++++

cols_bt_perf <- c("s0_1_s0_4", "s0_2_s0_4", "s0_3_s0_4", 
                  "eta_2_eta_5", "eta_3_eta_5", "eta_4_eta_5", 
                  "eta_2_eta_6", "eta_3_eta_6", "eta_4_eta_6", 
                  "eta_2_eta_7", "eta_3_eta_7", "eta_4_eta_7",
                  "eta_2_eta_8", "eta_3_eta_8", "eta_4_eta_8") 

real_values_bt <- data.frame(parameter = factor(cols_bt_perf, levels = cols_bt_perf),
                             real_value = c(plogis(beta_s0[1]+beta_s0[4]), plogis(beta_s0[2]+beta_s0[4]), plogis(beta_s0[3]+beta_s0[4]),
                                            plogis(beta_eta[2]+beta_eta[5]), plogis(beta_eta[3]+beta_eta[5]), plogis(beta_eta[4]+beta_eta[5]),
                                            plogis(beta_eta[2]+beta_eta[6]), plogis(beta_eta[3]+beta_eta[6]), plogis(beta_eta[4]+beta_eta[6]),
                                            plogis(beta_eta[2]+beta_eta[7]), plogis(beta_eta[3]+beta_eta[7]), plogis(beta_eta[4]+beta_eta[7]),
                                            plogis(beta_eta[2]+beta_eta[8]), plogis(beta_eta[3]+beta_eta[8]), plogis(beta_eta[4]+beta_eta[8])))

all_output_processed <- as.data.frame(cbind(do.call(rbind, outputs_processed), 
                                            rep(1:n_sim, each = 2*nrow(get(substr(file_list[l], 1, nchar(file_list[l])-6))$chain1)),
                                            rep(1:nrow(get(substr(file_list[l], 1, nchar(file_list[l])-6))$chain1), times = n_sim))) %>%
  rename(simulation = V32, iteration = V33) %>%
  janitor::clean_names() %>%
  mutate(s0_1_s0_4 = plogis(beta_s0[1] + beta_s0_4),
         s0_2_s0_4 = plogis(beta_s0[2] + beta_s0_4),
         s0_3_s0_4 = plogis(beta_s0[3] + beta_s0_4),
         eta_2_eta_5 = plogis(beta_eta[2] + beta_eta_5),
         eta_3_eta_5 = plogis(beta_eta[3] + beta_eta_5),
         eta_4_eta_5 = plogis(beta_eta[4] + beta_eta_5),
         eta_2_eta_6 = plogis(beta_eta[2] + beta_eta_6),
         eta_3_eta_6 = plogis(beta_eta[3] + beta_eta_6),
         eta_4_eta_6 = plogis(beta_eta[4] + beta_eta_6),
         eta_2_eta_7 = plogis(beta_eta[2] + beta_eta_7),
         eta_3_eta_7 = plogis(beta_eta[3] + beta_eta_7),
         eta_4_eta_7 = plogis(beta_eta[4] + beta_eta_7),
         eta_2_eta_8 = plogis(beta_eta[2] + beta_eta_8),
         eta_3_eta_8 = plogis(beta_eta[3] + beta_eta_8),
         eta_4_eta_8 = plogis(beta_eta[4] + beta_eta_8)) %>%
  dplyr::select(all_of(c(c("simulation", "iteration"),
                         cols_bt_perf))) %>%
  pivot_longer(cols = all_of(cols_bt_perf),
               names_to = "parameter")


df_performance_temp <- all_output_processed %>%
  group_by(simulation, parameter) %>%
  summarize(mean = mean(value),
            median = median(value)) %>%
  left_join(x = .,
            y = real_values_bt,
            by = "parameter") %>%
  mutate(diff = median - real_value,  
         RB = 100 * diff/real_value, # Calculate RB of each posterior for each scenario
         RMSE = sqrt(diff^2)) %>% # Calculate RMSE of each posterior for each scenario (= abs(diff))
  group_by(parameter) %>%
  summarize(RB = round( mean(RB), 1 ),
            RMSE = round( mean(RMSE), 3 ),
            abs_RB = abs(RB)) %>% 
  # Add coverage
  left_join(x = .,
            y = all_output_processed %>%
              group_by(simulation, parameter) %>%
              summarize(lower = quantile(value, 0.025),
                        upper = quantile(value, 0.975)) %>%
              left_join(x = .,
                        y = real_values_bt,
                        by = "parameter") %>%
              mutate(contains = ifelse(real_value > lower & real_value < upper, 1, 0)) %>%
              ungroup() %>%
              count(parameter, contains) %>%
              group_by(parameter) %>%
              mutate(coverage = round(100*n / sum(n))) %>%    # Calculate the coverage
              filter(contains == 1),
            by = "parameter") %>%
  mutate(parameter_comb = parameter,
         parameter = ifelse(parameter %in% c("s0_1_s0_4", "s0_2_s0_4", "s0_3_s0_4"), "beta_s0_4",
                            ifelse(parameter %in% c("eta_2_eta_5", "eta_3_eta_5", "eta_4_eta_5"), "beta_eta_5",
                                   ifelse(parameter %in% c("eta_2_eta_6", "eta_3_eta_6", "eta_4_eta_6"), "beta_eta_6",
                                          ifelse(parameter %in% c("eta_2_eta_7", "eta_3_eta_7", "eta_4_eta_7"), "beta_eta_7",
                                                 "beta_eta_8")))))


# Now let's create a character string with the two most extreme RB/RMSE/coverage values
# As the interval
df_performance_intervals <- df_performance_temp %>%
  filter(parameter_comb %in% c("s0_1_s0_4", "s0_2_s0_4", "s0_3_s0_4", 
                               "eta_2_eta_5", "eta_3_eta_5", "eta_4_eta_5", 
                               "eta_2_eta_6", "eta_3_eta_6", "eta_4_eta_6", 
                               "eta_2_eta_7", "eta_3_eta_7", "eta_4_eta_7",
                               "eta_2_eta_8", "eta_3_eta_8", "eta_4_eta_8")) %>%
  mutate(parameter = ifelse(parameter_comb %in% c("s0_1_s0_4", "s0_2_s0_4", "s0_3_s0_4"), "beta_s0_4",
                            ifelse(parameter_comb %in% c("eta_2_eta_5", "eta_3_eta_5", "eta_4_eta_5"), "beta_eta_5",
                                   ifelse(parameter_comb %in% c("eta_2_eta_6", "eta_3_eta_6", "eta_4_eta_6"), "beta_eta_6",
                                          ifelse(parameter_comb %in% c("eta_2_eta_7", "eta_3_eta_7", "eta_4_eta_7"), "beta_eta_7",
                                                 "beta_eta_8"))))) %>%
  group_by(parameter) %>%
  summarize(max_RB = max(RB),
            min_RB = min(RB),
            max_RMSE = max(RMSE),
            min_RMSE = min(RMSE),
            max_coverage = max(coverage),
            min_coverage = min(coverage))  %>%
  mutate(RB_final = paste0(min_RB,  "; ", max_RB),
         RMSE_final = paste0(min_RMSE,  "; ", max_RMSE),
         coverage_final = paste0(min_coverage,  "; ", max_coverage)) %>%
  dplyr::select(parameter, RB_final, RMSE_final, coverage_final)


df_performance_final <- df_performance_raw %>%
  left_join(x = .,
            y = df_performance_intervals,
            by = "parameter") %>%
  mutate(RB_final = ifelse(is.na(RB_final), RB, RB_final),
         RMSE_final = ifelse(is.na(RMSE_final), RMSE, RMSE_final)) %>%
  filter(parameter %in% cols_bt) 


# Add RB and RMSE values to the dataframe for the plots
df_plot_additive <- df_caterpillar %>%
  left_join(x = .,
            y = df_performance_final,
            by = "parameter")



# ~~~ b. Slopes --------------------------------------------------------------------

cols_bt <- c("beta_phi_4", "beta_eta_9", 
             "beta_beta_4","beta_gamma_4",
             "beta_s0_5", "beta_s1_4") 

real_values_bt <- data.frame(parameter = factor(cols_bt, levels = cols_bt),
                             real_value = c(beta_phi[4], beta_eta[9], 
                                            beta_beta[4], beta_gamma[4],
                                            beta_s0[5], beta_s1[4]))

df_caterpillar <- as.data.frame(cbind(do.call(rbind, outputs_processed), 
                                      rep(1:n_sim, each = 2*nrow(get(substr(file_list[l], 1, nchar(file_list[l])-6))$chain1)),
                                      rep(1:nrow(get(substr(file_list[l], 1, nchar(file_list[l])-6))$chain1), times = n_sim))) %>%
  rename(simulation = V32, iteration = V33) %>%
  janitor::clean_names() %>%
  dplyr::select(all_of(c(c("simulation", "iteration"),
                         cols_bt))) %>%
  pivot_longer(cols = all_of(cols_bt),
               names_to = "parameter") %>%
  group_by(simulation, parameter) %>%
  left_join(x = .,
            y = real_values_bt,
            by = "parameter") %>%
  mutate(parameter = factor(parameter, levels = cols_bt))


# Process the data for the RB, RMSE and coverage values ++++++++++++++++++++++++
df_performance_slopes <- df_caterpillar %>%
  # Calculate RB and RMSE
  group_by(simulation, parameter) %>%
  summarize(mean = mean(value),
            median = median(value)) %>%
  left_join(x = .,
            y = real_values_bt,
            by = "parameter") %>%
  mutate(diff = median - real_value,  
         RB = sign(real_value) * 100 * diff/real_value, # Calculate RB of each posterior for each scenario
         RMSE = sqrt(diff^2)) %>% # Calculate RMSE of each posterior for each scenario (= abs(diff))
  group_by(parameter) %>%
  summarize(RB = round( mean(RB), 1 ),
            RMSE = round( mean(RMSE), 3 ),
            abs_RB = abs(RB)) %>% 
  # Add coverage
  left_join(x = .,
            y = df_caterpillar %>%
              mutate(to_remove = ifelse(iteration %% 2 == 0, 1, 0)) %>%
              filter(to_remove == 0) %>%
              group_by(simulation, parameter) %>%
              summarize(lower = quantile(value, 0.025),
                        upper = quantile(value, 0.975)) %>%
              left_join(x = .,
                        y = real_values_bt,
                        by = "parameter") %>%
              mutate(contains = ifelse(real_value > lower & real_value < upper, 1, 0)) %>%
              ungroup() %>%
              count(parameter, contains) %>%
              group_by(parameter) %>%
              mutate(coverage = round(100*n / sum(n))) %>%    # Calculate the coverage
              filter(contains == 1),
            by = "parameter")

df_plot_slopes <- df_caterpillar %>%
  left_join(x = .,
            y = df_performance_slopes,
            by = "parameter") %>%
  dplyr::select(-abs_RB, -contains, -n) %>%
  mutate(RB = as.character(RB),
         RMSE = as.character(RMSE),
         coverage = as.character(coverage))

df_plot_title <- c(expression(paste(italic("\u03C6")[sea~ice])),
                   # Eta
                   expression(paste(italic("\u03B7")[A[FD]])),
                   expression(paste(italic("\u03B7")[A['0-']])),
                   expression(paste(italic("\u03B7")[A['1-']])),
                   expression(paste(italic("\u03B7")[A[S]])),
                   expression(paste(italic("\u03B7")[sea~ice])),
                   # Beta
                   expression(paste(italic("\u03B2")[sea~ice])),
                   # Gamma
                   expression(paste(italic("\u03B3")[sea~ice])),
                   # s0
                   expression(paste(italic("s")['0 singleton'])),
                   expression(paste(italic("s")['0 sea ice'])),
                   # s1
                   expression(paste(italic("s")['1 sea ice'])))

cols_bt <- c("beta_phi_4", 
             "beta_eta_5", "beta_eta_6", "beta_eta_7", "beta_eta_8", "beta_eta_9", 
             "beta_beta_4","beta_gamma_4",
             "beta_s0_4", "beta_s0_5", "beta_s1_4") 
additive <- c(0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0)

for (k in 1:length(cols_bt)) {
  if (additive[k] == 1) {
    df_plot_k <- df_plot_additive %>%
      filter(parameter == cols_bt[k])
    df_performance_final_k <- df_plot_k %>%
      distinct(RB, RB_final, RMSE, RMSE_final, coverage, coverage_final)
    
    plot_title_k <- parse(text = df_plot_title[k])[[1]]
    
    assign(x = paste0("plot_", k),
           value =  ggplot(data = df_plot_k,
                           aes(x = simulation , y = value)) +
             geom_hline(aes(yintercept = real_value), data = df_plot_k) +
             stat_summary(fun.data = get_median_and_CI,
                          fun.args = list(lower = 0.025, upper = 0.975),
                          geom = "pointrange", size = 0.05, linewidth = 0.2, alpha = 0.75) +
             theme_bw() +
             theme(axis.title = element_blank(), axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   plot.title = element_text(hjust = 0.5, size = 10),
                   plot.subtitle = element_text(hjust = 0.5, size = 7)) +
             coord_flip() +
             labs(title = bquote(.(plot_title_k)),
                  subtitle = paste0("RB = [", df_performance_final_k$RB_final,  "]% (", df_performance_final_k$RB, "%)", 
                                    # "RB = ", df_performance_final_k$RB_final,  "% (", df_performance_final_k$RB, "%)", 
                                    "\nRMSE = [",  df_performance_final_k$RMSE_final, "] (", df_performance_final_k$RMSE, ")",
                                    "\ncoverage = ",  df_performance_final_k$coverage, "%"))
    )
  } else {
    df_plot_k <- df_plot_slopes %>%
      filter(parameter == cols_bt[k])
    df_performance_final_k <- df_plot_k %>%
      distinct(RB, RMSE, coverage)
    
    plot_title_k <- parse(text = df_plot_title[k])[[1]]
    
    assign(x = paste0("plot_", k),
           value =  ggplot(data = df_plot_k,
                           aes(x = simulation , y = value)) +
             geom_hline(aes(yintercept = real_value), data = df_plot_k) +
             stat_summary(fun.data = get_median_and_CI,
                          fun.args = list(lower = 0.025, upper = 0.975),
                          geom = "pointrange", size = 0.05, linewidth = 0.2, alpha = 0.75) +
             theme_bw() +
             theme(axis.title = element_blank(), axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   plot.title = element_text(hjust = 0.5, size = 10),
                   plot.subtitle = element_text(hjust = 0.5, size = 7)) +
             coord_flip() +
             labs(title = bquote(.(plot_title_k)),
                  subtitle = paste0("RB = ", df_performance_final_k$RB, "%",
                                    "\nRMSE = ",  df_performance_final_k$RMSE,
                                    "\ncoverage = ",  df_performance_final_k$coverage, "%"))
    )
  }
  
}

full_plot <- plot_1 + plot_2 + plot_3 + plot_4 + plot_5 + plot_6 + plot_7 + 
  plot_8 + plot_9 + plot_10 + plot_11

ggsave(full_plot,
       filename = "03_outputs/Figure S9.png",
       width = 18, height = 16, units = "cm", dpi = 600)



# ~ Fig. S10. Posteriors model with biologging ---------------------------------

load("03_outputs/case_study/fit_with_biologging_1.RData")
load("03_outputs/case_study/fit_with_biologging_2.RData")
res <- rbind(fit_with_biologging_1, fit_with_biologging_2) 

name <- c("beta_phi[1]", "beta_phi[2]", "beta_phi[3]", "beta_phi[4]", "beta_phi[5]",
          "beta_eta[1]", "beta_eta[2]", "beta_eta[3]", "beta_eta[4]", "beta_eta[5]", "beta_eta[6]", "beta_eta[7]", "beta_eta[8]", "beta_eta[9]",
          "beta_beta[1]", "beta_beta[2]", "beta_beta[3]", "beta_beta[4]", "beta_beta[5]",
          "beta_gamma[1]", "beta_gamma[2]", "beta_gamma[3]", "beta_gamma[4]", "beta_gamma[5]", "beta_gamma[6]", 
          "beta_s0[1]", "beta_s0[2]", "beta_s0[3]", "beta_s0[4]", "beta_s0[5]", "beta_s0[6]", 
          "beta_s1[1]", "beta_s1[2]", "beta_s1[3]", "beta_s1[4]",
          "beta_p[1]", "beta_p[2]", "beta_p[3]", "beta_p[4]", "sigma_omega")

model_df <- data.frame(name = name,
                       new_name = c("a[1]", "a[2]", "a[3]", "a[4]", "a[5]",
                                    "b[1]", "b[2]", "b[3]", "b[4]", "b[5]", "b[6]", "b[7]", "b[8]", "b[9]",
                                    "c[1]", "c[2]", "c[3]", "c[4]", "c[5]",
                                    "d[1]", "d[2]", "d[3]", "d[4]", "d[5]", "d[6]",
                                    "e[1]", "e[2]", "e[3]", "e[4]", "e[5]", "e[6]",
                                    "f[1]", "f[2]", "f[3]", "f[4]",
                                    "g[1]", "g[2]", "g[3]", "g[4]", "sigma[omega]"),
                       covariate = c("survival 2-4", "survival 5-15", "survival 16+", "space-use unknown-pelagic", "sea ice", 
                                     "denning prob 4yr", "denning prob", "denning prob", "denning prob", "FDA", "A0-", "A1-", "SA", "sea ice",
                                     "early litter survival", "early litter survival", "early litter survival", "DOY", "sea ice",
                                     "twinning", "twinning", "twinning", "DOY", "space-use unknown-pelagic", "sea ice",
                                     "CoY survival", "CoY survival", "CoY survival", "singleton", "space-use unknown-pelagic", "sea ice",
                                     "yrl survival", "yrl survival", "yrl survival", "sea ice",
                                     "p resident", "p unk", "p pelagic", "p time pel", "sigma omega"),
                       parameter = c(rep("phi", times = 5),
                                     rep("eta", times = 9),
                                     rep("beta", times = 5),
                                     rep("gamma", times = 6),
                                     rep("s0", times = 6),
                                     rep("s1", times = 4),
                                     rep("p", times = 5)),
                       parameter2 = c(rep("survival", times = 5),
                                      rep("denning\nprobability", times = 9),
                                      rep("early litter\nsurvival", times = 5),
                                      rep("twinning\nprobability", times = 6),
                                      rep("CoY\nsurvival", times = 6),
                                      rep("yearling\nsurvival", times = 4),
                                      rep("recapture\nprobability", times = 5)))

cols <- model_df$name
caterpillar <- as.data.frame(res) %>%
  pivot_longer(cols = all_of(cols)) %>%
  left_join(x = .,
            y = model_df,
            by = "name") %>%
  mutate(parameter2 = factor(parameter2, levels = c("survival",
                                                    "denning\nprobability", 
                                                    "early litter\nsurvival", 
                                                    "twinning\nprobability", 
                                                    "CoY\nsurvival",
                                                    "yearling\nsurvival", 
                                                    "recapture\nprobability")),
         new_name = factor(new_name, levels = rev(c("a[1]", "a[2]", "a[3]", "a[4]", "a[5]",
                                                    "b[1]", "b[2]", "b[3]", "b[4]", "b[5]", "b[6]", "b[7]", "b[8]", "b[9]",
                                                    "c[1]", "c[2]", "c[3]", "c[4]", "c[5]",
                                                    "d[1]", "d[2]", "d[3]", "d[4]", "d[5]", "d[6]",
                                                    "e[1]", "e[2]", "e[3]", "e[4]", "e[5]", "e[6]",
                                                    "f[1]", "f[2]", "f[3]", "f[4]",
                                                    "g[1]", "g[2]", "g[3]", "g[4]", "sigma[omega]")))) %>%
  dplyr::select(new_name, name, value, covariate, parameter, parameter2) 

ggplot(caterpillar, 
       aes(x = value, y = new_name)) +
  geom_vline(xintercept = 0) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975),
               geom = "pointrange", orientation = "y", 
               size = 0.35) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975),
               geom = "pointrange", orientation = "y", 
               size = 0.35) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.25, upper = 0.75),
               geom = "pointrange", orientation = "y", 
               size = 0.35, linewidth = 1) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4, 6)) +
  scale_y_discrete(labels = function(labels) {
    valid_labels <- rev(model_df$new_name[model_df$new_name %in% labels])
    parse(text = valid_labels)
  }) +
  scale_color_brewer(palette = "Dark2", na.value = "grey40") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic")) +
  labs(x = "posterior distribution") +
  facet_grid(parameter2~., 
             space = "free",
             scales = "free")

ggsave(full_plot,
       filename = "03_outputs/Figure S10.png",
       width = 18, height = 20, units = "cm", dpi = 600)


# ~ Fig. S11. Age VS reproductive rates ---------------------------------------------

custom_palette_4 <- c("#FCBF49", "#F77F00", "#D62828", "#003049")
custom_palette_3 <- c("#F99418", "#D62828", "#003049")

load("03_outputs/case_study/fit_with_biologging_1.RData")
load("03_outputs/case_study/fit_with_biologging_2.RData")
res <- rbind(fit_with_biologging_1, fit_with_biologging_2) 

# Denning probability ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Age4 <- c(1, 0, 0, 0) ; AgeY <- c(0, 1, 0, 0) ; AgeM <- c(0, 0, 1, 0) ; AgeO <- c(0, 0, 0, 1)
probability <- matrix(NA, nrow = dim(res)[1], ncol = 4,
                      dimnames = list(1:dim(res)[1], c("4yr", "5-8yr", "9-14yr", "\u226515yr")))
for (i in 1:4) {              # age
  for (j in 1:dim(res)[1]) {       # MCMC iteration
    probability[j, i] <- plogis(res[j, "beta_eta[1]"] * Age4[i] +
                                  res[j, "beta_eta[2]"] * AgeY[i] + 
                                  res[j, "beta_eta[3]"] * AgeM[i] +
                                  res[j, "beta_eta[4]"] * AgeO[i])
    
  }
}

df_plot <- as.data.frame.table(probability) %>%
  rename(iteration = Var1, age_group = Var2, value = Freq)

plot_denning <- ggplot(data = df_plot,
                       aes(x = age_group, y = value, fill = age_group)) +
  geom_violin(alpha = 0.8) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975),
               geom = "pointrange", size = 0.4, color = "grey", shape = 21) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.25, upper = 0.75),
               geom = "pointrange", size = 0.5, stroke = 0.4, linewidth = 1, shape = 21, 
               color = "grey") +  
  theme_bw() +
  theme(legend.position = "none") +
  lims(y = c(0, 1)) +
  # scale_fill_brewer(palette = "Dark2") +
  scale_fill_manual(values = custom_palette_4) +
  labs(x = "age at mating", 
       y = expression(paste("denning probability (", italic("\u03B7"), ")")))



# Early cub survival +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

AgeY <- c(1, 0, 0) ; AgeM <- c(0, 1, 0) ; AgeO <- c(0, 0, 1)
probability <- matrix(NA, nrow = dim(res)[1], ncol = 3,
                      dimnames = list(1:dim(res)[1], c("4-8yr", "9-14yr", "\u226515yr")))
for (i in 1:3) {              # age
  for (j in 1:dim(res)[1]) {       # MCMC iteration
    probability[j, i] <- plogis(res[j, "beta_beta[1]"] * AgeY[i] + 
                                  res[j, "beta_beta[2]"] * AgeM[i] +
                                  res[j, "beta_beta[3]"] * AgeO[i])
    
  }
}

df_plot <- as.data.frame.table(probability) %>%
  rename(iteration = Var1, age_group = Var2, value = Freq)

plot_early_survival <- ggplot(data = df_plot,
                              aes(x = age_group, y = value, fill = age_group)) +
  geom_violin(alpha = 0.8) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975),
               geom = "pointrange", size = 0.4, color = "grey", shape = 21) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.25, upper = 0.75),
               geom = "pointrange", size = 0.5, stroke = 0.4, linewidth = 1, shape = 21, 
               color = "grey") + 
  theme_bw() +
  theme(legend.position = "none") +
  lims(y = c(0, 1)) +
  scale_fill_manual(values = custom_palette_3) +
  labs(x = "age at mating", y = expression(paste("early litter survival (", italic("\u03B2"), ")")))


get_probability_direction(res[, "beta_beta[3]"] - res[, "beta_beta[2]"])
get_probability_direction(res[, "beta_beta[1]"] - res[, "beta_beta[2]"])



# Litter size ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

AgeY <- c(1, 0, 0) ; AgeM <- c(0, 1, 0) ; AgeO <- c(0, 0, 1)
probability <- matrix(NA, nrow = dim(res)[1], ncol = 3,
                      dimnames = list(1:dim(res)[1], c("4-8yr", "9-14yr", "\u226515yr")))
for (i in 1:3) {              # age
  for (j in 1:dim(res)[1]) {       # MCMC iteration
    probability[j, i] <- plogis(res[j, "beta_gamma[1]"] * AgeY[i] + 
                                  res[j, "beta_gamma[2]"] * AgeM[i] +
                                  res[j, "beta_gamma[3]"] * AgeO[i])
    
  }
}

df_plot <- as.data.frame.table(probability) %>%
  rename(iteration = Var1, age_group = Var2, value = Freq)

plot_litter_size <- ggplot(data = df_plot,
                           aes(x = age_group, y = value, fill = age_group)) +
  geom_violin(alpha = 0.8) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975),
               geom = "pointrange", size = 0.4, color = "grey", shape = 21) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.25, upper = 0.75),
               geom = "pointrange", size = 0.5, stroke = 0.4, linewidth = 1, shape = 21, 
               color = "grey") + 
  theme_bw() +
  theme(legend.position = "none") +
  lims(y = c(0, 1)) +
  scale_fill_manual(values = custom_palette_3) +
  labs(x = "age at mating", y = expression(paste("twinning probability (", italic("\u03B3"), ")")))


# CoY survival +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

AgeY <- c(1, 0, 0) ; AgeM <- c(0, 1, 0) ; AgeO <- c(0, 0, 1)
probability <- matrix(NA, nrow = dim(res)[1], ncol = 3,
                      dimnames = list(1:dim(res)[1], c("4-8yr", "9-14yr", "\u226515yr")))
for (i in 1:3) {              # age
  for (j in 1:dim(res)[1]) {       # MCMC iteration
    probability[j, i] <- plogis(res[j, "beta_s0[1]"] * AgeY[i] + 
                                  res[j, "beta_s0[2]"] * AgeM[i] +
                                  res[j, "beta_s0[3]"] * AgeO[i])
    
  }
}

df_plot <- as.data.frame.table(probability) %>%
  rename(iteration = Var1, age_group = Var2, value = Freq)

plot_s0 <- ggplot(data = df_plot,
                  aes(x = age_group, y = value, fill = age_group)) +
  geom_violin(alpha = 0.8) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975),
               geom = "pointrange", size = 0.4, color = "grey", shape = 21) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.25, upper = 0.75),
               geom = "pointrange", size = 0.5, stroke = 0.4, linewidth = 1, shape = 21, 
               color = "grey") + 
  theme_bw() +
  theme(legend.position = "none") +
  lims(y = c(0, 1)) +
  scale_fill_manual(values = custom_palette_3) +
  labs(x = "age at mating", y = parse(text = "CoY~survival~(italic(s)[0])"))


# Yearling survival ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

AgeY <- c(1, 0, 0) ; AgeM <- c(0, 1, 0) ; AgeO <- c(0, 0, 1)
probability <- matrix(NA, nrow = dim(res)[1], ncol = 3,
                      dimnames = list(1:dim(res)[1], c("4-8yr", "9-14yr", "\u226515yr")))
# dimnames = list(1:dim(res)[1], c("young", "middle-aged", "old")))
for (i in 1:3) {              # age
  for (j in 1:dim(res)[1]) {       # MCMC iteration
    probability[j, i] <- plogis(res[j, "beta_s1[1]"] * AgeY[i] + 
                                  res[j, "beta_s1[2]"] * AgeM[i] +
                                  res[j, "beta_s1[3]"] * AgeO[i])
    
  }
}

df_plot <- as.data.frame.table(probability) %>%
  rename(iteration = Var1, age_group = Var2, value = Freq)

plot_s1 <- ggplot(data = df_plot,
                  aes(x = age_group, y = value, fill = age_group)) +
  geom_violin(alpha = 0.8) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975),
               geom = "pointrange", size = 0.4, color = "grey", shape = 21) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.25, upper = 0.75),
               geom = "pointrange", size = 0.5, stroke = 0.4, linewidth = 1, shape = 21, 
               color = "grey") + 
  theme_bw() +
  theme(legend.position = "none") +
  lims(y = c(0, 1)) +
  scale_fill_manual(values = custom_palette_3) +
  labs(x = "age at mating", y = parse(text = "yearling~survival~(italic(s)[1])"))


layout <- "
AAAAAABBBBBFFFF
AAAAAABBBBBFFFF
AAAAAABBBBBFFFF
AAAAAABBBBBFFFF
CCCCCDDDDDEEEEE
CCCCCDDDDDEEEEE
CCCCCDDDDDEEEEE
CCCCCDDDDDEEEEE
"

(plot_denning + plot_early_survival) / 
  (plot_litter_size + plot_s0 + plot_s1) +
  plot_annotation(tag_levels = 'A') 

ggsave(full_plot,
       filename = "03_outputs/Figure S11.png",
       width = 18, height = 14, units = "cm", dpi = 600)



# ~ Fig. S12. Cost of reproduction -----------------------------------------

load("03_outputs/case_study/fit_with_biologging_1.RData")
load("03_outputs/case_study/fit_with_biologging_2.RData")
res <- rbind(fit_with_biologging_1, fit_with_biologging_2) 


FDA <- c(1, 0, 0, 0) ; A0 <- c(0, 1, 0, 0)  ; A1 <- c(0, 0, 1, 0)  ; AS <- c(0, 0, 0, 1) 
effect <- matrix(NA, nrow = dim(res)[1], ncol = 4,
                      dimnames = list(1:dim(res)[1], c("AFD", "A0", "A1", "AS")))
for (i in 1:4) {              # reproductive status
  for (j in 1:dim(res)[1]) {       # MCMC iteration
    effect[j, i] <- res[j, "beta_eta[5]"] * FDA[i] +
                                  res[j, "beta_eta[6]"] * A0[i] +
                                  res[j, "beta_eta[7]"] * A1[i] +
                                  res[j, "beta_eta[8]"] * AS[i]
    
  }
}

df_plot <- as.data.frame.table(effect) %>%
  rename(iteration = Var1, age_group = Var2, value = Freq)

(plot_additive <- ggplot(data = df_plot,
                        aes(x = age_group, y = value, fill = age_group)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_violin(alpha = 0.8) +
    stat_summary(fun.data = get_median_and_CI,
                 fun.args = list(lower = 0.025, upper = 0.975),
                 geom = "pointrange", size = 0.4, color = "grey", shape = 21) +
    stat_summary(fun.data = get_median_and_CI,
                 fun.args = list(lower = 0.25, upper = 0.75),
                 geom = "pointrange", size = 0.5, stroke = 0.4, linewidth = 1, shape = 21, 
                 color = "grey") +  
    theme_bw() +
    theme(legend.position = "none") +
    scale_x_discrete(labels = c(bquote(A['FD']),
                                bquote(A['0']), bquote(A['1']), bquote(A['S']))) +
    # lims(y = c(0, 1)) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(5, "Set1")[c(2:5)]) +
    labs(x = "state", 
         y = expression(paste("    additive effect on\ndenning probability (", italic("\u03B7"), ")"))))
# y = expression(paste("additive effect on denning probability (", italic("\u03B7"), ")"))))



FDA <- c(0, 1, 0, 0, 0) ; A0 <- c(0, 0, 1, 0, 0)  ; A1 <- c(0, 0, 0, 1, 0)  ; AS <- c(0, 0, 0, 0, 1) 
probability <- matrix(NA, nrow = dim(res)[1], ncol = 5,
                      dimnames = list(1:dim(res)[1], c("ALND", "AFD", "A0", "A1", "AS")))
for (i in 1:5) {              # reproductive status
  for (j in 1:dim(res)[1]) {       # MCMC iteration
    probability[j, i] <- plogis(res[j, "beta_eta[3]"] +
                                  res[j, "beta_eta[5]"] * FDA[i] +
                                  res[j, "beta_eta[6]"] * A0[i] +
                                  res[j, "beta_eta[7]"] * A1[i] +
                                  res[j, "beta_eta[8]"] * AS[i])
    
  }
}

df_plot <- as.data.frame.table(probability) %>%
  rename(iteration = Var1, age_group = Var2, value = Freq)

(plot_denning <- ggplot(data = df_plot,
                        aes(x = age_group, y = value, fill = age_group)) +
    geom_violin(alpha = 0.8) +
    stat_summary(fun.data = get_median_and_CI,
                 fun.args = list(lower = 0.025, upper = 0.975),
                 geom = "pointrange", size = 0.4, color = "grey", shape = 21) +
    stat_summary(fun.data = get_median_and_CI,
                 fun.args = list(lower = 0.25, upper = 0.75),
                 geom = "pointrange", size = 0.5, stroke = 0.4, linewidth = 1, shape = 21, 
                 color = "grey") +  
    theme_bw() +
    theme(legend.position = "none") +
    scale_x_discrete(labels = c(bquote(A['LND']), bquote(A['FD']),
                                bquote(A['0']), bquote(A['1']), bquote(A['S']))) +
    lims(y = c(0, 1)) +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "state", 
         y = expression(paste("denning probability (", italic("\u03B7"), ")"))))


plot_additive + plot_denning +
  plot_annotation(tag_levels = list(c("A", "B"))) 
  
ggsave(full_plot,
       filename = "03_outputs/Figure S14.png",
       width = 18, height = 7.5, units = "cm", dpi = 600)


# ~ Fig. S13. Posteriors model without biologging ------------------------------------

load("03_outputs/case_study/fit_without_biologging_1.RData")
load("03_outputs/case_study/fit_without_biologging_2.RData")
res <- rbind(fit_without_biologging_1, fit_without_biologging_2) 


name <- c("beta_phi[1]", "beta_phi[2]", "beta_phi[3]", "beta_phi[4]", "beta_phi[5]",
          "beta_beta[1]", "beta_beta[2]", "beta_beta[3]", "beta_beta[4]", "beta_beta[5]", "beta_beta[6]", "beta_beta[7]", "beta_beta[8]", 
          "beta_gamma[1]", "beta_gamma[2]", "beta_gamma[3]", "beta_gamma[4]", "beta_gamma[5]", 
          "beta_s0[1]", "beta_s0[2]", "beta_s0[3]", "beta_s0[4]", "beta_s0[5]", 
          "beta_s1[1]", "beta_s1[2]", "beta_s1[3]", "beta_s1[4]",
          "beta_p[1]", "beta_p[2]", "beta_p[3]", "beta_p[4]", "sigma_omega")

model_df <- data.frame(name = name,
                       new_name = c("a[1]", "a[2]", "a[3]", "a[4]", "a[5]",
                                    "b[1]", "b[2]", "b[3]", "b[4]", "b[5]", "b[6]", "b[7]", "b[8]", 
                                    "d[1]", "d[2]", "d[3]", "d[4]", "d[5]",
                                    "e[1]", "e[2]", "e[3]", "e[4]", "e[5]",
                                    "f[1]", "f[2]", "f[3]", "f[4]",
                                    "g[1]", "g[2]", "g[3]", "g[4]", "sigma[omega]"),
                       covariate = c("survival 2-4", "survival 5-15", "survival 16+", "space-use unknown-pelagic", "sea ice", 
                                     "breeding prob 4yr", "breeding prob", "breeding prob", "breeding prob", "A0-", "A1-", "SA", "sea ice",
                                     "twinning", "twinning", "twinning", "DOY", "sea ice",
                                     "CoY survival", "CoY survival", "CoY survival", "singleton", "sea ice",
                                     "yrl survival", "yrl survival", "yrl survival", "sea ice",
                                     "p resident", "p unk", "p pelagic", "p time pel", "sigma omega"),
                       parameter = c(rep("phi", times = 5),
                                     rep("delta", times = 8),
                                     rep("gamma", times = 5),
                                     rep("s0", times = 5),
                                     rep("s1", times = 4),
                                     rep("p", times = 5)),
                       parameter2 = c(rep("survival", times = 5),
                                      rep("breeding\nprobability", times = 8),
                                      rep("twinning\nprobability", times = 5),
                                      rep("CoY\nsurvival", times = 5),
                                      rep("yearling\nsurvival", times = 4),
                                      rep("recapture\nprobability", times = 5)))


cols <- model_df$name
caterpillar <- as.data.frame(res) %>%
  pivot_longer(cols = all_of(cols)) %>%
  left_join(x = .,
            y = model_df,
            by = "name") %>%
  mutate(parameter2 = factor(parameter2, levels = c("survival",
                                                    "breeding\nprobability", 
                                                    "twinning\nprobability", 
                                                    "CoY\nsurvival",
                                                    "yearling\nsurvival", 
                                                    "recapture\nprobability")),
         new_name = factor(new_name, levels = rev(c("a[1]", "a[2]", "a[3]", "a[4]", "a[5]",
                                                    "b[1]", "b[2]", "b[3]", "b[4]", "b[5]", "b[6]", "b[7]", "b[8]", 
                                                    "d[1]", "d[2]", "d[3]", "d[4]", "d[5]",
                                                    "e[1]", "e[2]", "e[3]", "e[4]", "e[5]",
                                                    "f[1]", "f[2]", "f[3]", "f[4]",
                                                    "g[1]", "g[2]", "g[3]", "g[4]", "sigma[omega]")))) %>%
  dplyr::select(new_name, name, value, covariate, parameter, parameter2) 


ggplot(caterpillar, 
       aes(x = value, y = new_name)) +
  geom_vline(xintercept = 0) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975),
               geom = "pointrange", orientation = "y", 
               size = 0.35) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975),
               geom = "pointrange", orientation = "y", 
               size = 0.35) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.25, upper = 0.75),
               geom = "pointrange", orientation = "y", 
               size = 0.35, linewidth = 1) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4, 6)) +
  scale_y_discrete(labels = function(labels) {
    # Filter out only existing labels in each facet to avoid NA issues
    # valid_labels <- labels_greek_and_not_greek[labels_greek_and_not_greek %in% labels]
    valid_labels <- rev(model_df$new_name[model_df$new_name %in% labels])
    parse(text = valid_labels)
  }) +
  scale_color_brewer(palette = "Dark2", na.value = "grey40") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic")) +
  labs(x = "posterior distribution") +
  facet_grid(parameter2~., 
             space = "free",
             scales = "free")

ggsave(full_plot,
       filename = "03_outputs/Figure S13.png",
       width = 18, height = 16, units = "cm", dpi = 600)



# ~ Fig. S14-15. Improvement estimates roe deer with error -------------------------

# Real value
{beta <- eta <- delta <- NULL
a1_phi <- logit(0.4)  # Survival (Fawn)   
a2_phi <- logit(0.85)  # Survival 
a3_phi <- -0.4  # time effect   
time <- as.vector(scale(1:8))
beta[1] <- 0.3  # Breeding probability (Y) 
beta[2] <- 0.4  # Breeding probability (NB/FB) 
beta[3] <- 0.75  # Breeding probability (SB) 
eta[1] <- 0.3  # Early offspring survival probability (Y) 
eta[2] <- 0.4  # Early offspring survival probability (NB/FB) 
eta[3] <- 0.75  # Early offspring survival probability (SB) 
gamma <- 0.4 # Probability of large litter/brood
delta[1] <- beta[1] * eta[1]
delta[2] <- beta[2] * eta[2]
delta[3] <- beta[3] * eta[3]
gamma <- 0.4 # Probability of large litter/brood
}


# ~~~ a. Load and process data -------------------------------------------------

path <- c("03_outputs/simulations/model_outputs/main_example/")

list_files_all <- list.files(path)

p <- c(0.1, 0.25)
pGPS <- c(0, 0.1, 0.25, 0.5, 1) 
pError <- c(0, 0.1, 0.2)

cols <- c("phiJ", "phiA", "a3_phi", 
          "delta_1", "delta_2", "delta_3", 
          "gamma") 

real_values <- data.frame(parameter = factor(cols, levels = cols),
                          real_value = c(plogis(a1_phi), plogis(a2_phi),  a3_phi, 
                                         delta[1], delta[2], delta[3], 
                                         gamma))

# Need to add the loop for pError !!!!!!!!!!!!!!!!!!!!!!
df_all <- data.frame()
for (n in 1:length(p)) {
  for (m in 1:length(pGPS)) {
    
    #++++++++++++++++++++++++++++++ Without GPS ++++++++++++++++++++++++++++++++#
    if (pGPS[m] == 0) {
      
      print(paste0("p = ", p[n]))
      print(paste0("pGPS = ", pGPS[m]))
      print(paste0("pError = ", 0))
      
      scenario <- paste0("p_", p[n], "_pGPS_", pGPS[m], "_")
      file_list <- list_files_all[grep(list_files_all, pattern = scenario)]
      
      model_fits <- list() 
      for (l in 1:length(file_list)) {
        load(paste0(path, file_list[l]))
        
        
        model_fits[[l]] <- rbind(get(substr(file_list[l], 1, nchar(file_list[l])-6))[[1]],
                                 get(substr(file_list[l], 1, nchar(file_list[l])-6))[[2]])
        
        rm(list = substr(file_list[l], 1, nchar(file_list[l])-6))
      }
      
      
      
      df <- as.data.frame(do.call(rbind, model_fits)) %>%
        mutate(simulation = rep(1:length(file_list), 
                                each = nrow(model_fits[[1]]))) %>%
        janitor::clean_names() %>%
        mutate(phiJ = plogis(a1_phi),
               phiA = plogis(a2_phi)) %>%
        dplyr::select(all_of(cols), simulation) %>%
        pivot_longer(cols = all_of(cols),
                     names_to = "parameter") %>%
        group_by(parameter,
                 simulation) %>%
        summarize(median = median(value), sd = sd(value)) %>%
        group_by(parameter) %>%
        summarize(upper_95 = quantile(median, probs = 0.975),
                  lower_95 = quantile(median, probs = 0.025),
                  upper_50 = quantile(median, probs = 0.75),
                  lower_50 = quantile(median, probs = 0.25),
                  mean = mean(median)) %>%
        mutate(p = p[n],
               p_GPS = pGPS[m],
               p_Error = 0) 
      
      df_all <- rbind(df_all, df)
      
    } else { 
      
      #+++++++++++++++++++++++++++++ With GPS ++++++++++++++++++++++++++++++++#
      
      for (k in 1:length(pError)) {
        
        print(paste0("p = ", p[n]))
        print(paste0("pGPS = ", pGPS[m]))
        print(paste0("pError = ", pError[k]))
        
        scenario <- paste0("p_", p[n], "_pGPS_", pGPS[m], "_pError_", pError[k])
        file_list <- list_files_all[grep(list_files_all, pattern = scenario)]
        
        model_fits <- list() 
        for (l in 1:length(file_list)) {
          load(paste0(path, file_list[l]))
          
          
          model_fits[[l]] <- rbind(get(substr(file_list[l], 1, nchar(file_list[l])-6))[[1]],
                                   get(substr(file_list[l], 1, nchar(file_list[l])-6))[[2]])
          
          rm(list = substr(file_list[l], 1, nchar(file_list[l])-6))
        }
        
        
        df <- as.data.frame(do.call(rbind, model_fits)) %>%
          mutate(simulation = rep(1:length(file_list), 
                                  each = nrow(model_fits[[1]]))) %>%
          janitor::clean_names() %>%
          mutate(phiJ = plogis(a1_phi),
                 phiA = plogis(a2_phi),
                 delta_1 = beta_1 * eta_1,
                 delta_2 = beta_2 * eta_2,
                 delta_3 = beta_3 * eta_3) %>%
          dplyr::select(all_of(cols), simulation) %>%
          pivot_longer(cols = all_of(cols),
                       names_to = "parameter") %>%
          group_by(parameter,
                   simulation) %>%
          summarize(median = median(value), sd = sd(value)) %>%
          group_by(parameter) %>%
          summarize(upper_95 = quantile(median, probs = 0.975),
                    lower_95 = quantile(median, probs = 0.025),
                    upper_50 = quantile(median, probs = 0.75),
                    lower_50 = quantile(median, probs = 0.25),
                    mean = mean(median)) %>%
          mutate(p = p[n],
                 p_GPS = pGPS[m],
                 p_Error = pError[k]) 
        
        df_all <- rbind(df_all, df)
      }
    }
  }
}

df_plot <- df_all %>%
  mutate(p = as.factor(p),
         p_GPS = as.factor(p_GPS),
         p_Error = as.factor(p_Error))

df_plot_title <- c(expression(paste(italic("\u03C6")['1'])),
                   expression(paste(italic("\u03C6")['2'])),
                   expression(paste(italic("\u03C6")['time'])),
                   
                   # delta
                   expression(paste(italic("\u03B4")['1'])),
                   expression(paste(italic("\u03B4")['2'])),
                   expression(paste(italic("\u03B4")['3'])),
                   
                   # gamma
                   expression(paste(italic("\u03B3"))))
ymin <- c(0, 0, -1.5, 0, 0, 0, 0)
ymax <- c(1, 1, 1.5, 1, 1, 1, 1)


# ~~~ b. Plot ------------------------------------------------------------------

# p = 0.10
df_plot_1 <- df_plot %>%
  filter(p == 0.1)
for (k in 1:length(cols)) {
  df_plot_k <- df_plot_1 %>%
    filter(parameter == cols[k])
  
  plot_title_k <- parse(text = df_plot_title[k])[[1]]
  
  assign(x = paste0("plot_", k),
         value =  ggplot(data = df_plot_k,
                         aes(x = p_Error , y = mean, color = p_Error, alpha = p_GPS)) + 
           geom_pointrange(aes(ymin = lower_95, ymax = upper_95),
                           position = position_dodge(width = 0.4)) +
           geom_pointrange(aes(ymin = lower_50, ymax = upper_50),
                           position = position_dodge(width = 0.4), 
                           linewidth = 1) +
           geom_hline(yintercept = real_values$real_value[k]) +
           
           scale_alpha_discrete(range = c(0.1, 1)) +
           # scale_alpha_continuous(range = c(0.2, 0.9), breaks = c(0, 0.1, 0.25, 0.5, 1)) +
           scale_color_manual(values = c("black", "#ff3300", "darkred")) + 
           # scale_color_grey(start = 0.8, end = 0.2) +
           scale_y_continuous(limits = c(ymin[k], ymax[k])) +
           theme_bw() +
           theme(plot.title = element_text(hjust = 0.5, size = 10),
                 legend.title = element_text(hjust = 0.5),
                 axis.title.y = element_blank()) +
           labs(x = bquote(.(expression(paste(italic("\u03B5"))))), y = bquote(.(plot_title_k)), 
                title = bquote(.(plot_title_k)),
                alpha = expression(kappa),
                color = "false neg\nrate")
  )
}
full_plot <- plot_1 + plot_2 + plot_3 + plot_4 + plot_5 + plot_6 + plot_7 +
  plot_layout(guides = 'collect')

ggsave(full_plot,
       filename = "03_outputs/Figure S14.png",
       width = 18, height = 20, units = "cm", dpi = 600)

# p = 0.25
df_plot_1 <- df_plot %>%
  filter(p == 0.25)
for (k in 1:length(cols)) {
  df_plot_k <- df_plot_1 %>%
    filter(parameter == cols[k])
  
  plot_title_k <- parse(text = df_plot_title[k])[[1]]
  
  assign(x = paste0("plot_", k),
         value =  ggplot(data = df_plot_k,
                         aes(x = p_Error , y = mean, color = p_Error, alpha = p_GPS)) + 
           geom_pointrange(aes(ymin = lower_95, ymax = upper_95),
                           position = position_dodge(width = 0.4)) +
           geom_pointrange(aes(ymin = lower_50, ymax = upper_50),
                           position = position_dodge(width = 0.4), 
                           linewidth = 1) +
           geom_hline(yintercept = real_values$real_value[k]) +
           
           scale_alpha_discrete(range = c(0.1, 1)) +
           # scale_alpha_continuous(range = c(0.2, 0.9), breaks = c(0, 0.1, 0.25, 0.5, 1)) +
           scale_color_manual(values = c("black", "#ff3300", "darkred")) + 
           # scale_color_grey(start = 0.8, end = 0.2) +
           scale_y_continuous(limits = c(ymin[k], ymax[k])) +
           theme_bw() +
           theme(plot.title = element_text(hjust = 0.5, size = 10),
                 legend.title = element_text(hjust = 0.5),
                 axis.title.y = element_blank()) +
           labs(x = bquote(.(expression(paste(italic("\u03B5"))))), y = bquote(.(plot_title_k)), 
                title = bquote(.(plot_title_k)),
                alpha = expression(kappa),
                color = "false neg\nrate")
  )
}
full_plot <- plot_1 + plot_2 + plot_3 + plot_4 + plot_5 + plot_6 + plot_7 +
  plot_layout(guides = 'collect')

ggsave(full_plot,
       filename = "03_outputs/Figure S15.png",
       width = 18, height = 20, units = "cm", dpi = 600)




# __________________________________________________________________________----




