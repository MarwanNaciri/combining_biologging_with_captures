#==============================================================================#
#                                                                              #
#                        Plot main figures for Chapter 3                       #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(ggdist)
library(patchwork)
library(nimble) # for the logit function

get_median_and_CI <- function(x, lower, upper) {
  output <- data.frame(y = median(x), 
                       ymin = quantile(x, probs = lower),
                       ymax = quantile(x, probs = upper))
  return(output)
}


# ~ Fig. 2:  Improvement estimates roe deer ---------------------------------

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

p <- c(0.05, 0.1, 0.25, 0.5)
p_GPS <- c(0, 0.1, 0.25, 0.5, 1) 
p_Error <- 0

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


# ~~~ b. Plot -----------------------------------------------------------------

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
       filename = "03_outputs/Figure 1.png",
       units = 'cm', width = 18, height = 16, dpi = 600)

# ~ Fig. 3. Improvements increased information content -----------------------------------------------

# Here I use p = 0.25

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

p <- 0.25 # c(0.05, 0.1, 0.25, 0.5)
p_GPS <- c(0, 0.1, 0.25, 0.5, 1) 
p_Error <- 0

cols <- c("phiJ", "phiA", "a3_phi", 
          "delta_1", "delta_2", "delta_3", 
          "gamma") 

real_values <- data.frame(parameter = factor(cols, levels = cols),
                          real_value = c(plogis(a1_phi), plogis(a2_phi),  a3_phi, 
                                         delta[1], delta[2], delta[3], 
                                         gamma))
scenario <- c("", "_offspring_survival", "_litter_size")
df_all <- data.frame()
for (n in 1:length(scenario)) {
  
  print(paste0("scenario = ", scenario[n]))
  
  path <- paste0("07_results/model_outputs/CR_GPS_GLS/simulations_revisions/main_example",
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

# p = 0.25
df_plot_2 <- df_plot %>%
  filter(p == 0.25)

for (k in 1:length(cols)) {
  df_plot_k <- df_plot_2 %>%
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

ggsave(plot = full_plot, 
       filename = "03_outputs/Figure 1.png",
       units = 'cm', width = 18, height = 16, dpi = 600)

# ~ Fig. 4: Overall reproductive success rates ------------------------------

custom_palette_4 <- c("#006399", "#003049","#D62828", "#F99418")

# median expected number of weaned cubs
load("03_outputs/case_study/fit_with_biologging_1.RData")
load("03_outputs/case_study/fit_with_biologging_2.RData")
res <- rbind(fit_with_biologging_1, fit_with_biologging_2) 

mAge4 <- c(1, 0, 0, 0) ; mAgeY <- c(0, 1, 0, 0) ; mAgeM <- c(0, 0, 1, 0) ; mAgeO <- c(0, 0, 0, 1)
probability <- repr_success <- eta <- beta <- gamma <- s01 <- s02 <- s1 <- 
  matrix(NA, nrow = dim(res)[1], ncol = 4,
         dimnames = list(1:dim(res)[1],  c("4yr", "5-8yr", "9-14yr", "\u226515yr")))
for (i in 1:4) {              # maternal age group
  for (j in 1:dim(res)[1]) {      # MCMC iteration
    eta[j, i] <- plogis(res[j, "beta_eta[1]"] * mAge4[i] +
                          res[j, "beta_eta[2]"] * mAgeY[i] +
                          res[j, "beta_eta[3]"] * mAgeM[i] + 
                          res[j, "beta_eta[4]"] * mAgeO[i])
    
    beta[j, i] <- plogis(res[j, "beta_beta[1]"] * (mAge4[i] + mAgeY[i]) +
                           res[j, "beta_beta[2]"] * mAgeM[i] + 
                           res[j, "beta_beta[3]"] * mAgeO[i])
    gamma[j, i] <- plogis(res[j, "beta_gamma[1]"] * (mAge4[i] + mAgeY[i]) +
                            res[j, "beta_gamma[2]"] * mAgeM[i] + 
                            res[j, "beta_gamma[3]"] * mAgeO[i])
    s01[j, i] <- plogis(res[j, "beta_s0[1]"] * (mAge4[i] + mAgeY[i]) +
                          res[j, "beta_s0[2]"] * mAgeM[i] + 
                          res[j, "beta_s0[3]"] * mAgeO[i] +
                          res[j, "beta_s0[4]"])
    s02[j, i] <- plogis(res[j, "beta_s0[1]"] * (mAge4[i] + mAgeY[i]) +
                          res[j, "beta_s0[2]"] * mAgeM[i] + 
                          res[j, "beta_s0[3]"] * mAgeO[i])
    s1[j, i] <- plogis(res[j, "beta_s1[1]"] * (mAge4[i] + mAgeY[i]) +
                         res[j, "beta_s1[2]"] * mAgeM[i] + 
                         res[j, "beta_s1[3]"] * mAgeO[i])
    
    probability[j, i] <- eta[j, i] * beta[j, i] * 
      # A single cub (who survives)
      ((1-gamma[j, i])*s01[j, i]*s1[j, i] +
         # Two cubs who both survive   
         gamma[j, i]*s02[j, i]^2*s1[j, i]^2 +
         # Two cubs and only 1 survives
         gamma[j, i]*2*s02[j, i]*(1-s02[j, i])*s1[j, i] +
         gamma[j, i]*s02[j, i]^2*2*s1[j, i]*(1-s1[j, i]) 
      )
    
    repr_success[j, i] <- eta[j, i] * beta[j, i] * 
      # A single cub (who survives)
      ((1-gamma[j, i])*s01[j, i]*s1[j, i] +
         # Two cubs who both survive   
         2 * gamma[j, i]*s02[j, i]^2*s1[j, i]^2 +
         # Two cubs and only 1 survives
         gamma[j, i]*2*s02[j, i]*(1-s02[j, i])*s1[j, i] +
         gamma[j, i]*s02[j, i]^2*2*s1[j, i]*(1-s1[j, i]) 
      )
    
  }
}

# Probability of successful reproduction
apply(probability, 2, median)
get_probability_direction(probability[, 3] - probability[, 4])
df_plot <- as.data.frame.table(probability) %>%
  rename(iteration = Var1, age_group = Var2, value = Freq)



plot_prob <- ggplot(data = df_plot,
                    aes(x = age_group, y = value, fill = age_group)) +
  geom_violin(alpha = 0.8) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975),
               geom = "pointrange", size = 0.5, stroke = 0.5,
               color = "grey", shape = 21) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.25, upper = 0.75),
               geom = "pointrange", size = 0.5, stroke = 0.5,
               linewidth = 1, shape = 21, color = "grey") +
  theme_bw() +
  lims(y = c(0, 1)) +
  theme(legend.position = "none") +
  # scale_fill_brewer(palette = "Dark2") +
  scale_fill_manual(values = custom_palette_4) +
  labs(x = "age at mating", y = "probability of\nsuccessful reproduction")
plot_prob

ggsave("03_outputs/Figure 3.png",
       width = 9, height = 7.5, units = "cm", dpi = 600)




# ~ Fig. 5: Effect of sea-ice availability -------------------------------------

get_probability_direction <- function(posterior) {
  pd <- ifelse(median(posterior) > 0, sum(posterior > 0) / length(posterior),
               sum(posterior < 0) / length(posterior))
  return(pd)
}

load("03_outputs/case_study/fit_with_biologging_1.RData")
load("03_outputs/case_study/fit_with_biologging_2.RData")
res <- rbind(fit_with_biologging_1, fit_with_biologging_2) 

sea_ice_data <- read_csv("01_inputs/sea_ice_metrics.csv", show_col_types = F)
sea_ice <- NULL
years <- 1986:2024
for (t in 1:length(years)) {
  sea_ice[t] <- sea_ice_data$ice_free_days[which(sea_ice_data$year == years[t])] 
}
sea_ice_s <- as.vector(scale(sea_ice))

set.seed(1)
ice_free_days <- sea_ice_data %>% filter(year > 1986)
offset <- rnorm(n = nrow(ice_free_days), mean = 0, sd = 0.75)

lengthgrid <- 100
grid <- seq(min(sea_ice), max(sea_ice), length = lengthgrid) 
grid_scaled <- seq(min(sea_ice_s), max(sea_ice_s), length = lengthgrid) 



# Subadult survival ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subadult_survival  <- matrix(data = NA, nrow = dim(res)[1], ncol = lengthgrid,
                             dimnames = list(1:dim(res)[1], grid))

for (t in 1:lengthgrid) {                 # year
  for (k in 1:dim(res)[1]) {         # MCMC iteration
    subadult_survival[k, t] <- plogis(res[k, "beta_phi[1]"] +
                                        res[k, "beta_phi[5]"] * grid_scaled[t]) 
  }                               
}

df_plot_1 <- as.data.frame.table(subadult_survival) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  group_by(var) %>%
  summarize(median = median(value)) %>%
  mutate(var = as.numeric(as.character(var)))


n_iterations <- 100
set.seed(1)
iterations <- sample(1:dim(res)[1], size = n_iterations)

subadult_survival <- matrix(data = NA, nrow = n_iterations, ncol = lengthgrid,
                            dimnames = list(1:n_iterations, grid))

for (t in 1:lengthgrid) {              # date of emergence
  for (k in 1:n_iterations) {                # MCMC iteration
    subadult_survival[k, t] <- plogis(res[iterations[k], "beta_phi[1]"] +
                                        res[iterations[k], "beta_phi[5]"] * grid_scaled[t]) 
  }                               
}

df_plot_2 <- as.data.frame.table(subadult_survival) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  mutate(var = as.numeric(as.character(var)))

plot_subadult_survival <- ggplot() +
  geom_line(data = df_plot_2, aes(x = var, y = value, group = iteration), 
            color = "grey35", linewidth = 0.3, alpha = 0.2) +
  geom_line(data = df_plot_1, aes(x = var, y = median), linetype = "solid") +
  geom_rug(data = ice_free_days, aes(x = ice_free_days + offset, color = year), 
           alpha = 0.90, linewidth = 0.5) +
  scale_colour_stepsn(colours = c("#022C6B", "#1B72B6", "#FB6A49", "#CF1418"),
                      values = NULL,
                      space = "Lab",
                      na.value = "grey50",
                      guide = "colourbar",
                      aesthetics = "colour") +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "ice-free days", y = expression(paste("subadult survival (", italic("\u03C6"), ")")))


caterpillar <- as.data.frame(res) %>%
  janitor::clean_names() %>%
  mutate(sea_ice = beta_phi_5) %>%
  dplyr::select(sea_ice) %>%
  mutate(parameter = "sea_ice")

label_pd <- paste0("p[d] == ", round(get_probability_direction(caterpillar$sea_ice), 2))

caterpillar_plot <- ggplot(data = caterpillar,
                           aes(x = sea_ice, y = parameter)) +
  geom_vline(xintercept = 0, color = "black", linetype = "dotted") +
  stat_halfeye(aes(x = sea_ice, y = parameter, fill = after_stat(x < 0)), 
               color = NA, alpha = 1) + 
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975), orientation = "y",
               geom = "pointrange", size = 0.2, linewidth = 0.25, position = position_dodge(0.5)) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.25, upper = 0.75), orientation = "y",
               geom = "pointrange", size = 0.3, linewidth = 1, 
               position = position_dodge(0.5)) +
  annotate("text", x = median(caterpillar$sea_ice), y = 1.3, 
           label = label_pd, parse = TRUE, size = 3) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1)) +
  coord_cartesian(xlim = c(-0.85, 1.2)) + 
  scale_fill_manual(values = c("grey90", "grey75")) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_blank(),
        legend.position = "none") +
  labs(x = "", y = "") 

inset <- inset_element(
  caterpillar_plot,
  left = -0.06,
  bottom = -0.06,
  right = 0.50,
  top = 0.52
)

plot_subadult_survival <- plot_subadult_survival + inset



# Denning probability ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
denning_prob  <- matrix(data = NA, nrow = dim(res)[1], ncol = lengthgrid,
                        dimnames = list(1:dim(res)[1], grid))

for (t in 1:lengthgrid) {                 # year
  for (k in 1:dim(res)[1]) {         # MCMC iteration
    denning_prob[k, t] <- plogis(res[k, "beta_eta[3]"] +
                                   res[k, "beta_eta[9]"] * grid_scaled[t]) 
  }                               
}

df_plot_1 <- as.data.frame.table(denning_prob) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  group_by(var) %>%
  summarize(median = median(value)) %>%
  mutate(var = as.numeric(as.character(var)))


n_iterations <- 100
set.seed(1)
iterations <- sample(1:dim(res)[1], size = n_iterations)

denning_prob <- matrix(data = NA, nrow = n_iterations, ncol = lengthgrid,
                       dimnames = list(1:n_iterations, grid))

for (t in 1:lengthgrid) {              # date of emergence
  for (k in 1:n_iterations) {                # MCMC iteration
    denning_prob[k, t] <- plogis(res[iterations[k], "beta_eta[3]"] +
                                   res[iterations[k], "beta_eta[9]"] * grid_scaled[t]) 
  }                               
}

df_plot_2 <- as.data.frame.table(denning_prob) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  mutate(var = as.numeric(as.character(var)))


plot_denning_prob <- ggplot() +
  geom_line(data = df_plot_2, aes(x = var, y = value, group = iteration), 
            color = "grey35", linewidth = 0.3, alpha = 0.2) +
  geom_line(data = df_plot_1, aes(x = var, y = median), linetype = "solid") +
  geom_rug(data = ice_free_days, aes(x = ice_free_days + offset, color = year), 
           alpha = 0.90, linewidth = 0.5) +
  scale_colour_stepsn(colours = c("#022C6B", "#1B72B6", "#FB6A49", "#CF1418"),
                      values = NULL,
                      space = "Lab",
                      na.value = "grey50",
                      guide = "colourbar",
                      aesthetics = "colour") +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "ice-free days", y = expression(paste("denning probability (", italic("\u03B7"), ")")))


caterpillar <- as.data.frame(res) %>%
  janitor::clean_names() %>%
  mutate(sea_ice = beta_eta_9) %>%
  dplyr::select(sea_ice) %>%
  mutate(parameter = "sea_ice")

label_pd <- paste0("p[d] == ", round(get_probability_direction(caterpillar$sea_ice), 2))

caterpillar_plot <- ggplot(data = caterpillar,
                           aes(x = sea_ice, y = parameter)) +
  geom_vline(xintercept = 0, color = "black", linetype = "dotted") +
  stat_halfeye(aes(x = sea_ice, y = parameter, fill = after_stat(x < 0)), 
               color = NA, alpha = 1) + 
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975), orientation = "y",
               geom = "pointrange", size = 0.2, linewidth = 0.25, position = position_dodge(0.5)) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.25, upper = 0.75), orientation = "y",
               geom = "pointrange", size = 0.3, linewidth = 1, 
               position = position_dodge(0.5)) +
  annotate("text", x = median(caterpillar$sea_ice), y = 1.3, 
           label = label_pd, parse = TRUE, size = 3) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1)) +
  coord_cartesian(xlim = c(-0.85, 1.2)) + 
  scale_fill_manual(values = c("grey90", "grey75")) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_blank(),
        legend.position = "none") +
  labs(x = "", y = "") 

inset <- inset_element(
  caterpillar_plot,
  left = -0.06,
  bottom = -0.06,
  right = 0.50,
  top = 0.52
)

plot_denning_prob <- plot_denning_prob + inset


# Early litter survival ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
early_litter_survival  <- matrix(data = NA, nrow = dim(res)[1], ncol = lengthgrid,
                                 dimnames = list(1:dim(res)[1], grid))

for (t in 1:lengthgrid) {                 # year
  for (k in 1:dim(res)[1]) {         # MCMC iteration
    early_litter_survival[k, t] <- plogis(res[k, "beta_beta[2]"] +
                                            res[k, "beta_beta[5]"] * grid_scaled[t]) 
  }                               
}

df_plot_1 <- as.data.frame.table(early_litter_survival) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  group_by(var) %>%
  summarize(median = median(value)) %>%
  mutate(var = as.numeric(as.character(var)))


n_iterations <- 100
set.seed(1)
iterations <- sample(1:dim(res)[1], size = n_iterations)

early_litter_survival <- matrix(data = NA, nrow = n_iterations, ncol = lengthgrid,
                                dimnames = list(1:n_iterations, grid))

for (t in 1:lengthgrid) {              # date of emergence
  for (k in 1:n_iterations) {                # MCMC iteration
    early_litter_survival[k, t] <- plogis(res[iterations[k], "beta_beta[2]"] +
                                            res[iterations[k], "beta_beta[5]"] * grid_scaled[t]) 
  }                               
}

df_plot_2 <- as.data.frame.table(early_litter_survival) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  mutate(var = as.numeric(as.character(var)))


plot_early_litter_survival <- ggplot() +
  geom_line(data = df_plot_2, aes(x = var, y = value, group = iteration), 
            color = "grey35", linewidth = 0.3, alpha = 0.2) +
  geom_line(data = df_plot_1, aes(x = var, y = median), linetype = "solid") +
  geom_rug(data = ice_free_days, aes(x = ice_free_days + offset, color = year), 
           alpha = 0.90, linewidth = 0.5) +
  scale_colour_stepsn(colours = c("#022C6B", "#1B72B6", "#FB6A49", "#CF1418"),
                      values = NULL,
                      space = "Lab",
                      na.value = "grey50",
                      guide = "colourbar",
                      aesthetics = "colour") +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "ice-free days", y = expression(paste("early litter survival (", italic("\u03B2"), ")")))


caterpillar <- as.data.frame(res) %>%
  janitor::clean_names() %>%
  mutate(sea_ice = beta_beta_5) %>%
  dplyr::select(sea_ice) %>%
  mutate(parameter = "sea_ice")

label_pd <- paste0("p[d] == ", round(get_probability_direction(caterpillar$sea_ice), 2))

caterpillar_plot <- ggplot(data = caterpillar,
                           aes(x = sea_ice, y = parameter)) +
  geom_vline(xintercept = 0, color = "black", linetype = "dotted") +
  stat_halfeye(aes(x = sea_ice, y = parameter, fill = after_stat(x < 0)), 
               color = NA, alpha = 1) + 
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975), orientation = "y",
               geom = "pointrange", size = 0.2, linewidth = 0.25, position = position_dodge(0.5)) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.25, upper = 0.75), orientation = "y",
               geom = "pointrange", size = 0.3, linewidth = 1, 
               position = position_dodge(0.5)) +
  annotate("text", x = median(caterpillar$sea_ice), y = 1.3, 
           label = label_pd, parse = TRUE, size = 3) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1)) +
  coord_cartesian(xlim = c(-0.85, 1.2)) + 
  scale_fill_manual(values = c("grey90", "grey75")) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_blank(),
        legend.position = "none") +
  labs(x = "", y = "")  

inset <- inset_element(
  caterpillar_plot,
  left = -0.06,
  bottom = -0.06,
  right = 0.50,
  top = 0.52
)

plot_early_litter_survival <- plot_early_litter_survival + inset


# Twinning prob ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
twinning_prob  <- matrix(data = NA, nrow = dim(res)[1], ncol = lengthgrid,
                         dimnames = list(1:dim(res)[1], grid))

for (t in 1:lengthgrid) {                 # year
  for (k in 1:dim(res)[1]) {         # MCMC iteration
    twinning_prob[k, t] <- plogis(res[k, "beta_gamma[2]"] +
                                    res[k, "beta_gamma[6]"] * grid_scaled[t]) 
  }                               
}

df_plot_1 <- as.data.frame.table(twinning_prob) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  group_by(var) %>%
  summarize(median = median(value)) %>%
  mutate(var = as.numeric(as.character(var)))


n_iterations <- 100
set.seed(1)
iterations <- sample(1:dim(res)[1], size = n_iterations)

twinning_prob <- matrix(data = NA, nrow = n_iterations, ncol = lengthgrid,
                        dimnames = list(1:n_iterations, grid))

for (t in 1:lengthgrid) {              # date of emergence
  for (k in 1:n_iterations) {                # MCMC iteration
    twinning_prob[k, t] <- plogis(res[iterations[k], "beta_gamma[2]"] +
                                    res[iterations[k], "beta_gamma[6]"] * grid_scaled[t]) 
  }                               
}

df_plot_2 <- as.data.frame.table(twinning_prob) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  mutate(var = as.numeric(as.character(var)))



set.seed(1)
ice_free_days <- sea_ice_data %>% filter(year > 1986) %>%
  mutate(year = year + 1) %>%
  left_join(x = .,
            y = data_CR_events,
            by = "year")


plot_twinning_prob <- ggplot() +
  geom_line(data = df_plot_2, aes(x = var, y = value, group = iteration), 
            color = "grey35", linewidth = 0.3, alpha = 0.2) +
  geom_line(data = df_plot_1, aes(x = var, y = median), linetype = "solid") +
  geom_rug(data = ice_free_days, aes(x = ice_free_days + offset, color = year), 
           alpha = 0.90, linewidth = 0.5) +
  scale_colour_stepsn(colours = c("#022C6B", "#1B72B6", "#FB6A49", "#CF1418"),
                      values = NULL,
                      space = "Lab",
                      na.value = "grey50",
                      guide = "colourbar",
                      aesthetics = "colour") +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "ice-free days", y = expression(paste("twinning probability (", italic("\u03B3"), ")")))

caterpillar <- as.data.frame(res) %>%
  janitor::clean_names() %>%
  mutate(sea_ice = beta_gamma_6) %>%
  dplyr::select(sea_ice) %>%
  mutate(parameter = "sea_ice")

label_pd <- paste0("p[d] == ", round(get_probability_direction(caterpillar$sea_ice), 2))

caterpillar_plot <- ggplot(data = caterpillar,
                           aes(x = sea_ice, y = parameter)) +
  stat_halfeye(aes(x = sea_ice, y = parameter, fill = after_stat(x < 0)), 
               color = NA, alpha = 1) + 
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975), orientation = "y",
               geom = "pointrange", size = 0.2, linewidth = 0.25, position = position_dodge(0.5)) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.25, upper = 0.75), orientation = "y",
               geom = "pointrange", size = 0.3, linewidth = 1, 
               position = position_dodge(0.5)) +
  geom_vline(xintercept = 0, color = "black", linetype = "dotted") +
  annotate("text", x = median(caterpillar$sea_ice), y = 1.3, 
           label = label_pd, parse = TRUE, size = 3) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1)) +
  coord_cartesian(xlim = c(-0.85, 1.2)) + 
  scale_fill_manual(values = c("grey90", "grey75")) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_blank(),
        legend.position = "none") +
  labs(x = "", y = "") 

inset <- inset_element(
  caterpillar_plot,
  left = -0.06,
  bottom = -0.06,
  right = 0.50,
  top = 0.52
)

plot_twinning_prob <- plot_twinning_prob + inset


# CoY survival  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
coy_survival  <- matrix(data = NA, nrow = dim(res)[1], ncol = lengthgrid,
                        dimnames = list(1:dim(res)[1], grid))

for (t in 1:lengthgrid) {                 # year
  for (k in 1:dim(res)[1]) {         # MCMC iteration
    coy_survival[k, t] <- plogis(res[k, "beta_s0[2]"] +
                                   res[k, "beta_s0[6]"] * grid_scaled[t]) 
  }                               
}

df_plot_1 <- as.data.frame.table(coy_survival) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  group_by(var) %>%
  summarize(median = median(value)) %>%
  mutate(var = as.numeric(as.character(var)))


n_iterations <- 100
set.seed(1)
iterations <- sample(1:dim(res)[1], size = n_iterations)

coy_survival <- matrix(data = NA, nrow = n_iterations, ncol = lengthgrid,
                       dimnames = list(1:n_iterations, grid))

for (t in 1:lengthgrid) {              # date of emergence
  for (k in 1:n_iterations) {                # MCMC iteration
    coy_survival[k, t] <- plogis(res[iterations[k], "beta_s0[2]"] +
                                   res[iterations[k], "beta_s0[6]"] * grid_scaled[t]) 
  }                               
}

df_plot_2 <- as.data.frame.table(coy_survival) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  mutate(var = as.numeric(as.character(var)))


plot_coy_survival <- ggplot() +
  geom_line(data = df_plot_2, aes(x = var, y = value, group = iteration), 
            color = "grey35", linewidth = 0.3, alpha = 0.2) +
  geom_line(data = df_plot_1, aes(x = var, y = median), linetype = "solid") +
  geom_rug(data = ice_free_days, aes(x = ice_free_days + offset, color = year), 
           alpha = 0.90, linewidth = 0.5) +
  scale_colour_stepsn(colours = c("#022C6B", "#1B72B6", "#FB6A49", "#CF1418"),
                      values = NULL,
                      space = "Lab",
                      na.value = "grey50",
                      guide = "colourbar",
                      aesthetics = "colour") +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "ice-free days", y = parse(text = "CoY~survival~(italic(s)[0])"))


caterpillar <- as.data.frame(res) %>%
  janitor::clean_names() %>%
  mutate(sea_ice = beta_s0_6) %>%
  dplyr::select(sea_ice) %>%
  mutate(parameter = "sea_ice")

label_pd <- paste0("p[d] == ", round(get_probability_direction(caterpillar$sea_ice), 2))

caterpillar_plot <- ggplot(data = caterpillar,
                           aes(x = sea_ice, y = parameter)) +
  geom_vline(xintercept = 0, color = "black", linetype = "dotted") +
  # stat_halfeye(aes(x = sea_ice, y = parameter), 
  #              color = NA, alpha = 0.5) + 
  stat_halfeye(aes(x = sea_ice, y = parameter, fill = after_stat(x < 0)), 
               color = NA, alpha = 1) + 
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975), orientation = "y",
               geom = "pointrange", size = 0.2, linewidth = 0.25, position = position_dodge(0.5)) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.25, upper = 0.75), orientation = "y",
               geom = "pointrange", size = 0.3, linewidth = 1, 
               position = position_dodge(0.5)) +
  annotate("text", x = median(caterpillar$sea_ice), y = 1.3, 
           label = label_pd, parse = TRUE, size = 3) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1)) +
  coord_cartesian(xlim = c(-0.85, 1.2)) + 
  scale_fill_manual(values = c("grey90", "grey75")) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_blank(),
        legend.position = "none") +
  labs(x = "", y = "") 

inset <- inset_element(
  caterpillar_plot,
  left = -0.06,
  bottom = -0.06,
  right = 0.50,
  top = 0.52
)

plot_coy_survival <- plot_coy_survival + inset


# yearling survival  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
yearling_survival  <- matrix(data = NA, nrow = dim(res)[1], ncol = lengthgrid,
                             dimnames = list(1:dim(res)[1], grid))

for (t in 1:lengthgrid) {                 # year
  for (k in 1:dim(res)[1]) {         # MCMC iteration
    yearling_survival[k, t] <- plogis(res[k, "beta_s1[2]"] +
                                        res[k, "beta_s1[4]"] * grid_scaled[t]) 
  }                               
}

df_plot_1 <- as.data.frame.table(yearling_survival) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  group_by(var) %>%
  summarize(median = median(value)) %>%
  mutate(var = as.numeric(as.character(var)))


n_iterations <- 100
set.seed(1)
iterations <- sample(1:dim(res)[1], size = n_iterations)

yearling_survival <- matrix(data = NA, nrow = n_iterations, ncol = lengthgrid,
                            dimnames = list(1:n_iterations, grid))

for (t in 1:lengthgrid) {              # date of emergence
  for (k in 1:n_iterations) {                # MCMC iteration
    yearling_survival[k, t] <- plogis(res[iterations[k], "beta_s1[2]"] +
                                        res[iterations[k], "beta_s1[4]"] * grid_scaled[t]) 
  }                               
}

df_plot_2 <- as.data.frame.table(yearling_survival) %>%
  rename(iteration = Var1, var = Var2, value = Freq) %>%
  mutate(var = as.numeric(as.character(var)))


# set.seed(1)
# ice_free_days <- sea_ice_data %>% filter(year > 1986)
# offset <- rnorm(n = nrow(ice_free_days), mean = 0, sd = 0.75)

plot_yearling_survival <- ggplot() +
  geom_line(data = df_plot_2, aes(x = var, y = value, group = iteration), 
            color = "grey35", linewidth = 0.3, alpha = 0.2) +
  geom_line(data = df_plot_1, aes(x = var, y = median), linetype = "solid") +
  geom_rug(data = ice_free_days, aes(x = ice_free_days + offset, color = year), 
           alpha = 0.90, linewidth = 0.5) +
  scale_colour_stepsn(colours = c("#022C6B", "#1B72B6", "#FB6A49", "#CF1418"),
                      # colours = c("#92c5de","#f7f7f7", "#e88465", "#8d0622"),
                      values = NULL,
                      space = "Lab",
                      na.value = "grey50",
                      guide = "colourbar",
                      aesthetics = "colour") +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "ice-free days", y = parse(text = "yearling~survival~(italic(s)[1])"))

caterpillar <- as.data.frame(res) %>%
  janitor::clean_names() %>%
  mutate(sea_ice = beta_s1_4) %>%
  dplyr::select(sea_ice) %>%
  mutate(parameter = "sea_ice")

label_pd <- paste0("p[d] == ", round(get_probability_direction(caterpillar$sea_ice), 2))

caterpillar_plot <- ggplot(data = caterpillar,
                           aes(x = sea_ice, y = parameter)) +
  geom_vline(xintercept = 0, color = "black", linetype = "dotted") +
  # stat_halfeye(aes(x = sea_ice, y = parameter), 
  #              color = NA, alpha = 0.5) + 
  stat_halfeye(aes(x = sea_ice, y = parameter, fill = after_stat(x > 0)), 
               color = NA, alpha = 1) + 
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.025, upper = 0.975), orientation = "y",
               geom = "pointrange", size = 0.2, linewidth = 0.25, position = position_dodge(0.5)) +
  stat_summary(fun.data = get_median_and_CI,
               fun.args = list(lower = 0.25, upper = 0.75), orientation = "y",
               geom = "pointrange", size = 0.3, linewidth = 1, 
               position = position_dodge(0.5)) +
  annotate("text", x = median(caterpillar$sea_ice), y = 1.3, 
           label = label_pd, parse = TRUE, size = 3) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1)) +
  coord_cartesian(xlim = c(-0.85, 1.2)) + 
  scale_fill_manual(values = c("grey90", "grey75")) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_blank(),
        legend.position = "none") +
  labs(x = "", y = "") 

inset <- inset_element(
  caterpillar_plot,
  left = -0.06,
  bottom = -0.06,
  right = 0.50,
  top = 0.52
)

plot_yearling_survival <- plot_yearling_survival + inset


((plot_subadult_survival) + (plot_denning_prob)) /
  ((plot_early_litter_survival) + (plot_twinning_prob)) /
  ((plot_coy_survival) + (plot_yearling_survival)) +
  plot_annotation(tag_levels = list(c("A", "", "B", "", "C", "", 
                                      "D", "", "E", "", "F", ""))) 

ggsave("03_outputs/Figure 4.png",
       width = 18, height = 21, units = "cm", dpi = 600)

#___________________________________________________________________________----