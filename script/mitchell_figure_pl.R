################################################################################
# Preamble                                                                     #
################################################################################

library("here")
library("tidyverse")
library("janitor")

my_theme <- theme_classic(base_size = 10) 

theme_set(my_theme)

# Colours
col_hyp <- "#F8766D"
col_pow <- "#00BFC4"

################################################################################
# Mitchell - Data                                                              #
################################################################################

mitchell_phys <- read_csv(here("data","mitchell_data_phys.csv")) %>%
  clean_names() %>%
  rename(id = participant_id)

IDs <- mitchell_phys %>%
  pull(id) %>%
  unique() 

apo_tib <- mitchell_phys %>%
  select(id,apo) %>% rename(power = apo) %>%
  mutate(durations = 5,
         effort_id = 6)


mitchell_pow <- read_csv(here("data","mitchell_data_pow.csv")) %>%
  clean_names() %>%
  rename(durations = duration) %>%
  cbind(effort_id = rep(1:5, each = length(IDs))) %>%
  rbind(apo_tib) %>%
  filter(!(id == 8 & effort_id %in% c(3,4))) 

# Efforts 3 and 4 for athlete 8 do not appear to be maximal, hence the reason 
# the study re-tested 2 times.


################################################################################
# Mitchell - Fitting the power-law model to each athlete                       #
################################################################################

fit_pow_lr_pow <- function(data_pl,t){
  pl_lr <- lm(data = data_pl, I(log(power)) ~ I(log(durations)))
  a <- as.numeric(exp(coef(pl_lr)[1]))
  b <- as.numeric(-coef(pl_lr)[2])
  pl_power <- a/(t^b)
  return(list(S = a,
              E = 1-b,
              duration = t,
              power = pl_power))
}


aux_IDs <- mitchell_pow %>%
  pull(id) %>%
  unique() 

pl_tib <- tibble()

pl_tib_plot <- tibble()
t_p <- 1:(max(mitchell_pow$durations, na.rm = TRUE) * 1.1)

for (i in aux_IDs) {
  data_pl_i <- mitchell_pow %>%
    filter(id == i  & durations >= 120 & durations <= 1200) 
  pl_list_i <- fit_pow_lr_pow(data_pl = data_pl_i, t = t_p)
  pl_tib_aux <- tibble(id = i,
                       S = pl_list_i$S,
                       E = pl_list_i$E)
  pl_tib <- rbind(pl_tib,pl_tib_aux)
  
  pl_tib_plot_aux <- tibble(id = i,
                       durations = pl_list_i$duration,
                       power = pl_list_i$power)
  
  pl_tib_plot <- rbind(pl_tib_plot,pl_tib_plot_aux)
}

################################################################################
# Mitchell - Fitting the critical-power model to each athlete                  #
################################################################################

fit_cp_lr_pow <- function(data_cp,t){
  cp_lr <- lm(data = data_cp, power ~ I(1/durations))
  w_prime <- cp_lr$coefficients["I(1/durations)"]
  cp <- cp_lr$coefficients["(Intercept)"]
  cp_power <- w_prime/t + cp
  
  return(list(cp = cp,
              w_prime = w_prime,
              duration = t,
              power = cp_power))
}


aux_IDs <- mitchell_pow %>%
  pull(id) %>%
  unique() 

cp_tib <- tibble()

cp_tib_plot <- tibble()
t_p <- 1:(max(mitchell_pow$durations, na.rm = TRUE) * 1.1)

for (i in aux_IDs) {
  data_cp_i <- mitchell_pow %>%
    filter(id == i & durations >= 120 & durations <= 900) 
  cp_list_i <- fit_cp_lr_pow(data_cp = data_cp_i, t = t_p)
  cp_tib_aux <- tibble(id = i,
                       cp = cp_list_i$cp,
                       w_prime = cp_list_i$w_prime)
  cp_tib <- rbind(cp_tib,cp_tib_aux)
  
  cp_tib_plot_aux <- tibble(id = i,
                            durations = cp_list_i$duration,
                            power = cp_list_i$power)
  
  cp_tib_plot <- rbind(cp_tib_plot,cp_tib_plot_aux)
}

par_tib <- pl_tib %>% 
  left_join(cp_tib, by = "id")



################################################################################
# Mitchell - Correlation plot 30 minute power                                  #
################################################################################

corr_plot_tib <- mitchell_phys %>%
  left_join(pl_tib, by = "id") %>%
  mutate(pl_30 = S/((60*30)^(1-E))) 

corr_plot_tib_test <-  corr_plot_tib %>%
  select(id,pl_30,
         percent_type_1,percent_type_2, 
         mhc_i_csa,mhc_ii_csa,
         c_density_mm2,c_f,
         caf_mhc_i,caf_mhc_ii) %>%
  gather(key = "var", value = "value", -c("id","pl_30")) %>%
  group_by(var) %>%
  mutate(correlation = cor(value,pl_30,use="complete.obs"),
         p_val = cor.test(value,pl_30,use="complete.obs")$p.value,
         x_pos =  (max(value, na.rm = TRUE) - min(value, na.rm = TRUE)) * 0.8 + min(value, na.rm = TRUE),
         y_pos =  415, #max(pl_30, na.rm = TRUE) * 1.35, 
         lab = if_else(round(p_val,3) < 0.001,paste("r =", round(correlation, 3),"\np < 0.001") ,
                       paste("r =", round(correlation, 3),"\np =", round(p_val,3)))) %>%
  mutate(var = recode(var,
                      percent_type_1 = "Type I fibre %",
                      percent_type_2 = "Type II fibre %",
                      mhc_i_csa = "CSA type I (\u03BCm^-2)",
                      mhc_ii_csa = "CSA type II (\u03BCm^-2)",
                      c_density_mm2 = "Capillary density (cap.mm^-2)",
                      c_f = "Capillary to fibre ratio",
                      caf_mhc_i = "CC type I",
                      caf_mhc_ii = "CC type II"))

var_levels <- corr_plot_tib_test %>% 
  pull(var) %>%
  unique()

corr_plot_tib_test %>%
  ggplot(aes(x = value, y = pl_30)) +
  geom_point(aes()) +
  geom_text(mapping = aes(x = x_pos, y = y_pos, label = lab), size = 2.5) + 
  geom_smooth(method = "lm", formula = 'y~x', se = FALSE) +
  labs(x = "", y = "30-min. max. power (as predicted by the power-law model)") +
  scale_colour_manual(values = c("black", "red"), guide = "none")+
  facet_wrap(~ factor(var, levels = var_levels), scales = "free_x", ncol = 2, strip.position = "bottom") +
  coord_cartesian(ylim = c(150, 435)) +
  theme(strip.background = element_blank(), strip.placement = "outside") -> pl_30_cor

ggsave(here("figures","pl_30_par_cor_plot.png"), plot = pl_30_cor, dpi = 320, width = 0.5 * 8, height =  0.5 * 12)


################################################################################
# Mitchell - Correlation plot for E                                            #
################################################################################

corr_plot_tib_test <-  corr_plot_tib %>%
  select(id,E,
         percent_type_1,percent_type_2, 
         mhc_i_csa,mhc_ii_csa,
         c_density_mm2,c_f,
         caf_mhc_i,caf_mhc_ii) %>%
  gather(key = "var", value = "value", -c("id","E")) %>%
  group_by(var) %>%
  mutate(correlation = cor(value,E,use="complete.obs"),
         p_val = cor.test(value,E,use="complete.obs")$p.value,
         x_pos =  (max(value, na.rm = TRUE) - min(value, na.rm = TRUE)) * 0.85 + min(value, na.rm = TRUE),
         y_pos =  0.92,
         lab = if_else(round(p_val,3) < 0.001,paste("r =", round(correlation, 3),"\np < 0.001") ,
                       paste("r =", round(correlation, 3),"\np =", round(p_val,3)))) %>%
  mutate(var = recode(var,
                      percent_type_1 = "Type I fibre %",
                      percent_type_2 = "Type II fibre %",
                      mhc_i_csa = "CSA type I (\u03BCm^-2)",
                      mhc_ii_csa = "CSA type II (\u03BCm^-2)",
                      c_density_mm2 = "Capillary density (cap.mm^-2)",
                      c_f = "Capillary to fibre ratio",
                      caf_mhc_i = "CC type I",
                      caf_mhc_ii = "CC type II"))

var_levels <- corr_plot_tib_test %>% 
  pull(var) %>%
  unique()

corr_plot_tib_test  %>%
  ggplot(aes(x = value, y = E)) +
  geom_point() +
  geom_text(mapping = aes(x = x_pos, y = y_pos, label = lab), size = 2.5) +
  geom_smooth(method = "lm", formula = 'y~x', se = FALSE) +
  labs(x = "", y = "Endurance parameter, E, of the power-law model") +
  facet_wrap(~ factor(var, levels = var_levels), scales = "free_x", ncol = 2,strip.position = "bottom") +
  coord_cartesian(ylim = c(0.77,0.94)) +
  theme(strip.background = element_blank(), strip.placement = "outside") -> E_cor

ggsave(here("figures","Endurance_par_cor_plot.png"), plot = E_cor, dpi = 320, width = 0.5 * 8, height = 0.5 * 12)

