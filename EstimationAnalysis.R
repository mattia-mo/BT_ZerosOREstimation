# load required packages
if(!requireNamespace("tidyverse")) install.packages("tidyverse")
library(tidyverse)
if(!requireNamespace("ggplot2")) install.packages("ggplot2")
library(ggplot2)
if (!requireNamespace("naniar")) install.packages("naniar")
library(naniar)
if(!requireNamespace("ggpubr")) install.packages("ggpubr")
library(ggpubr)
if(!requireNamespace("latex2exp")) install.packages("latex2exp")
library(latex2exp)

# load generated data
# path <- ("insert path here")
# setwd(path)
estimates <- readRDS("estimates.rds")

# data wrangling
est_CI <- estimates %>%
  
  # breaking up the structure of one row per table to give each estimation of each method an own row
  pivot_longer(cols = starts_with("or_"), names_to = "method", values_to = "estimation_or") %>%
  mutate(method = as.factor(str_extract(method, "[^_]+$")),
         method = factor(method, labels = c("Fisher", "Midp", "Small", "Wald"))) %>% 
  group_by(rep,scenario) %>% 
  mutate(lower = case_when(
    method == "Wald" ~ lower_wald,
    method == "Fisher" ~ lower_fisher,
    method == "Small" ~ lower_small,
    method == "Midp" ~ lower_midp,
    TRUE ~ NA
  ),
  upper = case_when(
    method == "Wald" ~ upper_wald,
    method == "Fisher" ~ upper_fisher,
    method == "Small" ~ upper_small,
    method == "Midp" ~ upper_midp,
    TRUE ~ NA
  )) %>% 
  select(!(starts_with("upper_") | starts_with("lower_")))


# indicate all method failings (Inf, 0) as NA but keeping confidence intervals with infinite values
# also performing the log transformation for the estimation of the OR and the confidence intervals
# and indicate not appropriate results(NaN) as NA as well
est_CI_final <- est_CI %>% 
  mutate(across(starts_with("estimation_"), ~replace(., . == Inf | . == 0, NA)),
         log_or = round(log(or), 2), log_est = round(log(estimation_or), 2),
         log_lower = round(log(lower), 2), log_upper = round(log(upper), 2),
         across(c(lower, upper, log_lower, log_upper), ~replace_na(., NA)))


### handling of the zeros/method failing

# proceeding CCA to the estimations
# thereby estimations are only kept if all methods have been able to provide an outcome
est_CI_CCA <- est_CI_final %>% 
  group_by(table) %>% 
  mutate(estimation_or = if(any(is.na(estimation_or))) NA else estimation_or,
         lower = if(any(is.na(estimation_or))) NA else lower,
         upper = if(any(is.na(estimation_or))) NA else upper,
         log_est = if(any(is.na(estimation_or))) NA else log_est,
         log_lower = if(any(is.na(estimation_or))) NA else log_lower,
         log_upper = if(any(is.na(estimation_or))) NA else log_upper)
est_CI_CCA$handling <- "CCA"

# proving conducting CCA leads to same results as excluding tables with zeros
est_without_0 <- ungroup(est_CI_final %>% filter(zero == FALSE))
est_with_CCA <- ungroup(est_CI_CCA %>% filter(!is.na(estimation_or)) %>% select(!handling))

# Amount of different entries in the Dataset is 0
nrow(setdiff(est_without_0, est_with_CCA))

# ACA: as all available esimations are include in the data and all other not available estimations
# are marked as NA for ACA the data can stay as it is
est_CI_ACA <- est_CI_final
est_CI_ACA$handling <- "ACA"

# combine both handling approaches to one dataframe
est_handled <- rbind(est_CI_CCA, est_CI_ACA)


### getting an overview about zeros in tables

# compute the percentages of zeros per scenario
zero_perc <- est_handled %>% 
  group_by(scenario, rep) %>% 
  summarize(zero_in_table = any(zero), .groups = "drop") %>% 
  group_by(scenario) %>% 
  summarize(zero_perc = sum(zero_in_table)/max(rep))

# creating an overview of the scenarios including their parameters and the percentage of zeros
scenarios_overview <- estimates %>%
  select(scenario, or, n_obs, p_x) %>%
  unique() %>% 
  arrange(scenario) %>% 
  as.data.frame() %>% 
  merge(zero_perc, by = "scenario", all.x = T)

# extract number of scenarios that contain at least 1 table with zero
scen_zeros <- zero_perc %>% filter(zero_perc !=0) %>% 
  select(scenario) %>% 
  unlist() %>% 
  unname()

# getting a data frame of scenarios that contain at least 1 table with zero, 
# ordered by the percentage of zeros from lowest to highest
scen_zeros_ordered <- zero_perc %>% filter(zero_perc !=0) %>% 
  arrange(zero_perc) %>% 
  select(scenario, zero_perc)

# examine where the zeros are positioned in the tables
# thereby first for each scenario it is checked where zeros are occurring and how often in which position
# then the number of occurrences in the individual scenarios is added up to gain an insight 
# where the zeros occurred in the whole simulation
l <- estimates %>%
  group_by(scenario) %>%
  summarize(pos0 = list(Reduce("+", lapply(table, function(x) x == 0))))
zeros_positions <- Reduce("+", l$pos0)
zeros_positions

# proving that all fails to estimate an OR value are appearing in tables with zeros
est_failed <- est_CI_final %>% filter(is.na(estimation_or))
all(est_failed$zero == TRUE)

### Formatting functions for the plots
format_labels <- function(x) {
  TeX(x)
}
theta_labeller <- function(theta) {
  paste0("log(theta) == ", theta)
}
n_obs = n_labeller <- function(n) {
  paste0("n[obs] == ", n)
}
p_x = p_labeller <- function(p) {
  paste0("p[x] ==", p)
}


### Analysis of the generated Data

# creating the plot, that is displaying how high the percentage of zeros is in the scenarios
# as the parameters for the scenarios are 3, a new interaction variable needs to be created out of n_obs and p_x
n_rep <- max(estimates$rep)

estimates %>% 
  group_by(or, n_obs, p_x) %>% 
  summarize(zero_perc = sum(zero)/n_rep) %>% 
  mutate(n_obs = as.factor(as.numeric(n_obs)),
         or = as.factor(as.numeric(or)),
         p_x = as.factor(as.numeric(p_x)),
         n_and_p_x = interaction(p_x, n_obs, sep = ";"),
         n_and_p_x = factor(n_and_p_x, labels = paste0("$p_x = ", p_x, ", n_{obs} = ", n_obs, "$"))) %>% 
  ggplot() +
  geom_tile(aes(x = or, y = n_and_p_x, fill = zero_perc)) +
  scale_fill_continuous(low = "white", high = "black") +
  scale_y_discrete(labels = format_labels) +
  labs(x = TeX(r'(Odds Ratio $\theta$)'), y = "", fill = "Amount of \nZeros in %") +
  theme(legend.position = c(-0.1, -0.05), 
        legend.direction = "horizontal",
        legend.spacing.x = unit(0.5, "cm"),
        legend.background = element_rect(colour = "gray"),
        legend.title = element_text(size = 13),
        axis.text = element_text(size = 15), 
        axis.title.x = element_text(size = 15))

ggsave("plots/RatesofZeros_Scenarios.png", width = 10, height = 6)


### Visualization of estimates

# Swarmplot
# creating swarmplots for all scenarios, where x is the number of the table (rep) and y the OR estimation on log scale
# with the true OR also displayed by a horizontal line

est_CI_final %>%
  mutate(n_and_p_x = interaction(p_x, n_obs, sep = ";")) %>% 
  ggplot() +
  geom_point(aes(x=rep,y=log_est, color = method), alpha = 0.3) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  geom_hline(aes(yintercept = log_or, linetype = "true OR")) +
  facet_grid(log_or ~ n_obs + p_x, scales = "fixed", 
             labeller = labeller(.rows = as_labeller(theta_labeller, default = label_parsed), 
                                 n_obs = as_labeller(n_labeller, default = label_parsed),
                                 p_x = as_labeller(p_labeller, default = label_parsed))) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 17),
        legend.text = element_text(size = 17),
        legend.position = "bottom",
        legend.title = element_text(size = 17),
        strip.text = element_text(size = 14))+
  scale_linetype_manual(values = "solid", name = "")+
  labs(x="number of table", y=TeX(r'(\log(\hat{\theta}))'), color = "method")
ggsave("plots/swarmplot.png", width = 13, height = 9)

# displaying the OR estimations in a swarmplot with estimations despite a zero emphasized
# thereby it is differentiated between estimations of table without zeros (gray) and
# estimations of tables with zeros. In the latter a distinction is made between the estimation methods.

est_handled %>% 
  filter(handling == "ACA") %>% 
  mutate(zero_and_performed = ifelse(zero, paste0("yes_", method), "no")) %>% 
  ggplot(aes(rep, log_est, color = zero_and_performed)) +
  geom_point(data = . %>% filter(!startsWith(zero_and_performed, "yes_")), alpha = 0.7, position = "identity") +
  geom_point(data = . %>% filter(startsWith(zero_and_performed, "yes_")), alpha = 0.5) +
  geom_hline(aes(yintercept = log_or)) +
  scale_color_manual(values = c("no" = "gray", "yes_Fisher" = "green", "yes_Small" = "red", 
                                "yes_Wald" = "blue", "yes_Midp" = "orange")) +
  facet_grid(log_or~n_obs + p_x, labeller =labeller(.rows = as_labeller(theta_labeller, default = label_parsed), 
                                                    n_obs = as_labeller(n_labeller, default = label_parsed),
                                                    p_x = as_labeller(p_labeller, default = label_parsed))) +
  labs(title = "", x = "Number of table", y = expression(log(hat(theta))),
       color = "Table contains zero & \n method that performed") +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.position = "bottom",
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 15))
ggsave("plots/LogEst_allscens_0emphasized.png", width = 13, height = 10)

# EDA Histograms
# creating frequency histograms of the estimations for all scenarios and methods
# where the known OR value is displayed as well via a vertical line to get a better idea
# whether there is a tendency of over- or underestimation.

EDA_Plots <- est_CI %>% 
  group_by(scenario) %>% 
  group_split() %>% 
  map( ~ ggplot(.x) +
         geom_histogram(aes(estimation_or), bins = 40) +
         geom_vline(aes(xintercept = or)) +
         labs(x=TeX(r'(estimation $\hat{\theta}$)')) +
         facet_wrap(~method) +
         theme(strip.text = element_text(size = 15),
               axis.text = element_text(size = 15),
               axis.title = element_text(size = 17)))

lapply(seq_along(EDA_Plots), function(i){
  ggsave(paste0("plots/EDA_scenario_", i, ".png"), EDA_Plots[[i]], width = 10, height = 7)
})



# Plotting the estimations of tables with zeros, also missing values in the estimations (method failing)
# for both handling approaches
# missing values are plotted at the very bottom

est_handled %>%
  filter(zero) %>%
  ggplot() +
  geom_miss_point(aes(rep, estimation_or, color = method), alpha = 0.4) +
  geom_hline(aes(yintercept = or, linetype = "true OR")) +
  facet_wrap(scenario ~ handling, ncol = 6, labeller = partial(label_both, sep = " = ", 
                                                               multi_line = T)) +
  labs(title = "", y = "estimated odds ratio", x="number of table") +
  theme(strip.text = element_text(size = 10, lineheight = 0.3),
  axis.title = element_text(size = 15),
  legend.direction = "horizontal",
  legend.position = "bottom",
  legend.text = element_text(size = 13),
  legend.title = element_text(size = 15)
  ) +
  scale_linetype_manual(values = "solid", name = "")
ggsave("plots/zerotables_NAsincluded.png", width = 10, height = 8)


### Performance measurements #######################################################

# Coverage
# checking if known OR value of each scenario resides in between the borders of the confidence interval
# and then computing the mean of these checks that result either in 1 or 0.
# The coverage is computed for all scenarios individually

# Coverage can be computed including also incomplete borders such as Inf
est_handled %>% 
  # check if known value is in between the borders (log and normal scale)
  mutate(inside_CI = lower<or & upper>or,
         inside_CI_log = log_lower<log_or & log_upper>log_or) %>%
  # by grouping the coverage can be obtained for each scenario, method and handling approach individually
  group_by(scenario, method, handling) %>%
  # computing the mean of the checks if the known OR value is inside the interval
  summarize(n_performed = sum(!is.na(inside_CI)),
            coverage = sum(inside_CI, na.rm = TRUE)/n_performed,
            coverage_log = sum(inside_CI_log, na.rm = TRUE)/n_performed) %>% 
  merge(zero_perc, by = "scenario", all.x = T) %>% 
  ggplot() +
  geom_point(aes(zero_perc, coverage, color = handling), size = 4) +
  geom_line(aes(zero_perc, coverage, color = handling), size = 1) +
  facet_wrap(~method, ncol = 2)+
  labs(x = "scenario with zeros in %") +
  theme(legend.position = c(0.93,0.1),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.key.size = unit(2, "lines"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))
ggsave("plots/PM_coverage_allBorders.png", width = 13, height = 9)


# Coverage can be computed including only borders that are complete (upper and lower being finite and not NA)
est_handled %>% 
  mutate(across(ends_with("er"), ~replace(., . == Inf | . == -Inf, NA))) %>% 
  mutate(inside_CI = lower<or & upper>or,
         inside_CI_log = log_lower<log_or & log_upper>log_or) %>%
  # by grouping the coverage can be obtained for each scenario, method and handling approach individually
  group_by(scenario, method, handling) %>%
  # computing the mean of the checks if the known OR value is inside the interval
  summarize(n_performed = sum(!is.na(inside_CI)),
            coverage = sum(inside_CI, na.rm = TRUE)/n_performed,
            coverage_log = sum(inside_CI_log, na.rm = TRUE)/n_performed) %>%
  # add the percentages of tables with zeros to the scenarios to be able to order them by the percentage
  merge(zero_perc, by = "scenario", all.x = T) %>% 
  ggplot() +
  geom_point(aes(zero_perc, coverage, color = handling), size = 4) +
  geom_line(aes(zero_perc, coverage, color = handling), size = 1) +
  facet_wrap(~method, ncol = 2)+
  labs(x = "scenario with zeros in %") +
  theme(legend.position = c(0.93,0.1),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.key.size = unit(2, "lines"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))
ggsave("plots/PM_coverage_by_method.png", width = 13, height = 9)


# Closest estimation

# function to find the method providing the closest estimation out of the outcome of several methods 
# Thereby the absolute distance between the estimations and the known log(OR) is computed 
# and the method with the smallest difference determined
# 
# arguments: 
# x:        vector of different estimations on log scale
# log_or:   true, known Odds Ratio on log scale
# method:   vector of names of the methods
# 
# output: a character containing the name of the method that provided the closest estimation to log_or

find_closest_estimation <- function(x, log_or, method) {
  a <- x - log_or
  if(all(is.na(a))){
    closest_method <- NA
  } else {
    closest_index <- which.min(abs(a))
    closest_method <- as.character(method[closest_index])
  }
  return(closest_method)
}

# the closest estimations are checked for all scenarios individually 
# only scenarios that contain a table with zero are included

closest_in_scen_zeros <- est_handled %>% 
  filter(scenario %in% scen_zeros) %>% 
  group_by(scenario, rep, zero, handling) %>% 
  # closest estimation is determined for each table (ACA and CCA)
  summarize(best = find_closest_estimation(log_est, log_or, method), .groups = "drop") %>% 
  filter(!is.na(best)) %>%
  mutate(best=factor(best, levels = c("Fisher", "Midp", "Wald", "Small"))) %>% 
  # scenarios are ordered by percentage of zeros
  merge(scen_zeros_ordered, by = "scenario", all.x = T) %>% 
  mutate(scenario = as.character(scenario),
         scenario = factor(scenario, levels = as.character(scen_zeros_ordered$scenario), 
                           labels = paste(scen_zeros_ordered$scenario, scen_zeros_ordered$zero_perc, sep = ":  ")))
closest_in_scen_zeros %>% 
  ggplot() +
  geom_bar(aes(x =  scenario, fill = best), position = "fill") +
  facet_wrap(~handling) +
  labs(x = "scenario with zeros in %", y = expression(paste("clostest log(", hat(theta), ") in %"))) +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 15),
        legend.key.size = unit(1.5, "lines"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        strip.text = element_text(size = 15)) +
  scale_fill_discrete(name = "method")
ggsave("plots/PM_closest_OnlyZero.png", width = 10, height = 6)

# percentages of methods providing the closest estimation 
# in scenario with highest perc of zeros (nr.19)
closest_in_scen_zeros %>% 
  filter(str_detect(scenario, "^19")) %>% 
  group_by(handling, best) %>% 
  count() %>% 
  group_by(handling) %>% 
  mutate(n_handling = sum(n),
         perc_being_closest = n/n_handling)


# Bias
# computing the bias in each scenario is done by taking the ratio of the estimation to the true value
# and take the mean of this ratio (log scale)
est_handled %>% 
  mutate(ratio_ors = log(estimation_or/or),
         across(starts_with("ratio"), ~replace(., is.infinite(.), NA))) %>% 
  group_by(scenario, method, handling) %>%
  # take the mean of these ratios (mean error)
  summarize(n_performed = sum(!is.na(estimation_or)),
            me = sum(ratio_ors, na.rm = TRUE)/n_performed, .groups = "drop") %>% 
  # add the percentages of tables with zeros to the scenarios to be able to order them by the percentage
  merge(zero_perc, by = "scenario", all.x = T) %>% 
  ggplot() +
  geom_point(aes(zero_perc, me, color = handling), size = 4) +
  geom_line(aes(zero_perc, me, color = handling), size = 1) +
  facet_wrap(~method, ncol = 2) +
  theme(legend.position = c(0.93,0.1),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.key.size = unit(2, "lines"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15)) +
  labs(x = "scenario with zeros in %", y = expression(paste("mean of  ", log(frac(hat(theta),theta)))))
ggsave("plots/PM_bias.png", width = 13, height = 9) 


# MSE
# computing the mean squared error 
est_handled %>% 
  # computing the squared error of each estimation (normal and log scale)
  mutate(dist_squared = (estimation_or - or)^2,
         dist_squared_log = (log_est - log_or)^2) %>% 
  group_by(scenario, method, handling) %>%
  # compute the mean of the squared errors
  summarize(n_performed = sum(!is.na(estimation_or)),
            mse = sum(dist_squared, na.rm = TRUE)/n_performed,
            mse_log = sum(dist_squared_log, na.rm = TRUE)/n_performed) %>%
  # add the percentages of tables with zeros to the scenarios to be able to order them by the percentage
  merge(zero_perc, by = "scenario", all.x = T) %>% 
  ggplot() +
  geom_point(aes(zero_perc, mse_log, color = handling), size = 4) +
  geom_line(aes(zero_perc, mse_log, color = handling), size = 1) +
  facet_wrap(~method, ncol = 2) +
  labs(y = expression(paste("mean of MSE with log(", hat(theta), ")")),
       x = "scenario with zeros in %") +
  theme(legend.position = c(0.06,0.9),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.key.size = unit(2, "lines"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))
ggsave("plots/PM_mseLOG.png" , width = 13, height = 9)
