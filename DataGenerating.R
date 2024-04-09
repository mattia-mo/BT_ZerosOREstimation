# load required packages
if (!requireNamespace("epitools")) install.packages("epitools")
if (!requireNamespace("tidyverse")) install.packages("tidyverse")
if (!requireNamespace("utils")) install.packages("utils")
if (!requireNamespace("checkmate")) install.packages("checkmate")
library(epitools)
library(tidyverse)
library(utils)
library(checkmate)

# setting the seed to create reproducible results
set.seed(567)

# N as simulation sample size
N <-  1

# sample sizes of studies
n_obs <- c(50,100,500)

# true OR values
or <- c(1, 2, 10, 30)

# proportion exposed vs non-exposed
p_x <- c(0.35, 0.5)

# create all possible scenarios based on the predefined parameters
scenarios <- expand_grid(or, n_obs, p_x) %>% 
  mutate(number_scen = row_number())
n_scenarios <- nrow(scenarios)

# choice of estimation methods
methods <- c("midp", "fisher", "wald", "small")
names(methods) <- c("Mid-p", "Fisher", "Wald", "Small")
n_methods <- length(methods)


# function to create one fourfold table, based on the sample size, the true Odds Ratio value and 
# the proportion of individuals being exposed
# 
# arguments: 
# n:    sample size of the observed study (n_obs)
# OR:   true, known Odds Ratio (theta)
# p_x:  proportion of exposure of all individuals 
# 
# output: a table with dimension (2,2) of exposure variable X and outcome variable Y

simuldata <-  function(n, OR, p_x){
  
  # function checks
  assert_int(n, lower = 1)
  assert_number(OR, finite = TRUE)
  assert_number(p_x, lower = 0, upper = 1)
  
  # creating variable X based on proportion of exposures p_x
  times1 <- as.integer(n*p_x)
  times0 <- n - times1
  x <- c(rep(0, times0), rep(1, times1))
  if(n%%2 == 1) x <- c(x,0)
  
  # drawing Y based on the binomial distribution with the inverse logit function
  y <- rbinom(n, 1, plogis(x*log(OR)))
  
  # creating a 2x2 table of X and Y
  df <- data.frame(x = x, y = y)
  tab <- with(df, table(x,y))
  
  # making sure in Y being only 0s or only 1s the 2x2 structure is kept
  if(all(y == 0)){
    tab <- cbind(tab, "1" = c(0,0))
  } else if(all(y == 1)){
    tab <- cbind("0" = c(0,0), tab)
  }
  
  # flipping table to the same appearance as in the thesis 
  new_order <- c("1","0")
  tab <- tab[new_order, new_order]
  
  # output of a 2x2 table
  return(tab)
}


# function to catch errors of expressions
# used for oddsratio.midp() as this function is throwing an error, if there is a zero included, and not only returning NA
# 
# arguments: 
# expr:    expression of a function (with inputs) that needs to be checked for errors
# 
# output: a list of length 2 containing the value of the expression and the error object if an error occurs

catchError <- function(expr) {
  err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }))
  list(value=value, error=err)
}

# function to create one table and the process of estimations for this table 
# A 2x2 table is simulated based on the parameters of the given scenario. This table is checked for zeros. 
# Then all chosen estimation methods are applied to the table. For further combining all tables 
# the results (estimation of OR + confidence intervals) are stored in one line of a dataframe.
# 
# arguments: 
# rep:        number the table should get
# n_obs:      sample size of the observed study
# p_x:        proportion of exposure of all individuals, default is 0.5
# oddsratio:  true, known Odds Ratio (theta)
# 
# output: a dataframe with dimension (1,18) containing the parameters of the scenario

onerep <- function(rep, n_obs, p_x = 0.5, oddsratio) {
  
  # function check
  assert_number(rep, lower = 1)
  
  # simulate the table
  table_of_rep <- simuldata(n_obs, oddsratio, p_x)
  
  # check if zero in table
  zero <- any(table_of_rep == 0)
  
  # applying the estimation methods on table
  
  result_wald <- oddsratio.wald(table_of_rep)
  or_wald <- result_wald$measure[2,1]
  lower_wald <- result_wald$measure[2,2]
  upper_wald <- result_wald$measure[2,3]
  
  result_fisher <- oddsratio.fisher(table_of_rep)
  or_fisher <- result_fisher$measure[2,1]
  lower_fisher <-  result_fisher$measure[2,2]
  upper_fisher <-  result_fisher$measure[2,3]
  
  
  result_small <- suppressWarnings(oddsratio.small(table_of_rep))
  # suppression of warning as a check for the warning of not using 1/OR probably is implemented 
  # incorrectly in oddsratio.small (see https://rdrr.io/cran/epitools/src/R/oddsratio.small.R)
  or_small <- result_small$measure[2,1]
  lower_small <- result_small$measure[2,2]
  upper_small <- result_small$measure[2,3]
  
  
  # as the midp function is giving an error immediately, 
  # the error needs to be caught and denoted in the data as NA
  err <- catchError(oddsratio.midp(table_of_rep))
  if(is.null(err$error)) {
    result_midp <- oddsratio.midp(table_of_rep)
    or_midp <- result_midp$measure[2,1]
    lower_midp <- result_midp$measure[2,2]
    upper_midp <- result_midp$measure[2,3]
  } else {
    or_midp <- NA
    lower_midp <- NA
    upper_midp <- NA
  }
  
  # prepare output of the parameters, OR and CI
  out <- data.frame(
    rep = rep,
    or = oddsratio,
    n_obs = n_obs,
    p_x = p_x,
    table = NA,
    zero = zero,
    or_wald = or_wald, or_fisher = or_fisher, or_small = or_small, or_midp = or_midp,
    lower_wald = lower_wald, lower_fisher = lower_fisher, lower_small = lower_small, lower_midp = lower_midp,
    upper_wald = upper_wald, upper_fisher = upper_fisher, upper_small = upper_small, upper_midp = upper_midp
  )
  # also adding the table in its original form to the dataframe
  out$table <- list(table_of_rep)
  
  # returning 1 row with all information (table + estimations)
  return(out)
}


# function to conduct multiple repetitions of one scenario inclusive storing the random number states
# This function is simulating tables + their OR estimations of a whole scenario based on the parameters of the given scenario. 
# Thereby it is conducting multiple times single rows containing the a simulated table plus its OR estimations.
# For each of these single rows the random number states are stored as Morris(2019) is suggesting.
# 
# arguments: 
# n_obs:      sample size of the scenario
# or:         true, known Odds Ratio (theta)
# p_x:        proportion of exposure of all individuals, default is 0.5
# n_rep:      number of repetitions (equal to how many tables should be simulated)
# number_scen:number of the scenario that the tables should have
# 
# output: a list containing the states and all estimations of the scenario

multreps <- function(n_obs, or, p_x, n_rep, number_scen){
  
  # function check
  assert_number(number_scen, lower = 1)
  
  # predefine the output
  estimates <- data.frame(matrix(ncol = 6 + (n_methods*3), nrow = n_rep))
  x <- c("rep", "or", "n_obs", "p_x", "table", "zero", 
         "or_wald", "or_fisher", "or_small", "or_midp", 
         "lower_wald", "lower_fisher", "lower_small", "lower_midp", 
         "upper_wald", "upper_fisher", "upper_small", "upper_midp")
  colnames(estimates) <- x
  states <- matrix(ncol = 627, nrow = n_rep)
  
  # indicator of how far in the simulation process the function is
  cat("Current scenario:", number_scen, "\n")
  
  # Run all repetions of the scenario
  for (r in 1:n_rep) {
    states[r, ] <- c(rep = r, .Random.seed)
    estimates[r, ] <- onerep(rep = r, n_obs = n_obs, oddsratio = or, p_x = p_x)
  }
  
  return(list(estimates, states))
}

# creating multiple repetitions on all predefined scenarios and store the results (states inclusive)
# perform N estimations for all methods for each scenario with the parameters the scenario is holding
k <- scenarios %>% 
  apply(.,1, function(scen){
    multreps(n_obs = scen[2], or = scen[1], p_x = scen[3], n_rep = N, number_scen = scen[4])
  }) 

# get estimations out of the apply-result and store them with the allocated scenario number
estimates_all <- lapply(k, `[[`, 1) %>% 
  bind_rows() %>% 
  group_by(or, n_obs, p_x) %>% 
  mutate(scenario = cur_group_id(), .before = "rep")

# get states out of the apply-result and store them with the allocated scenario number
states_all <- lapply(k, `[[`, 2) %>%
  do.call(rbind,.) %>% 
  as.data.frame() %>% 
  rename(rep = V1) %>% 
  mutate(scenario = rep(1:n_scenarios, each = N), .before = "rep")


# divide states data frame into smaller data frames 
# so that one data frame contains all states of one scenario to avoid storing large files
states_list <- split(states_all, states_all$scenario)
for(i in seq_len(n_scenarios)) {
  filename <- paste0("states/states_", i, ".rds")
  states_of_scenario <- states_list[[i]]
  saveRDS(states_of_scenario, file = filename)
}

# save estimates file
saveRDS(estimates_all, file = "estimates.rds")



### how to create a specific table based on the states

# example: recreating table number i of scenario k
i <- 873
k <- 7

row_of_example <- estimates_all %>% filter(scenario == k & rep == i)
# table that should be recreated
row_of_example$table[[1]]

# recreating the table by using the stored states
param_of_example <- as.vector(scenarios[k,])
.Random.seed <- as.integer(states_list[[k]][i,-c(1, 2)])
recreated_row <- onerep(i, param_of_example$n_obs, param_of_example$p_x, param_of_example$or)
recreated_table <- recreated_row$table[[1]]
recreated_table

