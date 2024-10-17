
# Setup -------------------------------------------------------------------

# use parent directory ".\\Code_and_Data"
# to be able to run the code

#setwd(".\\Code_and_Data")

# Packages ----------------------------------------------------------------

library(tidyverse)

# Explanation of variables ------------------------------------------------

# props: set probabilities for categories and groups to create the multinomial data
# phi: set dispersion parameter
# m: cluster size
# b: sample size of clusters
# comp: Tukey or Dunnet comparisons with glht
# modeltype: method used for model fitting and dispersion calculation
# dfu: whether degrees of freedom are specified for singlep adjustment (either mvt- or mvn-distribution)
# nsim: total number of simulations per simulation run (not used for sum(sim))
# sim: number of simulations for each row (always 1)
# C: number of categories
# G: number of groups
# minpm: minimum of props*m
# minp: minimum of props
# maxp: maximum of props
# HAinc: whether an alternative hypothesis is included or not
# cH0: how much null hypotheses are included
# cHA: how much alternative hypotheses are included
# cmfz:  count how often less coefficients were estimated than should be estimated (dropped columns)
# cFWER: count how many times p < 0.05 even though hypothesis is under null
# cPower: count how many times p < 0.05 when hypothesis is under altenrative
# cgPower: count how many times p < 0.05 for all hypotheses
# cConf: count how many times trueprob lies outside interval
# disphat: estimated dispersion parameter
# cerrnas: count error: "NAs found in the working weights variable 'wz'"
# cerrfin: count error: "Some elements in the working weights variable 'wz' are not finite"
# cwarwz: count warning: "diagonal elements of the working weights variable 'wz' have been replaced by 1.819e-12"
# cwarcon: count warning: "convergence not obtained in 30 IRLS iterations"
# cwarany: count warning: any other warning
# cwarzc: count warning: check if there are any deleted columns due to zero counts  "Deleted 2 columns of the response matrix due to zero counts"
# cabseps: in glht function: Completion with error
# cglhtany: n glht function: count any other warning/error
# biasprob: calculated bias and then mean of biases of probs in simulation run
# biasabsprob: calculated absolute bias and the mean ob abs biases of probs in simulation run
# estprobs: vector with tracked estimated probabilities (ratios between groups)
# trueprobs: vector with true probabilities (ratios between groups)


# Main Simulation ---------------------------------------------------------

setwd(".\\Results\\main")

v_reslist <- list.files(pattern = ".rds")

l_dat <- list()

for (i in 1:length(v_reslist)) {
  
  print(v_reslist[i])
  df_temp <- readRDS(v_reslist[i])
  
  df_temp_NAs_cmfz <- df_temp %>% group_by(props,
                                           phi,
                                           b,
                                           m,
                                           comp,
                                           modeltype,
                                           C,
                                           G,
                                           minpm,
                                           minp,
                                           maxp,
                                           HAinc,
                                           cH0,
                                           cHA) %>%
    summarise(
      cNA = as.integer(across(cFWER, ~ sum(is.na(
        .
      )))),
      cmfz = sum(cmfz, na.rm = TRUE),
      totsim = sum(sim, na.rm = TRUE),
      cerrnas = sum(cerrnas, na.rm = TRUE),
      cerrfin = sum(cerrfin, na.rm = TRUE),
      cwarwz = sum(cwarwz, na.rm = TRUE),
      cwarcon = sum(cwarcon, na.rm = TRUE),
      cwarany = sum(cwarany, na.rm = TRUE),
      cwarzc = sum(cwarzc, na.rm = TRUE),
      cabseps = sum(cabseps, na.rm = TRUE),
      cglhtany = sum(cglhtany, na.rm = TRUE)
    )
  
  df_temp_propbias <- df_temp %>% group_by(props,
                                           phi,
                                           b,
                                           m,
                                           comp,
                                           modeltype,
                                           C,
                                           G,
                                           minpm,
                                           minp,
                                           maxp,
                                           HAinc,
                                           cH0,
                                           cHA) %>% drop_na(cFWER) %>%
    separate_longer_delim(c(estprobs, trueprobs), delim = ",") %>%
    mutate(
      estprobs = as.numeric(estprobs),
      trueprobs = as.numeric(trueprobs),
      estprobsmean = mean(estprobs),
      estbias = estprobs - trueprobs,
      estbiasSE = (estprobs - estprobsmean) ^ 2,
      estbiasMSE = (estprobs - trueprobs) ^ 2,
      phib = disphat - phi
    ) %>%
    summarise(
      probbias = sum(estbias),
      probbiasSE = sum(estbiasSE),
      probbiasMSE = sum(estbiasMSE),
      phibiasalt = sum(phib)
    )
  
  df_temp_sum <- df_temp %>% group_by(props,
                                      phi,
                                      b,
                                      m,
                                      comp,
                                      modeltype,
                                      C,
                                      G,
                                      minpm,
                                      minp,
                                      maxp,
                                      HAinc,
                                      cH0,
                                      cHA) %>% drop_na(cFWER) %>%
    mutate(
      phimean = mean(disphat),
      estphibias = disphat - phi,
      estphibiasSE = (disphat - phimean) ^ 2,
      estphibiasMSE = (disphat - phi) ^ 2
    ) %>%
    summarise(
      sim = sum(sim, na.rm = TRUE),
      cFWER = sum(cFWER, na.rm = TRUE),
      cPower = sum(cPower, na.rm = TRUE),
      cgPower = sum(cgPower, na.rm = TRUE),
      cConf = sum(cConf, na.rm = TRUE),
      phibias = (1 / sim) * sum(estphibias, na.rm = TRUE),
      phibiasSE = sqrt((1 / (sim - 1)) * sum(estphibiasSE, na.rm = TRUE)),
      phibiasMSE = (1 / sim) * sum(estphibiasMSE, na.rm = TRUE),
      
    )
  
  df_temp_fin <- cbind(df_temp_propbias[, -c(1:14)], df_temp_sum)
  df_temp_fin <-
    merge(
      df_temp_fin,
      df_temp_NAs_cmfz,
      by = c(
        "props",
        "phi",
        "b",
        "m",
        "comp",
        "modeltype",
        "C",
        "G",
        "minpm",
        "minp",
        "maxp",
        "HAinc",
        "cH0",
        "cHA"
      )
    )
  
  df_temp_fin <- df_temp_fin %>% group_by(props,
                                          phi,
                                          b,
                                          m,
                                          comp,
                                          modeltype,
                                          C,
                                          G,
                                          minpm,
                                          minp,
                                          maxp,
                                          HAinc,
                                          cH0,
                                          cHA) %>%
    mutate(
      probbias = (1 / sim) * probbias,
      probbiasSE = (1 / sim) * probbiasSE,
      probbiasMSE = (1 / sim) * probbiasMSE,
      FWER_woNA = cFWER / sim,
      FWER_wiNA = cFWER / totsim,
      Power_woNA = cPower / sim,
      Power_wiNA = cPower / totsim,
      gPower_woNA = cgPower / sim,
      gPower_wiNA = cgPower / totsim,
      covprob_woNA = cConf / sim,
      covprob_wiNA = cConf / totsim,
    )
  
  l_dat <- c(l_dat, list(df_temp_fin))
  
}


l_dat_main <- bind_rows(l_dat)

saveRDS(l_dat_main, "mult_sim_main.rds")


## Zero handling Data -----------------------------------------------------

#setwd(".\\Code_and_Data")
setwd(".\\Results\\zero_handling")

v_reslist <- list.files(pattern = ".rds")

l_dat <- list()

for (i in 1:length(v_reslist)) {
  
  print(v_reslist[i])
  df_temp <- readRDS(v_reslist[i])
  df_temp <- df_temp[,-19]
  
  df_temp_NAs_cmfz <- df_temp %>% group_by(props,
                                           phi,
                                           b,
                                           m,
                                           comp,
                                           modeltype,
                                           dfu,
                                           dataset,
                                           C,
                                           G,
                                           minpm,
                                           minp,
                                           maxp,
                                           HAinc,
                                           cH0,
                                           cHA) %>%
    summarise(
      cNA = as.integer(across(cFWER, ~ sum(is.na(
        .
      )))),
      cmfz = sum(cmfz, na.rm = TRUE),
      totsim = sum(sim, na.rm = TRUE),
      cerrnas = sum(cerrnas, na.rm = TRUE),
      cerrfin = sum(cerrfin, na.rm = TRUE),
      cwarwz = sum(cwarwz, na.rm = TRUE),
      cwarcon = sum(cwarcon, na.rm = TRUE),
      cwarany = sum(cwarany, na.rm = TRUE),
      cwarzc = sum(cwarzc, na.rm = TRUE),
      cabseps = sum(cabseps, na.rm = TRUE),
      cglhtany = sum(cglhtany, na.rm = TRUE)
    )
  
  df_temp_propbias <- df_temp %>% group_by(props,
                                           phi,
                                           b,
                                           m,
                                           comp,
                                           modeltype,
                                           dfu,
                                           dataset,
                                           C,
                                           G,
                                           minpm,
                                           minp,
                                           maxp,
                                           HAinc,
                                           cH0,
                                           cHA) %>% drop_na(cFWER) %>%
    separate_longer_delim(c(estprobs, trueprobs), delim = ",") %>%
    mutate(
      estprobs = as.numeric(estprobs),
      trueprobs = as.numeric(trueprobs),
      estprobsmean = mean(estprobs),
      estbias = estprobs - trueprobs,
      estbiasSE = (estprobs - estprobsmean) ^ 2,
      estbiasMSE = (estprobs - trueprobs) ^ 2,
      phib = disphat - phi
    ) %>%
    summarise(
      probbias = sum(estbias),
      probbiasSE = sum(estbiasSE),
      probbiasMSE = sum(estbiasMSE),
      phibiasalt = sum(phib)
    )
  
  df_temp_sum <- df_temp %>% group_by(props,
                                      phi,
                                      b,
                                      m,
                                      comp,
                                      modeltype,
                                      dfu,
                                      dataset,
                                      C,
                                      G,
                                      minpm,
                                      minp,
                                      maxp,
                                      HAinc,
                                      cH0,
                                      cHA) %>% drop_na(cFWER) %>%
    mutate(
      phimean = mean(disphat),
      estphibias = disphat - phi,
      estphibiasSE = (disphat - phimean) ^ 2,
      estphibiasMSE = (disphat - phi) ^ 2
    ) %>%
    summarise(
      sim = sum(sim, na.rm = TRUE),
      cFWER = sum(cFWER, na.rm = TRUE),
      cPower = sum(cPower, na.rm = TRUE),
      cgPower = sum(cgPower, na.rm = TRUE),
      cConf = sum(cConf, na.rm = TRUE),
      phibias = (1 / sim) * sum(estphibias, na.rm = TRUE),
      phibiasSE = sqrt((1 / (sim - 1)) * sum(estphibiasSE, na.rm = TRUE)),
      phibiasMSE = (1 / sim) * sum(estphibiasMSE, na.rm = TRUE),
      
    )
  
  df_temp_fin <- cbind(df_temp_propbias[, -c(1:16)], df_temp_sum)
  df_temp_fin <-
    merge(
      df_temp_fin,
      df_temp_NAs_cmfz,
      by = c(
        "props",
        "phi",
        "b",
        "m",
        "comp",
        "modeltype",
        "dfu",
        "dataset",
        "C",
        "G",
        "minpm",
        "minp",
        "maxp",
        "HAinc",
        "cH0",
        "cHA"
      )
    )
  
  df_temp_fin <- df_temp_fin %>% group_by(props,
                                          phi,
                                          b,
                                          m,
                                          comp,
                                          modeltype,
                                          dfu,
                                          dataset,
                                          C,
                                          G,
                                          minpm,
                                          minp,
                                          maxp,
                                          HAinc,
                                          cH0,
                                          cHA) %>%
    mutate(
      probbias = (1 / sim) * probbias,
      probbiasSE = (1 / sim) * probbiasSE,
      probbiasMSE = (1 / sim) * probbiasMSE,
      FWER_woNA = cFWER / sim,
      FWER_wiNA = cFWER / totsim,
      Power_woNA = cPower / sim,
      Power_wiNA = cPower / totsim,
      gPower_woNA = cgPower / sim,
      gPower_wiNA = cgPower / totsim,
      covprob_woNA = cConf / sim,
      covprob_wiNA = cConf / totsim,
    )
  
  l_dat <- c(l_dat, list(df_temp_fin))
  
}

l_dat_zero <- bind_rows(l_dat)

saveRDS(l_dat_zero, "mult_sim_zero.rds")

