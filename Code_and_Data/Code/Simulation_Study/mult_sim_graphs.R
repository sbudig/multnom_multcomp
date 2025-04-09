
# Setup -------------------------------------------------------------------

# use parent directory ".\\Code_and_Data"
# to be able to run the code

#setwd(".\\Code_and_Data")

# Packages ----------------------------------------------------------------

library(tidyverse)
library(patchwork)
#library(viridis)

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

# Data --------------------------------------------------------------------

## Main Simulation --------------------------------------------------------

l_dat_main_og <- readRDS(".\\Results\\main\\mult_sim_main.rds")

l_dat_main <- l_dat_main_og %>% mutate(
  modeltypef = fct_recode(
    modeltype,
    "Multinomial" = "multinomial",
    "Pearson" = "pearson",
    "Afroz" = "afroz",
    "Farrington" = "farrington",
    "Deviance" = "deviance",
    "DM_VGAM" = "DM_VGAM",
    "DM_MGLM" = "DM_MGLM"
  ),
  minpmb = minpm*b
)

l_dat_main$modeltypef <- factor(l_dat_main$modeltypef,
                                levels = c("Multinomial", "Deviance", "Pearson",
                                           "Farrington", "Afroz", "DM_VGAM", "DM_MGLM"))


l_dat_trans_main <- l_dat_main %>%
  ungroup() %>% select(-modeltype) %>%
  group_by(modeltypef) %>%
  pivot_wider(
    names_from = modeltypef,
    values_from = c(probbias, probbiasSE, probbiasMSE,
                    FWER_woNA, FWER_wiNA,
                    Power_woNA,Power_wiNA,
                    gPower_woNA,gPower_wiNA,
                    cpPower_woNA,cpPower_wiNA,
                    covprob_woNA,covprob_wiNA,
                    sim, totsim,
                    cFWER, cPower, cgPower, cpPower, cConf,
                    cNA, cmfz, cerrnas, cerrfin, cwarwz, cwarcon, cwarany, cwarzc, cabseps, cglhtany)
  ) %>% 
  type.convert(as.is=TRUE)

l_dat_trans_main_power <- l_dat_trans_main %>% filter(cHA > 0)


## Zero handling Data -----------------------------------------------------

l_dat_zero <- readRDS(".\\Results\\zero_handling\\mult_sim_zero.rds")

l_dat_trans_zero <- l_dat_zero %>% 
  group_by(modeltype) %>%
  pivot_wider(
    names_from = modeltype:dataset,
    values_from = c(probbias, probbiasSE, probbiasMSE,
                    FWER_woNA, FWER_wiNA,
                    Power_woNA,Power_wiNA,
                    gPower_woNA,gPower_wiNA,
                    cpPower_woNA,cpPower_wiNA,
                    covprob_woNA,covprob_wiNA,
                    sim, totsim,
                    cFWER, cPower, cgPower, cpPower, cConf,
                    cNA, cmfz, cerrnas, cerrfin, cwarwz, cwarcon, cwarany, cwarzc, cabseps, cglhtany)
  ) %>% 
  type.convert(as.is=TRUE)

l_dat_trans_zero_power <- l_dat_trans_zero %>% filter(cHA > 0)

# Graphs ------------------------------------------------------------------

disppalette <- c('#a7e3d7','#58b4b9','#2880a0','#1d4c81','#1a1662') # picked sequential cerulean

# Figure 1 - main - FWER --------------------------------------------------

plot_main_fwer_overall <- ggplot(l_dat_main[order(l_dat_main$props),], aes(x = log10(minpmb), y = FWER_wiNA))+ 
  geom_point(aes(colour = factor(phi), shape = comp), size = 2.5)+ 
  facet_grid(~modeltypef)+
  geom_hline(yintercept=0.05)+
  geom_hline(yintercept=0.03733, linetype="dashed") +
  geom_hline(yintercept=0.06539 , linetype="dashed") +
  theme_bw()+ 
  theme_bw(base_size = 18)+ 
  theme(legend.position="bottom")+
  xlab(expression(paste(log10(min( (m[g]*pi[gc]) )))))+
  ylab("Family-wise error rate")+ 
 guides(shape=guide_legend(title="Group comparison contrasts"))+
  #        fill = guide_legend(override.aes = list(shape = 21, color = "black"))
  #        )+
  labs( colour =  expression(paste("Dispersion ", (phi))))+
  #scale_fill_viridis_d(begin = 1, end = 0)+
  #scale_shape_manual(values = c(21,22)) +
  scale_colour_manual(values = disppalette) +
  ylim(0,0.4)
  
plot_main_fwer_overall

ggsave(
  ".\\Figures\\Budig_SimInfMult_Figure_1.eps",
  plot = plot_main_fwer_overall ,
  width = 14,
  height = 7,
  device = "eps",
  dpi = 900
)

# Figure 2 - main - Power comparison --------------------------------------

plot_main_power_pearson_afroz <-
  ggplot(l_dat_trans_main_power[order(l_dat_trans_main_power$props),],
         aes(x = cpPower_wiNA_Afroz, y = cpPower_wiNA_Pearson)) +
  geom_point(aes(colour = factor(phi)), size = 2.5) +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw(base_size = 18) +
  labs(x = "Power - Afroz",
       y = "Power - Pearson",
       colour = expression(paste("Dispersion ", (phi))))+
  #scale_fill_viridis_d(begin = 1, end = 0)+
  scale_colour_manual(values = disppalette)

plot_main_power_pearson_farrington <-
  ggplot(l_dat_trans_main_power[order(l_dat_trans_main_power$props),],
         aes(x = cpPower_wiNA_Afroz, y = cpPower_wiNA_Farrington)) +
  geom_point(aes(colour = factor(phi)), size = 2.5) +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw(base_size = 18)  +
  labs(x = "Power - Afroz",
       y = "Power - Farrington",
       colour = expression(paste("Dispersion ", (phi))))+
  #scale_fill_viridis_d(begin = 1, end = 0)+
  scale_colour_manual(values = disppalette)

plot_main_power_afroz_DM_MGLM <-
  ggplot(l_dat_trans_main_power[order(l_dat_trans_main_power$props),],
         aes(x = cpPower_wiNA_Afroz, y = cpPower_wiNA_DM_MGLM)) +
  geom_point(aes(colour = factor(phi)), size = 2.5) +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw(base_size = 18) +
  labs(x = "Power - Afroz",
       y = "Power - MGLM_DM",
       colour = expression(paste("Dispersion ", (phi))))+
  #scale_fill_viridis_d(begin = 1, end = 0)+
  scale_colour_manual(values = disppalette)

plot_main_power_pearson_afroz_farrington_dm_mglm <-
  plot_main_power_pearson_afroz + 
  plot_main_power_pearson_farrington + 
  plot_main_power_afroz_DM_MGLM+ 
  plot_layout(guides = "collect") &
  theme(legend.position='bottom')

plot_main_power_pearson_afroz_farrington_dm_mglm

ggsave(
  ".\\Figures\\Budig_SimInfMult_Figure_2.eps",
  plot = plot_main_power_pearson_afroz_farrington_dm_mglm ,
  width = 14,
  height = 5.49,
  device = "eps",
  dpi = 900
)


# Figure 5 - zero handling - FWER comparison ------------------------------

plot_zero_fwer_afroz <- 
  ggplot(l_dat_trans_zero[order(l_dat_trans_zero$props),],
         aes(x = FWER_wiNA_afroz_modeldf_original,
             y = FWER_wiNA_afroz_modeldf_onerow))+ 
  #geom_point(aes(colour = factor(phi), shape = factor(m)), size = 1.9)+ 
  geom_point(aes(colour = factor(phi)), size = 2.5)+ 
  geom_abline(slope=1, intercept=0) +
  geom_hline(yintercept=0.05, linetype="solid") +
  geom_hline(yintercept=0.03733, linetype="dashed") +
  geom_hline(yintercept=0.06539 , linetype="dashed") +
  geom_vline(xintercept=0.05, linetype="solid") +
  geom_vline(xintercept=0.03733, linetype="dashed") +
  geom_vline(xintercept=0.06539 , linetype="dashed")+
  xlim(0,0.15) + ylim(0,0.15)+
  labs(x = "FWER - Afroz - OG",
       y = "FWER - Afroz - AO",
       colour = expression(paste("Dispersion ", (phi))))+
  #shape = expression(paste("Cluster size ", (m[gb]))))+
  theme_bw(base_size = 18)+ 
  scale_colour_manual(values = disppalette)
  #scale_fill_viridis_d(begin = 1, end = 0)

plot_zero_fwer_pearson <- 
  ggplot(l_dat_trans_zero[order(l_dat_trans_zero$props),],
         aes(x = FWER_wiNA_pearson_modeldf_original,
             y = FWER_wiNA_pearson_modeldf_onerow))+ 
  #geom_point(aes(colour = factor(phi), shape = factor(m)), size = 1.9)+ 
  geom_point(aes(colour = factor(phi)), size = 2.5)+ 
  geom_abline(slope=1, intercept=0) +
  geom_hline(yintercept=0.05, linetype="solid") +
  geom_hline(yintercept=0.03733, linetype="dashed") +
  geom_hline(yintercept=0.06539 , linetype="dashed") +
  geom_vline(xintercept=0.05, linetype="solid") +
  geom_vline(xintercept=0.03733, linetype="dashed") +
  geom_vline(xintercept=0.06539 , linetype="dashed")+
  xlim(0,0.15) + ylim(0,0.15)+
  labs(x = "FWER - Pearson - OG",
       y = "FWER - Pearson - AO",
       colour = expression(paste("Dispersion ", (phi))))+
  #shape = expression(paste("Cluster size ", (m[gb]))))+
  theme_bw(base_size = 18)+ 
  scale_colour_manual(values = disppalette)
  #scale_fill_viridis_d(begin = 1, end = 0)

plot_zero_fwer_dirmult <- 
  ggplot(l_dat_trans_zero[order(l_dat_trans_zero$props),],
         aes(x = FWER_wiNA_DM_MGLM_modeldf_original,
             y = FWER_wiNA_DM_MGLM_modeldf_onerow))+ 
  #geom_point(aes(colour = factor(phi), shape = factor(m)), size = 1.9)+ 
  geom_point(aes(colour = factor(phi)), size = 2.5)+ 
  geom_abline(slope=1, intercept=0) +
  geom_hline(yintercept=0.05, linetype="solid") +
  geom_hline(yintercept=0.03733, linetype="dashed") +
  geom_hline(yintercept=0.06539 , linetype="dashed") +
  geom_vline(xintercept=0.05, linetype="solid") +
  geom_vline(xintercept=0.03733, linetype="dashed") +
  geom_vline(xintercept=0.06539 , linetype="dashed")+
  xlim(0,0.15) + ylim(0,0.15)+
  labs(x = "FWER - DM_MGLM - OG",
       y = "FWER - DM_MGLM - AO",
       colour = expression(paste("Dispersion ", (phi))))+
  #shape = expression(paste("Cluster size ", (m[gb]))))+
  theme_bw(base_size = 18)+ 
  scale_colour_manual(values = disppalette)
  #scale_fill_viridis_d(begin = 1, end = 0)

plot_zero_fwer_afr_pear_dirm <- plot_zero_fwer_pearson + 
  plot_zero_fwer_afroz + 
  plot_zero_fwer_dirmult+ 
  plot_layout(guides = "collect")&
  theme(legend.position='bottom')

plot_zero_fwer_afr_pear_dirm

ggsave(
  ".\\Figures\\Budig_SimInfMult_Figure_3.eps",
  plot = plot_zero_fwer_afr_pear_dirm ,
  width = 14,
  height = 5.49,
  device = "eps",
  dpi = 900
)

# Figure 6 - zero handling - Power comparison -----------------------------

plot_zero_power_pearson <- ggplot(
  l_dat_trans_zero_power[order(l_dat_trans_zero_power$props),],
  aes(x = cpPower_wiNA_pearson_modeldf_original,
      y = cpPower_wiNA_pearson_modeldf_onerow)) +
  geom_point(aes(colour = factor(phi)), size = 2.5) +
  geom_abline(slope = 1, intercept = 0) +
  xlim(0, 1) + ylim(0, 1) +
  labs(x = "Power - Pearson - OG",
       y = "Power - Pearson - AO",
       colour = expression(paste("Dispersion ", (phi))))+
  theme_bw(base_size = 18)+ 
  scale_colour_manual(values = disppalette)
  #scale_fill_viridis_d(begin = 1, end = 0)

plot_zero_power_afroz <- ggplot(
  l_dat_trans_zero_power[order(l_dat_trans_zero_power$props),],
  aes(x = cpPower_wiNA_afroz_modeldf_original,
      y = cpPower_wiNA_afroz_modeldf_onerow)) +
  geom_point(aes(colour = factor(phi)), size = 2.5) +
  geom_abline(slope = 1, intercept = 0) +
  xlim(0, 1) + ylim(0, 1) +
  labs(x = "Power - Afroz - OG",
       y = "Power - Afroz - AO",
       colour = expression(paste("Dispersion ", (phi))))+
  theme_bw(base_size = 18)+ 
  scale_colour_manual(values = disppalette)
  #scale_fill_viridis_d(begin = 1, end = 0)

plot_zero_power_dirmult <- ggplot(
  l_dat_trans_zero_power[order(l_dat_trans_zero_power$props),],
  aes(x = cpPower_wiNA_DM_MGLM_modeldf_original,
      y = cpPower_wiNA_DM_MGLM_modeldf_onerow)) +
  geom_point(aes(colour = factor(phi)), size = 2.5) +
  geom_abline(slope = 1, intercept = 0) +
  xlim(0, 1) + ylim(0, 1) +
  labs(x = "Power - DM_MGLM - OG",
       y = "Power - DM_MGLM - AO",
       colour = expression(paste("Dispersion ", (phi))))+
  theme_bw(base_size = 18)+ 
  scale_colour_manual(values = disppalette)
  #scale_fill_viridis_d(begin = 1, end = 0)

plot_zero_power_afr_pear_dirm <-
  plot_zero_power_pearson +
  plot_zero_power_afroz +
  plot_zero_power_dirmult +
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom')

plot_zero_power_afr_pear_dirm

ggsave(
  ".\\Figures\\Budig_SimInfMult_Figure_4.eps",
  plot = plot_zero_power_afr_pear_dirm ,
  width = 14,
  height = 5.49,
  device = "eps",
  dpi = 900
)


# Additional Graphs -------------------------------------------------------

## main - FWER comparison -------------------------------------------------

plot_main_fwer_pearson_afroz <- 
  ggplot(l_dat_trans_main, aes(x = FWER_wiNA_Afroz, y = FWER_wiNA_Pearson))+ 
  geom_point(aes(colour = factor(phi)),  size = 2.5)+ 
  geom_abline(slope=1, intercept=0) +
  geom_hline(yintercept=0.05, linetype="solid") +
  geom_hline(yintercept=0.03733, linetype="dashed") +
  geom_hline(yintercept=0.06539 , linetype="dashed") +
  geom_vline(xintercept=0.05, linetype="solid") +
  geom_vline(xintercept=0.03733, linetype="dashed") +
  geom_vline(xintercept=0.06539 , linetype="dashed")+
  xlim(0,0.15) + ylim(0,0.15)+
  labs(x = "FWER - Afroz",
       y = "FWER - Pearson",
       colour = expression(paste("Dispersion ", (phi))))+ 
  theme_bw(base_size = 18)+
  scale_colour_manual(values = disppalette) 
#scale_fill_viridis_d(begin = 1, end = 0) 

plot_main_fwer_pearson_farrington <-
  ggplot(l_dat_trans_main, aes(x = FWER_wiNA_Farrington, y = FWER_wiNA_Pearson))+ 
  geom_point(aes(colour = factor(phi)), size = 2.5)+ 
  geom_abline(slope=1, intercept=0) +
  geom_hline(yintercept=0.05, linetype="solid") +
  geom_hline(yintercept=0.03733, linetype="dashed") +
  geom_hline(yintercept=0.06539 , linetype="dashed") +
  geom_vline(xintercept=0.05, linetype="solid") +
  geom_vline(xintercept=0.03733, linetype="dashed") +
  geom_vline(xintercept=0.06539 , linetype="dashed")+
  xlim(0,0.15) + ylim(0,0.15)+
  labs(x = "FWER - Farrington",
       y = "FWER - Pearson",
       colour = expression(paste("Dispersion ", (phi))))+ 
  theme_bw(base_size = 18)+
  #scale_fill_viridis_d(begin = 1, end = 0) +
  scale_colour_manual(values = disppalette)

plot_main_fwer_afroz_farrington <-
  ggplot(l_dat_trans_main, aes(x = FWER_wiNA_Farrington, y = FWER_wiNA_Afroz))+ 
  geom_point(aes(colour = factor(phi)),  size = 2.5)+ 
  geom_abline(slope=1, intercept=0) +
  geom_hline(yintercept=0.05, linetype="solid") +
  geom_hline(yintercept=0.03733, linetype="dashed") +
  geom_hline(yintercept=0.06539 , linetype="dashed") +
  geom_vline(xintercept=0.05, linetype="solid") +
  geom_vline(xintercept=0.03733, linetype="dashed") +
  geom_vline(xintercept=0.06539 , linetype="dashed")+
  xlim(0,0.15) + ylim(0,0.15)+
  labs(x = "FWER - Farrington",
       y = "FWER - Afroz",
       colour = expression(paste("Dispersion ", (phi))))+ 
  theme_bw(base_size = 18)+
  #scale_fill_viridis_d(begin = 1, end = 0) +
  scale_colour_manual(values = disppalette)

plot_main_fwer_pearson_afroz_fletcher <- plot_main_fwer_pearson_afroz + 
  plot_main_fwer_pearson_farrington + 
  plot_main_fwer_afroz_farrington + 
  plot_layout(guides = "collect") &
  theme(legend.position='bottom')

plot_main_fwer_pearson_afroz_fletcher


## main - Coverage probability comparison ---------------------------------

# Pearson vs Afroz
plot_main_conf_pearson_afroz <- ggplot(l_dat_trans_main[order(l_dat_trans_main$props),], 
                                       aes(x = 1-covprob_wiNA_Afroz, 
                                           y = 1-covprob_wiNA_Pearson))+ 
  geom_point(aes(colour = factor(phi)))+ 
  geom_hline(yintercept=0.95)+
  theme_bw(base_size = 18)+
  geom_hline(yintercept=0.95, linetype="solid") +
  geom_hline(yintercept=0.9346095, linetype="dashed") +
  geom_hline(yintercept=0.9626646, linetype="dashed") +
  geom_vline(xintercept=0.95, linetype="solid") +
  geom_vline(xintercept=0.9346095, linetype="dashed") +
  geom_vline(xintercept=0.9626646, linetype="dashed")+
  labs(x = "Coverage Prob. - Afroz",
       y = "Coverage Prob. - Pearson",
       color = expression(paste("Dispersion ", (phi))))+ 
  scale_colour_manual(values = disppalette)

# Pearson vs farrington 
plot_main_conf_pearson_farrington  <- ggplot(l_dat_trans_main[order(l_dat_trans_main$props),], 
                                          aes(x = 1-covprob_wiNA_Farrington, 
                                              y = 1-covprob_wiNA_Pearson))+ 
  geom_point(aes(colour = factor(phi)))+ 
  geom_hline(yintercept=0.95)+
  theme_bw(base_size = 18)+
  geom_hline(yintercept=0.95, linetype="solid") +
  geom_hline(yintercept=0.9346095, linetype="dashed") +
  geom_hline(yintercept=0.9626646, linetype="dashed") +
  geom_vline(xintercept=0.95, linetype="solid") +
  geom_vline(xintercept=0.9346095, linetype="dashed") +
  geom_vline(xintercept=0.9626646, linetype="dashed")+
  labs(x = "Coverage Prob. - Farrington",
       y = "Coverage Prob. - Pearson",
       color = expression(paste("Dispersion ", (phi))))+ 
  scale_colour_manual(values = disppalette)

# Afroz vs farrington 
plot_main_conf_afroz_farrington <- ggplot(l_dat_trans_main[order(l_dat_trans_main$props),], 
                                        aes(x = 1-covprob_wiNA_Farrington, 
                                            y = 1-covprob_wiNA_Pearson))+ 
  geom_point(aes(colour = factor(phi)))+ 
  geom_hline(yintercept=0.95)+
  theme_bw(base_size = 18)+
  geom_hline(yintercept=0.95, linetype="solid") +
  geom_hline(yintercept=0.9346095, linetype="dashed") +
  geom_hline(yintercept=0.9626646, linetype="dashed") +
  geom_vline(xintercept=0.95, linetype="solid") +
  geom_vline(xintercept=0.9346095, linetype="dashed") +
  geom_vline(xintercept=0.9626646, linetype="dashed")+
  labs(x = "Coverage Prob. - Farrington",
       y = "Coverage Prob. - Afroz",
       color = expression(paste("Dispersion ", (phi))))+ 
  scale_colour_manual(values = disppalette)

plot_main_conf_pearson_afroz_farrington <- 
  plot_main_conf_pearson_afroz + 
  plot_main_conf_pearson_farrington + 
  plot_main_conf_afroz_farrington + 
  plot_layout(guides = "collect") &
  theme(legend.position='bottom')

plot_main_conf_pearson_afroz_farrington


## zero handling - Coverage probability comparison ------------------------

plot_zero_covprob_pearson <- ggplot(l_dat_trans_zero[order(l_dat_trans_zero$props),], 
                                    aes(x = 1-covprob_wiNA_pearson_modeldf_original, 
                                        y = 1-covprob_wiNA_pearson_modeldf_onerow))+  
  geom_point(aes(colour = factor(phi)))+ 
  geom_hline(yintercept=0.95, linetype="solid") +
  geom_hline(yintercept=0.9346095, linetype="dashed") +
  geom_hline(yintercept=0.9626646, linetype="dashed") +
  geom_vline(xintercept=0.95, linetype="solid") +
  geom_vline(xintercept=0.9346095, linetype="dashed") +
  geom_vline(xintercept=0.9626646, linetype="dashed")+
  geom_abline(slope = 1, intercept = 0)+
  labs(x = "Coverage Prob. - Pearson - OG",
       y = "Coverage Prob. - Pearson - AO",
       color = expression(paste("Dispersion ", (phi))))+
  theme_bw(base_size = 18)+
  xlim(0.8972741, 1) + ylim(0.8972741, 1)+ 
  scale_colour_manual(values = disppalette)

plot_zero_covprob_afroz <- ggplot(l_dat_trans_zero[order(l_dat_trans_zero$props),], 
                                  aes(x = 1-covprob_wiNA_afroz_modeldf_original, 
                                      y = 1-covprob_wiNA_afroz_modeldf_onerow))+ 
  geom_point(aes(colour = factor(phi)))+ 
  geom_hline(yintercept=0.95, linetype="solid") +
  geom_hline(yintercept=0.9346095, linetype="dashed") +
  geom_hline(yintercept=0.9626646, linetype="dashed") +
  geom_vline(xintercept=0.95, linetype="solid") +
  geom_vline(xintercept=0.9346095, linetype="dashed") +
  geom_vline(xintercept=0.9626646, linetype="dashed")+
  geom_abline(slope = 1, intercept = 0)+
  labs(x = "Coverage Prob. - Afroz - OG",
       y = "Coverage Prob. - Afroz - AO",
       color = expression(paste("Dispersion ", (phi))))+
  theme_bw(base_size = 18)+
  xlim(0.8972741, 1) + ylim(0.8972741, 1)+ 
  scale_colour_manual(values = disppalette)

plot_zero_covprob_DM_MGLM <- ggplot(l_dat_trans_zero[order(l_dat_trans_zero$props),], 
                                    aes(x = 1-covprob_wiNA_DM_MGLM_modeldf_original, 
                                        y = 1-covprob_wiNA_DM_MGLM_modeldf_onerow))+ 
  geom_point(aes(colour = factor(phi)))+ 
  geom_hline(yintercept=0.95, linetype="solid") +
  geom_hline(yintercept=0.9346095, linetype="dashed") +
  geom_hline(yintercept=0.9626646, linetype="dashed") +
  geom_vline(xintercept=0.95, linetype="solid") +
  geom_vline(xintercept=0.9346095, linetype="dashed") +
  geom_vline(xintercept=0.9626646, linetype="dashed")+
  geom_abline(slope = 1, intercept = 0)+
  labs(x = "Coverage Prob. - Dir-mult - OG",
       y = "Coverage Prob. - Dir-mult - AO",
       color = expression(paste("Dispersion ", (phi))))+
  theme_bw(base_size = 18)+
  xlim(0.8972741, 1) + ylim(0.8972741, 1)+
  scale_colour_manual(values = disppalette)

plot_zero_covprob_afr_pear_dirm <-
  plot_zero_covprob_pearson +
  plot_zero_covprob_afroz +
  plot_zero_covprob_DM_MGLM +
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom')

plot_zero_covprob_afr_pear_dirm
