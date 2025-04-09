# Setup -------------------------------------------------------------------

library(VGAM)
library(MGLM)
library(multcomp)
library(tidyverse)
library(forcats)
library(Matrix)
library(xtable)

#set working directory
#setwd(".\\Code_and_Data")

## Helper functions -------------------------------------------------------

### backtransform logit in p(event)
expit <- function(x){
  exp(x)/(1 + exp(x))
}

# Dispersion parameter 'phi' to scaling parameter 'asum'
# asum: single positive number, the precision parameter of the dirichlet distribution
# asum close to 0 means extreme overdispersion, asum -> Inf means multinomial data

phi2asum <- function(phi, m){
  (phi-m)/(1-phi)
}

# phi2asum(phi=5, size=c(10,50,100,500,1000))

## Functions to prepare model parameters for use in glht() ----------------

model.frame.vglm <- function(model){
  model.frame <- model.framevlm(model)
}

model.vglm <- function(model, ...){
  df <- (dim(model@y)[2]-1)*model@misc$n - length(coef(model))
  multcomp:::modelparm.default(model, coef. = coef, vcov. = vcov, df = df)
}

mcp2matrix <- function(model, linfct){
  fc <- multcomp:::factor_contrasts(model)
  contrasts <- fc$contrasts
  factors <- fc$factors
  intercept <- fc$intercept
  mf <- fc$mf
  mm <- fc$mm
  alternative <- NULL
  if (!is.list(linfct) || is.null(names(linfct)))
    stop(sQuote("linfct"), "is not a named list")
  nhypo <- names(linfct)
  checknm <- nhypo %in% rownames(factors)
  if (!all(checknm))
    stop("Variable(s) ", sQuote(nhypo[!checknm]), " have been specified in ",
         sQuote("linfct"), " but cannot be found in ", sQuote("model"),
         "! ")
  if (any(checknm)) {
    checknm <- sapply(mf[nhypo[checknm]], is.factor)
    if (!all(checknm))
      stop("Variable(s) ", sQuote(paste(nhypo[!checknm],
                                        collapse = ", ")), " of class ",
           sQuote(paste(sapply(mf[nhypo[!checknm]], class), collapse = ", ")),
           " is/are not contained as a factor in ", sQuote("model"), ".")
  }
  m <- c()
  ctype <- c()
  for (nm in nhypo) {
    if (is.character(linfct[[nm]])) {
      Kchr <- function(kch) {
        types <- eval(formals(contrMat)$type)
        pm <- pmatch(kch, types)
        if (!is.na(pm)) {
          tmpK <- contrMat(table(mf[[nm]]), type = types[pm])
          ctype <<- c(ctype, types[pm])
        }
        else {
          tmp <- chrlinfct2matrix(kch, levels(mf[[nm]]))
          tmpK <- tmp$K
          m <<- c(m, tmp$m)
          if (is.null(alternative)) {
            alternative <<- tmp$alternative
          }
          else {
            if (tmp$alternative != alternative)
              stop("mix of alternatives currently not implemented")
          }
        }
        if (is.null(rownames(tmpK)))
          rownames(tmpK) <- paste(kch, 1:nrow(tmpK),
                                  sep = "_")
        if (length(nhypo) > 1)
          rownames(tmpK) <- paste(nm, rownames(tmpK),
                                  sep = ": ")
        list(K = tmpK)
      }
      tmp <- lapply(linfct[[nm]], Kchr)
      linfct[[nm]] <- do.call("rbind", lapply(tmp, function(x) x$K))
    }
  }
  for (nm in nhypo) {
    if (is.character(contrasts[[nm]])) {
      C <- do.call(contrasts[[nm]], list(n = nlevels(mf[[nm]])))
    }
    else {
      C <- contrasts[[nm]]
    }
    if (intercept || (!intercept && nm != colnames(factors)[1])) {
      Kstar <- linfct[[nm]] %*% C
    }
    else {
      Kstar <- linfct[[nm]]
    }
    pos <- factors[nm, ] == 1
    if (sum(pos) > 1)
      warning("covariate interactions found -- ",
              "default contrast might be inappropriate")
    attr(mm,"assign") <- unlist(attr(mm,"assign"))
    newlist <- attr(mm,"vassign")
    flist <- newlist[grep(nm, names(newlist))]
    level <- lapply(flist, function(x) attr(mm,"assign") %in% x)
    hypo <- vector(mode = "list", length = length(level))
    for(dp in seq_along(level)){
      hypo[[dp]] <- list(K = Kstar, where = level[[dp]])
    }
  }
  Ktotal <- matrix(0, nrow = sum(sapply(hypo, function(x) nrow(x$K))),
                   ncol = ncol(mm))
  colnames(Ktotal) <- colnames(mm)
  count <- 1
  for (h in hypo) {
    Ktotal[count:(count + nrow(h$K) - 1), h$where] <- h$K
    count <- count + nrow(h$K)
  }
  if (!is.matrix(Ktotal))
    Ktotal <- matrix(Ktotal, nrow = 1)
  nlist <- lapply(hypo, function(x) rownames(x$K))
  # hinzugefügt!!
  if(is.null(model@extra$use.refLevel)){
    refkatnr <- model@misc$M
  } else {
    refkatnr <- model@extra$use.refLevel
  }
  # bis hier
  ratios <- paste(rep(dimnames(model@y)[[2]][-refkatnr],
                      each=length(nlist[[1]])),
                  dimnames(model@y)[[2]][refkatnr], sep="/")
  rnames <- paste(ratios,
                  unlist(nlist), sep=": ")
  rownames(Ktotal) <- rnames
  if (is.null(ctype))
    ctype <- "User-defined"
  ctype <- paste(unique(ctype), collapse = ", ")
  attr(Ktotal, "type") <- ctype
  if (length(m) == 0)
    m <- 0
  list(K = Ktotal, m = m, alternative = alternative, type = ctype)
}


multin2mcp <-
  function(object,
           dispersion = c("none",
                          "pearson",
                          "afroz",
                          "farrington",
                          "deviance",
                          "dirmultinomial")) {
    disptype <- match.arg(dispersion)
    
    switch(
      disptype,
      none = {
        vcov <- vcov(object)
        df = object@df.residual
      },
      pearson = {
        vcov <- vcov.disp(object, "pearson") * vcov(object)
        df = object@df.residual
      },
      afroz = {
        vcov <- vcov.disp(object, "afroz") * vcov(object)
        df = object@df.residual
      },
      farrington = {
        vcov <- vcov.disp(object, "farrington") * vcov(object)
        df = object@df.residual
      },
      deviance = {
        vcov <- vcov.disp(object, "deviance") * vcov(object)
        df = object@df.residual
      },
      dirmultinomial = {
        vcov <- vcov(object)
        df = object@misc$n * (object@misc$M - 1) - object@misc$ncol.X.vlm
      }
    )
    parm.object <- parm(coef(object), vcov = vcov, df = df)
    return(parm.object)
  }

# Estimate Dispersion Parameter
vcov.disp <- function(model, dispersion=c("multinomial","pearson","afroz","farrington","dirmultinomial","deviance")){
  
  disptype <- match.arg(dispersion)
  
  if (is.null(model)) {
    return(NA)
  }
  
  if (disptype == "dirmultinomial"){
    icc <- expit(model@coefficients[length(model@misc$predictors.names)][[1]])
    phi <- 1 + (unique(model@extra$n2)-1)*icc
    return(phi)
  }
  
  y_vglm <- model@y
  y_vglm <- matrix(y_vglm, ncol = ncol(y_vglm))
  N_vglm <- length(y_vglm)
  w_vglm <- as.vector(model@prior.weights)
  
  pi_vglm <- model@fitted.values
  pi_vglm <- as.matrix(pi_vglm, ncol = ncol(pi_vglm))
  
  n_vglm <- model@misc$n
  
  s.bar_vglm <- sum((y_vglm - pi_vglm)/pi_vglm/(N_vglm - n_vglm))
  
  X2_vglm <- sum((y_vglm*w_vglm - pi_vglm*w_vglm) ^ 2 / (pi_vglm*w_vglm))
  
  phi.Pearson <- X2_vglm / model@df.residual
  
  switch (disptype,
          multinomial = {phi <- 1},
          pearson = {phi <- phi.Pearson},
          afroz = {phi <- (phi.Pearson/(1 + s.bar_vglm))},
          farrington = {phi <- (phi.Pearson - (N_vglm - n_vglm)*s.bar_vglm/(model@df.residual))},
          deviance = {phi <- model@criterion$deviance/model@df.residual}
  )
  
  return(phi)
}

# Generate Contrast Matrix for MGLM ---------------------------------------

f_generate_contr_mat <- function(num_cats, num_groups) {
  
  # Number of contrasts needed: (categories - 1) for intercepts + (categories - 1) for each group
  num_contrasts <- (num_cats - 1) + (num_groups - 1) * (num_cats - 1)
  
  # Initialize contrast matrix
  mat_T <- matrix(0, nrow = num_contrasts, ncol = num_groups * num_cats)
  
  contrast_index <- 1
  
  # Intercept contrasts: V2 - V1, V3 - V1, etc.
  for (j in 1:(num_cats-1)) {
    mat_T[contrast_index, 1] <- -1       # Beta for intercept V1
    mat_T[contrast_index, j*num_groups+1] <- 1        # Beta for intercept Vj
    contrast_index <- contrast_index + 1
  }
  
  # Group-specific contrasts: V2 - V1, V3 - V1, etc., for group2, group3, etc.
  for (g in 2:num_groups) {
    for (j in 1:(num_cats-1)) {
      start_index <- g  # Start of the current group's coefficients
      mat_T[contrast_index, start_index] <- -1  # Beta for group V1
      mat_T[contrast_index, start_index + num_groups*j] <- 1   # Beta for group Vj
      contrast_index <- contrast_index + 1
    }
  }
  
  return(mat_T)
}

# Differential Blood Count in Rats ----------------------------------------
# Page 54 Statistics in Toxicology - 2016 Hothorn

load(file = ".\\Code\\Example_Analysis\\Differential_Blood\\dif.rda")

dbb <- cbind(dif[, c(1:3)], dif[,5:10]*2)
colnames(dbb) <- c("sex", "animal", "Group", "Eos", "Baso",
                   "Stab", "Seg", "Mono", "Ly")
dbb$Group <- factor(dbb$Group,
                    levels = c("control", "low dose", "mid dose", "high dose"))
dbb$factorcomb <- dbb$sex:dbb$Group


## Contrast matrices ------------------------------------------------------

# define and name all pairwise logits between
# the 4 categories, i.e. matrix A
dbb_mat_A <- rbind(
  "Seg/Lym" = c(-1, 1, 0, 0),
  "Mon/Lym" = c(-1, 0, 1, 0),
  "Eos/Lym" = c(-1, 0, 0, 1),
  "Mon/Seg" = c( 0,-1, 1, 0),
  "Eos/Seg" = c( 0,-1, 0, 1),
  "Eos/Mon" = c( 0, 0,-1, 1))

# # define and name the comparisons between dose groups
# # within each sex, i.e. matrix B
dbb_mat_B <- rbind(
  "L/C, in females" = c(-1, 1, 0, 0, 0, 0, 0, 0),
  "M/C, in females" = c(-1, 0, 1, 0, 0, 0, 0, 0),
  "H/C, in females" = c(-1, 0, 0, 1, 0, 0, 0, 0),
  "L/C, in males" = c( 0, 0, 0, 0, -1, 1, 0, 0),
  "M/C, in males" = c( 0, 0, 0, 0, -1, 0, 1, 0),
  "H/C, in males" = c( 0, 0, 0, 0, -1, 0, 0, 1))

# extract contrast names
dbb_v_mat_A_names <- attributes(dbb_mat_A)$dimnames[[1]]
dbb_v_mat_B_names <- attributes(dbb_mat_B)$dimnames[[1]]

# Create final contrast names
dbb_v_cnames <- as.vector(outer(dbb_v_mat_A_names, dbb_v_mat_B_names, paste, sep = ": "))
dbb_v_cnames

# Calculate kronecker product of matrix A and B and assign respective contrast names
dbb_mat_A <- dbb_mat_A[,-1]
dbb_mat_kronAB <- kronecker(dbb_mat_B, dbb_mat_A)
dbb_mat_kronAB[ ,c(1:3)] <- 0
rownames(dbb_mat_kronAB) <- dbb_v_cnames

## Multinomial Model ------------------------------------------------------

# Fit a model with multinomial distribution assumption
fit_dbc_mult <- vglm(cbind(Seg, Mono, Eos, Ly) ~ factorcomb,
                     family = multinomial, data = dbb)

# calculate dispersion
vcov.disp(fit_dbc_mult, dispersion = "pearson")
vcov.disp(fit_dbc_mult, dispersion = "afroz")

summary(fit_dbc_mult)

# set seed for reproducibility
set.seed(123456)

# Here we use a custom multin2mcp function, which is a slightly modified function
# from the multcomp package, so we can directly extract the model parameters
# the covariance matrix (adjusted by phi) and the degrees of freedom from the model

# Multivariate t Distribution
s_fit_dbc_mult_df <- summary(glht(
  model = multin2mcp(fit_dbc_mult, dispersion = "none"),
  linfct = dbb_mat_kronAB,
  df = multin2mcp(fit_dbc_mult, dispersion = "none")$df
))
s_fit_dbc_mult_df

# Confidence Intervals Multivariate t Distribution
c_fit_dbc_mult_df <- confint(glht(
  model = multin2mcp(fit_dbc_mult, dispersion = "none"),
  linfct = dbb_mat_kronAB,
  df = multin2mcp(fit_dbc_mult, dispersion = "none")$df
))
c_fit_dbc_mult_df

## Afroz ------------------------------------------------------------------

# Multivariate t Distribution
s_fit_dbc_mult_afroz_df <- summary(glht(
  model = multin2mcp(fit_dbc_mult, dispersion = "afroz"),
  linfct = dbb_mat_kronAB,
  df = multin2mcp(fit_dbc_mult, dispersion = "none")$df
))
s_fit_dbc_mult_afroz_df

# Multivariate t Distribution
c_fit_dbc_mult_afroz_df <- confint(glht(
  model = multin2mcp(fit_dbc_mult, dispersion = "afroz"),
  linfct = dbb_mat_kronAB,
  df = multin2mcp(fit_dbc_mult, dispersion = "none")$df
))
c_fit_dbc_mult_afroz_df

## Pearson ----------------------------------------------------------------

# Multivariate t Distribution
# s_fit_dbc_mult_pearson_df <- summary(glht(
#   model = multin2mcp(fit_dbc_mult, dispersion = "pearson"),
#   linfct = dbb_mat_kronAB,
#   df = multin2mcp(fit_dbc_mult, dispersion = "none")$df
# ))
# s_fit_dbc_mult_pearson_df
# 
# # Confidence Intervals Multivariate t Distribution
# c_fit_dbc_mult_pearson_df <- confint(glht(
#   model = multin2mcp(fit_dbc_mult, dispersion = "pearson"),
#   linfct = dbb_mat_kronAB,
#   df = multin2mcp(fit_dbc_mult, dispersion = "none")$df
# ))
# c_fit_dbc_mult_pearson_df


## Farrington -------------------------------------------------------------

# Multivariate t Distribution
# s_fit_dbc_mult_farrington_df <- summary(glht(
#   model = multin2mcp(fit_dbc_mult, dispersion = "farrington"),
#   linfct = dbb_mat_kronAB,
#   df = multin2mcp(fit_dbc_mult, dispersion = "none")$df
# ))
# s_fit_dbc_mult_farrington_df
# 
# # Confidence Intevals Multivariate t Distribution
# c_fit_dbc_mult_farrington_df <- confint(glht(
#   model = multin2mcp(fit_dbc_mult, dispersion = "farrington"),
#   linfct = dbb_mat_kronAB,
#   df = multin2mcp(fit_dbc_mult, dispersion = "none")$df
# ))
# c_fit_dbc_mult_farrington_df

## Deviance ---------------------------------------------------------------

# Multivariate t Distribution
# s_fit_dbc_mult_deviance_df <- summary(glht(
#   model = multin2mcp(fit_dbc_mult, dispersion = "deviance"),
#   linfct = dbb_mat_kronAB,
#   df = multin2mcp(fit_dbc_mult, dispersion = "none")$df
# ))
# s_fit_dbc_mult_deviance_df
# 
# # Confidence Intervals Multivariate t Distribution
# c_fit_dbc_mult_deviance_df <- confint(glht(
#   model = multin2mcp(fit_dbc_mult, dispersion = "deviance"),
#   linfct = dbb_mat_kronAB,
#   df = multin2mcp(fit_dbc_mult, dispersion = "none")$df
# ))
# c_fit_dbc_mult_deviance_df

## Dirichlet-multinomial VGAM implementation ------------------------------

# Fit a model with Dirichlet-multinomial distributional assumption (VGAM implementation)

# fit_dbc_dirmult <- vglm(cbind(Seg, Mono, Eos, Ly) ~ factorcomb,
#                         family = dirmultinomial, data = dbb)
# 
# # zwei Werte da in einer Gruppe nur 199 samples
# # vcov.disp(fit_dbc_dirmult, dispersion = "dirmultinomial")
# # summary(fit_dbc_dirmult)
# 
# # add one column at the beginning of the Matrix Kstar, 
# # as we have one more Intercept (logitlink(phi))
# dbb_mat_kronAB_dir <- cbind(numeric(nrow(dbb_mat_kronAB)), dbb_mat_kronAB)
# 
# # Multivariate t Distribution
# s_fit_dbc_dirmult_df <- summary(glht(
#   model = multin2mcp(fit_dbc_dirmult, dispersion = "dirmultinomial"),
#   linfct = dbb_mat_kronAB_dir,
#   df = multin2mcp(fit_dbc_dirmult, dispersion = "dirmultinomial")$df
# ))
# s_fit_dbc_dirmult_df
# 
# # Confidence Intervals Multivariate t Distribution
# c_fit_dbc_dirmult_df <- confint(glht(
#   model = multin2mcp(fit_dbc_dirmult, dispersion = "dirmultinomial"),
#   linfct = dbb_mat_kronAB_dir,
#   df = multin2mcp(fit_dbc_dirmult, dispersion = "dirmultinomial")$df
# ))
# c_fit_dbc_dirmult_df

## Dirichlet-multinomial MGLM implementation -------------------------------

# Fit a model with Dirichlet-multinomial distributional assumption (MGLM implementation)

fit_dbc_MGLM_DM <- MGLMreg(
  formula = cbind(Ly, Seg, Mono, Eos) ~ factorcomb,
  dist = "DM",
  data = dbb
)

## Matrices ---------------------------------------------------------------

# For use with glht we need to transform the MGLM coefficients into the same form as 
# the VGAM coefficients
m_dbc_contr_mat <- f_generate_contr_mat(num_cats = ncol(fit_dbc_MGLM_DM@coefficients), 
                                    num_groups = nrow(fit_dbc_MGLM_DM@coefficients))
m_dbc_mglm_fit_coefs <- matrix(fit_dbc_MGLM_DM@coefficients)

m_dbc_mglm_lors <- m_dbc_contr_mat %*% m_dbc_mglm_fit_coefs

m_dbc_vcov_mglm_lors <- m_dbc_contr_mat %*% solve(-fit_dbc_MGLM_DM@Hessian) %*% t(m_dbc_contr_mat)

dbc_mglm_df <- (ncol(fit_dbc_MGLM_DM@data$Y-1))*nrow(fit_dbc_MGLM_DM@data$Y)-nrow(fit_dbc_MGLM_DM@data$Y)-length(m_dbc_mglm_fit_coefs)

dbc_modelparms_mglm <- list(coef = m_dbc_mglm_lors, vcov = m_dbc_vcov_mglm_lors, df = dbc_mglm_df)
attr(dbc_modelparms_mglm, "class") <- "parm"

#

s_fit_dbc_MGLM_DM_df <- summary(glht(
  model = dbc_modelparms_mglm,
  linfct = dbb_mat_kronAB,
  df = dbc_modelparms_mglm$df
))

c_fit_dbc_MGLM_DM_df <- confint(glht(
  model = dbc_modelparms_mglm,
  linfct = dbb_mat_kronAB,
  df = dbc_modelparms_mglm$df
))


## Table ------------------------------------------------------------------

# Creation of the table (Table: 2) as shown in the Example 2: Differential Blood Count 
# in the manuscript

df_dbc_mult_afroz <- data.frame(linear_hypotheses = attributes(c_fit_dbc_mult_df$linfct)$dimnames[[1]],
                                est_mult = c_fit_dbc_mult_df$confint[,1],
                                mul_SE = s_fit_dbc_mult_df$test$sigma,
                                afroz_SE = s_fit_dbc_mult_afroz_df$test$sigma,
                                mul_p = s_fit_dbc_mult_df$test$pvalues,
                                afr_p = s_fit_dbc_mult_afroz_df$test$pvalues)

# All comparisons are shown here, and only a selection of these comparisons 
# was used in the manuscript

print(xtable(df_dbc_mult_afroz, digits = c(0,3,3,3,3,3,3)), include.rownames=FALSE)


# data frame which is used for the generation of the plot with the confidence intervals
df_dbc_confints_df <- data.frame(linear_hypotheses = attributes(c_fit_dbc_mult_df$linfct)$dimnames[[1]],
                                 est_multinomial = c_fit_dbc_mult_df$confint[,1],
                                 lwr_multinomial = c_fit_dbc_mult_df$confint[,2],
                                 upr_multinomial = c_fit_dbc_mult_df$confint[,3],
                                 est_afroz = c_fit_dbc_mult_df$confint[,1],
                                 lwr_afroz = c_fit_dbc_mult_afroz_df$confint[,2],
                                 upr_afroz = c_fit_dbc_mult_afroz_df$confint[,3],
                                 est_mglm_dm = c_fit_dbc_MGLM_DM_df$confint[,1],
                                 lwr_mglm_dm = c_fit_dbc_MGLM_DM_df$confint[,2],
                                 upr_mglm_dm = c_fit_dbc_MGLM_DM_df$confint[,3])

row.names(df_dbc_confints_df) <- NULL

## Plots ------------------------------------------------------------------
names(dbb)[names(dbb) == "Ly"] <- "Lym"
dbb <- dbb %>%  arrange(Group, sex, desc(Lym)) %>% group_by(Group, sex) %>% mutate(id = row_number())
dbb_long <- dbb %>% 
  pivot_longer(
    cols = c(Eos, Baso, Stab, Seg, Mono, Lym), 
    names_to = "wbc"
  ) %>%
  uncount(value) %>%
  filter(wbc != "Stab")

dbb_long <- dbb_long %>%
  mutate(Group = case_when(
    Group == "control" ~ "Control",
    Group == "low dose" ~ "Low Dose",
    Group == "mid dose" ~ "Mid Dose",
    Group == "high dose" ~ "High Dose"
  )) %>%
  mutate(Group = factor(Group)) %>%
  mutate(Group = fct_relevel(Group, "Control", "Low Dose", "Mid Dose", "High Dose"))

dbb_long$wbc <- factor(dbb_long$wbc, levels = c("Eos","Seg","Mono","Lym"))
dbb_long$animal <- as.factor(dbb_long$animal)
dbb_long$id <- as.factor(dbb_long$id)

plot_dbc <- ggplot(dbb_long, aes(x = id, fill = wbc))+
  geom_bar(position = "stack", width = 0.8)+
  facet_grid(sex ~ Group, scale = "free_x")+ 
  scale_fill_manual(values=c('#7e84fa','#f68511','#23bbb5','#72e06a'))+
  theme_bw(base_size = 18)+
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())+
  labs(y = "Count", x = "ID", fill = "Category") +
  theme(legend.position="bottom")

plot_dbc

ggsave(
  ".\\Figures\\Budig_SimInfMult_Figure_5.eps",
  plot = plot_dbc ,
  width = 14,
  height = 6,
  device = "eps",
  dpi = 900
)

df_dbc_confints_df_long <- reshape(
  data = df_dbc_confints_df,
  idvar = "linear_hypotheses", # Die Variable, die die Beobachtungen identifiziert
  varying = 2:10, # Die Spalten, die von breit nach lang umgeformt werden sollen
  v.names = c("est", "lwr", "upr"), # Die Namen der neuen Spalten für die Werte
  times = c("multinomial", "afroz", "mglm_dm"), # Die Namen der neuen Spalte für die Methoden
  timevar = "methode", # Der Name der neuen Spalte für die Methoden
  direction = "long" # Die Richtung der Umformung
)

df_dbc_confints_df_long$methode <- factor(df_dbc_confints_df_long$methode, 
                                          levels = c("multinomial","afroz","mglm_dm"))

plot_dbc_cis <-
  df_dbc_confints_df_long %>% filter(methode == "multinomial" |
                                       methode == "mglm_dm" |
                                       methode == "afroz") %>% 
  mutate(methode = recode(methode,
                          "multinomial" = "Multinomial",
                          "afroz" = "Afroz",
                          "mglm_dm" = "Dirichlet-multinomial (MGLM)")) %>%
  ggplot(aes(x = est, y = linear_hypotheses)) +
  geom_vline(xintercept = 0,
             size = 0.75) +
  geom_errorbarh(
    aes(xmin = lwr, xmax = upr),
    size = 0.75,
    alpha = 1,
    height = 0.5
  ) +
  geom_point(size = 2.5) +
  labs(x = " log odds ratio",
       y = "Comparison") +
  facet_grid(~ methode) +
  theme_bw(base_size = 17)

plot_dbc_cis

ggsave(
  ".\\Figures\\Budig_SimInfMult_Figure_6.eps",
  plot = plot_dbc_cis,
  width = 12,
  height = 7.5,
  device = "eps",
  dpi = 800
)

# DYME Reproduction Example -----------------------------------------------

load(file = ".\\Code\\Example_Analysis\\DYME_Toxicity\\bivar.rda")

levels(bivar$DEFECT_TYPE) <- c("alive","malformed","dead")
bivar$DOSE <- as.factor(bivar$DOSE)
# Reformat data
bivar$alive <- ifelse(bivar$DEFECT_TYPE == "alive", 1, 0)
bivar$malformed <- ifelse(bivar$DEFECT_TYPE == "malformed", 1, 0)
bivar$dead <- ifelse(bivar$DEFECT_TYPE == "dead", 1, 0)
bivar$DAM_ID <- as.factor(bivar$DAM_ID)

bivar.re <- bivar %>% group_by(DAM_ID, DOSE) %>% summarize(alive = sum(alive),
                                                           malformed = sum(malformed),
                                                           dead = sum(dead))

bivar.re <- bivar.re[order(bivar.re$DOSE),]
bivar.re[39,4] <- 1

## Multinomial Model ------------------------------------------------------

# set seed for same p-values and confidence intervals as in the manuscript
set.seed(1234)

# Fit the model with the assumption of a multinomial distribution. 
# alive is set as the reference category
fit_mult <-
  vglm(cbind(malformed, dead, alive) ~ DOSE, multinomial, data = bivar.re)

summary(fit_mult)

# Display of the various dispersion parameters
vcov.disp(fit_mult, dispersion = "pearson")
vcov.disp(fit_mult, dispersion = "afroz")
vcov.disp(fit_mult, dispersion = "farrington")
vcov.disp(fit_mult, dispersion = "deviance")

## Afroz ------------------------------------------------------------------

# Use of the Afroz dispersion estimator as in Example 1: Developmental Toxicity 
# in the manuscript

# Multivariate t Distribution
s_fit_mult_afroz_df <- summary(glht(
  model = multin2mcp(fit_mult, dispersion = "afroz"),
  linfct = mcp2matrix(fit_mult, linfct = mcp(DOSE = "Dunnet"))$K,
  df = multin2mcp(fit_mult, dispersion = "none")$df
))
s_fit_mult_afroz_df

# Multivariate t Distribution
c_fit_mult_afroz_df <- confint(glht(
  model = multin2mcp(fit_mult, dispersion = "afroz"),
  linfct = mcp2matrix(fit_mult, linfct = mcp(DOSE = "Dunnet"))$K,
  df = multin2mcp(fit_mult, dispersion = "none")$df
))
c_fit_mult_afroz_df

## Plain Multinomial model ------------------------------------------------

# Multivariate t Distribution
# s_fit_mult_df <- summary(glht(
#   model = multin2mcp(fit_mult, dispersion = "none"),
#   linfct = mcp2matrix(fit_mult, linfct = mcp(DOSE = "Dunnet"))$K,
#   df = multin2mcp(fit_mult, dispersion = "none")$df
# ))
# s_fit_mult_df
# 
# # Confidence Interval Multivariate t Distribution
# c_fit_mult_df <- confint(glht(
#   model = multin2mcp(fit_mult, dispersion = "none"),
#   linfct = mcp2matrix(fit_mult, linfct = mcp(DOSE = "Dunnet"))$K,
#   df = multin2mcp(fit_mult, dispersion = "none")$df
# ))
# c_fit_mult_df

## Pearson ----------------------------------------------------------------

# Multivariate t Distribution
# s_fit_mult_pearson_df <- summary(glht(
#   model = multin2mcp(fit_mult, dispersion = "pearson"),
#   linfct = mcp2matrix(fit_mult, linfct = mcp(DOSE = "Dunnet"))$K,
#   df = multin2mcp(fit_mult, dispersion = "none")$df
# ))
# s_fit_mult_pearson_df
# 
# # Confidence Interval Multivariate t Distribution
# c_fit_mult_pearson_df <- confint(glht(
#   model = multin2mcp(fit_mult, dispersion = "pearson"),
#   linfct = mcp2matrix(fit_mult, linfct = mcp(DOSE = "Dunnet"))$K,
#   df = multin2mcp(fit_mult, dispersion = "none")$df
# ))
# c_fit_mult_pearson_df

## Farrington -------------------------------------------------------------

# Multivariate t Distribution
# s_fit_mult_farrington_df <- summary(glht(
#   model = multin2mcp(fit_mult, dispersion = "farrington"),
#   linfct = mcp2matrix(fit_mult, linfct = mcp(DOSE = "Dunnet"))$K,
#   df = multin2mcp(fit_mult, dispersion = "none")$df
# ))
# s_fit_mult_farrington_df
# 
# # Confidence Interval Multivariate t Distribution
# c_fit_mult_farrington_df <- confint(glht(
#   model = multin2mcp(fit_mult, dispersion = "farrington"),
#   linfct = mcp2matrix(fit_mult, linfct = mcp(DOSE = "Dunnet"))$K,
#   df = multin2mcp(fit_mult, dispersion = "none")$df
# ))
# c_fit_mult_farrington_df

## Deviance ---------------------------------------------------------------

# Multivariate t Distribution
# s_fit_mult_deviance_df <- summary(glht(
#   model = multin2mcp(fit_mult, dispersion = "deviance"),
#   linfct = mcp2matrix(fit_mult, linfct = mcp(DOSE = "Dunnet"))$K,
#   df = multin2mcp(fit_mult, dispersion = "none")$df
# ))
# s_fit_mult_deviance_df
# 
# # Confidence Interval Multivariate t Distribution
# c_fit_mult_deviance_df <- confint(glht(
#   model = multin2mcp(fit_mult, dispersion = "deviance"),
#   linfct = mcp2matrix(fit_mult, linfct = mcp(DOSE = "Dunnet"))$K,
#   df = multin2mcp(fit_mult, dispersion = "none")$df
# ))
# c_fit_mult_deviance_df

## Dirichlet-mutinomial VGAM implementation -------------------------------

# fit_dirmult <-
#   vglm(cbind(malformed, dead, alive) ~ DOSE, dirmultinomial, data = bivar.re)
# 
# summary(fit_dirmult)
# vcov.disp(fit_dirmult, dispersion = "dirmultinomial")
# 
# # Multivariate t Distribution
# s_fit_dirmult_df <- summary(glht(
#   model = multin2mcp(fit_dirmult, dispersion = "dirmultinomial"),
#   linfct = mcp2matrix(fit_dirmult, linfct = mcp(DOSE = "Dunnet"))$K,
#   df = multin2mcp(fit_dirmult, dispersion = "dirmultinomial")$df
# ))
# s_fit_dirmult_df
# 
# # Confidence Interval Multivariate t Distribution
# c_fit_dirmult_df <- confint(glht(
#   model = multin2mcp(fit_dirmult, dispersion = "dirmultinomial"),
#   linfct = mcp2matrix(fit_dirmult, linfct = mcp(DOSE = "Dunnet"))$K,
#   df = multin2mcp(fit_dirmult, dispersion = "dirmultinomial")$df
# ))
# c_fit_dirmult_df

## Dirichlet-multinomial MGLM implementation ------------------------------

# fit_bivar_MGLM_DM <- MGLMreg(
#   formula = cbind(alive, malformed, dead) ~ DOSE,
#   dist = "DM",
#   data = bivar.re
# )
# 
# # For use with glht
# m_bivar_contr_mat <- f_generate_contr_mat(
#   num_cats = ncol(fit_bivar_MGLM_DM@coefficients),
#   num_groups = nrow(fit_bivar_MGLM_DM@coefficients)
# )
# m_bivar_mglm_fit_coefs <- matrix(fit_bivar_MGLM_DM@coefficients)
# 
# m_bivar_mglm_lors <- m_bivar_contr_mat %*% m_bivar_mglm_fit_coefs
# 
# m_bivar_vcov_mglm_lors <- m_bivar_contr_mat %*% solve(-fit_bivar_MGLM_DM@Hessian) %*% t(m_bivar_contr_mat)
# 
# bivar_mglm_df <- (ncol(fit_bivar_MGLM_DM@data$Y - 1)) * nrow(fit_bivar_MGLM_DM@data$Y) -
#   nrow(fit_bivar_MGLM_DM@data$Y) - length(m_mglm_fit_coefs)
# 
# bivar_modelparms_mglm <- list(coef = m_bivar_mglm_lors, vcov = m_bivar_vcov_mglm_lors, df = bivar_mglm_df)
# attr(bivar_modelparms_mglm, "class") <- "parm"
# 
# 
# # Multivariate t Distribution
# s_fit_MGLM_DM_df <- summary(
#   glht(
#     model = bivar_modelparms_mglm,
#     linfct = mcp2matrix(fit_mult, linfct = mcp(DOSE = "Dunnet"))$K,
#     df = bivar_modelparms_mglm$df
#   )
# )
# s_fit_MGLM_DM_df
# 
# # Multivariate t Distribution
# c_fit_MGLM_DM__df <- confint(
#   glht(
#     model = bivar_modelparms_mglm,
#     linfct = mcp2matrix(fit_mult, linfct = mcp(DOSE = "Dunnet"))$K,
#     df = bivar_modelparms_mglm$df
#   )
# )
# c_fit_MGLM_DM__df

## Table ------------------------------------------------------------------

# creation of the table as it can also be found in the example in the manuscript

df_dyme_afroz <- data.frame(linear_hypotheses = attributes(s_fit_mult_afroz_df$linfct)$dimnames[[1]],
                            est_afroz = c_fit_mult_afroz_df$confint[,1],
                            afr_SE = s_fit_mult_afroz_df$test$sigma,
                            afr_p = s_fit_mult_afroz_df$test$pvalues,
                            exp_est_afroz = exp(c_fit_mult_afroz_df$confint[,1]),
                            exp_lwr_afroz = exp(c_fit_mult_afroz_df$confint[,2]),
                            exp_upr_afroz = exp(c_fit_mult_afroz_df$confint[,3]))

row.names(df_dyme_afroz) <- NULL

print(xtable(df_dyme_afroz), include.rownames=FALSE)

## Plots ------------------------------------------------------------------

names(bivar)[names(bivar) == "DOSE"] <- "Dose"

plot_DYME <- ggplot(bivar, aes(x= forcats::fct_infreq(DAM_ID), fill = DEFECT_TYPE))+
  geom_bar(position = "stack", width = 0.8)+
  facet_grid(~Dose,scale="free_x", labeller=label_both)+
  scale_fill_manual(values=c('#72e06a','#23bbb5','#f68511'))+
  theme_bw(base_size = 18)+
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())+
  labs(y = "Count", x = "ID", fill = "Category") + 
  theme(legend.position="bottom")


plot_DYME

ggsave(
  ".\\Figures\\Budig_SimInfMult_Figure_7.eps",
  plot = plot_DYME ,
  width = 14,
  height = 6,
  device = "eps",
  dpi = 900
)

# Flowcytometer -----------------------------------------------------------

df_flow <- readRDS(".\\Code\\Example_Analysis\\Flow_Cytometrie\\dflowcytrw.rds")

df_flow <- df_flow %>% filter(med != "m2s")
df_flow$cultmed <- as.factor(df_flow$cultmed)
df_flow$donor <- as.factor(df_flow$donor)

## Contrast matrices ------------------------------------------------------

# define and name all pairwise logits between
# the 4 categories, i.e. matrix A
m_flow_oddsap <- rbind(
  "Ab/ab" = c(-1, 1, 0, 0),
  "aB/ab" = c(-1, 0, 1, 0),
  "AB/ab" = c(-1, 0, 0, 1))
# 
# # define and name the comparisons between dose groups
# # within each sex, i.e. matrix B
m_flow_trtvscon <- rbind(
  "med2-med1" = c(-0.5,0.5, -0.5, 0.5),
  "mix-single" = c(0.5, 0.5, -0.5, -0.5),
  "interaction" = c(1,-1,-1,1))

v_flow_oddsap_names <- attributes(m_flow_oddsap)$dimnames[[1]]
v_flow_trtvscon_names <- attributes(m_flow_trtvscon)$dimnames[[1]]

v_flow_cnames <- as.vector(outer(v_flow_oddsap_names, v_flow_trtvscon_names, paste, sep = ": "))
v_flow_cnames

m_flow_oddsap <- m_flow_oddsap[,-1]
m_flow_Kstar <- kronecker(m_flow_trtvscon, m_flow_oddsap)
m_flow_Kstar[ ,c(1:3)] <- 0
rownames(m_flow_Kstar) <- v_flow_cnames


m_flow_block <- matrix(
  0,
  ncol = nlevels(df_flow$donor) * ncol(m_flow_oddsap) - ncol(m_flow_oddsap),
  nrow = nrow(m_flow_Kstar)
)

m_flow_Kstar <- cbind(m_flow_block, m_flow_Kstar)

## Multinomial Model ------------------------------------------------------

fit_flow_mult <-
  vglm(cbind(Ab, aB, AB, ab) ~ donor + cultmed, 
       family = multinomial, 
       data = df_flow)

summary(fit_flow_mult)

vcov.disp(fit_flow_mult, dispersion = "pearson")
vcov.disp(fit_flow_mult, dispersion = "afroz")
vcov.disp(fit_flow_mult, dispersion = "farrington")
vcov.disp(fit_flow_mult, dispersion = "deviance")

## Afroz ------------------------------------------------------------------

set.seed(24234)

# Multivariate t Distribution
s_flow_fit_mult_afroz_df <- summary(glht(
  model = multin2mcp(fit_flow_mult, dispersion = "afroz"),
  linfct = m_flow_Kstar,
  df = multin2mcp(fit_flow_mult, dispersion = "afroz")$df
))
s_flow_fit_mult_afroz_df

# Confidence Intervals Multivariate t Distribution
c_flow_fit_mult_afroz_df <- confint(glht(
  model = multin2mcp(fit_flow_mult, dispersion = "afroz"),
  linfct = m_flow_Kstar,
  df = multin2mcp(fit_flow_mult, dispersion = "afroz")$df
))
c_flow_fit_mult_afroz_df


## Plain multinomial model ------------------------------------------------

# Multivariate t Distribution
# s_flow_fit_mult_df <- summary(glht(
#   model = multin2mcp(fit_flow_mult, dispersion = "none"),
#   linfct = m_flow_Kstar,
#   df = multin2mcp(fit_flow_mult, dispersion = "none")$df
# ))
# s_flow_fit_mult_df
# 
# # Confidence Intervals Multivariate t Distribution
# c_flow_fit_mult_df <- confint(glht(
#   model = multin2mcp(fit_flow_mult, dispersion = "none"),
#   linfct = m_flow_Kstar,
#   df = multin2mcp(fit_flow_mult, dispersion = "none")$df
# ))
# c_flow_fit_mult_df

## Pearson ----------------------------------------------------------------

# Multivariate t Distribution
# s_flow_fit_mult_pearson_df <- summary(glht(
#   model = multin2mcp(fit_flow_mult, dispersion = "pearson"),
#   linfct = m_flow_Kstar,
#   df = multin2mcp(fit_flow_mult, dispersion = "pearson")$df
# ))
# s_flow_fit_mult_pearson_df
# 
# # Confidence Intervals Multivariate t Distribution
# c_flow_fit_mult_pearson_df <- confint(glht(
#   model = multin2mcp(fit_flow_mult, dispersion = "pearson"),
#   linfct = m_flow_Kstar,
#   df = multin2mcp(fit_flow_mult, dispersion = "pearson")$df
# ))
# c_flow_fit_mult_pearson_df

## Farrington -------------------------------------------------------------

# Multivariate t Distribution
# s_flow_fit_mult_farrington_df <- summary(glht(
#   model = multin2mcp(fit_flow_mult, dispersion = "farrington"),
#   linfct = m_flow_Kstar,
#   df = multin2mcp(fit_flow_mult, dispersion = "farrington")$df
# ))
# s_flow_fit_mult_farrington_df
# 
# # Confidence Intervals Multivariate t Distribution
# c_flow_fit_mult_farrington_df <- confint(glht(
#   model = multin2mcp(fit_flow_mult, dispersion = "farrington"),
#   linfct = m_flow_Kstar,
#   df = multin2mcp(fit_flow_mult, dispersion = "farrington")$df
# ))
# c_flow_fit_mult_farrington_df

## Deviance ---------------------------------------------------------------

# Multivariate t Distribution
# s_flow_fit_mult_deviance_df <- summary(glht(
#   model = multin2mcp(fit_flow_mult, dispersion = "deviance"),
#   linfct = m_flow_Kstar,
#   df = multin2mcp(fit_flow_mult, dispersion = "deviance")$df
# ))
# s_flow_fit_mult_deviance_df
# 
# # Confidence Intervals Multivariate t Distribution
# c_flow_fit_mult_deviance_df <- confint(glht(
#   model = multin2mcp(fit_flow_mult, dispersion = "deviance"),
#   linfct = m_flow_Kstar,
#   df = multin2mcp(fit_flow_mult, dispersion = "deviance")$df
# ))
# c_flow_fit_mult_deviance_df

## Dirichlet-multinomial VGAM implementation ------------------------------

# fit_flow_dirmult <-
#   vglm(cbind(Ab, aB, AB, ab) ~ donor + cultmed,
#        family = dirmultinomial,
#        data = df_flow)
# 
# vcov.disp(fit_flow_dirmult, dispersion = "dirmultinomial")
# 
# #summary(fit_flow_dirmult)
# 
# m_flow_Kstar_dir <- cbind(numeric(nrow(m_flow_Kstar)), m_flow_Kstar)
# 
# # Multivariate t Distribution
# s_flow_fit_dirmult_df <- summary(glht(
#   model = multin2mcp(fit_flow_dirmult, dispersion = "dirmultinomial"),
#   linfct = m_flow_Kstar_dir,
#   df = multin2mcp(fit_flow_dirmult, dispersion = "dirmultinomial")$df
# ))
# s_flow_fit_dirmult_df
# 
# # Confidence Intervals Multivariate t Distribution
# c_flow_fit_dirmult_df <- confint(glht(
#   model = multin2mcp(fit_flow_dirmult, dispersion = "dirmultinomial"),
#   linfct = m_flow_Kstar_dir,
#   df = multin2mcp(fit_flow_dirmult, dispersion = "dirmultinomial")$df
# ))
# c_flow_fit_dirmult_df

## Dirichlet-multinomial MGLM implementation ------------------------------

fit_flow_MGLM_DM <- MGLMreg(
  formula = cbind(ab, Ab, aB, AB) ~ donor + cultmed,
  dist = "DM",
  data = df_flow
)
str(df_flow)


# Model warnings are obtained. In general this model yields not reliable estimates/results
# Shouldnt be used for interpretation (cannot be used for p-value adjustment via glht)

# m_contr_A <- matrix(c(-1, 1, 0, 0,
#                       -1, 0, 1, 0,
#                       -1, 0, 0, 1
#                       ), ncol = 4, byrow = TRUE)
# 
# mat_I <- diag(7)
# m_contr_beta <- mat_I %x% m_contr_A
# 
# # to obtain coefficients like from the VGAM fit
# m_flow_mglm_fit_coefs <- m_contr_beta %*% matrix(t(fit_flow_MGLM_DM@coefficients))

# For use with glht
m_flow_contr_mat <- f_generate_contr_mat(
  num_cats = ncol(fit_flow_MGLM_DM@coefficients),
  num_groups = nrow(fit_flow_MGLM_DM@coefficients)
)
m_flow_mglm_fit_coefs <- matrix(fit_flow_MGLM_DM@coefficients)

m_flow_mglm_lors <- m_flow_contr_mat %*% matrix((fit_flow_MGLM_DM@coefficients))

m_flow_vcov_mglm_lors <- m_flow_contr_mat %*% solve(-fit_flow_MGLM_DM@Hessian) %*% t(m_flow_contr_mat)

flow_mglm_df <- (ncol(fit_flow_MGLM_DM@data$Y - 1)) * nrow(fit_flow_MGLM_DM@data$Y) -
  nrow(fit_flow_MGLM_DM@data$Y) - length(m_flow_mglm_fit_coefs)

flow_modelparms_mglm <- list(coef = m_flow_mglm_lors, vcov = m_flow_vcov_mglm_lors, df = flow_mglm_df)
attr(flow_modelparms_mglm, "class") <- "parm"

# negative entries in the diagonal of the variance covariance matrix. 
# Indicates unreliable modelfit. Cant be used with glht 

# # Multivariate t Distribution
# s_flow_fit_MGLM_DM_df <- summary(
#   glht(
#     model = flow_modelparms_mglm,
#     linfct = m_flow_Kstar,
#     df = flow_modelparms_mglm$df
#   )
# )
# s_fit_MGLM_DM_df
# 
# # Multivariate t Distribution
# c_flow_fit_MGLM_DM_df <- confint(
#   glht(
#     model = flow_modelparms_mglm,
#     linfct = m_flow_Kstar,
#     df = flow_modelparms_mglm$df
#   )
# )


## Table ------------------------------------------------------------------

df_flow_confints_afroz <- data.frame(linear_hypotheses = attributes(c_flow_fit_mult_afroz_df$linfct)$dimnames[[1]],
                                  est_afroz = c_flow_fit_mult_afroz_df$confint[,1],
                                  afroz_SE = s_flow_fit_mult_afroz_df$test$sigma,
                                  afroz_pvalue = s_flow_fit_mult_afroz_df$test$pvalues,
                                  exp_est_afroz = exp(c_flow_fit_mult_afroz_df$confint[,1]),
                                  exp_lwr_afroz = exp(c_flow_fit_mult_afroz_df$confint[,2]),
                                  exp_upr_afroz = exp(c_flow_fit_mult_afroz_df$confint[,3])
                                  )

row.names(df_flow_confints_afroz) <- NULL

print(xtable(df_flow_confints_afroz), include.rownames=FALSE)


## Plots ------------------------------------------------------------------

df_flow$cultmed <- droplevels(df_flow$cultmed)
df_flow$cultmed <- factor(df_flow$cultmed, levels = c("mix:m1","mix:m2","single:m1","single:m2"))
names(df_flow)[names(df_flow) == "donor"] <- "Donor"

df_flow_long <- df_flow %>% 
  pivot_longer(
    cols = c(ab, Ab, AB, aB), 
    names_to = "catc"
  ) %>%
  uncount(value)

plot_flow <- ggplot(df_flow_long, aes(x = cultmed, fill = catc))+
  geom_bar(position = "stack", width = 0.8)+
  facet_grid( ~ Donor, scale = "free_x", labeller = label_both)+ 
  scale_fill_manual(values=c('#7e84fa','#f68511','#23bbb5','#72e06a'))+
  theme_bw(base_size = 18)+
  labs(y = "Count", x = "Culture:Medium", fill = "Cell type") +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plot_flow

ggsave(
  ".\\Figures\\Budig_SimInfMult_Figure_8.eps",
  plot = plot_flow,
  width = 14,
  height = 6.5,
  device = "eps",
  dpi = 900
)
