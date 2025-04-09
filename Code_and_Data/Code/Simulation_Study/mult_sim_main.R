# Setup -------------------------------------------------------------------

# use parent directory ".\\Code_and_Data"

#setwd(".\\Results\\main")

# Packages ----------------------------------------------------------------

library(multcomp)
library(MCMCpack)
library(VGAM)
library(tidyverse)
library(MGLM)

# Explanation of variables ------------------------------------------------

# props: set probabilities for categories and groups to create the multinomial data
# phi: set dispersion parameter
# m: cluster size
# b: sample size of clusters
# comp: Tukey or Dunnett comparisons with glht
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

# Functions ---------------------------------------------------------------

## Helper functions -------------------------------------------------------

### backtransform logit in p(event)
expit <- function(x){
  exp(x)/(1 + exp(x))
}

## Dispersion parameter 'phi' to scaling parameter 'asum' -----------------

# asum: single positive number, the precision parameter of the dirichlet distribution
# asum close to 0 means extreme overdispersion, asum -> Inf means multinomial data

phi2asum <- function(phi, m){
  (phi-m)/(1-phi)
}

# phi2asum(phi=5, size=c(10,50,100,500,1000))

# add one to all rows of groups that have only zeroes
f_add1allrows <- function(df_dat){
  for (col in colnames(df_dat)[- ncol(df_dat)]) {
    if(any(ave(df_dat[,col], df_dat$group) == 0)) {
      df_dat[,col][df_dat[,col]==0] <- 1
    }
  }
  return(df_dat)
}

f_add1onerow <- function(df_dat){
  df_sumzero <- df_dat %>% group_by(group)  %>%
    summarise(across(everything(), mean)) %>% select(-group)
  indices <- which(df_sumzero == 0, arr.ind = TRUE)
  
  samplesize <- nrow(df_dat)/nlevels(df_dat$group)
  if(any(indices)){
    for (i in 1:nrow(indices)) {
      df_dat[indices[i,][1]*samplesize,indices[i,2][1]] <- 1
    }
  }
  return(df_dat)
}

# count the number of zeroes per group and category
f_countzeroes <- function(df_dat){
  df_dat %>% group_by(group)  %>%
    summarise(across(everything(), mean))  %>% 
    summarise(countzero=sum(.==0))
}

# count minimum of colSums
f_colmin <- function(df_dat){
  df_dat %>% select(-group) %>%
    summarise(across(everything(), ~sum(.)))  %>% 
    summarise(Min = min(across(everything())))
}


# Dirichlet multinomial data sampling Parameters --------------------------

# m_prop <- matrix(c(0.33,0.33,0.33, #1 H0
#                    0.33,0.33,0.33,
#                    0.33,0.33,0.33,
#                    0.33,0.33,0.33),
#                  ncol = 3,
#                  byrow = TRUE)
# 
# v_prop <- m_prop[1,]
# b <- 5
# m <- 10
# phi <- 1.5

## Dirichlet multinomial data sampling ------------------------------------

# b: a single integer; number of multinomial-dirichlet vectors to draw

# m: single integer or vector of integers, 
#  the multinomial sample size for the individual multinomial vectors to generate
#  if vector: recycled to have length n

# prop: vector of proportions, length (prop) must be = no. of categories
#  prop should sum to 1, prop is rescaled such that it really sums to 1 

# catnam: name of categories, should equal length prop vector

rdirmultinom <- function(v_prop, b, m, phi, catnam = NULL) {
  
  v_size <- rep(m, length.out = b)
  v_props <- v_prop / sum(v_prop)
  asum <- phi2asum(phi = phi, m = m)
  v_alpha <- asum * v_props
  m_propdir <- rdirichlet(n = b, alpha = v_alpha)
  
  m_dat <- t(apply(
    X = cbind(v_size, m_propdir),
    MARGIN = 1,
    FUN = function(x) {
      rmultinom(n = 1,
                size = x[1],
                prob = x[-1])
    }
  ))
  
  if (!is.null(catnam)) {
    colnames(m_dat) <- catnam
  }
  
  return(m_dat)
}

# rdirmultinom(n=50, size=10, phi=2, v_prop=c(0.8,0.15,0.05), 
#              catnam=c("alive", "dead", "malf"))

## Updated multcomp functions ---------------------------------------------

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
  # ab hier von Vogel hinzugefügt/geändert
  nlist <- lapply(hypo, function(x) rownames(x$K))
  ## hier von Budig hinzugefügt
  if(is.null(model@extra$use.refLevel)){
    refkatnr <- model@misc$M
  } else {
    refkatnr <- model@extra$use.refLevel
  }
  # bist hier von Budig
  ratios <- paste(rep(dimnames(model@y)[[2]][-refkatnr],
                      each=length(nlist[[1]])),
                  dimnames(model@y)[[2]][refkatnr], sep="/")
  rnames <- paste(ratios,
                  unlist(nlist), sep=": ")
  rownames(Ktotal) <- rnames
  # bis her von Vogel hinzugefügt/geändert
  if (is.null(ctype))
    ctype <- "User-defined"
  ctype <- paste(unique(ctype), collapse = ", ")
  attr(Ktotal, "type") <- ctype
  if (length(m) == 0)
    m <- 0
  list(K = Ktotal, m = m, alternative = alternative, type = ctype)
}


## Functions to prepare model parameters for use in glht() ----------------

multin2mcp <- function(object, dispersion=c("none","pearson","afroz","farrington","deviance")){
  
  disptype <- match.arg(dispersion)
  
  switch(disptype,
         none = {
           vcov <- vcov(object)
           df = object@df.residual
         },
         pearson = {
           vcov <- vcov.disp(object, "pearson")*vcov(object)
           df = object@df.residual
         },
         afroz = {
           vcov <- vcov.disp(object, "afroz")*vcov(object)
           df = object@df.residual
         },
         farrington = {
           vcov <-vcov.disp(object, "farrington")*vcov(object)
           df = object@df.residual
         },
         deviance = {
           vcov <-vcov.disp(object, "deviance")*vcov(object)
           df = object@df.residual
         }
  )
  parm.object <- parm(coef(object), vcov = vcov, df = df)
  return(parm.object)
}

vcov.disp <- function(model, dispersion=c("multinomial","pearson","afroz","farrington","DM_VGAM","deviance")){
  
  disptype <- match.arg(dispersion)
  
  if (is.null(model)) {
    return(NA)
  }
  
  if (disptype == "DM_VGAM"){
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


## Fit multinomial model function -----------------------------------------

f_fit_mult <-
  function(df_i,
           modeltype = c("multinomial",
                         "pearson",
                         "afroz",
                         "farrington",
                         "deviance",
                         "DM_VGAM",
                         "DM_MGLM")) {
    
    modeltype <- match.arg(modeltype)
    
    t_formula <- paste0("cbind(",
                        paste0(colnames(df_i[2:(length(df_i) - 1)]),
                               collapse = ","),",V1) ~ group")
    #assign("last.warning", NULL, envir = globalenv())
    warning_text <- NULL
    
    switch(
      modeltype,
      multinomial = ,
      pearson = ,
      afroz = ,
      farrington = ,
      deviance = {
        modeloutput <-tryCatch( withCallingHandlers(
          {
            vglm(
              formula = t_formula,
              family = multinomial,
              model = TRUE,
              data = df_i
            )},
          warning = function(w){
            warning_text <<- c(warning_text, conditionMessage(w))
            invokeRestart("muffleWarning")
          }),
          error = function(e) {
            warning_text <<-  trimws(e)
            modeloutput <- NULL
          }
        )
        
        
      },
      DM_VGAM = {
        modeloutput <-tryCatch(withCallingHandlers({
          vglm(
            formula = t_formula,
            family = dirmultinomial,
            model = TRUE,
            data = df_i,
            maxit = 100
          )},
          warning = function(w){
            warning_text <<- c(warning_text, conditionMessage(w))
            invokeRestart("muffleWarning")
          }),
          error = function(e) {
            warning_text <<-   trimws(e)
            modeloutput <- NULL
          }
        )
        
        
      },
      DM_MGLM = {
        
        t_formula <- paste0("cbind(",
                            paste0(colnames(df_i[1:(length(df_i) - 1)]),
                                   collapse = ","),") ~ group")
        
        modeloutput <-tryCatch(withCallingHandlers({
          MGLMreg(
            formula = t_formula,
            dist = "DM",
            data = df_i
          )},
          warning = function(w){
            warning_text <<- c(warning_text, conditionMessage(w))
            invokeRestart("muffleWarning")
          }),
          error = function(e) {
            warning_text <<-   trimws(e)
            modeloutput <- NULL
          }
        )
      }
    )
    
    return(list(modeloutput,unique(warning_text)))
  }



## Apply single step glht function to model -------------------------------

f_sumglht <- function(model,
                      mat_mcp,
                      modeltype = c("multinomial",
                                    "pearson",
                                    "afroz",
                                    "farrington",
                                    "deviance",
                                    "DM_VGAM",
                                    "DM_MGLM"),
                      dfu = c("zero", "modeldf")) {
  
  modeltype <- match.arg(modeltype)
  dfu <- match.arg(dfu)
  
  switch(
    modeltype,
    multinomial = {
      modelparms <- multin2mcp(model, dispersion = "none")
    },
    pearson = {
      modelparms <- multin2mcp(model, dispersion = "pearson")
    },
    afroz = {
      modelparms <- multin2mcp(model, dispersion = "afroz")
    },
    farrington = {
      modelparms <- multin2mcp(model, dispersion = "farrington")
    },
    deviance = {
      modelparms <- multin2mcp(model, dispersion = "deviance")
    },
    DM_VGAM = {
      modelparms <- multin2mcp(model, dispersion = "none")
    },
    DM_MGLM = {
      m_contr_mat <- f_generate_contr_mat(num_cats = ncol(model@coefficients), 
                                          num_groups = nrow(model@coefficients))
      m_mglm_fit_coefs <- matrix(model@coefficients)
      
      m_lors <- m_contr_mat %*% m_mglm_fit_coefs
      
      m_vcov_mglm_lors <- m_contr_mat %*% solve(-model@Hessian, tol = 1e-20) %*% t(m_contr_mat)
      
      mglm_df <- (ncol(model@data$Y-1))*nrow(model@data$Y)-nrow(model@data$Y)-length(m_mglm_fit_coefs)
      
      modelparms <- list(coef = m_lors, vcov = m_vcov_mglm_lors, df = mglm_df)
      attr(modelparms, "class") <- "parm"
    }
  )
  
  warning_text <- NULL
  
  
  switch(dfu,
         zero = {
           glhtout <- tryCatch(withCallingHandlers({summary(glht(
             model = modelparms,
             linfct = mat_mcp
           ))},
           warning = function(w){
             warning_text <<- c(warning_text, conditionMessage(w))
             invokeRestart("muffleWarning")
           }),
           error = function(e) {
             warning_text <<-   trimws(e)
             glhtout <- NULL
           }
           )},
         modeldf = {
           glhtout <- tryCatch(withCallingHandlers({summary(glht(
             model = modelparms,
             linfct = mat_mcp,
             df = modelparms$df
           ))},
           warning = function(w){
             warning_text <<- c(warning_text, conditionMessage(w))
             invokeRestart("muffleWarning")
           }),
           error = function(e) {
             warning_text <<-   trimws(e)
             glhtout <- NULL
           })})
  
  return(list(glhtout,unique(warning_text)))
  
}


## Apply glht function to model to obtain simultaneous confints -----------

f_confglht <- function(model,
                       mat_mcp,
                       modeltype = c("multinomial",
                                     "pearson",
                                     "afroz",
                                     "farrington",
                                     "deviance",
                                     "DM_VGAM",
                                     "DM_MGLM"),
                       dfu = c("zero", "modeldf")) {
  
  modeltype <- match.arg(modeltype)
  dfu <- match.arg(dfu)
  
  switch(
    modeltype,
    multinomial = {
      modelparms <- multin2mcp(model, dispersion = "none")
    },
    pearson = {
      modelparms <- multin2mcp(model, dispersion = "pearson")
    },
    afroz = {
      modelparms <- multin2mcp(model, dispersion = "afroz")
    },
    farrington = {
      modelparms <- multin2mcp(model, dispersion = "farrington")
    },
    deviance = {
      modelparms <- multin2mcp(model, dispersion = "deviance")
    },
    DM_VGAM = {
      modelparms <- multin2mcp(model, dispersion = "none")
    },
    DM_MGLM = {
      m_contr_mat <- f_generate_contr_mat(num_cats = ncol(model@coefficients), 
                                          num_groups = nrow(model@coefficients))
      m_mglm_fit_coefs <- matrix(model@coefficients)
      
      m_lors <- m_contr_mat %*% m_mglm_fit_coefs
      
      m_vcov_mglm_lors <- m_contr_mat %*% solve(-model@Hessian, tol = 1e-20) %*% t(m_contr_mat)
      
      mglm_df <- (ncol(model@data$Y-1))*nrow(model@data$Y)-nrow(model@data$Y)-length(m_mglm_fit_coefs)
      
      modelparms <- list(coef = m_lors, vcov = m_vcov_mglm_lors, df = mglm_df)
      attr(modelparms, "class") <- "parm"
    }
  )
  
  switch(dfu,
         zero = {
           confintout <- tryCatch(confint(glht(
             model = modelparms,
             linfct = mat_mcp
           )),
           error = function(e) {
             return(NULL)
           },
           warning = function(w) {
             return(NULL)
           }
           )
         },
         modeldf = {
           confintout <- tryCatch(confint(glht(
             model = modelparms,
             linfct = mat_mcp,
             df = modelparms$df
           )),
           error = function(e) {
             return(NULL)
           },
           warning = function(w) {
             return(NULL)
           }
           )
         })
  
  return(confintout)
  
}


## Generate vector of LoRs for given probability matrix --------------------

f_LORov <- function(m_prop, comp = c("Tukey", "Dunnett")) {
  
  comptype <- match.arg(comp)
  
  C <- ncol(m_prop)
  G <- nrow(m_prop)
  
  m_cmCT <- contrMat(rep(1, times = C), type = "Dunnett")
  m_cmG <- contrMat(rep(1, times = G), type = comptype)
  
  m_delta <- matrix(apply(m_prop, 1, function(x) (m_cmCT %*% log(x))))
  
  # identitiy matrix of size I
  m_id <- diag(nrow(m_cmCT))
  
  theta <- (m_cmG %x% m_id) %*% m_delta
  
  if(comptype == "Dunnett"){
    v_CL <- rep((LETTERS[1:(C - 1)]), times = G - 1)
    v_GL <- rep(letters[1:(G - 1)], each = C - 1)
  }
  if(comptype == "Tukey"){
    v_CL <- rep((LETTERS[1:(C - 1)]), times = choose(G,2))
    v_GL <- rep(letters[1:(choose(G,2))], each = C-1)
  }
  
  df_LOR <- data.frame(LOR = theta, Comps = paste0(v_CL, v_GL))
  return(df_LOR[order(df_LOR$Comps), ]$LOR)
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

## Simulation Function ----------------------------------------------------

f_mult_sim <- function(m_prop, df_simdat){
 
  df_simdatu <- unique(df_simdat[,c("phi","b","m","comp","dfu")])
  
  # generate dirmult data for each row of df_simdat 
  l_dat <- lapply(apply(df_simdatu, 1, function(x)
    apply(m_prop, 1, function(y)
      rdirmultinom(
        phi = as.numeric(x[1]),
        b = as.numeric(x[2]),
        m = as.numeric(x[3]),
        v_prop = y
      ), simplify = FALSE)), function(z)
        cbind(as.data.frame(do.call(rbind, z)),
              "group" = as.factor(rep(
                seq(1:nrow(m_prop)),
                each = (nrow(as.data.frame(do.call(
                  rbind, z
                ))) / nrow(m_prop))
              ))))
  
  if (any(df_simdat$modeltype %in% c("multinomial" ,"pearson", "afroz","farrington", "deviance"))) {
    l_mods_mult <- lapply(seq_along(l_dat), function(i)
      f_fit_mult(
        df_i = l_dat[[i]],
        modeltype = "multinomial"
      ))
  }
  
  if (any(df_simdat$modeltype %in% c("DM_VGAM"))) {
    l_mods_dirmult <- lapply(seq_along(l_dat), function(i)
      f_fit_mult(
        df_i = l_dat[[i]],
        modeltype = "DM_VGAM"
      ))
  }
  
  if (any(df_simdat$modeltype %in% c("DM_MGLM"))) {
    l_mods_DM_MGLM <- lapply(seq_along(l_dat), function(i)
      f_fit_mult(
        df_i = l_dat[[i]],
        modeltype = "DM_MGLM"
      ))
  }
  
  l_mods <- list()
  for (i in levels(df_simdat$modeltype)) {
    if (i %in% c("multinomial", "pearson", "afroz", "farrington", "deviance")) {
      l_mods <- c(l_mods, l_mods_mult)
    } else if (i == "DM_VGAM") {
      l_mods <- c(l_mods, l_mods_dirmult)
    } else if (i == "DM_MGLM") {
      l_mods <- c(l_mods, l_mods_DM_MGLM)
    }
  }
  
  ### extract errors and warnings
  l_errwarn <- lapply(l_mods, `[[`,2)
  
  ### check for errors and warnings---
  # check if there is any error in model fitting process
  # potential error:
  
  # "NAs found in the working weights variable 'wz'"
  checkerrnas <- unlist(lapply(l_errwarn, function(i) {
    if(any(is.null(i))){0
    } else if (any(grepl("NAs", i, fixed = TRUE))) {1} else (0)
  }))
  
  # "Some elements in the working weights variable 'wz' are not finite"
  checkerrfin <- unlist(lapply(l_errwarn, function(i) {
    if(any(is.null(i))){0
    } else if (any(grepl("finite", i, fixed = TRUE))) {1} else (0)
  }))
  
  # check if there is a working wz warning
  # "diagonal elements of the working weights variable 'wz' have been replaced by 1.819e-12"
  checkwarwz <- unlist(lapply(l_errwarn, function(i) {
    if(any(is.null(i))){0
    } else if (
      any(grepl("diagonal", i, fixed = TRUE))) {1} else (0)
  }))
  
  # check if there is a convergence warning
  # "convergence not obtained in 30 IRLS iterations"
  checkwarcon <- unlist(lapply(l_errwarn, function(i) {
    if(any(is.null(i))){0
    } else if (
      any(grepl("convergence", i, fixed = TRUE))) {1} else {0}
  }))
  
  # check if there are any deleted columns due to zero counts
  # "Deleted 2 columns of the response matrix due to zero counts"
  checkwarzc <- unlist(lapply(l_errwarn, function(i) {
    if(any(is.null(i))){0
    } else if (
      any(grepl("Deleted", i, fixed = TRUE))) {1} else (0)
  }))
  
  # check if there is any other warning
  checkwarany <- unlist(lapply(l_errwarn, function(i) {
    if(any(is.null(i))){0
    } else if (
      any(grepl("diagonal", i, fixed = TRUE)) | 
      any(grepl("convergence", i, fixed = TRUE)) |
      any(grepl("Error", i, fixed = TRUE)) |
      any(grepl("Deleted", i, fixed = TRUE))) {0} else (1)
  }))
  #---
  
  # extract modelfit list
  l_mods <- lapply(l_mods, `[[`,1)
  
  # extract estimated dispersion parameters
  # can put out NA for missing model
  # not implemented for DM_MGLM
  # v_disphat <- unlist(lapply(seq_along(l_mods), function(j) 
  #   vcov.disp(
  #     model = l_mods[[j]],
  #     dispersion = as.character(df_simdat$modeltype[j])
  #   )))
  
  v_index_mult <- which(df_simdat$modeltype %in% c("multinomial", "pearson", "afroz", "farrington", "deviance"))
  
  l_mat_mcp <-  
    lapply(seq_along(l_mods), function(a){
      if (is.null(l_mods[[a]])) {
        NULL 
      } else if (df_simdat$modeltype[a] %in% c("multinomial", "pearson", "afroz", "farrington", "deviance", "DM_VGAM")) {
        
        mcp2matrix(l_mods[[a]], linfct = mcp(group = as.character(df_simdat$comp[a])))$K
        
      } else if (df_simdat$modeltype[a] %in% c("DM_MGLM")) {
        
        mcp2matrix(l_mods[[v_index_mult[1]]], linfct = mcp(group = as.character(df_simdat$comp[a])))$K
        
      }
    }
    )
  
  # apply glht function with single step adjustment to each model
  l_ressumglht <-
    lapply(seq_along(l_mods), function(i)
      if (is.null(l_mods[[i]])) {
        NULL
      } else {
        f_sumglht(
          model = l_mods[[i]],
          mat_mcp = l_mat_mcp[[i]],
          modeltype = as.character(df_simdat$modeltype[i]),
          dfu = as.character(df_simdat$dfu[i])
        )
      })
  
  l_errwarn_glht <- lapply(l_ressumglht, `[[`,2)
  l_ressumglht <- lapply(l_ressumglht, `[[`,1)
  
  checkabseps <- unlist(lapply(l_errwarn_glht, function(i) {
    if(any(is.null(i))){0
    } else if (any(grepl("Completion", i, fixed = TRUE))) {1} else (0)
  }))
  
  checkglhtany <- unlist(lapply(l_errwarn_glht, function(i) {
    if(any(is.null(i))){0
    } else if (
      any(grepl("Completion", i, fixed = TRUE))) {0} else (1)
  }))
  
  
  # apply glht function to model to obtain simultaneous confints
  l_resconfglht <-
    lapply(seq_along(l_mods), function(i)
      if (is.null(l_mods[[i]])) {
        NULL
      } else {
        f_confglht(
          model = l_mods[[i]],
          mat_mcp = l_mat_mcp[[i]],
          modeltype = as.character(df_simdat$modeltype[i]),
          dfu = as.character(df_simdat$dfu[i])
        )
      })
  
  l_LORov <- lapply(seq_along(l_mods), function(i)
    f_LORov(m_prop, comp = as.character(df_simdat$comp[i])))
  
  # estimated log odds as string
  estprobs <- sapply(seq_along(l_ressumglht), function(i)
    if (is.null(l_ressumglht[[i]])) {
      NA
    }
    else if (length(l_LORov[[i]]) != length(l_ressumglht[[i]]$test$coefficients)) {
      NA
    }
    else {paste(l_ressumglht[[i]]$test$coefficients, collapse = ",")}
  )
  
  # true log odds as string
  trueprobs <- sapply(seq_along(l_ressumglht), function(i)
  {paste(l_LORov[[i]], collapse = ",")}
  )
  
  # count if coefficients are missing
  # count 1 if less coefficients were estimated than should be estimated (dropped columns)
  countmfz <-
    sapply(seq_along(l_ressumglht), function(i)
      if (is.null(l_ressumglht[[i]])) {
        NA
      }
      else if (length(l_LORov[[i]]) != length(l_ressumglht[[i]]$test$coefficients)) {
        1
      }
      else {
        0
      })
  
  # count FWER
  # count 1 if pvalue < 0.05 for Nullhypothesis
  countFWER <-
    sapply(seq_along(l_ressumglht), function(i)
      if (is.null(l_ressumglht[[i]])) {
        NA
      } else if (length(l_LORov[[i]]) != length(l_ressumglht[[i]]$test$coefficients)) {
        NA
      } else {
        as.integer(any(l_ressumglht[[i]]$test$pvalues[!l_LORov[[i]]] < 0.05))
      })
  
  # count Power
  # count 1 if pvalue < 0.05 for Alternative hypothesis
  countPower <-
    sapply(seq_along(l_ressumglht), function(i)
      if (is.null(l_ressumglht[[i]])) {
        NA
      } else if (length(l_LORov[[i]]) != length(l_ressumglht[[i]]$test$coefficients)) {
        NA
      } else {
        as.integer(any(l_ressumglht[[i]]$test$pvalues[!!l_LORov[[i]]] < 0.05))
      })
  
  # count globalPower
  # count 1 if pvalue < 0.05 for any hypothesis
  countgPower <-
    sapply(seq_along(l_ressumglht), function(i)
      if (is.null(l_ressumglht[[i]])) {
        NA
      } else if (length(l_LORov[[i]]) != length(l_ressumglht[[i]]$test$coefficients)) {
        NA
      } else {
        as.integer(any(l_ressumglht[[i]]$test$pvalues < 0.05))
      })
  
  # count pPower
  # count percentage how many hypothesis are correctly rejected
  countpPower <-
    sapply(seq_along(l_ressumglht), function(i)
      if (is.null(l_ressumglht[[i]])) {
        NA
      } else if (length(l_LORov[[i]]) != length(l_ressumglht[[i]]$test$coefficients)) {
        NA
      } else {
        sum(l_ressumglht[[i]]$test$pvalues[!!l_LORov[[i]]] < 0.05)/length(l_ressumglht[[i]]$test$pvalues[!!l_LORov[[i]]])
      })
  
  # count coverage confint
  # count 1 if logitprob lies outside interval
  countConf <-
    sapply(seq_along(l_resconfglht), function(i)
      if (is.null(l_ressumglht[[i]])) {
        NA
      } else if (length(l_LORov[[i]]) != length(l_ressumglht[[i]]$test$coefficients)) {
        NA
      } else {
        as.integer(any(
          !(
            l_resconfglht[[i]]$confint[, 2] < l_LORov[[i]] &
              l_LORov[[i]] < l_resconfglht[[i]]$confint[, 3]
          )
        ))
      })
  
  return(cbind(
    df_simdat,
    data.frame(
      sim = 1,
      C = ncol(m_prop),
      G = nrow(m_prop),
      minpm = min(m_prop) * df_simdat$m,
      minp = min(m_prop),
      maxp = max(m_prop),
      HAinc = unlist(lapply(l_LORov, function(x)
        as.integer(any(
          !!x
        )))),
      cH0 = unlist(lapply(l_LORov, function(x)
        sum(as.integer((!x)
        )))),
      cHA = unlist(lapply(l_LORov, function(x)
        sum(as.integer((!!x)
        )))),
      cmfz = countmfz,
      cFWER = countFWER,
      cPower = countPower,
      cgPower = countgPower,
      cpPower = countpPower,
      cConf = countConf,
      #disphat = v_disphat #disabled for simulations with mglm
      cerrnas = checkerrnas,
      cerrfin = checkerrfin,
      cwarwz = checkwarwz,
      cwarcon = checkwarcon,
      cwarany = checkwarany,
      cwarzc = checkwarzc,
      cabseps = checkabseps,
      cglhtany = checkglhtany,
      estprobs = estprobs,
      trueprobs = trueprobs
    )
  ))
  
}


## Global Simulation Function ---------------------------------------------

f_global_sim <- function(m_iprop, df_simdat) {
  
  l_prop <- replicate(df_simdat$nsim[1], m_iprop, simplify = FALSE)
  
  print(m_iprop)
  
  #plan(sequential, gc = TRUE)

  df_multsim_res <-
    do.call(rbind,
            lapply(l_prop, f_mult_sim, df_simdat))

  # in principle parallelisable, though often runs into RAM issues
  #plan(sequential)
  
  # df_multsim_res <-
  #   do.call(rbind,
  #           parallel::mclapply(l_prop, f_mult_sim, df_simdat, mc.cores = availableCores()))
  # 
  print("dfdone")
  print(Sys.time())
  
  propsname <- paste(m_iprop, collapse = ".")
  
  df_multsim_res <- cbind(props = propsname, df_multsim_res)
  
  if (file.exists(paste("mult_sim_", propsname, '.rds', sep = ""))) {
    boolFalse <- F
    while(boolFalse == F){
      try({df_complete_simres <-
        readRDS(paste("mult_sim_", propsname, '.rds', sep = ""))
        boolFalse <- T}, silent = FALSE)
      
    }
    
    saveRDS(
      rbind(df_complete_simres, df_multsim_res),
      paste("mult_sim_", propsname, '.rds', sep = ""),
      compress = FALSE
    )
  } else {
    saveRDS(df_multsim_res,
            paste("mult_sim_", propsname, '.rds', sep = ""),
            compress = FALSE)
  }
  
  rm(l_prop)
  rm(df_multsim_res)
  rm(propsname)
  gc()
  
}

# Settings ----------------------------------------------------------------

l_props <- list(
  matrix(c(0.33,0.33,0.33, #1 H0
           0.33,0.33,0.33,
           0.33,0.33,0.33,
           0.33,0.33,0.33),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.5,0.3,0.2, #2
           0.5,0.3,0.2,
           0.5,0.3,0.2,
           0.5,0.3,0.2),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.5,0.4,0.1, #3
           0.5,0.4,0.1,
           0.5,0.4,0.1,
           0.5,0.4,0.1),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.5,0.45,0.05, #4
           0.5,0.45,0.05,
           0.5,0.45,0.05,
           0.5,0.45,0.05),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.46,0.45,0.09, #5
           0.46,0.45,0.09,
           0.46,0.45,0.09,
           0.46,0.45,0.09),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.8,0.1,0.1, #6
           0.8,0.1,0.1,
           0.8,0.1,0.1,
           0.8,0.1,0.1),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.8, 0.15, 0.05, #7
           0.8, 0.15, 0.05,
           0.8, 0.15, 0.05,
           0.8, 0.15, 0.05),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.80, 0.19, 0.01, #8
           0.80, 0.19, 0.01,
           0.80, 0.19, 0.01,
           0.80, 0.19, 0.01),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.90, 0.05, 0.05, #9
           0.90, 0.05, 0.05,
           0.90, 0.05, 0.05,
           0.90, 0.05, 0.05),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.90, 0.08, 0.02, #10
           0.90, 0.08, 0.02,
           0.90, 0.08, 0.02,
           0.90, 0.08, 0.02),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.90, 0.09, 0.01, #11
           0.90, 0.09, 0.01,
           0.90, 0.09, 0.01,
           0.90, 0.09, 0.01),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.20, 0.30, 0.50, #12
           0.20, 0.30, 0.50,
           0.20, 0.30, 0.50,
           0.20, 0.30, 0.50),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.10, 0.40, 0.50, #13
           0.10, 0.40, 0.50,
           0.10, 0.40, 0.50,
           0.10, 0.40, 0.50),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.05, 0.45, 0.50, #14
           0.05, 0.45, 0.50,
           0.05, 0.45, 0.50,
           0.05, 0.45, 0.50),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.08, 0.45, 0.46, #15
           0.08, 0.45, 0.46,
           0.08, 0.45, 0.46,
           0.08, 0.45, 0.46),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.10, 0.10, 0.80, #16
           0.10, 0.10, 0.80,
           0.10, 0.10, 0.80,
           0.10, 0.10, 0.80),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.05, 0.15, 0.80, #17
           0.05, 0.15, 0.80,
           0.05, 0.15, 0.80,
           0.05, 0.15, 0.80),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.01, 0.19, 0.80, #18
           0.01, 0.19, 0.80,
           0.01, 0.19, 0.80,
           0.01, 0.19, 0.80),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.05, 0.05, 0.90, #19
           0.05, 0.05, 0.90,
           0.05, 0.05, 0.90,
           0.05, 0.05, 0.90),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.02, 0.08, 0.90, #20
           0.02, 0.08, 0.90,
           0.02, 0.08, 0.90,
           0.02, 0.08, 0.90),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.01, 0.09, 0.90, #21
           0.01, 0.09, 0.90,
           0.01, 0.09, 0.90,
           0.01, 0.09, 0.90),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.33, 0.33, 0.33, #22 HA
           0.33, 0.33, 0.33,
           0.50, 0.30, 0.20,
           0.50, 0.30, 0.20),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.50, 0.30, 0.20, #23 HA
           0.50, 0.30, 0.20,
           0.50, 0.40, 0.10,
           0.50, 0.40, 0.10),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.50, 0.40, 0.10, #24 HA
           0.50, 0.40, 0.10,
           0.50, 0.45, 0.05,
           0.50, 0.45, 0.05),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.50, 0.45, 0.05, #25 HA
           0.50, 0.45, 0.05,
           0.46, 0.45, 0.09,
           0.80, 0.10, 0.10),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.46, 0.45, 0.09, #26 HA
           0.46, 0.45, 0.09,
           0.80, 0.10, 0.10,
           0.80, 0.15, 0.05),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.80, 0.10, 0.10, #27 HA
           0.80, 0.10, 0.10,
           0.80, 0.15, 0.05,
           0.80, 0.19, 0.01),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.80, 0.15, 0.05, #28 HA
           0.80, 0.15, 0.05,
           0.80, 0.19, 0.01,
           0.90, 0.05, 0.05),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.80, 0.19, 0.01, #29 HA
           0.80, 0.19, 0.01,
           0.80, 0.15, 0.05,
           0.80, 0.10, 0.10),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.90, 0.05, 0.05, #30 HA
           0.90, 0.05, 0.05,
           0.90, 0.08, 0.02,
           0.90, 0.09, 0.01),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.90, 0.08, 0.02, #31 HA
           0.90, 0.09, 0.01,
           0.90, 0.09, 0.01,
           0.90, 0.09, 0.01),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.20, 0.30, 0.50, #32 HA
           0.20, 0.30, 0.50,
           0.33, 0.33, 0.33,
           0.33, 0.33, 0.33),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.10, 0.40, 0.50, #33 HA
           0.10, 0.40, 0.50,
           0.20, 0.30, 0.50,
           0.33, 0.33, 0.33),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.05, 0.45, 0.50, #34 HA
           0.05, 0.45, 0.50,
           0.10, 0.40, 0.50,
           0.20, 0.30, 0.50),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.08, 0.45, 0.47, #35 HA
           0.08, 0.45, 0.47,
           0.05, 0.45, 0.50,
           0.08, 0.45, 0.47),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.10, 0.10, 0.80, #36 HA
           0.10, 0.10, 0.80,
           0.05, 0.45, 0.50,
           0.10, 0.40, 0.50),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.05, 0.15, 0.80, #37 HA
           0.05, 0.15, 0.80,
           0.10, 0.10, 0.80,
           0.08, 0.45, 0.47),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.01, 0.19, 0.80, #38 HA
           0.01, 0.19, 0.80,
           0.05, 0.15, 0.80,
           0.10, 0.10, 0.80),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.05, 0.05, 0.90, #39 HA
           0.05, 0.05, 0.90,
           0.05, 0.15, 0.80,
           0.08, 0.45, 0.46),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.02, 0.08, 0.90, #40 HA
           0.01, 0.09, 0.90,
           0.01, 0.09, 0.90,
           0.01, 0.09, 0.90),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.33, 0.33, 0.33, #41 HA
           0.33, 0.33, 0.33,
           0.20, 0.30, 0.50,
           0.20, 0.30, 0.50),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.50, 0.30, 0.20, #42 HA
           0.50, 0.30, 0.20,
           0.33, 0.33, 0.33,
           0.20, 0.30, 0.50),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.50, 0.40, 0.10, #43 HA
           0.50, 0.40, 0.10,
           0.50, 0.30, 0.20,
           0.33, 0.33, 0.33),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.50, 0.45, 0.05, #44 HA
           0.50, 0.40, 0.10,
           0.50, 0.45, 0.05,
           0.20, 0.30, 0.50),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.46, 0.45, 0.09, #45 HA
           0.46, 0.45, 0.09,
           0.50, 0.45, 0.05,
           0.20, 0.30, 0.50),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.80, 0.10, 0.10, #46 HA
           0.80, 0.10, 0.10,
           0.46, 0.45, 0.09,
           0.10, 0.40, 0.50),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.80, 0.15, 0.05, #47 HA
           0.80, 0.15, 0.05,
           0.80, 0.10, 0.10,
           0.46, 0.45, 0.09),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.80, 0.19, 0.01, #48 HA
           0.80, 0.19, 0.01,
           0.46, 0.45, 0.09,
           0.05, 0.45, 0.50),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.90, 0.05, 0.05, #49 HA
           0.90, 0.05, 0.05,
           0.80, 0.19, 0.01,
           0.20, 0.30, 0.05),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.90, 0.08, 0.02, #50 HA
           0.90, 0.05, 0.05,
           0.80, 0.19, 0.01,
           0.80, 0.15, 0.05),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.20, 0.30, 0.50, #51 HA
           0.20, 0.30, 0.50,
           0.10, 0.40, 0.50,
           0.05, 0.45, 0.50),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.10, 0.40, 0.50, #52 HA
           0.10, 0.40, 0.50,
           0.05, 0.45, 0.50,
           0.08, 0.45, 0.47),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.05, 0.45, 0.50, #53 HA
           0.05, 0.45, 0.50,
           0.08, 0.45, 0.47,
           0.10, 0.10, 0.80),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.08, 0.45, 0.46, #54 HA
           0.08, 0.45, 0.46,
           0.10, 0.10, 0.80,
           0.05, 0.15, 0.80),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.10, 0.10, 0.80, #55 HA
           0.10, 0.10, 0.80,
           0.05, 0.15, 0.80,
           0.01, 0.19, 0.80),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.05, 0.15, 0.80, #56 HA
           0.05, 0.15, 0.80,
           0.01, 0.19, 0.80,
           0.05, 0.05, 0.90),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.01, 0.19, 0.80, #57 HA
           0.01, 0.19, 0.80,
           0.05, 0.05, 0.90,
           0.02, 0.08, 0.90),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.05, 0.05, 0.90, #58 HA
           0.05, 0.05, 0.90,
           0.02, 0.08, 0.90,
           0.01, 0.09, 0.90),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.02, 0.08, 0.90, #59 HA
           0.05, 0.05, 0.90,
           0.01, 0.19, 0.80,
           0.05, 0.15, 0.80),
         ncol = 3,
         byrow = TRUE),
  matrix(c(0.90, 0.05, 0.05, #60 HA
           0.05, 0.90, 0.05,
           0.05, 0.05, 0.90,
           0.33, 0.33, 0.33),
         ncol = 3,
         byrow = TRUE)
)

# 50 Simulations per Iteration takes already quite long
# set nsim and sim_iterations depending on CPU core count
df_simdat <- expand.grid(
  phi = c(1.01, 1.5, 2, 5, 8),
  b = c(5, 10, 20, 50),
  m = c(10, 20, 50, 100),
  comp = c("Tukey","Dunnett"),
  modeltype = c("multinomial","pearson","afroz","farrington","deviance","DM_VGAM","DM_MGLM"),
  dfu = c("modeldf"), #only multivariate t-distribution was considered. Set to "zero" for M-normal distribution
  nsim = 20
)

# Split up into multiple Iterations to get intermediate results
sim_iterations <- 5
system.time(for (i in 1:sim_iterations) {
  print(paste("Iteration",i))
  system.time(for (i in 1:length(l_props)) {
    f_global_sim(l_props[[i]], df_simdat)
  })
}
)


