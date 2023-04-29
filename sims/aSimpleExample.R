# ******************************************************************************
#                             A SIMPLE EXAMPLE
# ******************************************************************************

# ------------------------------------------------------------------------------
# PREAMBLE
# ------------------------------------------------------------------------------
# load packages
pkgs <- c("dplyr", "Matrix", "lme4", "MASS", "expm", "ggplot2", "latex2exp", 
          "reshape2", "doParallel", "polynom")
for(pkg in pkgs){
  if(!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# setting the correct working directory
if("florianschwarb" %in% Sys.info()){
  wdir <- "/Users/florianschwarb/Desktop/Master-Thesis/Code/causalLMM/"
} else{
  wdir <- getwd()
}
setwd(paste0(wdir, "src"))

source("fitcausalLMM.R")
source("predictcausalLMM.R")

# set seed
set.seed(1)

# ------------------------------------------------------------------------------
# A SIMPLE EXAMPLE
# ------------------------------------------------------------------------------
# helper function to run the simulation
runSim_simple <- function(nenv = 50, ne = 200){
  # observational data
  eps_obs <- rnorm(nenv*ne)
  etax_obs <- rnorm(nenv*ne)
  etah_obs <- rnorm(nenv*ne)
  b_obs <- rnorm(nenv*1, sd = 0.5)
  S_obs <- rnorm(nenv*ne)
  Z_obs <- rnorm(nenv*ne)
  H_obs <- etah_obs
  X_obs <- S_obs + 2*H_obs + etax_obs
  Y_obs <- X_obs - 2*b_obs*Z_obs - 3*H_obs + eps_obs
  
  # shifted data
  eps_v <- rnorm(nenv*ne)
  etax_v <- rnorm(nenv*ne)
  etah_v <- rnorm(nenv*ne)
  b_v <- rnorm(nenv*1, sd = 0.5)
  Z_v <- rnorm(nenv*ne)
  H_v <- etah_v
  X_v <- 2 + 2*H_v + etax_v
  Y_v <- X_v + 2*b_v*Z_v - 3*H_v + eps_v
  
  # computing LMM estimator using package lme4
  ExpInd <- sapply(1:nenv, function(i){rep(i, ne)}) %>% as.vector()
  data_obs <- data.frame(Y = Y_obs, X = X_obs, Z = Z_obs)
  form <- Y ~ X + (Z | ExpInd)
  lmer.fit <- lmer(form, data = data_obs, REML = F) %>% 
    suppressMessages() %>%
    suppressWarnings()
  
  # computing the mse (over all envs)
  mse_env <- c()
  for(env in 1:nenv){
    ind <- which(ExpInd == env)
    mse_env[env] <- mean((Y_v[ind] - X_v[ind]*fixef(lmer.fit)[2])^2)
  }
  mse <- mean(mse_env)
  
  return(list(mse = mse, beta = fixef(lmer.fit)))
}

# MSE for varying beta
# ------------------------------------------------------------------------------
# tuning parameters
nsim <- 100
nenv <- 50
ne <- 200

# initializing the cluster
nCores <- detectCores()
cl <- makeCluster(nCores - 1)
clusterExport(cl, 'runSim_simple')

# running simulation in parallel
results <- foreach(sim = 1:nsim, .combine = rbind) %do% {
    runSim_simple(nenv = nenv, ne = ne)
}

# shutting down cluster
stopCluster(cl)

# averaging
summary <- data.frame(matrix(NA, ncol = 2, nrow = 3))
colnames(summary) <- c('beta', 'mse')
summary[3, 'mse'] <- mean(results[, 'mse'] %>% unlist() %>% as.vector())
tmp <- count <- 0
for(i in 1:length(results)){
  if(!is.na(results[[i]][2])) {
    tmp <- tmp + results[[i]][2] %>% as.double()
    count <- count + 1
  }
}
summary[3, 'beta'] <- tmp/count


# computing the quadratic through the 3 points using Lagrange polynomials
poly.calc(x = summary[, 'beta'], y = summary[, 'mse'])
parabola <- function(x){
  7.364637*x^2 - 4.70764*x + 8.332094
}

# plotting parabola through the 3 points
ind_PA <- which(results[, 'gamma'] == 'PA-LMM')
ind_LS <- which(results[, 'gamma'] == 'LS-LMM')
ind_IV <- which(results[, 'gamma'] == 'IV-LMM')
results <- results[, c('beta', 'mse')] %>% as.data.frame()
add_on <- (results[ind_PA, 'beta'] - results[ind_IV, 'beta'])/1000*20
beta <- seq(results[ind_PA, 'beta'] + add_on, results[ind_IV, 'beta'] - add_on, length.out = 1040)
mse <- sapply(beta, parabola) %>% as.vector()
data <- data.frame(beta = beta, mse = mse)
p_simple <- ggplot(data, aes(x = beta, y = mse)) +
  geom_line() + 
  geom_point(aes(x = results[ind_PA, 'beta'], y = results[ind_PA, 'mse']), size = 3) +
  annotate('text', x= results[ind_PA, 'beta'] - 4*add_on, y = results[ind_PA, 'mse'] + 0.3, label = 'PA-LMM') +
  geom_point(aes(x = results[ind_IV, 'beta'], y = results[ind_IV, 'mse']), size = 3) +
  annotate('text', x= results[ind_IV, 'beta'] + 4*add_on, y = results[ind_IV, 'mse'] + 0.3, label = 'IV-LMM') +
  geom_point(aes(x = results[ind_LS, 'beta'], y = results[ind_LS, 'mse']), size = 3) +
  annotate('text', x= results[ind_LS, 'beta'] - 4*add_on, y = results[ind_LS, 'mse'] + 0.3, label = 'LS-LMM') +
  theme(plot.title = element_text(face = 'bold'), 
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_rect(colour = "black", fill =NA, size = 1), 
        legend.position = 'none') +
  ylab('MSE') +
  xlab(TeX('$\\beta_{causalLMM}^{\\gamma}$')) +
  labs(title = 'Mean squared Error for shifted distribution')
p_simple

setwd("/Users/florianschwarb/Desktop/Master-Thesis/Code/causalLMM/fig")
ggsave('mse_a-simple-example.pdf', width = 5, height = 4)

# ------------------------------------------------------------------------------
# EMPIRICAL ANALYSIS
# ------------------------------------------------------------------------------
# helper function to run the simulation
runSim_advanced <- function(nenv = 10, ne = 200, gamma = 1, strong_shifts = F, 
                            random_slopes = F){
  # observational data
  eps_obs <- rnorm(nenv*ne)
  etax_obs <- mvrnorm(nenv*ne, mu = c(rep(0, 5)), Sigma = 0.5^2*diag(5))
  H_obs <- rnorm(nenv*ne)
  S_obs <- mvrnorm(nenv*ne, mu = c(0, 0), Sigma = diag(2))
  Z_obs <- mvrnorm(nenv*ne, mu = c(0, 0), Sigma = diag(2))
  X_obs <- S_obs[, 1] + 4*S_obs[, 2] + 2*H_obs + etax_obs
  if(!random_slopes){
    b_obs <- rnorm(nenv*1, sd = 1.5)
    Y_obs <- X_obs[, 1] + 5*X_obs[, 3] + b_obs*Z_obs[, 2] - H_obs + eps_obs
  } else{
    b_obs <- mvrnorm(nenv*1, mu = c(0, 0), Sigma = 1.5^2*diag(2))
    Y_obs <- (1 + b_obs[, 1])*X_obs[, 1] + (5 + b_obs[, 2])*X_obs[, 3] - H_obs + eps_obs
  }
  
  # shifted data
  eps_v <- rnorm(nenv*ne)
  etax_v <- mvrnorm(nenv*ne, mu = c(rep(0, 5)), Sigma = 0.5^2*diag(5))
  H_v <- rnorm(nenv*ne)
  if(!strong_shifts){
    mu <- mvrnorm(1, mu = c(2, 2), Sigma = 0.25^2*diag(2))
    t <- rnorm(1, mean = 1, sd = 0.25)
  } else{
    mu <- mvrnorm(1, mu = c(5, 5), Sigma = 0.5^2*diag(2))
    t <- rnorm(1, mean = 5, sd = 0.5)
  }
  delta_v <- mvrnorm(nenv*ne, mu = mu, Sigma = t^2*diag(2))
  Z_v <- mvrnorm(nenv*ne, mu = c(0, 0), Sigma = diag(2))
  X_v <- delta_v[, 1] + 4*delta_v[, 2] + 2*H_v + etax_v
  if(!random_slopes){
    b_v <- rnorm(nenv*1, sd = 1.5)
    Y_v <- X_v[, 1] + 5*X_v[, 3] + b_v*Z_v[, 2] - H_v + eps_v
  } else{
    b_v <- mvrnorm(nenv*1, mu = c(0, 0), Sigma = 1.5^2*diag(2))
    Y_v <- (1 + b_v[, 1])*X_v[, 1] + (5 + b_v[, 2])*X_v[, 3] - H_v + eps_v
  }
  
  # preparing for estimation
  ExpInd <- sapply(1:nenv, function(i){rep(i, ne)}) %>% as.vector()
  data_obs <- data.frame(Y = Y_obs, X1 = X_obs[, 1], X2 = X_obs[, 2], 
                         X3 = X_obs[, 3], X4 = X_obs[, 4], X5 = X_obs[, 5], 
                         Z1 = Z_obs[, 1], Z2 = Z_obs[, 2])
  data_v <- data.frame(Y = Y_v, X1 = X_v[, 1], X2 = X_v[, 2], X3 = X_v[, 3], 
                       X4 = X_v[, 4], X5 = X_v[, 5], Z1 = Z_v[, 1], Z2 = Z_v[, 2])
  
  # computing causalLME estimator
  causalLMM.fit <- fit_causalLMM(X_obs, Y_obs, Z_obs, S_obs, gamma, ExpInd)
  
  # computing LMM estimator using package lme4
  form <- Y ~ -1 + X1 + X2 + X3 + X4 + X5 + (Z1 + Z2 -1 | ExpInd)
  lmer.fit <- lmer(form, data = data_obs, REML = F) %>% 
    suppressMessages() %>%
    suppressWarnings()
  lmer_Notconverged <- function(m){
    df <- summary(m)
    !is.null(df$optinfo$conv$lme4$messages) && 
      grepl('failed to converge', df$optinfo$conv$lme4$messages)
  }
  
  # only use values if lmer converged
  if(!lmer_Notconverged(lmer.fit)){
    # predictions on shifted data
    Y_v_hat_causalLMM <- predict_causalLMM(causalLMM.fit, 
                                           newdata = list(X = X_v, Z = Z_v), 
                                           ExpInd)
    Y_v_hat_lmer <- predict(lmer.fit, newdata = data_v, re.form = NULL)
    
    # absolute prediction errors
    abs_err_causalLMM <- abs(Y_v - Y_v_hat_causalLMM)
    abs_err_lmer <- abs(Y_v - Y_v_hat_lmer) %>% as.vector()
    
    return(list(causalLMM = abs_err_causalLMM, lmer = abs_err_lmer))
  } else{
    return(list())
  }
}

# Moderate/Strong Shifts
# ------------------------------------------------------------------------------
# tuning parameters
nsim <- 25
nenv <- 50
ne <- 200
gamma <- 7
strong_shifts <- T

# initializing the cluster
nCores <- detectCores()
cl <- makeCluster(nCores - 1)
clusterExport(cl, 'runSim_advanced')

# running simulation in parallel
results <- foreach(sim = 1:nsim) %do% {
  sim_out <- runSim_advanced(nenv = nenv, ne = ne, gamma = gamma, 
                             strong_shifts = strong_shifts) 
  
  # compute empirical quantiles of absolute prediction errors
  if(length(sim_out) > 0){
    quant_causalLMM <- quantile(sim_out$causalLMM, probs = seq(0, 1, length.out = 10))
    quant_lmer <- quantile(sim_out$lmer, probs = seq(0, 1, length.out = 10))
    quants <- data.frame(alpha = seq(0, 1, length.out = 10), causalLMM = quant_causalLMM, 
                         LMM = quant_lmer)
    rownames(quants) <- 1:10
    quants
  }
}

# shutting down cluster
stopCluster(cl)

# average the quantiles over the simulation runs
quants <- data.frame(matrix(NA, ncol = 3, nrow = 10))
colnames(quants) <- c('alpha', 'causalLMM', 'LMM')
quants$alpha <- seq(0, 1, length.out = 10)
quants_causalLMM <- quants_LMM <- rep(0, 10)
for(sim in 1:nsim){
  if(length(results[[sim]] > 0)){
    quants_causalLMM <- quants_causalLMM + results[[sim]][, 'causalLMM']
    quants_LMM <- quants_LMM + results[[sim]][, 'LMM']
  }
}
quants$causalLMM <- 1/nsim*quants_causalLMM
quants$LMM <- 1/nsim*quants_LMM

# Plotting
if(strong_shifts){
  title <- 'Strong Shifts'
  filename <- 'empirical-analysis_strong-shifts.pdf'
} else{
  title <- 'Moderate Shifts'
  filename <- 'empirical-analysis_moderate-shifts.pdf'
}
quants <- melt(quants, id = 'alpha')
p_strong <- ggplot(quants, aes(x = alpha, y = value, color = variable)) +
  geom_point(shape = 1) + 
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_rect(colour = "black", fill =NA, size = 0.5), 
        legend.key=element_blank(), 
        legend.position = "bottom", 
        axis.title=element_text(size=8, face = "plain")) +
  scale_color_manual(values = c('red', 'blue')) +
  ylab('alpha-quantile of absolute prediction error') +
  labs(title = title) +
  guides(color = guide_legend(title = 'Legend'))
p_strong

setwd("/Users/florianschwarb/Desktop/Master-Thesis/Code/causalLMM/fig")
ggsave(filename, width = 5, height = 5)

# Model Missspecifications - Random Slopes
# ------------------------------------------------------------------------------
# tuning parameters
nsim <- 25
nenv <- 50
ne <- 200
gamma <- 7
strong_shifts <- TRUE

# initializing the cluster
nCores <- detectCores()
cl <- makeCluster(nCores - 1)
clusterExport(cl, 'runSim_advanced')

# running simulation in parallel
results <- foreach(sim = 1:nsim) %do% {
  sim_out <- runSim_advanced(nenv = nenv, ne = ne, gamma = gamma, 
                             strong_shifts = strong_shifts, random_slopes = T) 
  
  # compute empirical quantiles of absolute prediction errors
  if(length(sim_out) > 0){
    quant_causalLMM <- quantile(sim_out$causalLMM, probs = seq(0, 1, length.out = 10))
    quant_lmer <- quantile(sim_out$lmer, probs = seq(0, 1, length.out = 10))
    quants <- data.frame(alpha = seq(0, 1, length.out = 10), causalLMM = quant_causalLMM, 
                         LMM = quant_lmer)
    rownames(quants) <- 1:10
    quants
  }
}

# shutting down cluster
stopCluster(cl)

# average the quantiles over the simulation runs
quants <- data.frame(matrix(NA, ncol = 3, nrow = 10))
colnames(quants) <- c('alpha', 'causalLMM', 'LMM')
quants$alpha <- seq(0, 1, length.out = 10)
quants_causalLMM <- quants_LMM <- rep(0, 10)
for(sim in 1:nsim){
  if(length(results[[sim]] > 0)){
    quants_causalLMM <- quants_causalLMM + results[[sim]][, 'causalLMM']
    quants_LMM <- quants_LMM + results[[sim]][, 'LMM']
  }
}
quants$causalLMM <- 1/nsim*quants_causalLMM
quants$LMM <- 1/nsim*quants_LMM

# Plotting
quants <- melt(quants, id = 'alpha')
p_strong <- ggplot(quants, aes(x = alpha, y = value, color = variable)) +
  geom_point(shape = 1) + 
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_rect(colour = "black", fill =NA, size = 0.5), 
        legend.key=element_blank(), 
        legend.position = "bottom", 
        axis.title=element_text(size=8, face = "plain")) +
  scale_color_manual(values = c('red', 'blue')) +
  ylab('alpha-quantile of absolute prediction error') +
  guides(color = guide_legend(title = 'Legend'))
p_strong

setwd("/Users/florianschwarb/Desktop/Master-Thesis/Code/causalLMM/fig")
ggsave('empirical-analysis_model-missspecifications.pdf', width = 5, height = 5)

# Effect of gamma
# ------------------------------------------------------------------------------
# tuning parameters
nsim <- 25
nenv <- 50
ne <- 200
gammas <- c(0, 0.5, 1, 3, 7, 10, 16)
strong_shifts <- TRUE

# initializing the cluster
nCores <- detectCores()
cl <- makeCluster(nCores - 1)
clusterExport(cl, 'runSim_advanced')

# running simulation in parallel
for(gamma in gammas){
  results <- foreach(sim = 1:nsim) %do% {
    sim_out <- runSim_advanced(nenv = nenv, ne = ne, gamma = gamma, 
                               strong_shifts = strong_shifts) 
    
    # compute empirical quantiles of absolute prediction errors
    if(length(sim_out) > 0){
      quant_causalLMM <- quantile(sim_out$causalLMM, probs = seq(0, 1, length.out = 10))
      quant_lmer <- quantile(sim_out$lmer, probs = seq(0, 1, length.out = 10))
      quants <- data.frame(alpha = seq(0, 1, length.out = 10), causalLMM = quant_causalLMM, 
                           LMM = quant_lmer)
      rownames(quants) <- 1:10
      quants
    }
  }
  
  # average the quantiles over the simulation runs
  quants <- data.frame(matrix(NA, ncol = 3, nrow = 10))
  colnames(quants) <- c('alpha', 'causalLMM', 'LMM')
  quants$alpha <- seq(0, 1, length.out = 10)
  quants_causalLMM <- quants_LMM <- rep(0, 10)
  for(sim in 1:nsim){
    if(length(results[[sim]] > 0)){
      quants_causalLMM <- quants_causalLMM + results[[sim]][, 'causalLMM']
      quants_LMM <- quants_LMM + results[[sim]][, 'LMM']
    }
  }
  quants$causalLMM <- 1/nsim*quants_causalLMM
  quants$LMM <- 1/nsim*quants_LMM
  
  # merging with overall estimated empirical
  if(gamma == gammas[1]){
    quants_overall <- quants
    colnames(quants_overall) <- c('alpha', gamma, 'LMM')
  } else{
    quants_overall[, as.character(gamma)] <- quants$causalLMM
  }
}
quants_overall <- quants_overall[, c('alpha', as.character(gammas), 'LMM')]

# shutting down cluster
stopCluster(cl)

# Plotting
quants_overall_melted <- melt(quants_overall, id = 'alpha')
p_comp_gamma <- ggplot(quants_overall_melted, aes(x = alpha, y = value, 
                                                  color = variable, shape = variable)) +
  geom_point() + 
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_rect(colour = "black", fill =NA, size = 0.5), 
        legend.key=element_blank(), 
        legend.position = "bottom", 
        axis.title=element_text(size=8, face = "plain")) +
  ylab('alpha-quantile of absolute prediction error') +
  guides(color = guide_legend(title = 'gamma')) + 
  scale_colour_manual(name = 'gamma', 
                      labels = c('0', '0.5', '1', '3', '7', '10', '16', 'LMM'), 
                      values = c('yellow2', 'orange', 'mediumpurple1', 'green3', 
                                 'royalblue', 'blue', 'cyan1', 'red')) + 
  scale_shape_manual(name = 'gamma', 
                     labels = c('0', '0.5', '1', '3', '7', '10', '16', 'LMM'), 
                     values = c(1, 1, 1, 1, 1, 1, 1, 15))
p_comp_gamma

setwd("/Users/florianschwarb/Desktop/Master-Thesis/Code/causalLMM/fig")
ggsave('empirical-analysis_comparison-gamma.pdf', width = 5, height = 5)
