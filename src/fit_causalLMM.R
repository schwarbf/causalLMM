#' @title fit_causalLMM
#'
#' @description This function computes the causalLMM estimator for a prespecified
#'    value of gamma which corresponds to the amount of causal regularization/
#'    expected strength of shift perturbations.
#'
#' @param X data.frame: Contains the fixed-effects predictor
#' @param Y vector: Contains the response value for all environments
#' @param Z data.frame: Contains the random-effects predictors 
#' @param S data.frame: Contains the shift variables.
#' @param gamma float: Causal regularization tuning parameter.
#' @param ExpInd vector: Contains the indicators for the environments
#' @param tol float: Specifies the convergence tolerance for when the iterative 
#'    procedure to compute the causalLMM estimator should terminate. 
#' @param verbose boolean: If TRUE, then the alternating fitting steps for 
#'    stimating beta_causalLMM and Sigma_e are printed. 
#' 
#' @return Returns the model fit consisting of a list with 2 elements: 
#'    $fixed-effects: vector of estimated fixed-effects
#'    $random-effects: data.frame of estimated random-effects within each 
#'                     environment. 

fit_causalLMM <- function(X, Y, Z, S, gamma, ExpInd, tol = 1e-03, verbose = F){
  # extract tuning parameters
  m <- unique(ExpInd) %>% length()
  n <- length(ExpInd)
  
  # add an intercept to model if it is one-dimensional and vector operations 
  # don't work anymore
  if(is.vector(X)){X <- cbind(1, X)}
  if(is.vector(Z)){Z <- cbind(1, Z)}
  
  # extract parameters per environment 
  Xe <- Ye <- Ze <- Se <- PSe <- ne <- list()
  for(i in 1:m){
    ind <- which(ExpInd == i)
    ne[[i]] <- length(ind)
    Xe[[i]] <- X[ind, ]
    Ye[[i]] <- Y[ind]
    Ze[[i]] <- Z[ind, ]
    Se[[i]] <- if(is.vector(S)){S[ind]} else{S[ind, ]}
    PSe[[i]] <- Se[[i]] %*% solve(crossprod(Se[[i]])) %*% t(Se[[i]])
  }
  
  # helper function for profiled log likelihood
  estim_beta_causalLME <- function(Sigmae_hat){
    # compute Xe_tilde and Ye_tilde
    Xe_tilde <- Ye_tilde <- Sigmae_hat_inv <- list()
    for(i in 1:m){
      if(ncol(X) == 2){Xe[[i]] <- Xe[[i]][, 2]}
      Sigmae_hat_inv[[i]] <- solve(Sigmae_hat[[i]])
      Xe_tilde[[i]] <- (diag(ne[[i]]) + (sqrt(gamma) - 1)*PSe[[i]]) %*% Xe[[i]]
      Ye_tilde[[i]] <- (diag(ne[[i]]) + (sqrt(gamma) - 1)*PSe[[i]]) %*% Ye[[i]]
      if(ncol(X) == 2){Xe_tilde[[i]] <- cbind(1, Xe_tilde[[i]])}
    }
    
    # helper functions
    helper_part1 <- function(i){
      1/ne[[i]]*crossprod(Xe_tilde[[i]], Sigmae_hat_inv[[i]] %*% Xe_tilde[[i]])
    }
    helper_part2 <- function(i){
      1/ne[[i]]*crossprod(Xe_tilde[[i]], Sigmae_hat_inv[[i]] %*% Ye_tilde[[i]])
    }

    # computing beta_causalLME
    part1 <- lapply(1:m, helper_part1) %>% Reduce('+', .)
    part2 <- lapply(1:m, helper_part2) %>% Reduce('+', .)
    beta_causalLME_hat <- solve(part1, part2)
    
    return(beta_causalLME_hat)
  }
  
  # helper function to estimate Sigmae
  estim_Sigmae <- function(beta_hat){
    # estimating be
    Z_blockdiag <- bdiag(Ze) %>% as.matrix()
    b <- lm(Y - X %*% beta_hat ~ Z_blockdiag -1)$coefficients %>% as.vector()
    be_hat <- list()
    for(i in 0:(m-1)){
      ind <- (ncol(Z)*i + 1):(ncol(Z)*i + ncol(Z))
      be_hat[[i+1]] <- b[ind]
    }
    
    # estimating sigma^2 and tau^2
    helper_sigma2_hat <- function(i){
      return(norm(Ye[[i]] - Xe[[i]] %*% beta_hat - Ze[[i]] %*% be_hat[[i]], '2'))
    }
    helper_tau2_hat <- function(i, sigma2_hat){
      trZZt <- sum(diag(Ze[[i]] %*% t(Ze[[i]])))
      return(1/trZZt*(norm(Ye[[i]] - Xe[[i]] %*% beta_hat, '2') - ne[[i]]*sigma2_hat))
    }
    
    sigma2_hat <- 1/n*sapply(1:m, helper_sigma2_hat) %>% sum()
    tau2_hat <- sapply(1:m, helper_tau2_hat, sigma2_hat = sigma2_hat) %>% sum()
    Sigmae_hat <- lapply(1:m, function(i){tau2_hat*Ze[[i]] %*% t(Ze[[i]]) + sigma2_hat*diag(ne[[i]])})
    
    return(list(be_hat = be_hat, Sigmae_hat = Sigmae_hat, sigma2_hat = sigma2_hat, 
                tau2_hat = tau2_hat))
  }
  
  # alternating estimation of beta_hat and Sigmae_hat
  step <- 1
  not_converged <- TRUE
  while(not_converged){
    # optimization over beta
    if(step == 1){
      Sigmae_hat <- lapply(1:m, function(i){diag(ne[[i]])})
      beta_hat <- estim_beta_causalLME(Sigmae_hat)
      res <- matrix(list(beta_hat, 1, 1), ncol = 3)
      colnames(res) <- c('beta', 'sigma2', 'tau2')
    } else{
      beta_hat <- estim_beta_causalLME(Sigmae_hat)
    }

    # estimtation of Sigmae
    out <- estim_Sigmae(beta_hat)
    be_hat <- out$be_hat
    Sigmae_hat <- out$Sigmae_hat
    sigma2_hat <- out$sigma2_hat
    tau2_hat <- out$tau2_hat
    res_inter <- matrix(list(beta_hat, sigma2_hat, tau2_hat), ncol = 3)
    colnames(res_inter) <- c('beta', 'sigma2', 'tau2')
    res <- rbind(res, res_inter) %>% as.data.frame()
      
    # check for convergence
    if(step > 1 && norm(beta_hat - prev_beta_hat, '2') < tol){not_converged <- FALSE}
    if(not_converged){
      prev_beta_hat <- beta_hat
      step <- step + 1
    }
  }
  
  # concatenate results to list
  results <- list(fixed_effects = beta_hat, random_effects = be_hat, 
                  sigma2 = sigma2_hat, tau2 = tau2_hat)

  # output the intermediate results of the optimization
  if(verbose){print(res)}
  
  # return the last vector of fitted fixed effects
  return(results)
}