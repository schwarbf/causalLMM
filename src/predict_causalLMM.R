#' @title predict_causalLMM
#'
#' @description This function computes predictions on new test data using a 
#'    fitted causalLMM model. The fit here corresponds to the estimated
#'    fixed- and random effects. The random effects need to be estimated 
#'    seperately within each environment. 
#'
#' @param fit       list: Contains 2 elements: $fixed_effects: vector of estimated
#'                  fixed effects and $random_effects: data.frame of random
#'                  effects predictors for all environments. 
#' @param newdata   list: Contains 2 elements: $X data.frame of fixed-effects 
#'                  predictors and $Z data.frame of random-effects predictors. 
#' @param ExpInd    vector: Contains the indicators for the environments
#' 
#' @return vector: Predictions of response in each existing environment using
#'   fixed- and random effects. 

predict_causalLMM <- function(fit, newdata, ExpInd){
  # extract tuning parameters
  m <- unique(ExpInd) %>% length()
  
  # compute the prediction for existing environments including random effects
  Ye_hat <- list()
  for(i in 1:m){
    ind <- which(ExpInd == i)
    Ye_hat[[i]] <- newdata$X[ind, ] %*% fit$fixed_effects + newdata$Z[ind, ] %*% fit$random_effects[[i]]
  }
  Y_hat <- unlist(Ye_hat)
  
  return(Y_hat)
}