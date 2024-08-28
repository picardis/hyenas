predict_tp <- function (mod, newdata, which) {
  
  mm <- model.matrix(mod$conditions$formula, newdata)
  
  nStates <- length(mod$stateNames)
  
  start <- length(mod$CIbeta$step$est) + 
    length(mod$CIbeta$angle$est) + 1
  end <- start + ((nStates^2) - nStates) * ncol(mm) - 1
  
  est_mat <- matrix(nrow = ncol(mm), 
                    ncol = ((nStates^2) - nStates),
                    data = start:end)
  
  cnames <- c()
  
  for (i in 1:nStates) {
    for (j in 1:nStates) {
      if (i != j) {
        cnames <- c(cnames, paste0(i, "t", j))
      }
    }
  }
  
  colnames(est_mat) <- cnames
  
  # Are we estimating a transition to a different state or
  # the probability of staying in the current?
  states <- unlist(strsplit(which, split = "t"))
  
  # Transition to different state
  if (states[1] != states[2]) {
    
  B <- mod$mod$est[est_mat[, which]]
  
  pred <- mm %*% B
  
  Sigma <- mod$mod$Sigma[est_mat[, which], est_mat[, which]]
  
  var <- diag(mm %*% Sigma %*% t(mm))
  
  if (any(var < 0)){
    warning(paste("The model estimated at least one negative variance",
    "for transition", which))
    ## NOTE: this is probably a bad fix. Perhaps the function should stop if 
    # any(var < 0)
    #If variance is < 0, set to 0
    var[which(var < 0)] <- 0
  }
  
  #Standard error
  se <- sqrt(var)
  
  lwr <- pred - 1.96 * se
  upr <- pred + 1.96 * se
  
  res <- data.frame(lwr = lwr, est = pred, upr = upr)
  # res <- cbind(newdata, res)
  
  } else if (states[1] == states[2]) {
    # Remaining in the current state
    
    pred <- rep(0, nrow(newdata))
    lwr <- rep(0, nrow(newdata))
    upr <- rep(0, nrow(newdata))
  }
  
  res <- data.frame(lwr = lwr, est = pred, upr = upr, tp = which)
  res <- cbind(newdata, res) 
  return(res)
  
}

####

predict_hmm <- function (mod, newdata) {
  
  nStates <- length(mod$stateNames)
  
  trans <- c()
  
  link <- list()
  resp <- list()
  
  for (i in 1:nStates) {
    link[[i]] <- list()
    for (j in 1:nStates) {
      tt <- paste0(i, "t", j)
        link[[i]][[j]] <- predict_tp(mod = mod, newdata = newdata, which = tt)
      }
    }
    # Once all transitions have been calculated on link scale, backtransform
    # Link is multinomial logit. See `momentuHMM` vignette for reference to
    # Michelot et al. 2016 (moveHMM Methods in Ecology & Evolution article)
  for (i in 1:nStates){  
  explink <- lapply(link[[i]], function(x){
      x$lwr <- exp(x$lwr)
      x$est <- exp(x$est)
      x$upr <- exp(x$upr)
      return(x)
    })
    
    exp_vals <- lapply(explink, function(x){
      return(x[, c("lwr", "est", "upr")])
    })
    
    sum_exp_vals <- Reduce("+", exp_vals)
    
    resp_vals <- lapply(exp_vals, function(x){
      x/sum_exp_vals
    })
    
    bind_vals <- bind_rows(resp_vals)
    
    resp[[i]] <- bind_rows(link[[i]])
    resp[[i]]$lwr <- bind_vals$lwr
    resp[[i]]$est <- bind_vals$est
    resp[[i]]$upr <- bind_vals$upr
    rownames(resp[[i]]) <- NULL
  }
  
  tp_preds <- bind_rows(resp)
  return(tp_preds)
  
}


