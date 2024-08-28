# Write a function to compute shape and scale of Gamma given mean and sd
# Otherwise I'm just picking initial values with my eyes closed
gamma_pars <- function(mean = NULL, sd = NULL, 
                       scale = NULL, shape = NULL) {
  if(exists("mean") & exists("sd")) {
    
    scale_theta <- sqrt(sd^2/mean)
    shape_kappa <- mean/scale_theta
    
    pars <- c(scale_theta = scale_theta,
              shape_kappa = shape_kappa)
  } else if(exists("scale") & exists("shape")) {
    
    mean <- shape * scale
    sd <- sqrt(scale^2 * shape)
    
    pars <- c(mean = mean,
              sd = sd)
    
  }
  
  return(pars)
}