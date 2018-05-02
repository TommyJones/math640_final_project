################################################################################
# This script runs the control model (to be copied to "Math 640 Final Paper.Rmd")
################################################################################

### load libraries ----
library(textmineR)
library(gtools)
library(mcmcplots)
library(EnvStats)

### Load data and create data-related variables/constants ----
y <- colSums(nih_sample_dtm) # nih_sample_dtm from the textmineR package

y <- y[ order(y, decreasing = TRUE) ]

k <- length(y) # 20

y <- y[ 1:k ]

### set up sampler ----

# declare a function for the control sampler
control_sampler <- function(y, B, seed, theta0, alpha0) {
  
  # set up sampling functions for M-H
  f_alpha_k <- function(alpha_k, theta_k) {
    theta_k ^ alpha_k
  }
  
  j_alpha <- function(prob = FALSE, n = NULL, x = NULL, rate = 2) {
    # if prob is false, draw a sample, otherwise find p(x)
    if (! prob) {
      out <- rexp(n = n, rate = rate)
    } else {
      out <- dexp(x, rate)
    }
    
    out
  }
  
  # set up constants
  k <- length(y)
  
  theta <- matrix(0, nrow = B, ncol = k)
  
  theta[1,] <- theta0
  
  alpha <- matrix(0, nrow = B, ncol = k)
  
  alpha[1,] <- alpha0
  
  acc_alpha <- numeric(B)
  
  # run the sampler
  for (j in 2:B) {
    # sample theta
    theta[j,] <- rdirichlet(1, y + alpha[j-1,])
    
    # sample alpha
    alpha_star <- j_alpha(n = k)
    
    r <- (f_alpha_k(alpha_k = alpha_star, theta_k = theta[j-1,]) / f_alpha_k(alpha[j-1,],theta[j-1,])) /
      (j_alpha(prob = TRUE, x = alpha_star) / j_alpha(prob = TRUE, x = alpha[j-1,]))
    
    u <- runif(1)
    
    keep <- pmin(r,1) > u
    
    alpha[j,keep] <- alpha_star[keep]
    
    alpha[j,!keep] <- alpha[j-1,!keep]
    
    acc_alpha[j] <- sum(keep) / k
    
  }
  
  # return the result
  list(theta = theta, alpha = alpha, acc_alpha = acc_alpha, seed = seed, 
       theta0 = theta0, alpha0 = alpha0)
}

### test to see if it works ----
B <- 10000

result <- control_sampler(y = y, B = B, seed = 90210, 
                          theta0 = y / sum(y), alpha0 = rep(0.1,length(y)))

# acceptance rate
mean(result$acc_alpha)

# geweke diagnostic
g <- apply(result$alpha, 2, function(x) geweke.diag(x[ (B/4):B ])$z)

# by random chance, expect 5% to be greater than 1.96
mean(g >= 1.96)

