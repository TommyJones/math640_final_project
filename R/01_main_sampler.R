################################################################################
# This script runs the main model (to be called from "Math 640 Final Paper.Rmd")
################################################################################

### load libraries ----
library(textmineR)
library(gtools)
library(mcmcplots)
library(EnvStats)
library(statmod)

### set up sampler ----

main_sampler <- function(y, B, seed, theta0, alpha0, beta0, gamma) {
  
  # set up functions for M-H sampler for alpha
  
  # these functions are just for a single alpha_k
  f_alpha_k <- function(alpha_k, theta_k, beta) {
    alpha_k ^ (-beta - 1) * theta_k ^ alpha_k
  }
  
  j_alpha <- function(mu = 0.1, prob = FALSE, n = NULL, x = NULL, rate = 0.01) {
    # if prob is false, draw a sample, otherwise find p(x)
    if (! prob) {
      # out <- rexp(n = n, rate = rate)
      out <- rinvgauss(n = n, mean = mu, shape = rate)
      # out <- rpareto(n = n, location = mu, shape = rate)
    } else {
      # out <- dexp(x, rate)
      out <- dinvgauss(x, mean = mu, shape = rate)
      # out <- dpareto(x, location = mu, shape = rate)
    }
    
    out
  }
  
  # set up sampler
  k <- length(y)
  
  theta <- matrix(0, nrow = B, ncol = k)
  
  theta[ 1, ] <- theta0
  
  alpha <- matrix(0, nrow = B, ncol = k)
  
  alpha[ 1, ] <- alpha0
  
  acc_alpha <- numeric(ncol(alpha))
  
  beta <- numeric(B)
  
  beta[ 1 ] <- 2
  
  a <- 0; b <- 0 # reduces to non-informative prior for beta
  
  g <- gamma
  
  # run the sampler
  # run sampler
  for (j in 2:B) {
    
    # sample theta
    theta[j,] <- rdirichlet(1, y + alpha[j-1,])
    
    # sample beta
    beta[j] <- rgamma(1, k + a, sum(log(alpha[j-1,]) - k * log(g)) + b)
    
    # sample alpha
    alpha_star <- j_alpha(n = k)
    
    r <- (f_alpha_k(alpha_k = alpha_star, theta_k = theta[j-1,], beta = beta[j-1]) / 
            f_alpha_k(alpha[j-1,],theta[j-1,],beta[j-1])) /
      (j_alpha(prob = TRUE, x = alpha_star) / j_alpha(prob = TRUE, x = alpha[j-1,]))
    
    u <- runif(1)
    
    keep <- pmin(r,1) > u
    
    alpha[j,keep] <- alpha_star[keep]
    
    alpha[j,!keep] <- alpha[j-1,!keep]
    
    # acc_alpha[j] <- sum(keep) / k
    
    acc_alpha <- acc_alpha + keep
    
  }
  
  acc_alpha <- acc_alpha / B
  
  # return the result
  list(theta = theta, alpha = alpha, acc_alpha = acc_alpha, beta = beta,
       seed = seed, theta0 = theta0, alpha0 = alpha0)
  
}

### test to see if it works ----

# # Load data and create data-related variables/constants
# y <- colSums(nih_sample_dtm) # nih_sample_dtm from the textmineR package
# 
# y <- y[ order(y, decreasing = TRUE) ]
# 
# k <- length(y) # 20
# 
# y <- y[ 1:k ]
# 
# # run sampler
# B <- 10000
# 
# result <- main_sampler(y = y, B = B, seed = 8675309,
#                        theta0 = y / sum(y), alpha0 = rep(0.1,length(y)),
#                        beta0 = 2, gamma = 0.01)
# 
# # acceptance rate
# mean(result$acc_alpha)
# 
# summary(result$acc_alpha)
# 
# # geweke diagnostic
# g <- apply(result$alpha, 2, function(x) geweke.diag(x[ (B/4):B ])$z)
# 
# # by random chance, expect 5% to be greater than 1.96
# mean(abs(g) >= 1.96 | is.na(g))
# 
# sum(is.na(g))

