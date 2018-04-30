library(textmineR)
library(gtools)

k <- 20

set.seed(220)

# a <- colSums(nih_sample_dtm)
# 
# a <- a[ order(a, decreasing = TRUE) ]
# 
# a <- a[ 1:k ] / 50

y <- colSums(nih_sample_dtm)

y <- y[ order(y, decreasing = TRUE) ]

y <- y[ 1:k ]

set.seed(10)
# theta <- rdirichlet(1, a)


lBeta <- function(alpha) {
  lgamma(sum(alpha)) - sum(lgamma(alpha))
}

f_alpha <- function(alpha, theta, beta, log = TRUE){
  
  # out <- Beta(alpha) * prod(theta ^ (alpha - 1) * alpha ^ -(beta + 1))
  
  out <- lBeta(alpha) + sum((alpha - 1)  * log(theta) -(beta + 1) * log(alpha))
  
  if (! log)
    out <- exp(out)
  
  out
}

f_beta <- function(beta, alpha, log = TRUE) {
  k <- length(alpha)
  
  # out <- beta ^ k * prod(alpha) ^ -(beta + 1)
  
  out <- k * log(beta) - (beta + 1) * sum(log(alpha))
  
  if (! log)
    out <- exp(out)
  
  out
}

curve(f_beta(beta = x, alpha = rep(0.01, k)), from = 0, to = 1000)

# delcare sampling functions for alpha and beta
J_alpha <- function(beta, k) {
  
  EnvStats::rpareto(n = k, location = 1, s = beta)
  
}

J_alpha_a <- function(x, beta, k, log = TRUE) {
  
  out <- sum(log(EnvStats::dpareto(x, location = 1, s = beta)))
  
  if (! log)
    out <- exp(out)
  
  out
}

J_beta <- function(beta) {
  # rgamma(1, beta, 3)
  out <- rnorm(1, mean = beta, sd = 1)
  
  if (out <= 0) {
    out <- 0.001
  }
  
  out
}

J_beta_a <- function(x, beta) {
  # dgamma(x, beta, 3)
  dnorm(x, mean = beta, sd = 1)
}

# set up sampler
B <- 20000

theta <- matrix(0, nrow = B, ncol = k)

theta[ 1, ] <- y / sum(y)

alpha <- matrix(0, nrow = B, ncol = k)

alpha[ 1, ] <- y / 50

acc_alpha <- numeric(B)

beta <- numeric(B)

beta[ 1 ] <- 2

acc_beta <- numeric(B)

# run sampler

for (j in 2:B) {
  
  # sample theta
  theta[j,] <- rdirichlet(1, y + alpha[j-1,])
  
  # sample alpha
  alpha_star <- J_alpha(beta[j-1], ncol(alpha))
  
  r1 <- (f_alpha(alpha_star, theta[j-1,], beta[j-1]) - f_alpha(alpha[j-1,], theta[j-1,], beta[j-1])) - 
    (J_alpha_a(alpha_star, beta[j-1], ncol(alpha)) - J_alpha_a(alpha[j-1,], beta[j-1], ncol(alpha)))
  
  r1 <- exp(r1)
  
  u <- runif(1)
  
  if (u < min(r1, 1)) {
    alpha[j,] <- alpha_star
    acc_alpha[j] <- 1
  } else {
    alpha[j,] <- alpha[j-1,]
  }
  
  # sample beta
  beta_star <- J_beta(beta[j-1]) # beta[j-1]
  
  # r2 <- (f_beta(beta_star, alpha = alpha[j-1,]) - f_beta(beta[j-1], alpha = alpha[j-1])) -
  #   (log(J_beta_a(beta_star,beta[j-1])) - log(J_beta_a(beta[j-1], beta[j-1])))
  
  r2 <- (f_beta(beta_star, alpha = alpha[j-1,]) - f_beta(beta[j-1], alpha = alpha[j-1]))
  
  r2 <- exp(r2)
  
  u <- runif(1)
  
  if (u < min(r2, 1)) {
    beta[j] <- beta_star
    acc_beta[j] <- 1
  } else {
    beta[j] <- beta[j-1]
  }
  
  
}

mean(acc_alpha)

mean(acc_beta)

mcmcplots::mcmcplot1(matrix(beta,ncol = 1))
