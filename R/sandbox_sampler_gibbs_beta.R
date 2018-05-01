library(textmineR)
library(gtools)
library(mcmcplots)
library(EnvStats)

set.seed(90210)

y <- colSums(nih_sample_dtm)

y <- y[ order(y, decreasing = TRUE) ]

k <- length(y) # 20

y <- y[ 1:k ]

# declare functions to help with MH step for alpha

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

J_alpha <- function(beta, alpha, gamma) {
  k <- length(alpha)
  
  out <- EnvStats::rpareto(n = k, location = gamma, s = beta)
  
  out <- sort(out, decreasing = TRUE)
  # out <- out[ order(alpha) ]  
  
  out
  
}

J_alpha_a <- function(x, beta, k, gamma, log = TRUE) {
  
  out <- sum(log(EnvStats::dpareto(x, location = gamma, s = beta)))
  
  if (! log)
    out <- exp(out)
  
  out
}


# set up sampler
B <- 10000

theta <- matrix(0, nrow = B, ncol = k)

theta[ 1, ] <- y / sum(y)

alpha <- matrix(0, nrow = B, ncol = k)

alpha[ 1, ] <- y / 50

acc_alpha <- numeric(B)

beta <- numeric(B)

beta[ 1 ] <- 2

a <- 0; b <- 0 # reduces to non-informative prior for beta

gamma <- 0.01 # fixed parameter for pareto

# run sampler
for (j in 2:B) {
  
  # sample theta
  theta[j,] <- rdirichlet(1, y + alpha[j-1,])
  
  # sample beta
  beta[j] <- rgamma(k + a, sum(log(alpha[j-1] - k * log(gamma))) + b)
  
  # sample alpha
  alpha_star <- J_alpha(beta[j-1], ncol(alpha), gamma)
  
  r1 <- (f_alpha(alpha_star, theta[j-1,], beta[j-1]) - 
           f_alpha(alpha[j-1,], theta[j-1,], beta[j-1])) - 
    (J_alpha_a(alpha_star, beta[j-1], ncol(alpha), gamma) - 
       J_alpha_a(alpha[j-1,], beta[j-1], ncol(alpha), gamma))
  
  # r1 <- exp(r1)
  
  u <- runif(1)
  
  if (log(u) < min(r1, 0)) {
    alpha[j,] <- alpha_star
    acc_alpha[j] <- 1
  } else {
    alpha[j,] <- alpha[j-1,]
  }
  
}

mean(acc_alpha)

sum(acc_alpha)


