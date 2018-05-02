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

J_alpha <- function(beta, theta) {

  # out <- EnvStats::rpareto(n = length(alpha), location = g, s = beta)
  # 
  # out <- sort(out, decreasing = TRUE)
  
  out <- sapply(theta, function(x) rpareto(n=1, location=x, s=beta))

  out
  
}

J_alpha_a <- function(x, beta, theta, log = TRUE) {
  
  # out <- sum(log(EnvStats::dpareto(x, location = g, s = beta)))
  
  # out <- sum(sapply(theta, function(y) log(dpareto(x, location = y, s = beta))))
  
  out <- sum(log(dpareto(x, location = theta, s = beta)))
  
  if (! log)
    out <- exp(out)
  
  out
}

# these functions are just for a single alpha_k
f_alpha_k <- function(alpha_k, theta_k, beta) {
  alpha_k ^ (-beta - 1) * theta_k ^ alpha_k
}

j_alpha <- function(prob = FALSE, n = NULL, x = NULL, rate = 5) {
  # if prob is false, draw a sample, otherwise find p(x)
  if (! prob) {
    out <- rexp(n = n, rate = rate)
  } else {
    out <- dexp(x, rate)
  }
  
  out
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

a <- 0; b <- 0 # reduces to non-informative prior for beta

g <- 0.01 # fixed parameter for pareto

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
  
  acc_alpha[j] <- sum(keep) / k
  
}

mean(acc_alpha)




