library(textmineR)
library(gtools)

k <- 20

set.seed(220)

y <- colSums(nih_sample_dtm)

y <- y[ order(y, decreasing = TRUE) ]

y <- y[ 1:k ]

set.seed(10)
# theta <- rdirichlet(1, a)


lBeta <- function(alpha) {
  lgamma(sum(alpha)) - sum(lgamma(alpha))
}

f_alpha <- function(alpha, theta, log = TRUE){
  
  # out <- Beta(alpha) * prod(theta ^ (alpha - 1) * alpha ^ -(beta + 1))
  
  out <- lBeta(alpha) + sum((alpha - 1)  * log(theta))
  
  if (! log)
    out <- exp(out)
  
  out
}

# delcare sampling functions for alpha 
J_alpha <- function(beta, k) {
  
  EnvStats::rpareto(n = k, location = 1, s = beta)
  
}

J_alpha_a <- function(x, beta, k, log = TRUE) {
  
  out <- sum(log(EnvStats::dpareto(x, location = 1, s = beta)))
  
  if (! log)
    out <- exp(out)
  
  out
}

# set up sampler
B <- 20000

theta <- matrix(0, nrow = B, ncol = k)

theta[ 1, ] <- y / sum(y)

alpha <- matrix(0, nrow = B, ncol = k)

alpha[ 1, ] <- y / 50

acc_alpha <- numeric(B)

k <- length(y)

b <- 2.5

# run sampler

for (j in 2:B) {
  
  # sample theta
  theta[j,] <- rdirichlet(1, y + alpha[j-1,])
  
  # sample alpha
  alpha_star <- J_alpha(beta = b, k = k)
  
  r1 <- (f_alpha(alpha_star, theta[j-1,]) - f_alpha(alpha[j-1,], theta[j-1,])) - 
    (J_alpha_a(x = alpha_star, beta = b, k = k) - J_alpha_a(x = alpha[j-1,], beta = b, k = k))
  
  r1 <- exp(r1)
  
  u <- runif(1)
  
  if (u < min(r1, 1)) {
    alpha[j,] <- alpha_star
    acc_alpha[j] <- 1
  } else {
    alpha[j,] <- alpha[j-1,]
  }
  
}

mean(acc_alpha)

mcmcplots::mcmcplot1(matrix(alpha[,1], ncol = 1))

