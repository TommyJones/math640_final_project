library(textmineR)
library(gtools)

set.seed(90210)

y <- colSums(nih_sample_dtm)

y <- y[ order(y, decreasing = TRUE) ]

k <- length(y) # 20

y <- y[ 1:k ]

B <- 20000

theta <- matrix(0, nrow = B, ncol = k)

theta[ 1, ] <- y / sum(y)

alpha <- rep(0.1, k)

for (j in 2:B) {
  
  theta[j, ] <- rdirichlet(n = 1, alpha = alpha + y)
  
}



