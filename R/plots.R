# load libraries
library(textmineR)
library(gtools)
library(mcmcplots)
library(EnvStats)
library(coda)

load(".RData")

# check acceptance rates
main_alpha_acc <- sapply(main_chains, function(x) x$acc_alpha)

control_alpha_acc <- sapply(control_chains, function(x) x$acc_alpha)

plot(density(rowMeans(control_alpha_acc)), lwd = 4, col = "blue", lty = 2,
     main = expression(paste("Acceptance rates of ", alpha, " parameters in both models")),
     xlab = "Acceptance rate",
     xlim = c(0.25, 0.7))
lines(density(rowMeans(main_alpha_acc)), lwd = 4, lty = 1, col = "red")
legend("topright", legend = c("Main", "Control"), lwd = 4, lty = c(1,3),
       col = c("red", "blue"))

alpha_acc_table <- data.frame(Main = quantile(rowMeans(main_alpha_acc), c(0,0.25,0.5,0.75,1)),
                              Control = quantile(rowMeans(control_alpha_acc), c(0,0.25,0.5,0.75,1)),
                              stringsAsFactors = FALSE)

rownames(alpha_acc_table) <- c("min.", "25%", "50%", "75%", "max.")

knitr::kable(alpha_acc_table, digits = 2)

# check convergence
main_beta_conv
main_alpha_conv
control_alpha_conv
main_theta_conv
control_theta_conv

conv_table <- data.frame(main_theta_conv, main_alpha_conv,
                         control_theta_conv, control_alpha_conv,
                         stringsAsFactors = FALSE)

par(mfrow = c(2,1), mar = c(2.1,4.1,2.1,2.1))
plot(rowMeans(main_posterior$alpha), type = "l",
     main = "Main model", col = "red", xaxt = "n", xlab = "",
     ylab = expression(alpha))
abline(v = c(1000,2000,3000))
plot(rowMeans(main_posterior$alpha), type = "l",
     main = "Control model", col = "blue", xaxt = "n", xlab = "",
     ylab = expression(alpha))
abline(v = c(1000,2000,3000))


# get indices of max, median, 3rd quartile words
q <- quantile(y, c(1,0.95, 0.75))

ind <- sapply(q, function(x) which(y == x)[1])

# look at ACF/mcmcplot

capture <- lapply(ind, function(k){
  a <- matrix(main_posterior$alpha[,k], ncol = 1)
  colnames(a) <- "alpha"
  mcmcplot1(a, greek = TRUE, col = "red")
  
  a <- matrix(control_posterior$alpha[,k], ncol = 1)
  colnames(a) <- "alpha"
  mcmcplot1(a, greek = TRUE)
  
})

capture <- lapply(ind, function(k){
  a <- matrix(main_posterior$theta[,k], ncol = 1)
  colnames(a) <- "theta"
  mcmcplot1(a, greek = TRUE, col = "red")
  
  a <- matrix(control_posterior$theta[,k], ncol = 1)
  colnames(a) <- "theta"
  mcmcplot1(a, greek = TRUE)
})

b <- matrix(main_posterior$beta, ncol = 1)
colnames(b) <- "beta"
mcmcplot1(b, greek = TRUE, col = "red")

# log posterior comparison
main_p <- apply(main_posterior$theta,1,function(x) dmultinom(y, prob = x, log = TRUE))

control_p <- apply(control_posterior$theta,1,function(x) dmultinom(y, prob = x, log = TRUE))

d1 <- density(main_p)
d2 <- density(control_p)

plot(d1, col = "red", lwd = 4,
     main = "Posterior Log Likelihood",
     xlim = range(c(d1$x, d2$x)))
lines(d2, col = "blue", lwd = 4, lty = 2)


# calculate DIC for each model
calc_dic <- function(y, theta_mat) {
  
  llik <- apply(theta_mat, 1, function(x) dmultinom(y, prob = x, log = TRUE))
  
  pdic <- 2 * var(llik)
  
  dic <- -2 * dmultinom(y, prob = colMeans(theta_mat)) + 2 * pdic
  
  dic
  
}

main_dic <- calc_dic(y, main_posterior$theta)

control_dic <- calc_dic(y, control_posterior$theta)

barplot(c(Main = main_dic, Control = control_dic), col = c("red", "blue"),
        density = c(-1,25))
