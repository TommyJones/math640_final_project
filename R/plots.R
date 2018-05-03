# load libraries
library(textmineR)
library(gtools)
library(mcmcplots)
library(EnvStats)
library(coda)

load(".RData")

# combine posteriors
control_posterior <- list(theta = do.call(rbind, lapply(control_chains, function(x) x$theta)),
                          alpha = do.call(rbind, lapply(control_chains, function(x) x$alpha)))

main_posterior <- list(theta = do.call(rbind, lapply(main_chains, function(x) x$theta)),
                       alpha = do.call(rbind, lapply(main_chains, function(x) x$alpha)),
                       beta = do.call(c, lapply(main_chains, function(x) x$beta)))


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

conv <- parallel::mclapply(list(main = main_posterior, control = control_posterior), 
                           function(x){
                             apply(x$theta,1,function(z) dmultinom(y, prob = z, log = T))
                           })

par(mfrow = c(2,1), mar = c(2.1,4.1,2.1,2.1))
plot(conv$main, type = "l",
     main = "Main model", col = "red", xaxt = "n", xlab = "",
     ylab = expression(paste("log[P(y|",theta,")]")))
abline(v = c(1000,2000,3000))
plot(conv$control, type = "l",
     main = "Control model", col = "blue", xaxt = "n", xlab = "",
     ylab = expression(paste("log[P(y|",theta,")]")))
abline(v = c(1000,2000,3000))

g <- lapply(conv, geweke.diag)

# conv_table <- data.frame(main_theta_conv, main_alpha_conv,
#                          control_theta_conv, control_alpha_conv,
#                          stringsAsFactors = FALSE)
# 
# par(mfrow = c(2,1), mar = c(2.1,4.1,2.1,2.1))
# plot(rowMeans(main_posterior$alpha), type = "l",
#      main = "Main model", col = "red", xaxt = "n", xlab = "",
#      ylab = expression(alpha))
# abline(v = c(1000,2000,3000))
# plot(rowMeans(main_posterior$alpha), type = "l",
#      main = "Control model", col = "blue", xaxt = "n", xlab = "",
#      ylab = expression(alpha))
# abline(v = c(1000,2000,3000))


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
calc_dic <- function(y, theta_mat, B = NULL) {
  
  llik <- apply(theta_mat, 1, function(x) dmultinom(y, prob = x, log = TRUE))
  
  pdic <- 2 * var(llik)
  
  
  if (! is.null(B)) {
    dic <- sapply(seq_len(B), function(th){
      -2 * dmultinom(y, prob = theta_mat[ sample(seq_len(nrow(theta_mat)), 1) , ], log = TRUE)
    })
    
    dic <- dic + 2 * pdic
  } else {
    dic <- -2 * dmultinom(y, prob = colMeans(theta_mat), log = TRUE) + 2 * pdic
  }
  
  dic
  
}

main_dic <- calc_dic(y, main_posterior$theta, B = 10000)

control_dic <- calc_dic(y, control_posterior$theta, B = 10000)

d1 <- density(main_dic)

d2 <- density(control_dic)

plot(d1, col = "red", lwd = 4,
     main = "Posterior DIC",
     xlim = range(c(d1$x, d2$x)))
lines(d2, col = "blue", lwd = 4, lty = 2)
abline(v = quantile(main_dic, probs = c(0.025, 0.5, 0.975)), col = "red", lwd = 2)
abline(v = quantile(control_dic, probs = c(0.025, 0.5, 0.975)), col = "blue", lwd = 2, lty = 2)


# barplot(c(Main = main_dic, Control = control_dic), col = c("red", "blue"),
#         density = c(-1,25))

# log-log plots of alpha and theta vs. y

# alpha
plot(log10(seq_along(y)), log10(y), 
     lwd = 3, yaxt = "n", ylab = "log10(x)", xlab = "rank(x)", type = "l",
     lty = 3, main = expression(paste("Comparing ", alpha)))
par(new = TRUE)
plot(log10(seq_along(y)), log10(colMeans(main_posterior$alpha)), 
     lwd = 3, yaxt = "n", ylab = "", xaxt = "n", xlab = "", type = "l", lty = 1, 
     col = rgb(1,0,0,0.5))
par(new = TRUE)
plot(log10(seq_along(y)), log10(colMeans(control_posterior$alpha)), 
     lwd = 3, yaxt = "n", ylab = "", xaxt = "n", xlab = "", type = "l", lty = 2, 
     col = rgb(0,0,1,0.5))
legend("bottomleft", legend = c("y", "Main", "Control"),
       lwd = 3, lty = c(3,1,2), col = c("black", "red", "blue"))

# theta
plot(log10(seq_along(y)), log10(y), 
     lwd = 3, yaxt = "n", ylab = "log10(x)", xlab = "rank(x)", type = "l",
     lty = 3, main = expression(paste("Comparing ", theta)))
par(new = TRUE)
plot(log10(seq_along(y)), log10(colMeans(main_posterior$theta)), 
     lwd = 3, yaxt = "n", ylab = "", xaxt = "n", xlab = "", type = "l", lty = 1, 
     col = rgb(1,0,0,0.5))
par(new = TRUE)
plot(log10(seq_along(y)), log10(colMeans(control_posterior$theta)), 
     lwd = 3, yaxt = "n", ylab = "", xaxt = "n", xlab = "", type = "l", lty = 2, 
     col = rgb(0,0,1,0.5))
legend("bottomleft", legend = c("y", "Main", "Control"),
       lwd = 3, lty = c(3,1,2), col = c("black", "red", "blue"))


