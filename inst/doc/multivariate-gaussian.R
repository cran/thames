## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----setup--------------------------------------------------------------------
set.seed(2023)
library(thames)

## -----------------------------------------------------------------------------
s0      <- 1
mu_star <- 1
n       <- 20
d       <- 1

## -----------------------------------------------------------------------------
library(mvtnorm)

Y =  rmvnorm(n, mu_star, diag(d))
sn = 1/(n+1/s0)
mn = sum(Y)/(n + 1/s0)

## -----------------------------------------------------------------------------
mu_sample = rmvnorm(2000, mean=mn,  sn*diag(d))

## -----------------------------------------------------------------------------
reg_log_prior <- function(param,sig2) {
  d <- length(param)
  p <- dmvnorm(param, rep(0,d), (sig2)*diag(d), log = TRUE)
  return(p)
}

log_prior <- apply(mu_sample, 1, reg_log_prior, s0)

## -----------------------------------------------------------------------------
reg_log_likelihood <- function(param, X) {
  n <- nrow(X)
  d <- length(param)
  sum(dmvnorm(X, param, diag(d),log = TRUE))
}

log_likelihood <- apply(mu_sample, 1, reg_log_likelihood, Y)

## -----------------------------------------------------------------------------
log_post <- log_prior + log_likelihood

## -----------------------------------------------------------------------------
result <- thames(log_post,mu_sample)

## -----------------------------------------------------------------------------
-result$log_zhat_inv

## -----------------------------------------------------------------------------
-result$log_zhat_inv_L 
-result$log_zhat_inv_U

## -----------------------------------------------------------------------------
result_90 <- thames(log_post,mu_sample, p = 0.05)
-result_90$log_zhat_inv_L 
-result_90$log_zhat_inv_U

## -----------------------------------------------------------------------------
- n*d*log(2*pi)/2 - d*log(s0*n+1)/2 - sum(Y^2)/2 + sum(colSums(Y)^2)/(2*(n+1/s0))

## -----------------------------------------------------------------------------
s0      <- 1
n       <- 20
d       <- 2
mu_star <- rep(1,d)
Y =  rmvnorm(n, mu_star, diag(d))
sn = 1/(n+1/s0)
mn = apply(Y,2,sum)/(n + 1/s0)
mu_sample = rmvnorm(2000, mean=mn,  sn*diag(d))

mvg_log_post <- function(param, X, sig2){
  n <- nrow(X)
  d <- length(param)
  l <- sum(dmvnorm(X, param, diag(d),log = TRUE))
  p <- dmvnorm(param, rep(0,d), (sig2)*diag(d), log = TRUE)
  return(p + l)
}
log_post <- apply(mu_sample, 1,mvg_log_post,Y,s0)
result <- thames(log_post,mu_sample)
-result$log_zhat_inv

## -----------------------------------------------------------------------------
- n*d*log(2*pi)/2 - d*log(s0*n+1)/2 - sum(Y^2)/2 + sum(colSums(Y)^2)/(2*(n+1/s0))

