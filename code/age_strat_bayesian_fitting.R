## This code is a modified version of this project
# https://github.com/thimotei/covid_underreporting

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
#library(RColorBrewer)
library(rriskDistributions)
library(stringr)
library(greta)
library(greta.gp)

# fit a Bayesian GP regression model to critical and known cases
# data, integrating over uncertainty in the true CFR,
# and return the posterior mean and 95% CI
age_stratified_bayesian_model <- function (casesKnownStrat,
                                           outcomeVec,
                                           slopeMean = 0.0826,
                                           slopeCI = c(0.0754,0.0897),
                                           interceptMean = -4.520,
                                           interceptCI = c(-4.894, -4.146),
                                           n_inducing = 5, verbose = TRUE) {

  data <- casesKnownStrat
  cases_known <- dplyr::select(data, -date_num, -date) %>%
    rowSums(as.matrix(.))

  # only fit to time points where there are known cases
  data <- data[which(cases_known>0),]
  n <- nrow(data)
  #times <- seq(min(data$date_num), max(data$date_num))
  times <- as.integer(data$date_num)
  dates <- data$date

  data <- dplyr::select(data, -date_num, -date) %>%
    as.matrix(.)

  # GP parameters for squared-exponential kernel plus a bias term (intercept)
  # for reporting rate
  lengthscale <- greta::lognormal(4, 0.5)
  sigma <- greta::lognormal(-1, 1)
  temporal <- greta.gp::rbf(lengthscales = lengthscale,
                            variance = sigma ^ 2)
  intercept <- greta.gp::bias(1)
  reporting_kernel <- temporal + intercept
 
  # IID noise kernel for observation overdispersion (clumped death reports)
  sigma_obs <- greta::normal(0, 0.5, truncation = c(0, Inf))
  observation_kernel <- greta.gp::white(sigma_obs ^ 2)

  # combined kernel (marginalises a bunch of parameters for easier sampling)
  kernel <- reporting_kernel + observation_kernel
 
  # a set of inducing points at which to estimate the GPs; using a subset of
  # regressors approximation (put an inducing point at the last time, since we
  # care a lot about that estimate)
  #inducing_points <- seq(min(times), max(times), length.out = n_inducing + 1)[-1]
  inducing_inds <- round(seq(1, length(times), length.out = n_inducing + 1)[-1])
  inducing_points <- times[inducing_inds]
 
  # GP for the (probit-) reporting rate
  z <- greta.gp::gp(times, inducing = inducing_points, kernel)
 
  # convert to probabilities
  reporting_rate <- greta::iprobit(z)

  
  # distribution over plausible baseline CFR values from China study. The 95%
  # CIs are symmetric around the estimate, so we assume it's an approximately
  # Gaussian distribution, truncated to allowable values.
  slopeSD <- mean(abs(slopeCI - slopeMean)) / 1.96
  outcomeSlope <- greta::normal(slopeMean, slopeSD)
  interceptSD <- mean(abs(interceptCI - interceptMean)) / 1.96
  outcomeIntercept <- greta::normal(interceptMean, interceptSD)

  # get mean age vec
  ageBins <- colnames(data)
  meanAges <- mid_bin_age(ageBins)
  probOutcomeStrat <- exp(meanAges * outcomeSlope + outcomeIntercept)/100

  # compute the expected number of critical at each timepoint, given the true CFR,
  # number of reported cases with known outcomes, and reporting rate
  log_expected_outcome <-
    log(data %*% probOutcomeStrat) - log(reporting_rate) 
  expected_outcome <- exp(log_expected_outcome)

  # define sampling distribution
  greta::distribution(outcomeVec) <- greta::poisson(expected_outcome)
  
  # construct the model
  m <- greta::model(reporting_rate)
  
  n_chains <- 500
  
  # sample initial values for hyperparameters from within their priors
  inits <- replicate(
    n_chains,
    greta::initials(
      lengthscale = rlnorm(1, 4, 0.5),
      sigma = abs(rnorm(1, 0, 0.5)),
      sigma_obs = abs(rnorm(1, 0, 0.5)),
      outcomeSlope = rnorm(1, slopeMean, slopeSD),
      outcomeIntercept = rnorm(1, interceptMean, interceptSD)
    ),
    simplify = FALSE
  )
  
#  if (verbose) {
#    country <- data$country[1]
#    message("running model for ", country)
#  }
  

  # draw a bunch of mcmc samples
  draws <- greta::mcmc(
    m,
    sampler = greta::hmc(Lmin = 15, Lmax = 20),
    chains = n_chains,
    warmup = 1000,
    n_samples = 1000,
    initial_values = inits,
    one_by_one = TRUE,
    verbose = verbose
  ) 
  

  # extend the number of chains until convergence (or give up)
  for (i in 1:5) {
    r_hats <- coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]
    n_eff <- coda::effectiveSize(draws)
    decent_samples <- max(r_hats) <= 1.1 & min(n_eff) > 1000 
    if (!decent_samples) {
      if (verbose) {
        message("maximum R-hat: ", max(r_hats),
                "\nminimum n-eff: ", min(n_eff))
      }
      draws <- greta::extra_samples(draws, 2000, one_by_one = TRUE, verbose = verbose)
    }
  }
  
  # predict without IID noise (true reporting rate, without clumped death reporting)
  # could predict to more granular times here too
  z_smooth <- greta.gp::project(z, times, kernel = reporting_kernel)
  reporting_rate_smooth <- greta::iprobit(z_smooth)
  draws_pred <- greta::calculate(reporting_rate_smooth, values = draws)
  
  # get estimates
  draws_pred_mat <- as.matrix(draws_pred)
  
  # compute posterior mean and 95% credible interval and return
  subreportingMat <- tibble::tibble(
    date = dates,
    estimate = colMeans(draws_pred_mat),
    lower = apply(draws_pred_mat, 2, quantile, 0.025),
    upper = apply(draws_pred_mat, 2, quantile, 0.975)
  )
  rm(draws)
  
  return(subreportingMat)
}

