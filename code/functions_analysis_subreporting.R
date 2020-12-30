library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(rriskDistributions)
library(stringr)
source("./functions_auxiliary.R")


# receive as input the quartiles of the time from symptoms to outcome
# and return a function that returns the values of a gamma delay
# function with such quartiles
onset2Outcome <- function(onset2Outcomequartiles) { #
  pars <- get.gamma.par(p=c(0.25, 0.5, 0.75), q=onset2Outcomequartiles,
                        plot=FALSE, verbose=FALSE)
  delayFun <- function(x){
    delayDensity <- dgamma(x+1, pars[1], pars[2])
  }
  return(delayFun)
}

# same as above, but for the cumulative delay function
onset2Outcome_cumulative <- function(onset2OutcomeQuartiles) { #
  pars <- get.gamma.par(p=c(0.25, 0.5, 0.75), q=onset2OutcomeQuartiles,
                        plot=FALSE, verbose=FALSE)
  delayFun <- function(x){
    delayDensity <- pgamma(x+1, pars[1], pars[2])
  }
  return(delayFun)
}


# Calculates how many of past cases are expected to be resolved by
# day. Whether it calculates a cummulative or instantaneous
# resolving of cases depends on whether delay_fun is cummulative
# or instantaneous
# sum C^t * F_D(T-t)
cases_with_known_outcome <- function(newCases, delay_fun){
  pKnownOutcome <- rev(delay_fun(c(0:(length(newCases)-1))))
  casesKnownOutcome <- sum(pKnownOutcome * newCases)
  return(casesKnownOutcome)
}


# calculate poisson confint for each element in a vector
poisson_confint <- function(countVec) {
  upper <- NULL
  lower <- NULL
  for (l in c(1:length(countVec))) {
    poissonTest <- poisson.test(countVec[l])$conf.int
    lower <- c(lower, poissonTest[1])
    upper <- c(upper, poissonTest[2])
  }
  output <- list(lower = lower, upper = upper)
  return(output)
}


# Estimate porportion of cases resulting in a given outcome,
# corrected for delay. Returns confidence intervals computed
# from a poisson outcome model
proportion_outcome_delay <- function(newCases, newOutcome,
                                     delay_fun, cumulative=FALSE){
  if (cumulative) {
    outcomeIncidence <- cumsum(newOutcome)
  } else {
    outcomeIncidence <- newOutcome
  }
  cumulative_known_t <- NULL # cumulative cases with known outcome at time tt
  # Sum over cases up to time tt
  for(ii in 1:length(newCases)){
    # number of cases with known outcome at time ii
    known_i <- cases_with_known_outcome(newCases[1:ii], delay_fun) 
    cumulative_known_t <- c(cumulative_known_t,known_i) # Tally cumulative known
  }
  confintOutcome <- poisson_confint(outcomeIncidence)
  proportionOutcome <- outcomeIncidence/cumulative_known_t
  proportionOutcomeL <- confintOutcome$lower/cumulative_known_t
  proportionOutcomeU <- confintOutcome$upper/cumulative_known_t
  proportionDf <- data.frame(proportionOutcome = proportionOutcome,
                             proportionOutcomeL = proportionOutcomeL,
                             proportionOutcomeU = proportionOutcomeU)
  return(proportionDf)
}


# Function to estimate % of a given outcome in infected in data
# corrected by delay. It can calculate the instantaneous %
# or the cumulative %. For that, delay_fun must be appropiately
# instantaneous or cumulative, and the cumulative flag must
# be turned correspondingly.
scale_cfr_temporal <- function(newCases, newOutcome,
                               delay_fun, cumulative=FALSE){
  if (cumulative) {
    outcomeIncidence <- cumsum(newOutcome)
  } else {
    outcomeIncidence <- newOutcome
  }
  casesKnownT <- NULL # cases with known outcome at time tt
  # Sum over cases up to time tt
  for(ii in 1:length(newCases)){
    known_i <- cases_with_known_outcome(newCases[1:ii], delay_fun) # number of cases with known outcome at time ii
    casesKnownT <- c(casesKnownT, known_i) # Tally cumulative known
  }
  # naive outcome proportion
  b_tt <- sum(outcomeIncidence)/sum(newCases) 
  # corrected outcome proportion
  p_tt <- (outcomeIncidence/casesKnownT) %>% pmin(.,1)
  
#  cfrDf <- data.frame(nCFR = b_tt, cCFR = p_tt,
#                      total_outcomes = sum(outcomeIncidence),
#                      outcomes = outcomeIncidence,
#                      cum_known_t = round(casesKnownT),
#                      total_cases = sum(newCases))
  cfrDf <- data.frame(naiveOutcomeProp = b_tt,
                      correctedOutcomeProp = p_tt,
                      totalOutcomes = sum(outcomeIncidence),
                      outcomes = outcomeIncidence,
                      casesKnownT = round(casesKnownT),
                      total_cases = sum(newCases))
  return(cfrDf)
}


# Function to estimate subreporting from a given baseline outcome
# and the observed outcomes
estimate_subreporting_cum <- function(newCases, delay_fun, newOutcome,
                                      baselinePercentOutcome,
                                      confLevel = 0.95){
  propOutcome <- proportion_outcome_delay(newCases=newCases,
                                          newOutcome=newOutcome,
                                          delay_fun=delay_fun,
                                          cumulative=TRUE)
  baselineProp <- baselinePercentOutcome/100
  # compare corrected ICUR with baseline ICU rate
  underreportingEstimate <- baselineProp/propOutcome$proportionOutcome
  underreportingLower <- baselineProp/propOutcome$proportionOutcomeU
  underreportingUpper <- baselineProp/propOutcome$proportionOutcomeL
  # cap values at 1
  underreportingEstimate[which(underreportingEstimate>1)] <- 1
  underreportingLower[which(underreportingLower>1)] <- 1
  underreportingUpper[which(underreportingUpper>1)] <- 1
  # put together output df
  underestimateReport <- data.frame(cumCases = cumsum(newCases),
                              estimate = underreportingEstimate,
                              lower = underreportingLower,
                              upper = underreportingUpper,
                              cumOutcome = cumsum(newOutcome))
  return(underestimateReport)
}


# calculate the expected number of outcomes, and the corrected
# outcome proportion in time 
expected_proportion_in_time <- function(individualCaseDf, proportionOutcomeDf,
                                         delay_fun, maxDate=NA) {

  individualCaseDf$outcomeProb <- assign_bin_prop(individualCaseDf$Edad,
                                        proportionOutcomeDf$age,
                                        proportionOutcomeDf$outcome)
  meanOutcome <- mean(individualCaseDf$outcomeProb, na.rm = TRUE)
  individualCaseDf$outcomeProb[which(is.na(individualCaseDf$outcomeProb))] <- meanOutcome
  # sum all proportions of mortality for a given day
  dailyProb <- group_by(individualCaseDf, FechaInicioSintomas) %>%
    summarize(., outcomeSum = sum(outcomeProb), cases = n()) %>%
    dplyr::rename(., day = FechaInicioSintomas)
  # Put in a vector including all dates
  minDate <- min(date(dailyProb$day))
  if (is.na(maxDate)) {
    maxDate <- date(max(individualCaseDf$FechaResultado))
  }
  nDays <- lubridate::interval(minDate, maxDate) / lubridate::days(1)
  daysVec <- seq(from=0, to=nDays) + lubridate::ymd(minDate) 
  #allDays <- data.frame(day=as.character(daysVec))
  allDays <- data.frame(day=daysVec)
  # put together in a df with the dynamics
  outcomeSumDyn <- merge(dailyProb, allDays, all.y=TRUE) %>%
    dplyr::mutate(., outcomeInTodaysCases = ifelse(is.na(outcomeSum), 0, outcomeSum),
                  cases = ifelse(is.na(cases), 0, cases)) %>%
    dplyr::select(., -outcomeSum) %>%
    as_tibble(.)
  # calculate the N outcomes to be expected on each day accounting for the
  # delay function
  expectedOutcome <- NULL
  expectedResolvedCases <- NULL
  for (nr in c(1:nrow(outcomeSumDyn))) {
    delayFunT <- rev(delay_fun(c(0:(nr-1))))
    expectedOutcome[nr] <- sum(outcomeSumDyn$outcomeInTodaysCases[1:nr] *
                               delayFunT)
    expectedResolvedCases[nr] <- sum(outcomeSumDyn$cases[1:nr] *
                                     delayFunT)
  }
  correctedProbOutcome <- expectedOutcome/expectedResolvedCases
  outcomeSumDyn$expectedOutcomes <- expectedOutcome
  outcomeSumDyn$cases_known <- expectedResolvedCases
  outcomeSumDyn$correctedProbOutcome <- correctedProbOutcome
  return(outcomeSumDyn)
}


get_fitting_data <- function(outcomeDf, delay_fun, baselineOutcomeProp) {
  outcome_threshold_date <- outcomeDf %>% 
    dplyr::mutate(outcome_cum_sum = cumsum(newOutcome)) %>% 
    dplyr::filter(outcome_cum_sum >= 10) %>% 
    dplyr::pull(day) %>% 
    min()
  #return adjusted date and reporting_estimate
  fittingData <- scale_cfr_temporal(outcomeDf$newCases,
                                  outcomeDf$newOutcome,
                                  delay_fun=delay_fun,
                                  cumulative=FALSE) %>% 
    dplyr::as_tibble() %>% 
    dplyr::mutate(reporting_estimate = baselineOutcomeProp/correctedOutcomeProp) %>% 
    dplyr::mutate(reporting_estimate = pmin(reporting_estimate, 1),
           date = outcomeDf$day,
           date_num = as.numeric(outcomeDf$day),
           #date_num = date_num - min(date_num) + 1,
           date_num = date_num,
           newOutcome = outcomeDf$newOutcome,
           cases_known = casesKnownT) %>% 
    dplyr::filter(date >= outcome_threshold_date) %>% 
    dplyr::select(date, date_num, reporting_estimate, newOutcome, cases_known)
  return(fittingData)
}

get_fitting_data_strat <- function(outcomeDf, delay_fun) {
  casesKnownDf <- data.frame(row1 = numeric(nrow(outcomeDf)))
  for (strat in colnames(outcomeDf)) {
    newCasesStrat <- outcomeDf[[strat]]
    casesKnownT <- NULL
    for(ii in 1:length(newCasesStrat)){
      known_i <- cases_with_known_outcome(newCasesStrat[1:ii], delay_fun) # number of cases with known outcome at time ii
      casesKnownT <- c(casesKnownT, known_i)
    }
    stratDf <- data.frame(noName = casesKnownT)
    colnames(stratDf) <- strat
    casesKnownDf <- cbind(casesKnownDf, stratDf)
  }
  casesKnownDf <- dplyr::select(casesKnownDf, -1)
  return(casesKnownDf)
}


estimate_icu_occupancy <- function(previousICU, newICU, delay_fun) {
  indicatorVec <- c(rep("known", length(previousICU)),
  rep("predicted", length(newICU)))
  allNewICU <- c(previousICU, newICU)
  cumulativeICU <- cumsum(allNewICU)
  cumulativeRecoveries <- NULL
  for (day in c(1:length(cumulativeICU))) {
    recoveryDist <- rev(delay_fun(c(1:day)-1))
    cumulativeRecoveries[day] <- sum(allNewICU[1:day] * recoveryDist)
  }
  icuOccupancy <- cumulativeICU - cumulativeRecoveries
  return(icuOccupancy)
}
  
