# Simulate disease dynamics with underreporting and
# test the estimation of the underreporting in time

library(dplyr)
library(tidyr)
library(lubridate)
library(rriskDistributions)
library(greta)
library(greta.gp)
set.seed(2691)

source("./functions_auxiliary.R")
source("./functions_analysis_subreporting.R")
source("./bayesian_subreporting_fitting.R")

###########
# Simulation parameters
###########

onset2ICUquartiles <- c(6, 10, 12) # datos para ajustar curva de delay
percentageOutcome <- 3 # percentage of cases going critical 
percentageOutcomeRange <- c(2, 4) #uncertainty range to use in the fitting
simulationPoints <- 100
nRep <- 2
functionRun <- 2
f <- functionRun

# make the functions of underlying infectious dynamics and 
# underlying reporting
t <- seq(0, simulationPoints-1)/(simulationPoints-1)

# infectivity functions
funcInf <- list()
funcInf[[1]] <- 15 + 10*t - 5*t^2 + 6*t^3
funcInf[[2]] <- pmax(sin(t*2*pi)*30 + 25, 0)
funcInf[[3]] <- exp(-(t-0.5)^2/0.08^2)*120 + 10
funcInf[[4]] <- 10 + exp((2.5*t)^2)
funcInf[[5]] <- pmax(sin(t*15*pi)*10 + 25, 0)
# subreporting functions
funcSR <- list()
funcSR[[1]] <- function(inf){out = rep(0.8, length(inf))}
funcSR[[2]] <- function(inf){out = 1 - (inf / (inf+150))}
funcSR[[3]] <- function(inf){out = 0.3 +
  0.55/(1+exp(-0.3+4*-c(1:length(inf)/length(inf))))}

##########
# Define functions
##########
 
# simulate a stochastic sample of an underlying infectious
# dynamic given in infectionFun, of the daily reporting
# process, and of the appearence of critical patients
simulate_reporting_data <- function(delay_fun, infectionFun,
                                    reportingFun, percentageOutcome) {
  newCasesReal <- NULL
  newCasesReported <- NULL
  resolved <- NULL
  newOutcome <- NULL
  reportVals <- reportingFun(infectionFun)
  for (t in c(1:length(infectionFun))) {
    # sample new infected and reported per day
    newCasesReal[t] <- rpois(1, infectionFun[t])
    newCasesReported[t] <- rbinom(1, newCasesReal[t], reportVals[t])
    resolved[t] <- 0
    for (tt in c(0:(t-1))) {
      known_tt <- newCasesReal[t-tt] * delay_fun(tt)
      resolved[t] <- resolved[t] + known_tt
    }
    newOutcome[t] <- rbinom(1, round(resolved[t]),
                             percentageOutcome/100)
  }
  simulationDf <- data.frame(day = c(1:length(infectionFun)),
                             newCasesReal = newCasesReal,
                             newCasesReported = newCasesReported,
                             newOutcome = newOutcome,
                             infectionMean = infectionFun,
                             realReporting = reportVals) 
  return(simulationDf)
}


############
# Run simulation and analysis
############
delay_fun <- onset2Outcome(onset2ICUquartiles)
# make data offset to move predictions date back in time
xx <- c(0:30)
dateOffset <- round(sum(xx * delay_fun(xx)))
fullDf <- NULL

#for (f in c(1:length(funcInf))) {
  for (s in c(1:length(funcSR))) {
    for (r in c(1:nRep)) {

      # simulate the disease and reporting dynamics
      simulationDf <- simulate_reporting_data(delay_fun=delay_fun,
                          infectionFun=funcInf[[f]],
                          reportingFu=funcSR[[s]],
                          percentageOutcome=percentageOutcome) %>%
        as_tibble(.)

      # make a dataframe with the same structure as that used in the
      # actual data
      simDates <- simulationDf$day - max(simulationDf$day) + lubridate::today()
      simulationDf$day <- simDates
      criticalDf <- data.frame(day=simDates,
                               newCases=simulationDf$newCasesReported,
                               newOutcome=simulationDf$newOutcome) %>%
        as_tibble(.)

      # fit underreporting model to simulated data
      fittingData <- get_fitting_data(outcomeDf = criticalDf,
                                      delay_fun = delay_fun,
                                      baselineOutcomeProp = percentageOutcome/100)

      prediction <- run_bayesian_model(fittingData,
                                       percentageOutcome=percentageOutcome,
                                       percentageOutcomeRange=percentageOutcomeRange)

      predictionOffset <- dplyr::mutate(prediction, date = date - dateOffset)

      croppedSim <- dplyr::filter(simulationDf, day %in% predictionOffset$date)

      repDf <- cbind(croppedSim, predictionOffset) %>%
        dplyr::mutate(., funcInfInd = f, funcSRInd = s, nRep = r)

      fullDf <- rbind(fullDf, repDf)

      write.csv(fullDf, file = paste("./fit_test_simulations",
                                     as.character(f), ".csv", sep=""),
                sep = ",", row.names = FALSE)
      saveRDS(fullDf, file = paste("./fit_test_simulations",
              as.character(f), ".RDS", sep=""))
    }
  }
#}

# put together many simulations
allSims <- NULL
for (i in c(1:5)) {
  fileName <- paste("../datos_procesados/simulations_parts/fit_test_simulations",
                    i, ".RDS", sep="")
  tempSim <- readRDS(fileName)
  allSims <- rbind(allSims, tempSim)
}

saveRDS(allSims, "../datos_procesados/5_simulations.RDS")

allSims <- NULL
for (i in c(1:5)) {
  fileName <- paste("../datos_procesados/simulations_parts/5_fit_test_simulations",
                    i, "_shorter.RDS", sep="")
  tempSim <- readRDS(fileName)
  allSims <- rbind(allSims, tempSim)
}

saveRDS(allSims, "../datos_procesados/5_simulations_shorter.RDS")

