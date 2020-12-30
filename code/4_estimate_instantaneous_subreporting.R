library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(rriskDistributions)
library(stringr)
library(greta)
library(greta.gp)
source("./functions_auxiliary.R")
source("./functions_analysis_subreporting.R")
source("./bayesian_subreporting_fitting.R")

dataFile <- "../datos_procesados/1_outcomeDynamics.RDS"
outcomeDf <- readRDS(dataFile) %>%
  as_tibble(.)


# Get the parameters for the simulation
criticalUru <- read.csv("../datos_procesados/2_critical_proportion_uru.csv")
severeUru <- read.csv("../datos_procesados/2_severe_proportion_uru.csv")

criticalUru <- dplyr::filter(criticalUru, study %in% c("Salje", "Driscoll", "Brazeau",
                                                    "Levin"))
severeUru <- dplyr::filter(severeUru, study %in% c("Salje", "Driscoll", "Brazeau",
                                                    "Levin", "Verity"))
criticalMean <- mean(criticalUru$critical_cases)
criticalIC <- dplyr::filter(criticalUru, !(study %in% c("Salje", "Brazeau"))) %>%
  with(., c(min(criticalL_cases), max(criticalH_cases)))

severeMean <- mean(severeUru$severe_cases)
severeIC <- dplyr::filter(severeUru, study %in% c("Brazeau", "Levin")) %>%
  with(., c(min(severeL_cases), max(severeH_cases)))

onset2ICUquartiles <- c(4, 7, 10) # datos para ajustar curva de delay
onset2Hospquartiles <- c(3, 5, 9) # datos para ajustar curva de delay
fechaAnalisis <- Sys.Date()

# Fit to critical cases
delay_fun_crit <- onset2Outcome(onset2ICUquartiles)
outcomeDf$newOutcome <- outcomeDf$newCritical
#
fittingData <- get_fitting_data(outcomeDf, delay_fun_crit,
                                baselineOutcomeProp = criticalMean)
predictionCritical <- run_bayesian_model(fittingData,
                                 percentageOutcome=criticalMean,
                                 percentageOutcomeRange=criticalIC)

# Fit to severe cases
delay_fun_hosp <- onset2Outcome(onset2Hospquartiles)
outcomeDf$newOutcome <- outcomeDf$newHosp
#
fittingData <- get_fitting_data(outcomeDf, delay_fun_hosp,
                                baselineOutcomeProp = severeMean)
predictionSevere <- run_bayesian_model(fittingData,
                                 percentageOutcome=severeMean,
                                 percentageOutcomeRange=severeIC)

saveRDS(predictionCritical, "../datos_procesados/4_subreportingEstimate_critical.RDS")
saveRDS(predictionSevere, "../datos_procesados/4_subreportingEstimate_severe.RDS")

