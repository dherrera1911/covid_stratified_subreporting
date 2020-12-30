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
source("./age_strat_bayesian_fitting.R")

parameterFile <- "../datos_procesados/severity_stratified_calculation/1_model_summaries.csv"
modelParams <- read.csv(parameterFile, stringsAsFactors=FALSE)

dynamicsFile <- "../datos_procesados/1_outcomeDynamics.csv"
dynamicsDf <- read.csv(dynamicsFile, stringsAsFactors=FALSE) %>%
  dplyr::mutate(., day=lubridate::as_date(day))

individualCasesFile <- "../datos_procesados/1_resumen_casos_individuales.csv"
individualCases <- read.csv(individualCasesFile, stringsAsFactors=FALSE) %>%
  dplyr::filter(., !is.na(Edad)) %>%
  dplyr::mutate(., FechaInicioSintomas=lubridate::as_date(FechaInicioSintomas))


binsVec <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34",
            "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69",
            "70-74", "75-79", "80-84", "85+") 

stratInd <- bin_ages(individualCases$Edad, binsVec = binsVec)
stratInd[is.na(stratInd)] <- length(binsVec)+1
binsVec <- c(binsVec, NA)
individualCases$ageStrat <- factor(binsVec[stratInd], levels = binsVec)

casesByAge <- NULL
minDate <- min(dynamicsDf$day)
maxDate <- max(dynamicsDf$day)

for (strat in levels(individualCases$ageStrat)) {
  stratCases <- dplyr::filter(individualCases, ageStrat == strat)
  newCasesStrat <- events_by_day(datesVec=stratCases$FechaInicioSintomas,
                                 minDate=minDate, maxDate=maxDate)
  dynamicsDf[[strat]] <- newCasesStrat
}

outcomeDf <- dplyr::select(dynamicsDf, levels(individualCases$ageStrat))

onset2ICUDelays <- read.csv("../datos_procesados/2_critical_delays.csv",
                            stringsAsFactors=FALSE)
onset2ICUquartiles <- quantile(onset2ICUDelays[[1]], na.rm=TRUE)[2:4]

fechaAnalisis <- Sys.Date()

# Fit to critical cases
delay_fun_crit <- onset2Outcome(onset2ICUquartiles)
casesKnownStrat <- get_fitting_data_strat(outcomeDf, delay_fun_crit)
casesKnownStrat$date_num <- as.numeric(dynamicsDf$day)
casesKnownStrat$date <- dynamicsDf$day

criticalParams <- dplyr::filter(modelParams, Fit=="Critical")

slopeCrit <- criticalParams$Slope
interceptCrit <- criticalParams$Intercept
slopeCritCI <- with(criticalParams, c(SlopeL, SlopeH))
interceptCritCI <- with(criticalParams, c(InterceptL, InterceptH))

critVec <- dynamicsDf$newCritical

predictionCritical <- age_stratified_bayesian_model(casesKnownStrat,
                                                    outcomeVec=critVec,
                                                    slopeMean=slopeCrit,
                                                    slopeCI=slopeCritCI,
                                                    interceptMean=interceptCrit,
                                                    interceptCI=interceptCritCI)

write.csv(predictionCritical, "./8_estimate_subreporting_crit_stratified.csv",
          row.names=FALSE)

# Fit to severe cases
onset2HospDelays <- read.csv("../datos_procesados/2_hospital_delays.csv")
onset2Hospquartiles <- quantile(onset2HospDelays[[1]], na.rm=TRUE)[2:4]

delay_fun_hosp <- onset2Outcome(onset2Hospquartiles)
casesKnownStrat <- get_fitting_data_strat(outcomeDf, delay_fun_hosp)
casesKnownStrat$date_num <- as.numeric(dynamicsDf$day)
casesKnownStrat$date <- dynamicsDf$day

severeParams <- dplyr::filter(modelParams, Fit=="Severe")
slopeSevere <- severeParams$Slope
interceptSevere <- severeParams$Intercept
slopeSevereCI <- with(severeParams, c(SlopeL, SlopeH))
interceptSevereCI <- with(severeParams, c(InterceptL, InterceptH))

severeVec <- dynamicsDf$newHosp

predictionSevere <- age_stratified_bayesian_model(casesKnownStrat,
                                                    outcomeVec=severeVec,
                                                    slopeMean=slopeSevere,
                                                    slopeCI=slopeSevereCI,
                                                    interceptMean=interceptSevere,
                                                    interceptCI=interceptSevereCI)

write.csv(predictionSevere, "./8_estimate_subreporting_severe_stratified.csv",
          row.names=FALSE)

