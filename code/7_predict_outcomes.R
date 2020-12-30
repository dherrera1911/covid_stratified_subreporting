library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(RColorBrewer)
library(zoo)
library(ggpubr)
source("./functions_auxiliary.R")
source("./functions_analysis_subreporting.R")
source("./functions_fitting_predict.R")

forecastDate <- "2020-12-31"
daysUsed <- 60
daysFirstCum <- 8
cutLastDays <- 2
fecha <- "20201215"
fitModel <- "exponential"

individualCasesFile <- paste("../datos_procesados/1_resumen_casos_individuales_",
                             fecha, ".csv", sep="")
outcomeDynamicsFile <- paste("../datos_procesados/1_outcomeDynamics_", 
                             fecha, ".csv", sep="")
#subreportFile <- paste("../datos_procesados/8_estimate_subreporting_crit_stratified_",
#                       fecha, "_shorter.csv", sep="")
subreportFile <- paste("../datos_procesados/8_estimate_subreporting_crit_stratified_",
                       fecha, ".csv", sep="")
icuStayFile <- "../datos_procesados/1_estadias_cti.csv"

###################
###################
# load data
###################
###################

# literature estimates for different outcomes
criticalPercentageLit <-
  read.csv("../datos_procesados/severity_stratified_calculation/1_fitted_critical.csv",
           stringsAsFactors=FALSE) %>%
dplyr::filter(., study=="Fitted")

severePercentageLit <-
  read.csv("../datos_procesados/severity_stratified_calculation/1_fitted_severe.csv",
           stringsAsFactors=FALSE) %>%
dplyr::filter(., study=="Fitted")

ifrPercentageLit <-
  read.csv("../datos_procesados/severity_stratified_calculation/1_fitted_ifr.csv",
           stringsAsFactors=FALSE) %>%
dplyr::filter(., study=="Fitted")

# dataframe with individual cases and with summary dynamics
individualCasesDf <- read.csv(individualCasesFile, stringsAsFactors=FALSE) %>%
  dplyr::mutate(., FechaResultado=lubridate::as_date(FechaResultado),
                FechaInicioSintomas=lubridate::as_date(FechaInicioSintomas),
                FechaRec=lubridate::as_date(FechaRec),
                dateCritical=lubridate::as_date(dateCritical),
                dateHosp=lubridate::as_date(dateHosp))

outcomeDynamics <- read.csv(outcomeDynamicsFile, stringsAsFactors=FALSE) %>%
  dplyr::rename(., date=day) %>%
  dplyr::mutate(., date=lubridate::as_date(date))

# estimates of subreporting
subreportDf <- read.csv(subreportFile, stringsAsFactors=FALSE) %>%
  dplyr::mutate(., date=lubridate::as_date(date))

# data of different delays of cases in Uruguay
onset2ICUDelays <- read.csv("../datos_procesados/2_critical_delays.csv",
                            stringsAsFactors=FALSE)
onset2ICUquartiles <- quantile(onset2ICUDelays[[1]], na.rm=TRUE)[2:4]
delay_fun_crit <- onset2Outcome(onset2ICUquartiles)

onset2HospDelays <- read.csv("../datos_procesados/2_hospital_delays.csv",
                            stringsAsFactors=FALSE)
onset2Hospquartiles <- quantile(onset2HospDelays[[1]], na.rm=TRUE)[2:4]
delay_fun_hosp <- onset2Outcome(onset2Hospquartiles)

onset2ReportDelays <- read.csv("../datos_procesados/2_reporting_delay.csv",
                            stringsAsFactors=FALSE)

icuStay <- read.csv(icuStayFile) %>%
  dplyr::filter(., estadiaCTI < 40) %>%
  dplyr::mutate(., delay=estadiaCTI)
icuStayQuartiles <- quantile(icuStay$estadiaCTI)[2:4]
empiricalDist <- table(icuStay$estadiaCTI)/nrow(icuStay)

icu_stay_dist_cum <- function(x) {
  output <- rep(1, length(x))
  cumDist <- cumsum(empiricalDist)
  validInd <- which(x < length(cumDist))
  if (length(validInd>0)) {
    output[validInd] <- cumDist[x[validInd]+1]
  }
  return(output)
}

lastDay <- max(individualCasesDf$FechaInicioSintomas) - cutLastDays
forecastLength <- lubridate::interval(lastDay, forecastDate)/lubridate::days(1)

## public uruguay data
datosUruguayDir <- "https://raw.githubusercontent.com/GUIAD-COVID/datos-y-visualizaciones-GUIAD/master/datos/estadisticasUY.csv"
datosUruguay <- read.csv(url(datosUruguayDir), stringsAsFactors=FALSE) %>%
  dplyr::mutate(., cantTest=as.numeric(cantTest),
                porcentajePositivo=cantTestPositivos/cantTest*100,
                cantCasosNuevos=as.numeric(cantCasosNuevos),
                cantCTI=as.numeric(cantCTI),
                fecha=lubridate::dmy(fecha))

##########################
##########################
# Nowcast the number of cases with onset symptoms
# today, by calculating te proportion of cases that
# are expected to have been reported each day
# according to historical delay times
##########################
##########################

# correct by delay
correctedCases <- outcomeDynamics$newCases
multiplier <- dplyr::filter(onset2ReportDelays, delay>=0)$cumProb
multiplier <- rev(c(multiplier, rep(1, length(correctedCases) - length(multiplier))))
correctedCases <- correctedCases/multiplier
outcomeDynamics$correctedCases <- correctedCases

# correct by subreporting and cut last days indicated
correctedDynamics <- as_tibble(merge(outcomeDynamics, subreportDf)) %>%
  dplyr::mutate(., realCases = correctedCases/estimate) %>%
  dplyr::filter(., date <= lastDay)

# get initial cumulative value to use
initCum <- round(sum(correctedDynamics[c(-daysFirstCum:-1)-daysUsed+1+
                             nrow(correctedDynamics),][["realCases"]]))
correctedDynamics <- correctedDynamics[c(-daysUsed:-1)+1+nrow(correctedDynamics),]

##########################
##########################
# Forecast cases in next few days
##########################
##########################

# Fit subexponential curve to the current cases
casosDiarios <- round(correctedDynamics$realCases)
casosCumulativos <- cumsum(casosDiarios)
casosCumulativo <- c(initCum, casosDiarios[2:length(casosDiarios)])

if (fitModel=="subexponential") {
  #function to be used to predict
  prediction_fun_subexp <- function(modelCoefs) {
    generalizedGrowdthFun <- function(casosCumulativo) {
      dailyCases <- modelCoefs[2]*casosCumulativo^modelCoefs[1]
      return(dailyCases)
    }
    return(generalizedGrowdthFun)
  }
  reg <- modelo_sub_exponencial(casosDiarios, initCum)
  modelCoefs <- coef(reg)
  # fitted prediction function
  fitted_fun <- prediction_fun_subexp(modelCoefs)
  expectedNew <- fitted_fun(casosCumulativos)
  #lastCumulative <- sum(expectedNew)
  lastCumulative <- casosCumulativos[length(casosCumulativos)]
  extrapolatedNew <- round(propagate_sub_exponential(lastCumulative,
                                                 forecastLength, fitted_fun))
} else if (fitModel=="exponential") {
  #function to be used to predict
  prediction_fun_exp <- function(modelCoefs) {
    exponentialGrowthFun <- function(tiempo) {
      tiempo <- c(tiempo[1]-1, tiempo)
      inputCumulative <- modelCoefs[1]*exp(tiempo*modelCoefs[2])
      return(diff(inputCumulative))
    }
    return(exponentialGrowthFun)
  }
  reg <- modelo_exponencial(casosCumulativos)
  modelCoefs <- coef(reg)
  # fitted prediction function
  fitted_fun <- prediction_fun_exp(modelCoefs)
  expectedNew <- fitted_fun(c(1:length(casosDiarios)))
  extrapolatedNew <- round(fitted_fun(c(1:forecastLength)+length(casosDiarios)))
}

completeCases <- c(casosDiarios, extrapolatedNew)

#################################
#################################
# Estimate the expected percentage of
# outcome given the time-varying age
# distribution
#################################
#################################

criticalPercentageLit <- criticalPercentageLit %>%
  dplyr::mutate(., outcome=critical/100, outcomeL=criticalL/100,
                outcomeH=criticalH/100)

#ifrPercentageLit <- ifrPercentageLit %>%
#  dplyr::mutate(., outcome=IFR/100, outcomeL=IFR_L/100,
#                outcomeH=IFR_H/100)
#

critOutcomeDf <- expected_proportion_in_time(individualCasesDf,
                                             criticalPercentageLit,
                                             delay_fun_crit) %>%
  dplyr::mutate(., critProp=outcomeInTodaysCases/cases) %>%
  dplyr::filter(., day <= lastDay)

critOutcomeDf <- critOutcomeDf[c(-daysUsed:-1)+1+nrow(critOutcomeDf),]
recentProp <- mean(critOutcomeDf$critProp[c(-19:-1)+nrow(critOutcomeDf)])

completeProp <- c(critOutcomeDf$critProp, rep(recentProp, forecastLength))

critOutcomeDf$meanProp <- zoo::rollmean(critOutcomeDf$critProp, k=5, fill=NA)
critOutcomeDf$meanVal <- zoo::rollmean(critOutcomeDf$expectedOutcomes, k=5, fill=NA)
plotProp <- dplyr::filter(critOutcomeDf, day >= "2020-10-01") %>%
  ggplot(., aes(x=day, y=meanProp*100)) +
  geom_line() +
  ylab("Porcentaje criticos teorico") +
  ylim(c(0.5, 2)) +
  theme_bw()

#################################
#################################
# Put everything together to estimate
# critical cases
#################################
#################################
outcomeDynamics <- dplyr::filter(outcomeDynamics, date<=lastDay)
dates <- outcomeDynamics$date[c(-daysUsed:-1)+1+nrow(outcomeDynamics)]
extended <- max(dates)+c(1:forecastLength)
dates <- c(dates, extended)

indicatorVector <- c(rep("known", daysUsed), rep("predicted", forecastLength))
indList <- which(indicatorVector=="predicted")

futureOutcomes <- completeCases * completeProp

todaysOutcomes <- NULL
for (nr in c(1:length(futureOutcomes))) {
  delayFunT <- rev(delay_fun_crit(c(0:(nr-1))))
  todaysOutcomes[nr] <- sum(futureOutcomes[1:nr] * delayFunT)
}

predictedCrit <- todaysOutcomes[indList]

# compare predicted with recent outcomes
realOutcomes <- outcomeDynamics$newCritical[c(-daysUsed:-1)+1+nrow(outcomeDynamics)]
realOutcomes <- c(realOutcomes, rep(NA, forecastLength))
predictedCompare <- todaysOutcomes[-indList]

predictionDf <- data.frame(date=dates,
                           cases=completeCases,
                           criticalPredicted=todaysOutcomes,
                           criticalObserved=realOutcomes,
                           status=indicatorVector) %>%
  dplyr::mutate(., smoothObserved=zoo::rollmean(criticalObserved,
                                                k=3, fill=NA, na.rm=TRUE))

totalPredicted <- dplyr::filter(predictionDf, status=="predicted") %>%
  summarize(., totalNewPredicted=sum(criticalPredicted))

predictedICUNew <- dplyr::filter(predictionDf, status=="predicted")$criticalPredicted

# dias desde que tomar los casos criticos conocidos
diasTomar <- 30
critCases <- dplyr::filter(individualCasesDf, !is.na(IngresoCTI) &
                           IngresoCTI != "1-01-01"  &
                           IngresoCTI >= lastDay-diasTomar)

ctiRecientesFechas <- events_by_day(critCases$IngresoCTI)
ctiRecientesFechas <- dplyr::filter(datosUruguay, fecha==lastDay)$cantCTI

predictedICUOccupancy <- estimate_icu_occupancy(previousICU=ctiRecientesFechas,
                                                newICU=predictedICUNew,
                                                delay_fun=icu_stay_dist_cum)

predictedICUOccupancy <- predictedICUOccupancy[1-c(forecastLength:1)+length(predictedICUOccupancy)]

########################
# compare predicted data with what happened
# later
########################

##load data from future
#outcomeDynamicsFuture <- 
#  read.csv("../datos_procesados/1_outcomeDynamics.csv",
#                            stringsAsFactors=FALSE) %>%
#  dplyr::rename(., date=day) %>%
#  dplyr::mutate(., date=lubridate::as_date(date))
#
## crop data from future and put together with predictions
#outcomeDynamicsFuture <- dplyr::filter(outcomeDynamicsFuture,
#                                       date %in% predictionDf$date)
#croppedPredict <- dplyr::filter(predictionDf, date %in% outcomeDynamicsFuture$date)
#croppedPredict$futureCases <- outcomeDynamicsFuture$newCases
#croppedPredict$futureCrit <- outcomeDynamicsFuture$newCritical
#croppedPredict$futureCritSmooth <- zoo::rollmean(croppedPredict$futureCrit,
#                                                 k=3, fill=NA, na.rm=TRUE)
#
#futureCritPlot <- ggplot(croppedPredict, aes(x=date, y=criticalPredicted)) +
#  geom_line(color="red") +
#  geom_line(aes(y=futureCritSmooth), color="black")
#
#predictionSummary <- dplyr::filter(croppedPredict, status=="predicted") %>%
#  summarize(., totalPredicted=sum(criticalPredicted),
#            totalObserved=sum(futureCrit))
#
#futureCasePlot <- ggplot(croppedPredict, aes(x=date, y=criticalPredicted)) +
#  geom_line(color="red") +
#  geom_line(aes(y=futureCritSmooth), color="black")
#

############################
############################
# plot analysis
############################
############################

cutDay <- lastDay - daysUsed + 1
cutDynamics <- dplyr::filter(outcomeDynamics, date>=cutDay)

casosObservados <- cutDynamics %>%
  ggplot(., aes(x=date, y=newCases)) +
  geom_point() +
  #title("Casos observados") +
  theme_bw()
correccionDelay <- cutDynamics %>%
  ggplot(., aes(x=date, y=correctedCases)) +
  geom_point() +
  #title("Casos corregido delay reporte") +
  theme_bw()
correccionSubreport <- dplyr::filter(correctedDynamics, date>=cutDay) %>%
  ggplot(., aes(x=date, y=realCases)) +
  geom_point() +
  #title("Casos corregidos por subreporte") +
  theme_bw()
dinamicaCasos <- ggpubr::ggarrange(casosObservados, correccionDelay,
                          correccionSubreport, nrow=1)

naFillVec <- rep(NA, forecastLength)
fitDf <- data.frame(dia=cutDay-1+c(1:length(completeCases)),
                     observados=c(casosDiarios, naFillVec),
                     ajustados=c(expectedNew, extrapolatedNew),
                     proportionCrit=completeProp,
                     status=indicatorVector,
                     estimatedNewCrit=todaysOutcomes,
                     observedCTI=c(cutDynamics$newCritical, naFillVec))

fitPlot <- ggplot(fitDf, aes(x=dia, y=observados)) +
  geom_point(color="black") +
  geom_line(aes(y=ajustados), color="red") +
  theme_bw() +
  ylab("Casos reales estimados y ajustados")

critPropPlot <- ggplot(fitDf, aes(x=dia, y=proportionCrit, color=status,
                                  fill=status)) +
  geom_point() +
  theme_bw()

icuNewFit <- ggplot(fitDf, aes(x=dia)) +
  geom_line(aes(y=zoo::rollmean(observedCTI, k=3, fill=NA)), color="black") +
  geom_line(aes(y=estimatedNewCrit), color="red") +
  theme_bw() +
  ylab("Casos críticos nuevos") +
  xlab("Fecha")

# predicted ICU occupancy
occupancyDf <- data.frame(dia=lastDay+c(1:forecastLength),
                          icuOccupancy=predictedICUOccupancy)
occupancyPlot <- ggplot(occupancyDf, aes(x=dia, y=icuOccupancy)) +
  geom_line()

saturationDay <- dplyr::filter(occupancyDf, icuOccupancy>350)
saturationDay <- saturationDay$dia[1]

# compare predicted occupancy with observed occupancy
observedOccupancy <- dplyr::filter(datosUruguay,
                                   fecha %in% occupancyDf$dia) %>%
  dplyr::mutate(., dia=fecha)

compareOccupancyDf <- merge(occupancyDf, observedOccupancy, by="dia",
                            all=TRUE)
compareOccupancyPlot <- ggplot(compareOccupancyDf, aes(x=dia, y=icuOccupancy)) +
  geom_line() +
  geom_point(aes(y=cantCTI)) +
  xlab("Fecha") +
  ylab("Ocupación CTI (pacientes COVID-19)") +
  theme_bw()

ggsave("../graficas/7_prediccion_ocupacionCTI.png", compareOccupancyPlot,
       height=9, width=12, units="cm")

# compare predictions with new observations
realData <- read.csv("../datos_procesados/1_outcomeDynamics_20201220.csv",
                     stringsAsFactors=FALSE) %>%
  dplyr::filter(., day>(lastDay+cutLastDays))

predictionsFromToday <- dplyr::filter(fitDf, dia > (lastDay+cutLastDays) &
                                      dia <= (max(realData$day)))

sum(realData$newCritical)
sum(predictionsFromToday$estimatedNewCrit)

ggsave("../graficas/000_casos_observados.png", dinamicaCasos,
       height=8, width=22, units="cm")

ggsave("../graficas/000_casos_ajustados.png", fitPlot,
       height=9, width=12, units="cm")

ggsave("../graficas/000_proporcion_criticos.png", critPropPlot,
       height=9, width=12, units="cm")

ggsave("../graficas/000_criticos_ajustados.png", icuNewFit,
       height=9, width=12, units="cm")


