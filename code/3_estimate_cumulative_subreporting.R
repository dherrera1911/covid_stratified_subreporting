source("./functions_auxiliary.R")
source("./functions_analysis_subreporting.R")
source("./plotting_functions.R")

dataFile <- "../datos_procesados/1_resumen_casos_individuales.csv"
data <- read.csv(dataFile, stringsAsFactors=FALSE) %>%
  dplyr::mutate(., FechaInicioSintomas = lubridate::date(FechaInicioSintomas),
                dateCritical = lubridate::date(dateCritical),
                dateHosp = lubridate::date(dateHosp)) %>%
  as_tibble(.)

# Get the parameters for the simulation
criticalUru <- read.csv("../datos_procesados/severity_stratified_calculation/1_fitted_critical.csv")
severeUru <- read.csv("../datos_procesados/2_severe_proportion_uru.csv")

criticalUru <- dplyr::filter(criticalUru, study %in% c("Salje", "Driscoll", "Brazeau",
                                                    "Levin"))
severeUru <- dplyr::filter(severeUru, study %in% c("Salje", "Driscoll", "Brazeau",
                                                    "Levin", "Verity"))
criticalMean <- mean(criticalUru$critical_cases)
criticalIC <- dplyr::filter(criticalUru, !(study %in% c("Salje", "Brazeau"))) %>%
  with(., c(min(criticalL_cases), max(criticalH_cases)))

severeMean <- mean(severeUru$severe_cases)
severeIC <- dplyr::filter(severeUru, !(study %in% c("Brazeau", "Verity"))) %>%
  with(., c(min(severeL_cases), max(severeH_cases)))

percentCritical <- c(criticalIC[1], criticalMean, criticalIC[2])
percentHosp <- c(severeIC[1], severeMean, severeIC[2])

onset2ICUquartiles <- c(4, 7, 10) # datos para ajustar curva de delay
onset2Hospquartiles <- c(3, 5, 9) # datos para ajustar curva de delay



# get vectors of new ocurrences for cases and relvant outcomes
minDate <- date(min(data$FechaInicioSintomas, na.rm=TRUE))
maxDateVec <- with(data, c(FechaInicioSintomas, dateCritical, FechaResultado,
                            dateHosp))
maxDateVec <- maxDateVec[which(!is.na(as.character(maxDateVec)))]
maxDate <- date(max(maxDateVec))

newCriticalUru <- dplyr::filter(data, critico | fallecido) %>%
  with(., events_by_day(datesVec=dateCritical, minDate=minDate, maxDate=maxDate))
newHospUru <- dplyr::filter(data, critico | fallecido | hospitalizado) %>%
  with(., events_by_day(datesVec=.$dateHosp, minDate=minDate, maxDate=maxDate))
newCasesUru <- events_by_day(datesVec=data$FechaInicioSintomas,
                             minDate=minDate, maxDate=maxDate)

# fit the delay fun to empirical data
delay_fun_crit <- onset2Outcome_cumulative(onset2ICUquartiles)
delay_fun_hosp <- onset2Outcome_cumulative(onset2Hospquartiles)

# make column with dates
nDays <- lubridate::interval(minDate, maxDate) / lubridate::days(1)
daysColCrit <- seq(from=0, to=nDays) + lubridate::ymd(minDate) - onset2ICUquartiles[2]
daysColHosp <- seq(from=0, to=nDays) + lubridate::ymd(minDate) - round(onset2Hospquartiles[2])

#####################################################
# Estimate underreporting using the critical cases
#####################################################
subreportCrit <- list()
subreportCritPlot <- list()
subreportHosp <- list()
subreportHospPlot <- list()
for (pC in c(1:length(percentCritical))) {
  # subreporting using critical cases
  subreportCrit[[pC]] <- estimate_subreporting_cum(newCases = newCasesUru,
                                                 delay_fun = delay_fun_crit,
                                                 newOutcome = newCriticalUru,
                                                 baselinePercentOutcome = percentCritical[pC])
  subreportCrit[[pC]]$date <- daysColCrit
  subreportCritPlot[[pC]] <- plot_subreporting(subreportCrit[[pC]])
  # subreporting using hospitalised cases
  subreportHosp[[pC]] <- estimate_subreporting_cum(newCases = newCasesUru,
                                                 delay_fun = delay_fun_hosp,
                                                 newOutcome = newHospUru,
                                                 baselinePercentOutcome = percentHosp[pC])
  subreportHosp[[pC]]$date <- daysColHosp
  subreportHospPlot[[pC]] <- plot_subreporting(subreportHosp[[pC]])
}
names(subreportCrit) <- paste("rate_crit_", percentCritical, sep="")
names(subreportHosp) <- paste("rate_severe_", percentHosp, sep="")


#####################################################
# Estimate underreporting using the cases that started in the
# last N days
#####################################################
dateCut <- maxDate - 60

dataCut <- dplyr::filter(data, FechaInicioSintomas > dateCut)
newCriticalUruCut <- dplyr::filter(dataCut, critico | fallecido) %>%
  with(., events_by_day(datesVec=dateCritical, minDate=dateCut, maxDate=maxDate))
newHospUruCut <- dplyr::filter(dataCut, critico | fallecido | hospitalizado) %>%
  with(., events_by_day(datesVec=.$dateHosp, minDate=dateCut, maxDate=maxDate))
newCasesUruCut <- events_by_day(datesVec=dataCut$FechaInicioSintomas,
                             minDate=dateCut, maxDate=maxDate)
nDays <- lubridate::interval(dateCut, maxDate) / lubridate::days(1)
daysColCritCut <- lubridate::ymd(maxDate) - seq(from=nDays, to=0) + - onset2ICUquartiles[2]
daysColHospCut <- lubridate::ymd(maxDate) - seq(from=nDays, to=0) + - onset2Hospquartiles[2]

subreportCritCut <- list()
subreportCritCutPlot <- list()
subreportHospCut <- list()
subreportHospCutPlot <- list()
for (pC in c(1:length(percentCritical))) {
  # subreporting using critical cases
  subreportCritCut[[pC]] <- estimate_subreporting_cum(newCases = newCasesUruCut,
                                                 delay_fun = delay_fun_crit,
                                                 newOutcome = newCriticalUruCut,
                                                 baselinePercentOutcome = percentCritical[pC])
  subreportCritCut[[pC]]$date <- daysColCritCut
  subreportCritCutPlot[[pC]] <- plot_subreporting(subreportCritCut[[pC]])
  # subreporting using hospitalised cases
  subreportHospCut[[pC]] <- estimate_subreporting_cum(newCases = newCasesUruCut,
                                                 delay_fun = delay_fun_hosp,
                                                 newOutcome = newHospUruCut,
                                                 baselinePercentOutcome = percentHosp[pC])
  subreportHospCut[[pC]]$date <- daysColHospCut
  subreportHospCutPlot[[pC]] <- plot_subreporting(subreportHospCut[[pC]])
}
names(subreportCritCut) <- paste("rate_crit_", percentCritical, sep="")
names(subreportHospCut) <- paste("rate_severe_", percentHosp, sep="")

saveRDS(subreportCrit, "../datos_procesados/3_subreport_cumulative_crit.RDS")
saveRDS(subreportHosp, "../datos_procesados/3_subreport_cumulative_severe.RDS")
saveRDS(subreportCritCut, "../datos_procesados/3_subreport_cumulative_crit_recent.RDS")
saveRDS(subreportHospCut, "../datos_procesados/3_subreport_cumulative_severe_recent.RDS")

