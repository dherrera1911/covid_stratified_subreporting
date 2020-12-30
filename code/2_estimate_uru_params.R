library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(RColorBrewer)
source("./functions_auxiliary.R")
source("./functions_analysis_subreporting.R")

fecha <- "20201222"
dataFile <- paste("../datos_procesados/1_resumen_casos_individuales_",
                  fecha, ".csv", sep="")

dataDf <- read.csv(dataFile, stringsAsFactors=FALSE) %>%
  dplyr::mutate(FechaInicioSintomas = string2date(FechaInicioSintomas),
                         FechaResultado = string2date(FechaResultado),
                         dateCritical = string2date(dateCritical),
                         dateHosp = string2date(dateHosp),
                         FechaRec = string2date(FechaRec))

####################################
# Bin the ages for the parameters calculations
####################################

binsVec <- c("0-15", "16-30", "31-45", "46-60", "61-75", "75+")
stratInd <- bin_ages(dataDf$Edad, binsVec = binsVec)
stratInd[is.na(stratInd)] <- length(binsVec)+1
binsVec <- c(binsVec, NA)
dataDf$ageStrat <- factor(binsVec[stratInd], levels = binsVec)

####################################
# Calculate delay distributions
####################################

# extract information tidily
criticalDf <- dplyr::filter(dataDf, critico | fallecido) %>%
  dplyr::mutate(., symptom2Onset = round(lubridate::interval(FechaInicioSintomas,
                                                    dateCritical)/lubridate::days(1)),
                   icuStay = round(lubridate::interval(dateCritical,
                                                    FechaRec)/lubridate::days(1)),
                    res2Crit = round(lubridate::interval(FechaResultado,
                                                    dateCritical)/lubridate::days(1)))

hospDf <- dplyr::filter(dataDf, hospitalizado) %>%
  dplyr::mutate(., symptom2Hosp = round(lubridate::interval(FechaInicioSintomas,
                                                            dateHosp)/lubridate::days(1)))

# extract information tidily
deathDf <- dplyr::filter(dataDf, fallecido) %>%
  dplyr::mutate(., res2death = round(lubridate::interval(FechaResultado,
                                                    FechaRec)/lubridate::days(1)))


# put into vectors and extract quantiles
symptom2CritVec <- criticalDf$symptom2Onset
symptom2CritVec <- symptom2CritVec[symptom2CritVec > 0]
delayCritQuantiles <- quantile(symptom2CritVec, na.rm=TRUE)

# put into vectors and extract quantiles
res2CritVec <- criticalDf$res2Crit
res2CritVec <- res2CritVec[res2CritVec > 0]
res2CritQuantiles <- quantile(res2CritVec, na.rm=TRUE)

icuStayVec <- criticalDf$icuStay
icuStayQuantiles <- quantile(icuStayVec, na.rm=TRUE)

symptom2HospVec <- hospDf$symptom2Hosp
symptom2HospVec <- symptom2HospVec[symptom2HospVec > 0]
delayHospQuantiles <- quantile(symptom2HospVec, na.rm=TRUE)

write.csv(symptom2CritVec, "../datos_procesados/2_critical_delays.csv",
          row.names=FALSE)
write.csv(symptom2HospVec, "../datos_procesados/2_hospital_delays.csv",
          row.names=FALSE)

# reporting delay
cutDate <- "2020-11-01"
delayDf <- dplyr::filter(dataDf, FechaInicioSintomas >= cutDate) %>%
  dplyr::mutate(., symptom2Test = round(lubridate::interval(FechaInicioSintomas,
                                                    FechaResultado)/lubridate::days(1)))
cumDelayFun <- ecdf(delayDf$symptom2Test)
reportDelayDf <- data.frame(delay = c(-11:11), cumProb = cumDelayFun(c(-11:11)))
write.csv(reportDelayDf, "../datos_procesados/2_reporting_delay.csv",
          row.names=FALSE)

####################################
# Calculate outcome proportion by age strat
####################################
outcomeByAge <- group_by(dataDf, ageStrat) %>%
  summarize(., count = n(),
            critCount = sum(as.integer(critico | fallecido)),
            hospCount = sum(as.integer(hospitalizado | critico | fallecido)),
            deathCount = sum(as.integer(fallecido)), 
            critical = 100*mean(critico | fallecido),
            death = 100*mean(fallecido),
            hospitalized = 100*mean(hospitalizado),
            asymptomatic = 100*mean(asymptomatic, na.rm=TRUE)) %>%
  ungroup()
write.csv(outcomeByAge, "../datos_procesados/2_outcomes_uruguay.csv", row.names=FALSE)

letalityByAge <- group_by(dataDf, ageStrat) %>%
  summarize(., critCount = sum(as.integer(critico | fallecido)),
            hospCount = sum(as.integer(hospitalizado | critico | fallecido)),
            deathCount = sum(as.integer(fallecido)), 
            critLetality = deathCount/critCount*100,
            hospLetality = deathCount/hospCount*100) %>%
  ungroup()
write.csv(letalityByAge, "../datos_procesados/2_letality_uruguay.csv", row.names=FALSE)


# Global outcome not stratified 
globalOutcome <- summarize(dataDf, critical = 100*mean(critico | fallecido),
            death = 100*mean(fallecido),
            hospitalized = 100*mean(hospitalizado),
            asymptomatic = 100*mean(asymptomatic, na.rm=TRUE))
write.csv(letalityByAge, "../datos_procesados/2_desenlaces_uruguay.csv", row.names=FALSE)


#####################
##### Calculate population parameters for the Uruguay demographics ####
#####################
demographyVals <- read.csv("../datos_procesados/severity_stratified_calculation/2_uru_demographic_outcomes.csv",
                           stringsAsFactors=FALSE) %>%
  dplyr::filter(., outcome %in% c("Critical", "Severe") & study == "Fitted") %>%
  dplyr::select(., -X)

literatureCrit <- read.csv("../datos_procesados/severity_stratified_calculation/1_fitted_critical.csv",
                             stringsAsFactors=FALSE) %>%
  dplyr::filter(., study=="Fitted")
literatureSevere <- read.csv("../datos_procesados/severity_stratified_calculation/1_fitted_severe.csv",
                             stringsAsFactors=FALSE) %>%
  dplyr::filter(., study=="Fitted")

litAgeVec <- literatureCrit$age
stratInd <- bin_ages(dataDf$Edad, binsVec = litAgeVec)
ageCounts <- table(stratInd)
ageProps <- ageCounts/sum(ageCounts)

literatureCrit$covidAgeDist <- ageProps
literatureSevere$covidAgeDist <- ageProps

casesCritical <- summarize(literatureCrit, critical=sum(critical*covidAgeDist),
                           criticalL=sum(criticalL*covidAgeDist),
                           criticalH=sum(criticalH*covidAgeDist))
casesSevere <- summarize(literatureSevere, severe=sum(severe*covidAgeDist),
                           severeL=sum(severeL*covidAgeDist),
                           severeH=sum(severeH*covidAgeDist))

percentageVec <- c(casesCritical$critical,
                   dplyr::filter(demographyVals, outcome=="Critical")[["percentage"]],
                   casesSevere$severe,
                   dplyr::filter(demographyVals, outcome=="Severe")[["percentage"]])
percentageLVec <- c(casesCritical$criticalL,
                   dplyr::filter(demographyVals, outcome=="Critical")[["percentageL"]],
                   casesSevere$severe,
                   dplyr::filter(demographyVals, outcome=="Severe")[["percentageL"]])
percentageHVec <- c(casesCritical$criticalH,
                   dplyr::filter(demographyVals, outcome=="Critical")[["percentageH"]],
                   casesSevere$severe,
                   dplyr::filter(demographyVals, outcome=="Severe")[["percentageH"]])
outcomeVec <- c(rep("critical", 4), rep("severe", 4))
populationVec <- rep(c("cases", "demography"), 4)

averageOutcomeDf <- data.frame(outcome=outcomeVec, population=populationVec,
                               percentage=percentageVec,
                               percentageL=percentageLVec,
                               percentageH=percentageHVec)

# print table to copy into document
write.csv(averageOutcomeDf,
          "../datos_procesados/2_expected_outcome_uruguay.csv", row.names=FALSE)

