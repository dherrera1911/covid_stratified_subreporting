library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(RColorBrewer)
library(zoo)
source("./functions_auxiliary.R")
source("./functions_analysis_subreporting.R")

savingDir <- "../datos_procesados/"
fecha <- "20201228"
dataFile <- paste("../datos_crudos_GACH/gach", fecha, ".csv", sep="")

cutDate <- "2020-11-01"

##################################
# Load data, rename variables, some tidying
##################################

dataDf <- read.csv2(dataFile) %>%
  as_tibble(.)

names(dataDf) <- c("caso", "persona", "FechaNac", "edad",
                   "edadTipo", "sexo", "codDep", "Departamento",
                   "codLoc", "Localidad", "codBarrio", "Barrio",
                   "cod2Barrio", "Salud", "Bombero", "Policia",
                   "Militar", "Residencial", "Contacto", "FechaExp",
                   "FechaInicioSintomas", "Importado", "PaisImp",
                   "FechaRegreso", "Estudio", "FechaResultado",
                   "Resultado", "Internacion", "InternacionHosp",
                   "EgresoHosp", "CI", "IngresoCI", "EgresoCI",
                   "CTI", "IngresoCTI", "EgresoCTI", "ecEgreso",
                   "Fallece", "FFS", "FechaRec", "Sintomatico")

dataDf <- dplyr::select(dataDf, -edad, -edadTipo, -Salud, -Bombero,
                        -Policia, -Militar, -Residencial,
                        -PaisImp, -FechaRegreso)

dataDf <- dataDf %>%
  dplyr::mutate(., FechaNac = string2date(FechaNac),
                FechaInicioSintomas = string2date(FechaInicioSintomas),
                IngresoHosp = string2date(InternacionHosp),
                IngresoCI = string2date(IngresoCI),
                IngresoCTI = string2date(IngresoCTI),
                FechaRec = string2date(FechaRec),
                FechaResultado = string2date(FechaResultado),
                FechaExp = string2date(FechaExp),
                EgresoCTI = string2date(EgresoCTI),
                Edad = round(lubridate::interval(FechaNac, lubridate::today())/
                             lubridate::years(1)),
                hospitalizado = (Internacion=="S" | CI=="S" | CTI=="S"),
                critico = (CI=="S" | CTI=="S"),
                fallecido =  ecEgreso=="F")

# Filter out imported cases
dataDf <- dplyr::filter(dataDf, Importado != "S")
# Fix broken date
brokenDate <- with(dataDf, which(FechaInicioSintomas == "1-01-01" |
                                 is.na(FechaInicioSintomas)))
dataDf$FechaInicioSintomas[brokenDate] <- dataDf$FechaResultado[brokenDate]

#############################
# Do some more tidying
#############################

# Get the date in which critical cases become critical
criticalDf <- dataDf %>%
  dplyr::filter(., critico | fallecido) %>%
  dplyr::rowwise(.) %>%
  dplyr::mutate(., dateCritical = ifelse(critico & !(is.na(IngresoCI)
                                                     & is.na(IngresoCTI)),
                                         min(IngresoCI, IngresoCTI, na.rm=TRUE),
                                         ifelse(is.na(FechaRec),
                                                FechaResultado, FechaRec)),
                dateHosp = IngresoHosp, asymptomatic = FALSE)

# Get the date in which hospitalized cases were hospitalized
hospDf <- dataDf %>%
  dplyr::filter(., hospitalizado & !(critico | fallecido)) %>%
  dplyr::rowwise(.) %>%
  dplyr::mutate(., dateHosp = IngresoHosp, dateCritical = NA,
                asymptomatic = FALSE)

# merge hospitalization and critical dfs
expandedHospDf <- rbind(hospDf, criticalDf)

# Make a dataframe with summary of each individual case
nonCriticalDf <- dplyr::filter(dataDf, !(critico | fallecido | hospitalizado)) %>%
  dplyr::select(caso, sexo, Edad, FechaInicioSintomas, FechaResultado, FechaRec,
                IngresoCTI, EgresoCTI) %>%
  dplyr::mutate(., critico = FALSE, fallecido = FALSE, hospitalizado = FALSE,
                asymptomatic = FechaInicioSintomas == FechaResultado,
                dateCritical = NA, dateHosp = NA)

reducedCriticalDf <- dplyr::select(expandedHospDf, caso, sexo, Edad,
                                   FechaInicioSintomas, FechaResultado, critico,
                                   fallecido, dateCritical, hospitalizado,
                                   dateHosp, asymptomatic, FechaRec,
                                   IngresoCTI, EgresoCTI)

summarizedDf <- rbind(nonCriticalDf, reducedCriticalDf)
write.csv(summarizedDf, paste(savingDir, "1_resumen_casos_individuales_",
                              fecha, ".csv", sep=""), row.names=FALSE)

summarizedDfCut <- dplyr::filter(summarizedDf, FechaInicioSintomas>="2020-11-01")
write.csv(summarizedDfCut, paste(savingDir, "1_resumen_casos_individuales_",
                              fecha, "_cortado.csv", sep=""), row.names=FALSE)

#############################
# Make dynamics dataframe
#############################

# Do some dates tidying
minDate <- min(dataDf$FechaInicioSintomas, na.rm=TRUE)
maxDate <- with(dataDf, max(c(FechaInicioSintomas, FechaRec, FechaResultado),
                            na.rm=TRUE))
nDays <- lubridate::interval(minDate, maxDate) / lubridate::days(1)
allDays <- seq(from=0, to=nDays) + lubridate::ymd(minDate) 
# Get the number of each new event per day
newCritical <- events_by_day(datesVec=criticalDf$dateCritical,
                             minDate=minDate, maxDate=maxDate)
newHosp <- events_by_day(datesVec=expandedHospDf$dateHosp,
                             minDate=minDate, maxDate=maxDate)
newCases <- events_by_day(datesVec=dataDf$FechaInicioSintomas,
                             minDate=minDate, maxDate=maxDate)

outcomeDynamicsDf <- data.frame(day = allDays, newCases = newCases,
                         newCritical = newCritical, newHosp = newHosp)

write.csv(outcomeDynamicsDf, paste(savingDir, "1_outcomeDynamics_",
                                   fecha, ".csv", sep=""), row.names=FALSE)

# Get outcome dynamics from cut cases
minDate <- cutDate
criticalDfCut <- dplyr::filter(criticalDf, FechaInicioSintomas>="2020-11-01")
newCriticalCut <- events_by_day(datesVec=criticalDfCut$dateCritical,
                             minDate=minDate, maxDate=maxDate)
newHospCut <- events_by_day(datesVec=expandedHospDf$dateHosp,
                             minDate=minDate, maxDate=maxDate)
newCasesCut <- events_by_day(datesVec=dataDf$FechaInicioSintomas,
                             minDate=minDate, maxDate=maxDate)

nDays <- lubridate::interval(minDate, maxDate) / lubridate::days(1)
allDays <- seq(from=0, to=nDays) + lubridate::ymd(minDate) 
outcomeDynamicsCutDf <- data.frame(day=allDays, newCases=newCasesCut,
                         newCritical=newCriticalCut, newHosp=newHospCut)

write.csv(outcomeDynamicsCutDf, paste(savingDir, "1_outcomeDynamics_",
                                   fecha, "_cut.csv", sep=""), row.names=FALSE)

##########################
# Export ICU stay length vector
##########################
critDelayDf <- dplyr::filter(dataDf, CTI=="S") %>%
  dplyr::filter(., IngresoCTI<=lubridate::date("2020-11-01")) %>%
  dplyr::mutate(., estadiaCTI= round(lubridate::interval(IngresoCTI, EgresoCTI)/
                             lubridate::days(1))) %>%
  dplyr::filter(., !is.na(estadiaCTI)) %>%
  dplyr::select(., estadiaCTI)

write.csv(critDelayDf, paste("../datos_procesados/1_estadias_cti_",
                             fecha, ".csv", sep=""), row.names=FALSE)

