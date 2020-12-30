library(lubridate)
library(gridExtra)
library(ggpubr)
source("./plotting_functions.R")
source("./functions_analysis_subreporting.R")
source("./functions_auxiliary.R")

#######################################################
# 1 Plot histograms of symptoms to outcome delays
#######################################################
symptom2CritVec <- read.csv("../datos_procesados/2_critical_delays.csv")
symptom2HospVec <- read.csv("../datos_procesados/2_hospital_delays.csv")

critDelayDf <- data.frame(delay = symptom2CritVec[[1]]) %>%
  dplyr::filter(., delay <40)
severeDelayDf <- data.frame(delay = symptom2HospVec[[1]]) %>%
  dplyr::filter(., delay <40)

critDelayPlot <- plot_delay_histogram(critDelayDf, nBins=35,
                                      title="Tiempo en volverse crítico") +
  xlab("Días entre inicio de síntomas y crítico")

severeDelayPlot <- plot_delay_histogram(severeDelayDf, nBins=42,
                                      title="Tiempo en volverse severo") +
  xlab("Días entre inicio de síntomas y severo")
delaysPlot <- gridExtra::grid.arrange(critDelayPlot, severeDelayPlot,
                                          ncol = 2)

ggsave("../graficas/1_histograma_delays.png", delaysPlot,
       width = 16, height = 8, units = "cm")

#######################################################
# 3 Plot results from cumulative subreporting estimate
#######################################################
#subreportCritCum <- readRDS("../datos_procesados/3_subreport_cumulative_crit.RDS")
#subreportSevereCum <- readRDS("../datos_procesados/3_subreport_cumulative_severe.RDS")
#subreportCritCutCum <- readRDS("../datos_procesados/3_subreport_cumulative_crit_recent.RDS")
#subreportSevereCutCum <- readRDS("../datos_procesados/3_subreport_cumulative_severe_recent.RDS")
#
#subreportCritPlot <- list()
#subreportSeverePlot <- list()
#subreportCritCutPlot <- list()
#subreportSevereCutPlot <- list()
#critVals <- paste("% críticos referencia: ", c(1, 1.4, 2.5), "%", sep="")
#severeVals <- paste("% severos referencia: ", c(1.8, 5, 7.5), "%", sep="")
#for (val in c(1:length(subreportCritCum))) {
#  subreportCritPlot[[val]] <- plot_subreporting(subreportCritCum[[val]],
#                                         title = critVals[val])
#  subreportSeverePlot[[val]] <- plot_subreporting(subreportSevereCum[[val]],
#                                           title = severeVals[val])
#  subreportCritCutPlot[[val]] <- plot_subreporting(subreportCritCutCum[[val]],
#                                         title = critVals[val])
#  subreportSevereCutPlot[[val]] <- plot_subreporting(subreportSevereCutCum[[val]],
#                                           title = severeVals[val])
#}
#
#fullList <- c(subreportCritPlot, subreportSeverePlot)
#fullListCut <- c(subreportCritCutPlot, subreportSevereCutPlot)
#
#cumulativeSubreportPlot <- gridExtra::grid.arrange(grobs=fullList, ncol = 2,
#                                                   as.table=FALSE)
#cumulativeSubreportCutPlot <- gridExtra::grid.arrange(grobs=fullListCut, ncol = 2,
#                                                   as.table=FALSE)
#
#ggsave("../graficas/3_subreporte_cumulativo.png", cumulativeSubreportPlot,
#       width = 16, height = 16, units = "cm")
#ggsave("../graficas/3_subreporte_cumulativo_reciente.png", cumulativeSubreportCutPlot,
#       width = 16, height = 16, units = "cm")
#
#
########################################################
## 4 Plot results from instantaneous subreporting estimate
########################################################
#subreportCrit <- readRDS("../datos_procesados/4_subreportingEstimate_critical_20201127.RDS")
#subreportSevere <- readRDS("../datos_procesados/4_subreportingEstimate_severe_20201127.RDS")
#
#subreportCritPlot <- plot_subreporting(subreportCrit,
#                                       title = "Estimacion con criticos")
#subreportSeverePlot <- plot_subreporting(subreportSevere,
#                                         title = "Estimacion con severos")
#estimateSubreport <- gridExtra::grid.arrange(subreportCritPlot, subreportSeverePlot,
#                                          ncol = 2)
#
#ggsave("../graficas/4_estimado_instantaneo.png", estimateSubreport,
#       width = 16, height = 8, units = "cm")
#
#subreportCritShort <- readRDS("../datos_procesados/4_subreportingEstimate_critical_shorter.RDS")
#subreportSevereShort <- readRDS("../datos_procesados/4_subreportingEstimate_severe_shorter.RDS")
#subreportCritShortPlot <- plot_subreporting(subreportCritShort,
#                                       title = "Estimacion con criticos")
#subreportSevereShortPlot <- plot_subreporting(subreportSevereShort,
#                                         title = "Estimacion con severos")
#estimateSubreportShort <- gridExtra::grid.arrange(subreportCritShortPlot,
#                                                  subreportSevereShortPlot,
#                                                  ncol = 2)
#
#ggsave("../graficas/4_estimado_instantaneo_corto.png", estimateSubreportShort,
#       width = 16, height = 8, units = "cm")
#

#######################################################
# 4 Plot estimated real dynamics
#######################################################
#daysCut <- 6
#caseDynamics <- readRDS("../datos_procesados/1_outcomeDynamics.RDS") %>%
#  dplyr::rename(., date = day) %>%
#  dplyr::mutate(., date = date(date))
#
#realDynCrit <- merge(caseDynamics, subreportCrit, by = "date") %>%
#  dplyr::mutate(., caseEst = newCases/estimate, caseEstH = newCases/lower,
#                caseEstL = newCases/upper) %>%
#  as_tibble(.)
#realDynSevere <- merge(caseDynamics, subreportSevere, by = "date") %>%
#  dplyr::mutate(., caseEst = newCases/estimate, caseEstH = newCases/lower,
#                caseEstL = newCases/upper) %>%
#  as_tibble(.)
#
## remove last daysCut days, to prevent distortions from symptoms backlog
#realDynCrit <- realDynCrit[c(1:(nrow(realDynCrit)-daysCut)),]
#realDynSevere <- realDynSevere[c(1:(nrow(realDynSevere)-daysCut)),]
#
#realDynCritPlot <- plot_estimated_dyn(realDynCrit, meanWindow=3,
#                                      title = "Estimacion con criticos")
#realDynSeverePlot <- plot_estimated_dyn(realDynSevere, meanWindow=3,
#                                      title = "Estimacion con severos")
#absSubreportCritPlot <- plot_abs_subreport(realDynCrit, meanWindow=3,
#                                      title = "Estimacion con criticos")
#absSubreportSeverePlot <- plot_abs_subreport(realDynSevere, meanWindow=3,
#                                      title = "Estimacion con severos")
#
#estimatedCases <- gridExtra::grid.arrange(realDynCritPlot, realDynSeverePlot,
#                                          ncol = 1)
#absoluteSubreport <- gridExtra::grid.arrange(absSubreportCritPlot,
#                                             absSubreportSeverePlot,
#                                             ncol = 1)
#
#ggsave("../graficas/4_estimado_casos_reales.png", estimatedCases,
#       width = 16, height = 12, units = "cm")
#ggsave("../graficas/4_estimado_casos_no_detectados.png", absoluteSubreport,
#       width = 16, height = 12, units = "cm")
#
########################################################
## 5 Plot simulated dynamics
########################################################
#simsResults <- readRDS("../datos_procesados/5_simulations.RDS") %>%
#  as_tibble(.)
#simsResultsShort <- readRDS("../datos_procesados/5_simulations_shorter.RDS") %>%
#  as_tibble(.)
#
## plot data standard
#simsPlots <- list()
#infPlots <- list()
#estimatePlots <- list()
#plottingSamples <- list()
#estimatedDyn <- list()
#for (f in c(1:length(unique(simsResults$funcInfInd)))) {
#  simsPlots[[f]] <- list()
#  infDf <- dplyr::filter(simsResults, funcInfInd==f)
#  for (s in c(1:length(unique(infDf$funcSRInd)))) {
#    simsPlots[[f]][[s]] <- list()
#    infSRDf <- dplyr::filter(infDf, funcSRInd==s)
#    for (r in c(1:length(unique(infSRDf$nRep)))) {
#      infSRDf_r <- dplyr::filter(infSRDf, nRep==r)
#      simsPlots[[f]][[s]][[r]] <- plot_subreport_sim(infSRDf_r)
#      if (r == 1) {
#        plottingSamples[[length(plottingSamples)+1]] <- simsPlots[[f]][[s]][[r]]
#        infPlots[[length(infPlots)+1]] <- plot_infections_sim(infSRDf_r)
#        estimatePlots[[length(estimatePlots)+1]] <- plot_infections_sim(infSRDf_r,
#                                                                   plotEst=TRUE)
#      }
#    }
#  }
#}
#
#simsSamplesPlot <- ggpubr::ggarrange(plotlist=plottingSamples, ncol=5, nrow=3,
#                                     common.legend=TRUE)
#examplesInfections <- ggpubr::ggarrange(plotlist=infPlots, ncol=5, nrow=3,
#                                        common.legend=TRUE)
#infectionsEstimates <- ggpubr::ggarrange(plotlist=estimatePlots, ncol=5, nrow=3,
#                                        common.legend=TRUE)
#
#simsSamplesPlot <- ggpubr::annotate_figure(simsSamplesPlot,
#                top = text_grob("Estimacion subreporte simulaciones",
#                                color = "black", face = "bold", size = 14),
#                left = text_grob("% infecciones reportadas",
#                                 color = "black", rot = 90))
#examplesInfections <- ggpubr::annotate_figure(examplesInfections,
#                top = text_grob("Epidemias simuladas",
#                                color = "black", face = "bold", size = 14),
#                left = text_grob("Número de casos",
#                                 color = "black", rot = 90))
#infectionsEstimates <- ggpubr::annotate_figure(infectionsEstimates,
#                top = text_grob("Estimacion de infecciones en epidemias simuladas",
#                                color = "black", face = "bold", size = 14),
#                left = text_grob("Número de casos",
#                                 color = "black", rot = 90))
#
#ggsave("../graficas/5_subreporte_simulado.png", simsSamplesPlot,
#       width = 24, height = 17, units = "cm")
#ggsave("../graficas/5_dinamicas_simuladas.png", examplesInfections,
#       width = 24, height = 17, units = "cm")
#ggsave("../graficas/5_estimaciones_metodo_simulacion.png", examplesInfections,
#       width = 24, height = 17, units = "cm")
#
#
## plot data with shorter temporal correlation
#simsPlotsShort <- list()
#plottingSamplesShort <- list()
#infPlotsShort <- list()
#estimatePlotsShort <- list()
#infFuncs <- unique(simsResultsShort$funcInfInd)
#for (f in c(1:length(infFuncs))) {
#  simsPlotsShort[[f]] <- list()
#  ff <- infFuncs[f]
#  infDf <- dplyr::filter(simsResultsShort, funcInfInd==ff)
#  for (s in c(1:length(unique(infDf$funcSRInd)))) {
#    simsPlotsShort[[f]][[s]] <- list()
#    infSRDf <- dplyr::filter(infDf, funcSRInd==s)
#    for (r in c(1:length(unique(infSRDf$nRep)))) {
#      infSRDf_r <- dplyr::filter(infSRDf, nRep==r)
#      simsPlotsShort[[f]][[s]][[r]] <- plot_subreport_sim(infSRDf_r)
#      if (r == 1) {
#        plottingSamplesShort[[length(plottingSamplesShort)+1]] <- simsPlotsShort[[f]][[s]][[r]]
#        infPlotsShort[[length(infPlotsShort)+1]] <- plot_infections_sim(infSRDf_r)
#        estimatePlotsShort[[length(estimatePlotsShort)+1]] <- plot_infections_sim(infSRDf_r,
#                                                                   plotEst=TRUE)
#      }
#    }
#  }
#}
#
#simsSamplesShortPlot <- ggpubr::ggarrange(plotlist=plottingSamplesShort,
#                                          ncol=5, nrow=3, common.legend=TRUE)
#infectionsEstimatesShort <- ggpubr::ggarrange(plotlist=estimatePlotsShort,
#                                          ncol=5, nrow=3, common.legend=TRUE)
#
#simsSamplesShortPlot <- ggpubr::annotate_figure(simsSamplesShortPlot,
#                top = text_grob("Estimacion subreporte simulaciones",
#                                color = "black", face = "bold", size = 14),
#                left = text_grob("% infecciones reportadas",
#                                 color = "black", rot = 90))
#infectionsEstimatesShort <- ggpubr::annotate_figure(infectionsEstimatesShort,
#                top = text_grob("Estimacion de infecciones en epidemias simuladas",
#                                color = "black", face = "bold", size = 14),
#                left = text_grob("Número de casos",
#                                 color = "black", rot = 90))
#
#ggsave("../graficas/5_subreporte_simulado_short.png", simsSamplesShortPlot,
#       width = 24, height = 17, units = "cm")
#ggsave("../graficas/5_estimaciones_metodo_simulacion_short.png",
#       infectionsEstimatesShort, width = 24, height = 17, units = "cm")
#
## plot summary of estimates
#simPerformancePlot <- plot_simulation_fit_summary(simsResults)
#simPerformancePlotShort <- plot_simulation_fit_summary(simsResultsShort)
#ggsave("../graficas/5_rendimiento_simulaciones.png", simPerformancePlot,
#       width = 24, height = 17, units = "cm")
#
#
########################################################
## 7 Plot predicted outcomes vs actual outcomes
########################################################
#predictionDf <- readRDS("../datos_procesados/7_expected_outcomes.RDS")
#criticalPrediction <- predictionDf$critical %>%
#  dplyr::mutate(., smoothOutcomeRate = zoo::rollmean(outcomeInTodaysCases/cases,
#                                 k=5, fill=NA), newOutcome = newCritical)
#severePrediction <- predictionDf$severe %>%
#  dplyr::mutate(., smoothOutcomeRate = zoo::rollmean(outcomeInTodaysCases/cases,
#                                 k=5, fill=NA), newOutcome = newHosp)
#
#critPlot <- plot_estimated_outcomes(criticalPrediction,
#                                    title="Predicción de críticos") +
#  ylab("Casos críticos")
#severePlot <- plot_estimated_outcomes(severePrediction,
#                                      title="Predicción de severos") +
#  ylab("Casos severos")
#
#predictionPlot <- ggpubr::ggarrange(plotlist=list(critPlot, severePlot),
#                                          ncol=1, nrow=2, common.legend=TRUE)
#
#ggsave("../graficas/7_prediccion_hospitalaria.png", predictionPlot,
#       width=12, height=7)
#
#
#adjustedCriticalPred <- dplyr::rename(subreportCritStrat, day = date) %>%
#  merge(., criticalPrediction, by="day") %>%
#  dplyr::mutate(., expectedOutcomes = expectedOutcomes / estimate)
#adjustedSeverePred <- dplyr::rename(subreportSevereStrat, day = date) %>%
#  merge(., severePrediction, by="day") %>%
#  dplyr::mutate(., expectedOutcomes = expectedOutcomes / estimate)
#critPlotAdj <- plot_estimated_outcomes(adjustedCriticalPred,
#                                    title="Predicción de críticos") +
#  ylab("Casos críticos")
#severePlotAdj <- plot_estimated_outcomes(adjustedSeverePred,
#                                    title="Predicción de severos") +
#  ylab("Casos severos")
#
#predictionPlotAdj <- ggpubr::ggarrange(plotlist=list(critPlotAdj, severePlotAdj),
#                                          ncol=1, nrow=2, common.legend=TRUE)
#
#ggsave("../graficas/7_prediccion_hospitalaria_ajustada.png", predictionPlotAdj,
#       width=12, height=7)
#

#######################################################
# 8 Plot results from instantaneous subreporting with age stratified dynamics
#######################################################

fecha <- "20201222"
critSRFile <- paste("../datos_procesados/8_estimate_subreporting_crit_stratified_",
                    fecha, "_shorter.csv", sep="")
severeSRFile <- paste("../datos_procesados/8_estimate_subreporting_severe_stratified_",
                    fecha, "_shorter.csv", sep="")
subreportCritStrat <- read.csv(critSRFile, stringsAsFactors=FALSE)
subreportCritStrat$date <- lubridate::as_date(subreportCritStrat$date)

subreportSevereStrat <- read.csv(severeSRFile, stringsAsFactors=FALSE)
subreportSevereStrat$date <- lubridate::as_date(subreportSevereStrat$date)

medianCritDelay <- median(critDelayDf[[1]])
medianSevereDelay <- median(severeDelayDf[[1]])
subreportCritStrat$date <- subreportCritStrat$date - medianCritDelay
subreportSevereStrat$date <- subreportSevereStrat$date - medianSevereDelay

subreportCritStratPlot <- plot_subreporting(subreportCritStrat,
                                       title = "Estimacion con criticos")
subreportSevereStratPlot <- plot_subreporting(subreportSevereStrat,
                                         title = "Estimacion con severos")
estimateSubreportStrat <- gridExtra::grid.arrange(subreportCritStratPlot,
                                                  subreportSevereStratPlot,
                                                  ncol = 2)

plotSaveFile <- paste("../graficas/8_estimado_subreporte_estratificado",
                      fecha, ".png", sep="")
ggsave(plotSaveFile, estimateSubreportStrat, width=16, height=8, units="cm")


#######################################################
# 8 Plot estimated real dynamics
#######################################################
daysCut <- 7
caseDynamicsFile <- paste("../datos_procesados/1_outcomeDynamics_", 
                          fecha, "_cut.csv", sep="")
caseDynamics <- read.csv(caseDynamicsFile) %>%
  dplyr::rename(., date = day) %>%
  dplyr::mutate(., date = lubridate::as_date(date))

realDynCrit <- merge(caseDynamics, subreportCritStrat, by = "date") %>%
  dplyr::mutate(., caseEst = newCases/estimate, caseEstH = newCases/lower,
                caseEstL = newCases/upper) %>%
  as_tibble(.)

realDynSevere <- merge(caseDynamics, subreportSevereStrat, by = "date") %>%
  dplyr::mutate(., caseEst = newCases/estimate, caseEstH = newCases/lower,
                caseEstL = newCases/upper) %>%
  as_tibble(.)

# remove last 7 days, to prevent distortions from symptoms backlog
realDynCrit <- realDynCrit[c(1:(nrow(realDynCrit)-daysCut)),]
realDynSevere <- realDynSevere[c(1:(nrow(realDynSevere)-daysCut)),]

realDynCritPlot_strat <- plot_estimated_dyn(realDynCrit, meanWindow=3,
                                      title = "Estimacion con criticos")
absSubreportCritPlot_strat <- plot_abs_subreport(realDynCrit, meanWindow=3,
                                      title = "Estimacion con criticos")
realDynSeverePlot_strat <- plot_estimated_dyn(realDynSevere, meanWindow=3,
                                      title = "Estimacion con severos")
absSubreportSeverePlot_strat <- plot_abs_subreport(realDynSevere, meanWindow=3,
                                      title = "Estimacion con severos")

estimatedCasesStrat <- gridExtra::grid.arrange(realDynCritPlot_strat,
                                          realDynSeverePlot_strat,
                                          ncol = 1)
absoluteSubreportStrat <- gridExtra::grid.arrange(absSubreportCritPlot_strat,
                                             absSubreportSeverePlot_strat,
                                             ncol = 1)

ggsave("../graficas/8_estimado_casos_reales_estratificado.png",
       estimatedCasesStrat, width = 16, height = 12, units = "cm")
ggsave("../graficas/8_estimado_casos_no_detectados_estratificado.png",
       absoluteSubreportStrat, width = 16, height = 12, units = "cm")

