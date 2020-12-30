library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(zoo)
library(ggpubr)

plot_literature_outcome <- function(literatureDf, title=NA, escalaLog=TRUE,
                                    colorVec=NA) {
  outcomePlot <- ggplot(literatureDf, aes(x=meanAge, y=outcome,
                                     color=study, fill=study)) +
    geom_line(size=0.6, linetype = "dashed") +
    geom_ribbon(aes(ymin = outcomeL, ymax = outcomeH), alpha = 0.2, colour=NA,
                show.legend=FALSE) +
    xlab("Edad") +
    ylab("% infectados evento") +
    labs(color = "Estudio") +
    theme_bw()
  if (!is.na(colorVec)) {
    outcomePlot <- outcomePlot +
      scale_colour_manual(values=colorVec) +
      scale_fill_manual(values=colorVec)
  }
  if (escalaLog) {
    outcomePlot <- outcomePlot +
    scale_y_continuous(trans = 'log10', labels=scales::comma)
  }
  if (!is.na(title)) {
    outcomePlot <- outcomePlot +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, size=10))
  }
  return(outcomePlot)
}


plot_delay_histogram <- function(delayDf, nBins = 40, title=NA) {
  quantiles <- quantile(delayDf$delay)[2:4]
  delay_fun <- onset2Outcome(quantiles)
  delayProb <- data.frame(delay = c(0:30)) %>%
    dplyr::mutate(., prob = delay_fun(c(0:30)))
  delayPlot <- ggplot(data = delayDf, aes(x=delay)) +
    geom_histogram(aes(y=100*..count../sum(..count..)), bins=nBins) +
    geom_line(data=delayProb, aes(x=delay, y=prob*100)) +
    xlab("Días entre inicio de síntomas y evento") +
    ylab("Porcentaje") +
    theme_bw()
  if (!is.na(title)) {
    delayPlot <- delayPlot +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, size=10))
  }
  return(delayPlot)
}

plot_subreporting <- function(subreportDf, title = NA) {
  subreportPlot <- ggplot(subreportDf, aes(x = date, y = estimate*100)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower*100, ymax = upper*100),
                fill="grey70", alpha=0.6) +
    ylim(0,100) +
    xlab("Fecha") +
    scale_x_date(date_breaks = "months" , date_labels = "%b") +
    ylab("% infecciones reportadas") +
    theme_bw()
    if (!is.na(title)) {
      subreportPlot <- subreportPlot +
        ggtitle(title) +
        theme(plot.title = element_text(hjust = 0.5, size=10))
    }
  return(subreportPlot)
}


plot_estimated_dyn <- function(estimatedDynDf, meanWindow=NA,
                               title=NA) {
  if (!is.na(meanWindow)) {
    estimatedDynDf <- dplyr::mutate(estimatedDynDf,
                  newCases = zoo::rollmean(newCases, k=meanWindow, fill=NA),
                  caseEst = zoo::rollmean(caseEst, k=meanWindow, fill=NA),
                  caseEstH = zoo::rollmean(caseEstH, k=meanWindow, fill=NA),
                  caseEstL = zoo::rollmean(caseEstL, k=meanWindow, fill=NA))
  }
  cols <- c("Reportados" = "black", "Estimados" = "blue")
  estDynPlot <- ggplot(estimatedDynDf, aes(x = date)) +
    geom_line(aes(y=newCases, color="Reportados")) +
    geom_line(aes(y=caseEst, color="Estimados")) +
    geom_ribbon(aes(ymin = caseEstL, ymax = caseEstH),
                fill="grey70", alpha=0.6) +
    xlab("Fecha") +
    scale_x_date(date_breaks = "months" , date_labels = "%b") +
    ylab("Casos nuevos") +
    scale_colour_manual(name=element_blank(),values=cols) +
    theme_bw()
  if (!is.na(title)) {
    estDynPlot <- estDynPlot +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, size=10))
  }
  return(estDynPlot)
}


plot_abs_subreport <- function(estimatedDynDf, meanWindow=NA,
                               title=NA) {
  if (!is.na(meanWindow)) {
    estimatedDynDf <- dplyr::mutate(estimatedDynDf,
                  newCases = zoo::rollmean(newCases, k=meanWindow, fill=NA),
                  caseEst = zoo::rollmean(caseEst, k=meanWindow, fill=NA),
                  caseEstH = zoo::rollmean(caseEstH, k=meanWindow, fill=NA),
                  caseEstL = zoo::rollmean(caseEstL, k=meanWindow, fill=NA),
                  subreportDiff = caseEst-newCases,
                  subreportDiffH = caseEstH-newCases,
                  subreportDiffL = caseEstL-newCases)
  }
  estDynPlot <- ggplot(estimatedDynDf, aes(x = date)) +
    geom_line(aes(y=subreportDiff)) +
    geom_ribbon(aes(ymin = subreportDiffL, ymax =  subreportDiffH),
                fill="grey70", alpha=0.6) +
    xlab("Fecha") +
    scale_x_date(date_breaks = "months" , date_labels = "%b") +
    ylab("Casos no detectados") +
    theme_bw()
  if (!is.na(title)) {
    estDynPlot <- estDynPlot +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, size=10))
  }
  return(estDynPlot)
}


plot_subreport_sim <- function(simulationDf) {
  cols <- c("Subreporte real" = "black", "Estimado" = "blue")
  simulationDf <- dplyr::mutate(simulationDf,
                                dailyRep = newCasesReported/newCasesReal)
  simulationPlot <- ggplot(simulationDf, aes(x = day)) +
        geom_line(aes(y = estimate*100, color = "Estimado"), size=1) + 
        geom_point(aes(y = dailyRep*100, color = "Subreporte real"), size=0.8) + 
        geom_ribbon(aes(ymin = lower*100, ymax = upper*100),
                    fill = "dodgerblue", alpha = 0.3) +
        geom_line(aes(y = realReporting*100, color = "Subreporte real"), size=1) +
        scale_colour_manual(name=element_blank(),values=cols) +
        ylab(" ") +
        xlab("Fecha") +
        scale_x_date(date_breaks = "months" , date_labels = "%b") +
        theme_bw() +
        ylim(c(0,100))
  return(simulationPlot)
}

plot_infections_sim <- function(simulationDf, meanWindow=5,
                                plotEst=FALSE) {
  simulationDf <- dplyr::mutate(simulationDf,
                realSmooth = zoo::rollmean(newCasesReal, k=meanWindow, fill=NA),
                reportedSmooth = zoo::rollmean(newCasesReported,
                                               k=meanWindow, fill=NA),
                estimateSmooth = zoo::rollmean(newCasesReported/estimate,
                                              k=meanWindow, fill=NA))
  cols <- c("Inf reales" = "blue",
            "Inf reportadas" = "black",
            "Hospitalizados" = "#3CAEA3",
            "Estimacion metodo" = "#FF5733")
  simulationPlot <- ggplot(simulationDf, aes(x = day)) +
        geom_line(aes(y = realSmooth,  color = "Inf reales"), size=1) + 
        geom_line(aes(y = reportedSmooth,  color = "Inf reportadas"), size=1) + 
        geom_point(aes(y = newOutcome, color = "Hospitalizados"), size=1) +
        scale_colour_manual(name=element_blank(),values=cols) +
        ylab(" ") +
        xlab(" ") +
        scale_x_date(date_breaks = "months" , date_labels = "%b") +
        theme_bw()
  if (plotEst) {
    simulationPlot <- simulationPlot +
        geom_line(aes(y = estimateSmooth,  color = "Estimacion metodo"), size=1)
  }
  return(simulationPlot)
}

# get performance of plain reporting and of estimate
plot_simulation_fit_summary <- function(simulationDf) {
  simulationDf <- dplyr::mutate(simulationDf,
                                estimated = newCasesReported/estimate)
  result <- group_by(simulationDf, funcInfInd, funcSRInd, nRep) %>%
    summarize(., subreport = mean(newCasesReported/newCasesReal,
                                      na.rm=TRUE),
              subreportAdjusted = mean(estimated/newCasesReal,
                                      na.rm=TRUE),
              subreportCum = sum(newCasesReported)/sum(newCasesReal),
              subreportCumAdjusted = sum(estimated)/sum(newCasesReal),
              subreportAbs = sum(newCasesReported)-sum(newCasesReal),
              subreportAbsAdjusted = sum(estimated)-sum(newCasesReal)) %>%
    tidyr::pivot_longer(., c("subreport", "subreportAdjusted",
                             "subreportCum", "subreportCumAdjusted",
                             "subreportAbs", "subreportAbsAdjusted"), 
                        names_to = "Method",
                        values_to = "Performance")

  meanSRHistogram <- dplyr::filter(result, Method %in% c("subreport",
                                                         "subreportAdjusted")) %>%
    ggplot(., aes(x=Performance, fill=Method)) +
    geom_histogram(position="identity", alpha=0.6) +
    theme_bw() +
    xlab("Subreporte medio en el tiempo") +
    ylab(" ") +
    scale_fill_discrete(name = "Metodo", labels = c("Reportados", "Estimados")) +
    geom_vline(xintercept = 1, linetype="dotted")
  totalSRHistogram <- dplyr::filter(result, Method %in% c("subreportCum",
                                                         "subreportCumAdjusted")) %>%
    ggplot(., aes(x=Performance, fill=Method)) +
    geom_histogram(position="identity", alpha=0.6) +
    theme_bw() +
    xlab("Subreporte global") +
    ylab(" ") +
    scale_fill_discrete(name = "Metodo", labels = c("Reportados", "Estimados")) +
    geom_vline(xintercept = 1, linetype="dotted")
  absoluteSRHistogram <- dplyr::filter(result, Method %in% c("subreportAbs",
                                                         "subreportAbsAdjusted")) %>%
    ggplot(., aes(x=Performance, fill=Method)) +
    geom_histogram(position="identity", alpha=0.6) +
    theme_bw() +
    xlab("Subreporte absoluto") +
    ylab(" ") +
    scale_fill_discrete(name = "Metodo", labels = c("Reportados", "Estimados")) +
    geom_vline(xintercept = 0, linetype="dotted")

  allPlots <- ggpubr::ggarrange(plotlist=list(meanSRHistogram, totalSRHistogram,
                                             absoluteSRHistogram), ncol=1, nrow=3)
  return(allPlots)
}


plot_estimated_dyn_simulation <- function(estimatedDynDf, meanWindow=NA,
                               title=NA) {
  if (!is.na(meanWindow)) {
    estimatedDynDf <- dplyr::mutate(estimatedDynDf,
                  newCases = zoo::rollmean(newCases, k=meanWindow, fill=NA),
                  caseEst = zoo::rollmean(caseEst, k=meanWindow, fill=NA),
                  caseEstH = zoo::rollmean(caseEstH, k=meanWindow, fill=NA),
                  caseEstL = zoo::rollmean(caseEstL, k=meanWindow, fill=NA))
  }
  cols <- c("Reportados" = "black", "Estimados" = "blue")
  estDynPlot <- ggplot(estimatedDynDf, aes(x = date)) +
    geom_line(aes(y=newCases, color="Reportados")) +
    geom_line(aes(y=caseEst, color="Estimados")) +
    geom_ribbon(aes(ymin = caseEstL, ymax = caseEstH),
                fill="grey70", alpha=0.6) +
    xlab("Fecha") +
    scale_x_date(date_breaks = "months" , date_labels = "%b") +
    ylab("Casos nuevos") +
    scale_colour_manual(name=element_blank(),values=cols) +
    theme_bw()
  if (!is.na(title)) {
    estDynPlot <- estDynPlot +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, size=10))
  }
  return(estDynPlot)
}

plot_estimated_outcomes <- function(predictionDf, title=NA, meanWindow=3) {
  predictionDf <- dplyr::mutate(predictionDf,
                                smoothOutcome = zoo::rollmean(newOutcome,
                                                              k=meanWindow,
                                                              fill=NA))
  cols <- c("Observados" = "black", "Esperados" = "red")
  outcomePlot <- ggplot(predictionDf, aes(x=day)) +
    geom_line(aes(y=smoothOutcome, color="Observados")) +
    geom_line(aes(y=expectedOutcomes, color="Esperados")) +
    scale_colour_manual(name=element_blank(),values=cols) +
    scale_x_date(date_breaks = "months" , date_labels = "%b") +
    theme_bw() +
    xlab(" ")
  if (!is.na(title)) {
    outcomePlot <- outcomePlot +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, size=10))
  }
  return(outcomePlot)
}

